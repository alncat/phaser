//(c) 2000-2014 Cambridge University Technical Services Ltd
//All rights reserved
#include <phaser/include/Bin.h>
#include <phaser/src/DataMR.h>
#include <phaser/io/Errors.h>
#include <scitbx/constants.h>
#include <phaser/src/Ensemble.h>
#include <phaser/src/UnitCell.h>
#include <phaser/io/Output.h>
#include <phaser/lib/jiffy.h>
#include <phaser/lib/round.h>
#include <phaser/lib/xyzRotMatDeg.h>
#include <phaser/lib/sphericalY.h>
#include <phaser/lib/aniso.h>
#include <phaser/lib/mean.h>
#include <scitbx/math/erf.h>
#include <cmtzlib.h>
#include <scitbx/math/r3_rotation.h>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace phaser {

DataMR::DataMR(std::string SG_HALL,af::double6 UNIT_CELL,data_refl& REFLECTIONS,Bin& DATABINS,data_norm& SIGMAN,data_outl& OUTLIER,data_composition& COMPOSITION,floatType res1,floatType res2,data_resharp& RESHARP,data_tncs& PTNCS)
          : DataB(SG_HALL,UNIT_CELL,REFLECTIONS,DATABINS,SIGMAN,COMPOSITION,res1,res2,OUTLIER,RESHARP,PTNCS)
{
  debug_flag = false;
  RiceV = 0;
  frac_known_B.resize(NREFL,0);
  frac_search_B.resize(NREFL,0);
  EM_known.resize(NREFL,0);
  EM_search.resize(NREFL,0);
  sum_Esqr_search.resize(NREFL,0);
  totvar_known.resize(NREFL,0);
  totvar_search.resize(NREFL,0);
  totvar_move.resize(NREFL,0);

  //always allocate memory for RICE too
  max_Esqr_search.resize(NREFL);
}

void DataMR::calcKnownB(map_str_float mapEnsScat,std::string search_modlid,floatType searchBfac)
{
  //What is the total molecular weight of the models?
  float1D mapEnsScatKnown(models_known.size());
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    mapEnsScatKnown[k] = mapEnsScat[modlid];
  }
  floatType mapEnsScatSearch(search_modlid.size() ? mapEnsScat[search_modlid]: 0);
  for (unsigned r = 0; r < NREFL; r++)
  {
    frac_known_B[r] = 0;
    frac_search_B[r] = 0;
    floatType rssqr = ssqr(r);
    for (unsigned k = 0; k < models_known.size(); k++)
    {
      std::string modlid = models_known[k].getModlid();
      floatType frac = mapEnsScatKnown[k];
                frac /= models_known[k].getMult();
                frac *= -1 + std::exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
      frac_known_B[r] += frac;
    }
    if (search_modlid.size() && searchBfac)
    {
      floatType frac = mapEnsScatSearch;
                frac *= -1 + std::exp(-2.0*searchBfac*rssqr/4.0);
      frac_search_B[r] += frac; //frac_search_B
    }
  }
}

void DataMR::initSearchMR()
{
  for (unsigned r = 0; r < NREFL; r++)
  {
    totvar_search[r] = 0;
    sum_Esqr_search[r] = 0;
    max_Esqr_search[r] = 0;
  }
}

void DataMR::init()
{
  models_known.resize(0);
  models_perturbRot.resize(0);
  models_perturbTrans.resize(0);
  models_initBfac.resize(0);
  for (unsigned r = 0; r < NREFL; r++)
  {
    totvar_search[r] = 0;
    totvar_known[r] = 0;
    totvar_move[r] = 0;
    sum_Esqr_search[r] = 0;
    max_Esqr_search[r] = 0;
    EM_known[r] = 0;
    EM_search[r] = 0;
    frac_known_B[r] = 0;
    frac_search_B[r] = 0;
  }
}

void DataMR::calcSearchVAR(Ensemble& ens_modlid,bool halfR,floatType searchBfac)
{
  // This function is used for the interpolation of the variance
  // at the identity rotation. It assumes that the variances
  // in the model are spherically symmetric and therefore do not
  // need to be calculated for each rotation separately
  // Numerical differences in the interpolation mean that there
  // will be slightly different values even if the variances are
  // spherically symmetric
  floatType scatFactor = ens_modlid.ensScat(TOTAL_SCAT);
  //rotation only
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ens_modlid.set_sqrt_ssqr(sqrtrssqr,rssqr);
    floatType thisV = ens_modlid.InterpV(); // D^2 for ensemble
              thisV *= scatFactor;                      // becomes sigmaA^2
    if (searchBfac)
              thisV *= exp(-2.0*searchBfac*ssqr(r)/4.0);
              //search does not have MULT factor
    if (PTNCS.use_and_present())
    {
      // halfR=true (look for average orientation of tNCS pair)
      // when called from FRF or BRF, or from BTF or FTF before rotref,
      // false otherwise (BTF or FTF after rotref, RNP)
      if (PTNCS.NMOL > 2 || !halfR)
        thisV *= PTNCS.NMOL*G_Vterm[r];
      else //halfR
      {
        gfun.calcReflTerms(r);
        thisV *= PTNCS.NMOL*gfun.refl_Gsqr_Vterm;
      }
    }
    totvar_search[r] = thisV;
  }
}

void DataMR::calcSearchVAR(std::string MODLID,mrgyre1D& GYRE,MapEnsemble& ensemble,floatType searchBfac)
{
  for (unsigned r = 0; r < NREFL; r++)
    totvar_search[r] = 0;
  for (unsigned g = 0; g < GYRE.size(); g++)
  {
    std::string modlid = GYRE[g].modlid(MODLID);
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    floatType fracScat = ens_modlid->ensScat(TOTAL_SCAT);
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      floatType rssqr = ssqr(r);
      ensemble.set_ssqr(rssqr);
      floatType  thisV = ens_modlid->InterpV();
      thisV *= fracScat;
      if (PTNCS.use_and_present())
        thisV *= PTNCS.NMOL*G_Vterm[r];
      totvar_search[r] += thisV;
    }
  }
}

void DataMR::calcSearchROT(dmat33 ROT,Ensemble& ens_modlid,floatType searchBfac) //halfR default false
{
  const cmplxType cmplxZERO(0.,0.),cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
  dmat33 R = ROT*ens_modlid.Frac2Orth();
  dmat33 Q1 = UnitCell::doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();
  floatType scatFactor = ens_modlid.ensScat(TOTAL_SCAT)/NSYMP;
  //rotation only
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType repsn = epsn(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ens_modlid.set_sqrt_ssqr(sqrtrssqr,rssqr);
    sum_Esqr_search[r] = 0;
    max_Esqr_search[r] = 0;
    EM_search[r] = cmplxZERO;
    Miller1D  rhkl = rotMiller(r);
    if (PTNCS.use_and_present())
    {
      gfun.calcReflTerms(r);
    }
    for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
    {
      if (!duplicate(isym,rhkl))
      {
        dvect3    RotSymHKL = Q1tr*rhkl[isym];
        cmplxType thisE = ens_modlid.InterpE(RotSymHKL);
        if (searchBfac)
                  thisE *= std::exp(-searchBfac*rssqr/4.0);
              //search does not have MULT factor
        if (PTNCS.use_and_present())
        {
          floatType theta  = rhkl[isym][0]*PTNCS.TRA.VECTOR[0] +
                             rhkl[isym][1]*PTNCS.TRA.VECTOR[1] +
                             rhkl[isym][2]*PTNCS.TRA.VECTOR[2];
          cmplxType ptncs_scat = cmplxONE; //iMOL=0
          for (floatType iMOL = 1; iMOL < PTNCS.NMOL; iMOL++)
            ptncs_scat += std::exp(TWOPII*iMOL*theta);
          ptncs_scat *= gfun.refl_G[isym]; //FRF,BRF rotation unknown
          thisE *= ptncs_scat;
        }
        floatType thisEsqr = thisE.real()*thisE.real() + thisE.imag()*thisE.imag();
                  thisEsqr *= repsn;
                  thisEsqr *= scatFactor;
        sum_Esqr_search[r] += thisEsqr;
        max_Esqr_search[r]  = std::max(max_Esqr_search[r],thisEsqr);
      }
    }
  }
}

// Calculate the non-varying parameters for a translation function
// Variance values interpolated from a molecular transform are
// precalculated and passed to the function. The solvent content is
// calculated from the search model.
void DataMR::calcSearchTRA1(dmat33 ROT,Ensemble& ens_modlid)
{
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    sum_Esqr_search[r] = 0;
    max_Esqr_search[r] = 0;
  }
}

// Calculate the varying parameters for a translation function
// E values interpolated from a molecular transform are
// precalculated and passed to the function
void DataMR::calcSearchTRA2(dvect3 TRA,bool halfR,floatType searchBfac,cmplx1D* interpE,cmplx1D* interpE2)
{
  const cmplxType cmplxZERO(0.,0.),cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);

  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType sqrt_repsn = sqrt_epsn[r];
    floatType rssqr = ssqr(r);
    EM_search[r] = cmplxZERO;
    Miller1D rhkl = rotMiller(r);
    if (PTNCS.use_and_present())
    {
      gfun.calcReflTerms(r);
    }
    for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
    {
      if (!duplicate(isym,rhkl))
      {
        cmplxType thisE = (*interpE)[r*SpaceGroup::NSYMP+isym];
                  thisE *= sqrt_repsn;
                  thisE *= dphi(r,isym,rhkl[isym],TRA);
        if (searchBfac)
                  thisE *= std::exp(-searchBfac*rssqr/4.0);
        if (PTNCS.use_and_present())
        {
          floatType theta  = rhkl[isym][0]*PTNCS.TRA.VECTOR[0] +
                             rhkl[isym][1]*PTNCS.TRA.VECTOR[1] +
                             rhkl[isym][2]*PTNCS.TRA.VECTOR[2];
          if (PTNCS.NMOL > 2)
          {
            cmplxType ptncs_scat = cmplxONE; //iMOL=0
            for (floatType iMOL = 1; iMOL < PTNCS.NMOL; iMOL++)
              ptncs_scat += std::exp(TWOPII*iMOL*theta);
            thisE *= ptncs_scat;
            //Gfunction term unity since no rotation
          }
          else if (halfR) //second molecule comes from doubling at same rot
          {
            thisE *= cmplxONE + std::exp(TWOPII*theta);
            thisE *= gfun.refl_G[isym];
          }
          else //second molecule is refined position
          {
            cmplxType thisE2 = (*interpE2)[r*SpaceGroup::NSYMP+isym];
                      thisE2 *= sqrt_repsn;
                      thisE2 *= dphi(r,isym,rhkl[isym],TRA);
                      thisE2 *= std::exp(TWOPII*theta);
            thisE += thisE2;
          }
        }
        EM_search[r] += thisE;
      }
    }
  }
}

void DataMR::calcSearchTRA2(dvect3 FRAC,std::string MODLID,MapEnsemble& ensemble,mrgyre1D& GYRE,cmplx1D* interpE)
{
  dmat33 ROT,PRtr,Rperturb;
  dvect3 TRA,CMori,CMrot,CMperturb;
  for (unsigned g = 0; g < GYRE.size(); g++)
  {
    std::string modlid = GYRE[g].modlid(MODLID);
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    PRtr = ens_modlid->principalRtr();
    CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = euler2matrixDEG(GYRE[g].EULR);
    CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform
    //orientation (for refinement)
    Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,0)*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented ensemble)
    ROT = ROT*PRtr;

    TRA = GYRE[g].FRAC + FRAC;
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of ensemble)
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA = UnitCell::doOrth2Frac(eOrthT);

    floatType sqrt_scatFactor = std::sqrt(ens_modlid->ensScat(TOTAL_SCAT)/NSYMP);
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      Miller1D rhkl = rotMiller(r);
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          cmplxType thisE = (*interpE)[(r*SpaceGroup::NSYMP+isym)*GYRE.size()+g];
                    thisE *= dphi(r,isym,rhkl[isym],TRA);
          //no G function correction (gfun), models in separate (real) rotations
          EM_known[r] += thisE;
        }
      }
    }
  }
}

void DataMR::calcKnownVAR(MapEnsemble& ensemble)
{
  // This function is used for the interpolation of the variance
  // at the identity rotation. It assumes that the variances
  // in the model are spherically symmetric and therefore do not
  // need to be calculated for each rotation separately
  // Numerical differences in the interpolation mean that there will
  // be slightly different values even if the variances are
  // spherically symmetric
  dmat33 ROT(1,0,0,0,1,0,0,0,1);
  for (int r = 0; r < NREFL; r++)
    totvar_known[r] = 0;

  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType rssqr = ssqr(r);
    ensemble.set_ssqr(rssqr);
    //save this interpolation, so not doing multiple times for multiple copies
    map_str_float Vmap;
  //stored  if (PTNCS.use_and_present()) gfun.calcReflTerms(r);
    for (unsigned k = 0; k < models_known.size(); k++)
    {
      std::string modlid = models_known[k].getModlid();
      floatType thisV(0);
      if (Vmap.find(modlid) != Vmap.end())
      {
        thisV = Vmap[modlid];
      }
      else
      {
        Ensemble* ens_modlid = &ensemble.find(modlid)->second;
        thisV = ens_modlid->InterpV();
        thisV *= ens_modlid->ensScat(TOTAL_SCAT);
        Vmap[modlid] = thisV;
      }
      thisV *= exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
      thisV /= models_known[k].getMult();
      if (PTNCS.use_and_present())
                thisV *= G_Vterm[r]; //not correct, B-factor differences
      totvar_known[r] += thisV;
    }
  }
}

void DataMR::calcKnownVAR(MapEnsemble& ensemble, map_str_pdb& PDB)
{
  // This function is used for the interpolation of the variance
  // at the identity rotation with presumed bfactors. It assumes that the variances
  // in the model are spherically symmetric and therefore do not
  // need to be calculated for each rotation separately
  // Numerical differences in the interpolation mean that there will
  // be slightly different values even if the variances are
  // spherically symmetric
  dmat33 ROT(1,0,0,0,1,0,0,0,1);
  for (int r = 0; r < NREFL; r++)
    totvar_known[r] = 0;

  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType rssqr = ssqr(r);
    ensemble.set_ssqr(rssqr);
    //save this interpolation, so not doing multiple times for multiple copies
    map_str_float Vmap;
  //stored  if (PTNCS.use_and_present()) gfun.calcReflTerms(r);
    for (unsigned k = 0; k < models_known.size(); k++)
    {
      std::string modlid = models_known[k].getModlid();
      int mzero(0);
      Molecule& pdb_mol = PDB[modlid].Coords(mzero);
      floatType thisV(0);
      if (Vmap.find(modlid) != Vmap.end())
      {
        thisV = Vmap[modlid];
      }
      else
      {
        Ensemble* ens_modlid = &ensemble.find(modlid)->second;
        thisV = ens_modlid->InterpV();
        //thisV *= ens_modlid->ensScat(TOTAL_SCAT);
        Vmap[modlid] = thisV;
      }
      thisV *= exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
      thisV /= models_known[k].getMult();
      if (PTNCS.use_and_present())
                thisV *= G_Vterm[r]; //not correct, B-factor differences
      floatType tmp_var(0);
      for(unsigned i = 0; i < pdb_mol.nAtoms(0); i++)
      {
        cctbx::eltbx::tiny_pse::table elementi(pdb_mol.getElement(0,i));
        floatType B(pdb_mol.getIsoB(mzero,i));
        tmp_var += elementi.atomic_number()*elementi.atomic_number()*exp(-2.0*B*rssqr/4.0);
      }
      totvar_known[r] += tmp_var/TOTAL_SCAT*thisV;
    }
  }
}


void DataMR::addInterpV(float1D& interpV)
{
  PHASER_ASSERT(interpV.size() == NREFL);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    totvar_known[r] += interpV[r];
}

void DataMR::addInterpE(cmplx1D& interpE)
{
  PHASER_ASSERT(interpE.size() == NREFL);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    EM_known[r] += interpE[r];
}

void DataMR::setInterpV(float1D& interpV)
{
  PHASER_ASSERT(interpV.size() == NREFL);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    totvar_known[r] = interpV[r];
}

void DataMR::setInterpE(cmplx1D& interpE)
{
  PHASER_ASSERT(interpE.size() == NREFL);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    EM_known[r] = interpE[r];
}

cmplx1D DataMR::getInterpE()
{
  return EM_known;
}

float1D DataMR::getInterpV()
{
  return totvar_known;
}

// Calculate the parameters for MR
// E values are interpolated from a molecular transform of the
// ensemble of structures.
// If shift (perturbation) parameters are changed, also have to change
// RefineMR::applyShift and RefineMR::getRefinePars
void DataMR::calcKnown(MapEnsemble& ensemble)
{
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    EM_known[r] = cmplxType(0,0);
    sum_Esqr_search[r] = 0; //because there is no sum_Esqr_known
    max_Esqr_search[r] = 0;
  }

  dmat33 ROT,PRtr,Rperturb;
  dvect3 TRA,CMori,CMrot,CMperturb;

  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    PRtr = ens_modlid->principalRtr();
    CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = models_known[k].getR();
    CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform
    //orientation (for refinement)
    Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[k][dir])*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented ensemble)
    ROT = ROT*PRtr;

    TRA = models_known[k].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[k];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of ensemble)
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA = UnitCell::doOrth2Frac(eOrthT);

    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr = Q1.transpose();

    floatType sqrt_scatFactor = std::sqrt(ens_modlid->ensScat(TOTAL_SCAT)/NSYMP);
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      floatType repsn = epsn(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      Miller1D rhkl = rotMiller(r);
      floatType terms  = sqrt_scatFactor;
                terms *= sqrt_epsn[r];
                terms *= std::exp(-models_known[k].getBfac()*rssqr/4.0);
                terms /= models_known[k].getMult();
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3    RotSymHKL = Q1tr*rhkl[isym];
          cmplxType thisE =  ens_modlid->InterpE(RotSymHKL);
                    thisE *= terms;
                    thisE *= dphi(r,isym,rhkl[isym],TRA);
          //no G function correction (gfun), models in separate (real) rotations
          EM_known[r] += thisE;
        }
      }
    }
  }
}


floatType DataMR::Rice()
{
//  Compute likelihood score for all x and r.  If there are any
//  relative phase ambiguities,use Wilson or Rice function approximation.
//  For the Rice function approximation,the largest single contribution
//  is taken as a partial structure factor and the variance is incremented
//  by the remaining contributions
//
//  The likelihood over all the data sets is the product of the likelihood for
//  the individual data sets (sum of the logs)
//  This does not take into account the correlation between them

  floatType totalLL(0);

  bool first_r(true);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    totalLL += Rice_refl(r,first_r);
    first_r = false;
  } // r loop
  return totalLL;
}

floatType DataMR::Rice_refl(int r,bool checkV)
{
  floatType reflLL(0);
  floatType fracB = 1 + frac_known_B[r] + frac_search_B[r];
  cmplxType EM(EM_known[r]+EM_search[r]);
  floatType EM_real = EM.real();
  floatType EM_imag = EM.imag();
  floatType phasedEsqr = EM_real*EM_real+EM_imag*EM_imag;
  floatType sumEsqr = phasedEsqr + sum_Esqr_search[r];
  floatType maxEsqr = std::max(phasedEsqr,max_Esqr_search[r]);
  // Note that, after interpolation of molecular transform, EM is systematically underestimated
  // by a few percent

  floatType V(0);
  V = PTNCS.EPSFAC[r];
  V -= totvar_known[r];
  V -= totvar_search[r];
  V += sumEsqr;
  V -= maxEsqr;
  V /= fracB;

  bool rcent = cent(r);
  floatType E = F[r]/sqrt_epsnSigmaN[r];
  floatType Esqr = fn::pow2(E);
  floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
  if (!rcent) wtVarE *= 2.;
  V += wtVarE; //centric correction has already been included
  RiceV = V;
  if (checkV && V <= 0) { return 0; }

  floatType FCsqr = maxEsqr/fracB; //fc^2/fracB
  PHASER_ASSERT(FCsqr >= 0);
  floatType FC = std::sqrt(FCsqr);
  floatType X(0);
//#define PHASER_NEGVAR_DEBUG
#ifdef PHASER_NEGVAR_DEBUG
  if (V <= 0) goto debug_code;
#else
  PHASER_ASSERT(V > 0);
#endif
  X = 2.0*E*FC/V;
  reflLL = -(std::log(V)+(Esqr + FCsqr)/V);
  PHASER_ASSERT(X >= 0);
  if (rcent)
  {
    X /= 2.0;
    reflLL = reflLL/2.0 + m_alogch.getalogch(X);
  }
  else reflLL += m_alogchI0.getalogI0(X); //acentric
  //if (!debug_flag) return reflLL;
  return reflLL;

//debugging
debug_code:
#ifdef PHASER_NEGVAR_DEBUG
    {
      if (debug_flag && checkV)
      {
        std::cout << "test " << "reflection # " << r << std::endl;
        for (int k = 0; k < models_known.size(); k++)
        {
          std::cout << "test known #" << k <<  " bfac=" << models_known[k].getBfac() << std::endl;
        }
        std::cout << "test " << "reflLL " <<  reflLL << std::endl;
        std::cout << "test " << "reso " <<  reso(r) << std::endl;
        int s = rbin(r);
        std::cout << "test " << "bin #" <<  s << " of " << bin.numbins() <<  std::endl;
        if (PTNCS.use_and_present())
        {
          std::cout << "calc PTNCS " << std::endl;
          bool b(false);
          unsigned rr = r;
          gfun.calcReflTerms(rr);
          std::cout << "   test G Vterm: " << gfun.refl_G_Vterm << std::endl;
          std::cout << "   test stored G Vterm: " << G_Vterm[r] << std::endl;
          std::cout << "   test Gsqr Vterm: " << gfun.refl_Gsqr_Vterm << std::endl;
          for (unsigned isym = 0; isym < gfun.refl_G.size(); isym++)
            std::cout << "   test G isym:" << isym << " " << gfun.refl_G[isym] << std::endl;
        }
        std::cout << "test epsfac: " << PTNCS.EPSFAC[r] << std::endl;
        std::cout << "test nmol: " << PTNCS.NMOL << std::endl;
        std::cout << "test " << "epsn " <<  epsn(r) << std::endl;
        std::cout << "test " << "EM_known " << EM_known[r] << std::endl;
        std::cout << "test " << "EM_search " << EM_search[r] << std::endl;
        std::cout << "test " << "EM_known+EM_search " << EM_known[r] + EM_search[r] << std::endl;
        std::cout << "test " << "FC " << FC << " FCsqr " << FCsqr << std::endl;
        std::cout << "test " << "E " << E << " Esqr " << Esqr << std::endl;
        std::cout << "test " << "phasedEsqr " << phasedEsqr << std::endl;
        std::cout << "test " << "max_Esqr_search " <<max_Esqr_search[r] <<  std::endl;
        std::cout << "test " << "sum_Esqr_search " <<sum_Esqr_search[r] <<  std::endl;
        std::cout << "test " << "ptncs " << PTNCS.EPSFAC[r] <<std::endl;
        std::cout << "test " << "-totvar_known " << -totvar_known[r] <<std::endl;
        std::cout << "test " << "-totvar_search " << -totvar_search[r] <<  std::endl;
        std::cout << "test " << "-totvar_known-totvar_search " << -totvar_known[r]-totvar_search[r] <<  std::endl;
        std::cout << "test " << "frac_known_B " << frac_known_B[r] <<  std::endl;
        std::cout << "test " << "frac_search_B " << frac_search_B[r] <<  std::endl;
        std::cout << "test " << "frac_known_B+frac_search_B " << frac_known_B[r]+frac_search_B[r] <<  std::endl;
        std::cout << "test " << "fracB " << fracB <<  std::endl;
        std::cout << "test " << "sumEsqr " << sumEsqr << std::endl;
        std::cout << "test " << "-maxEsqr " << -maxEsqr << std::endl;
        std::cout << "test " << "V =" << V << std::endl;
        std::cout << std::endl;
      }
      PHASER_ASSERT(V > 0);
      PHASER_ASSERT(FCsqr >= 0);
    }//debug
#endif
  return reflLL;
}

pair_flt_flt DataMR::ScaleFC()
{
  // Get parameters to put everything but anisotropy back into SigmaN
  floatType isotropicB = AnisoBeta2IsoB(sigmaN.ANISO,A(),B(),C(),cosAlpha(),cosBeta(),cosGamma());
  dmat6 isoBasAnisoBeta = IsoB2AnisoBeta(isotropicB,aStar(),bStar(),cStar(),cosAlphaStar(),cosBetaStar(),cosGammaStar());

  af_float F_ISO = getCorrectedF();
  PHASER_ASSERT(F_ISO.size() == NREFL);

  floatType sum(0),sumx(0),sumx2(0),sumy(0),sumxy(0);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType sigmaA_Ecalc(std::abs(EM_known[r]+EM_search[r]));
    floatType FC(sigmaA_Ecalc*sqrt_epsnSigmaN[r]);
    floatType W = 1;
    if (FC > 0 && F_ISO[r] > 0)
    {
 //scale FC to F (scale applies to F2)
      floatType Y = std::log(F_ISO[r]/FC);
      floatType s2 = ssqr(r);
      sum += W;
      sumx += W*s2;
      sumx2 += W*s2*s2;
      sumy += W*Y;
      sumxy += W*s2*Y;
    }
  }
  PHASER_ASSERT(sum*sumx2-fn::pow2(sumx));
  floatType slope = (sum*sumxy-(sumx*sumy))/(sum*sumx2-fn::pow2(sumx));
  floatType intercept = (sumy - slope*sumx)/sum;
  floatType scaleK = exp(intercept);
  floatType scaleB = -4*slope;
  return std::pair<floatType,floatType>(scaleK,scaleB);
}

floatType DataMR::Rfactor()
{
  // Get parameters to put everything but anisotropy back into SigmaN
  floatType isotropicB = AnisoBeta2IsoB(sigmaN.ANISO,A(),B(),C(),cosAlpha(),cosBeta(),cosGamma());
  dmat6 isoBasAnisoBeta = IsoB2AnisoBeta(isotropicB,aStar(),bStar(),cStar(),cosAlphaStar(),cosBetaStar(),cosGammaStar());

  af_float F_ISO = getCorrectedF();

  std::pair<floatType,floatType> scale = ScaleFC();
  floatType scaleK(scale.first);
  floatType scaleB(scale.second);
  long double totalR_numerator(0),totalR_denominator(0);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType sigmaA_Ecalc(std::abs(EM_known[r]+EM_search[r]));
    floatType FC(sigmaA_Ecalc*sqrt_epsnSigmaN[r]);
    totalR_numerator += std::fabs(F_ISO[r]-FC*scaleK*exp(-scaleB*ssqr(r)/4));
    totalR_denominator += F_ISO[r];
  } // r loop
  return (totalR_denominator > 0.) ? 100*totalR_numerator/totalR_denominator : -1.0 ;
}

// Linear algebra for the interpolation routines:
// ----------------------------------------------
//     PR          is the rotation matrix applied to molecule from PDB to orient
//                 along principal axes before computing molecular transform (MT)
//     PRtr        transpose (inverse) of P (=MTrT)
//     PT          is the translation that centres PDB molecule on origin
//                 after rotating by PR
//     Xpdb        are the original Angstrom coordinates from PDB file(s).
//     Xmol        are the orthogonal Angstrom coordinates of molecule in orientation
//                 of the MT
//
//     Xmol = PR*Xpdb + PT,or
//     Xpdb = PRtr*Xmol - PRtr*PT,
//
//     xmol        are the fractional coordinates in MT cell of Xmol
//     Frac2OrthMT is the orthogonalising matrix for the MT cell
//
//     Xmol = Frac2OrthMT * xmol
//
//     R           the applied test rotation
//     T           the applied test translation in Angstroms
//     t           the applied test translation in fractional coordinates
//                 (this is the translation given,not T)
//     Orth2FracX  is the fractionalising matrix for the crystal cell
//
//     t = Orth2FracX*T
//
//     Xx          are orthogonal Angstrom coordinates in crystal,after rotating
//                 and translating PDB coordinates by R and T
//
//     Xx = R*Xpdb + T
//
//     xx          are fractional coordinates in crystal
//
//     To interpolate molecular transform,we need to relate xx to xmol
//
//     xx = Orth2FracX*Xx
//        = Orth2FracX*(R*Xpdb+T)
//        = Orth2FracX*R*Xpdb + Orth2FracX*T
//        = Orth2FracX*R*Xpdb + t
//        = Orth2FracX*R*(PRtr*Xmol-PRtr*PT) + t
//        = Orth2FracX*R*PRtr*Xmol - Orth2FracX*R*PRtr*PT + t
//        = Orth2FracX*R*PRtr*Frac2OrthMT*xmol - Orth2FracX*R*PRtr*PT + t
//
//     Interpolation
//     indexes into MT by multiplying crystal HKL by transpose of Q1
//       Q1 = Orth2FracX*R*PRtr*Frac2OrthMT
//
//     corrects phase by adding 2*Pi*(HKL.v1)
//       v1 = t - Orth2FracX*R*PRtr*PT
//

void DataMR::initKnownMR(mr_set MRSET,MapEnsemble& ensemble)
{
// this loads input set into the models arrays
  models_known.clear();
  models_perturbRot.clear();
  models_perturbTrans.clear();
  models_initBfac.clear();
  for (unsigned k = 0 ; k < MRSET.KNOWN.size(); k++)
  {
    models_known.push_back(MRSET.KNOWN[k]);
    models_perturbRot.push_back(dvect3(0,0,0));
    models_perturbTrans.push_back(dvect3(0,0,0));
    models_initBfac.push_back(MRSET.KNOWN[k].getBfac());
  }
  map_str_float mapEnsScat = ensemble.ensScat(TOTAL_SCAT);
  float1D mapEnsScatKnown(models_known.size());
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    mapEnsScatKnown[k] = mapEnsScat[modlid];
  }
//now convert the non-fractional coordinates to fractional
  for (unsigned k = 0; k < models_known.size(); k++)
    if (!models_known[k].isFrac() ) //order important
    {
      dvect3 tmp = models_known[k].getT();
      dvect3 tmp2 = UnitCell::doOrth2Frac(tmp);
      models_known[k].setFracT(tmp2);
    }
//----------- Bfactors
  //start by putting values in range
  floatType minB(0);
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    if (models_known[k].getBfac() < lowerB()+tolB())
      models_known[k].setBfac(lowerB()+tolB());
    if (models_known[k].getBfac() > upperB()-tolB())
      models_known[k].setBfac(upperB()-tolB());
    minB = std::min(minB,models_known[k].getBfac());
  }
  //the limitB does not take account of the correlations between the B-factors
  //final B-factor may be outside this range when more than one in ASU
  //correlations accounted for in getMaxDistSpecial
  floatType tot_var(0);
  floatType kssqr = -2*fn::pow2(1/fullHiRes())/4;
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    floatType thisV = 1; //variance max
              thisV *= mapEnsScatKnown[k];
              thisV *= std::exp(kssqr*models_known[k].getBfac());
    tot_var += thisV;
  }
  //1.0e-06 is too small in some cases (large number of copies)
  floatType tol = 1.0e-06; //tolerance, otherwise V can be fractionally negative
  floatType V(1);
  float1D Bfac(models_known.size(),0); //store initial B
  for (unsigned k = 0; k < models_known.size(); k++)
    Bfac[k] =  models_known[k].getBfac();
  do {
  floatType tot_var_tol = tot_var + tol;
  floatType Bdiff = std::log(1/tot_var_tol)/kssqr;
  if (Bdiff > 0)
  {
    for (unsigned k = 0; k < models_known.size(); k++)
      models_known[k].setBfac(Bfac[k]+Bdiff);
  }
  //put upper B back in range if it overshot
  for (unsigned k = 0; k < models_known.size(); k++)
    if (models_known[k].getBfac() > upperB()-tolB())
      models_known[k].setBfac(upperB()-tolB());
  //final check for limits
  V = 1;  //V = PTNCS.EPSFAC[r] - totvar_known;
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    floatType thisV = 1; //variance max InterpV
              thisV *= mapEnsScatKnown[k];
              thisV *= std::exp(kssqr*models_known[k].getBfac());
    V -= thisV;
  }
  tol += 1.0e-06; //ready for next time, if unsuccessful
  if (tol > 1) break;
  }
  while ( V <= 0 );
  if (models_known.size())
  {
    ensemble.setup_vrms(MRSET.DRMS);
    bool halfR(false);
    initGfunction(halfR);//init for Known
    float1D Bfac(models_known.size());
    for (unsigned k = 0; k < models_known.size(); k++)
      Bfac[k] = models_known[k].getBfac();
    int nloop(0);
//paranoia, I think if all the Bfactors are positive then the variance has to be OK
    int1D selr = get_refl_negative_variance();
    int ir(0);
    while (negative_variance(ensemble,selr,ir) && -minB < upperB())
    {
      nloop++;
      for (unsigned k = 0; k < models_known.size(); k++)
        models_known[k].setBfac(std::min(upperB(),Bfac[k]-minB));
      minB--;
    }
  }
}

void DataMR::deleteLastKnown()
{
  models_known.pop_back();
  models_perturbRot.pop_back();
  models_perturbTrans.pop_back();
  models_initBfac.pop_back();
}

void DataMR::logModls(outStream where,Output& output,std::string added_txt)
{
  output.logTab(1,where,added_txt + "Known MR solutions");
  if (!models_known.size())
    output.logTab(1,where,"(empty solution set - no components)");
  else
  {
    output.logTab(1,where,"SOLU SPAC " + SpaceGroup::spcgrpname());
    for (unsigned k = 0; k < models_known.size(); k++)
      output.logTab(0,where,models_known[k].logfile(1));
  }
  output.logBlank(where);
}

cmplx1D  DataMR::getInterpE(dmat33 ROT,Ensemble& ens_modlid)
{
  cmplx1D interpE;
  dmat33 R = ROT*ens_modlid.Frac2Orth();
  dmat33 Q1 = UnitCell::doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();
  interpE.resize(SpaceGroup::NSYMP*NREFL);
  floatType sqrt_scatFactor = std::sqrt(ens_modlid.ensScat(TOTAL_SCAT)/NSYMP);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ens_modlid.set_sqrt_ssqr(sqrtrssqr,rssqr);
    for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
    {
      cmplxType thisE = ens_modlid.InterpE(Q1tr*rotMiller(r,isym));
                thisE *= sqrt_scatFactor;
      //store InterpE*fraction_scattering
      interpE[r*SpaceGroup::NSYMP+isym] = thisE;
    }
  }
  return interpE;
}
//Calculate variance for protein with atoms of different b factors

void DataMR::calcSearchVAR(data_pdb& PDB)
{
  for(int r = 0; r < NREFL; r++)
    totvar_move[r] = 0;
  int mzero(0);
  Molecule& pdb_mol = PDB.Coords(mzero);
  for(int r = 0; r < NREFL; r++)
  {
    if(selected[r])
    {

      for(unsigned i = 0; i < pdb_mol.nAtoms(0); i++)
      {
        cctbx::eltbx::tiny_pse::table elementi(pdb_mol.getElement(0,i));
        floatType B(pdb_mol.getIsoB(mzero,i));
        totvar_move[r] += elementi.atomic_number()*elementi.atomic_number()*exp(-2.0*B*ssqr(r)/4.0);
      }
    }
  }
}

/*void DataMR::calcRotVAR(MapEnsemble& ensemble)
{
  for (int r = 0; r < NREFL; r++)
    totvar_move[r] = 0;
}*/

//the variance considering atomic b factor.
void DataMR::calcMoveVAR(data_pdb& PDB){
  std::vector<std::pair<int,int> > MODEL(0);
  int pdb_num = 0;
  for(int ipdb = 0; ipdb < PDB.nMols(); ipdb++)
  {
    Molecule& pdb_mol = PDB.Coords(ipdb);
    for(unsigned m = 0; m < pdb_mol.nModels(); m++)
    {
      MODEL.push_back(std::pair<int,int>(ipdb,m));
    }
  }
  int model_size = boost::numeric_cast<int>(MODEL.size());
  for(int iewald = 0; iewald < model_size; iewald++)
  {
    int ipdb = MODEL[iewald].first;
    int m = MODEL[iewald].second;
    Molecule& pdb_mol = PDB.Coords(ipdb);
    for(unsigned i = 0; i < pdb_mol.nAtoms(m); i++)
    {
      cctbx::eltbx::tiny_pse::table elementi(pdb_mol.getElement(m,i));
      floatType b_fact = pdb_mol.getIsoB(m,i);
      for(unsigned r = 0; r < NREFL; r++)
        if (selected[r])
        {
          totvar_move[r] += elementi.atomic_number()*elementi.atomic_number()*exp(-2.0*b_fact*ssqr(r)/4.0);
        }
    }
  }

  for (unsigned r = 0; r < NREFL; r++)
    if(selected[r])
    {
      totvar_move[r] /= model_size;
    }

}
//The expected structure factor with a known configuration and correctly oriented,
//but with unknown position.

void DataMR::calcMoveVAR(dmat33 ROT, data_pdb& PDB)
{
  //for(int r = 0; r < NREFL; r++)
  //  totvar_move[r] = 0;
  dmat33 Q1 = UnitCell::doOrth2Frac(ROT);
  dmat33 Q1tr = Q1.transpose();
  //float1D EFcalc(NREFL,0);
  std::vector<std::pair<int,int> > MODEL(0);
  int pdb_num = 0;
  for (int ipdb = 0; ipdb < PDB.nMols(); ipdb++)
  {
    Molecule& pdb_mol = PDB.Coords(ipdb);
    for (unsigned m = 0; m < pdb_mol.nModels(); m++)
    {
      MODEL.push_back(std::pair<int,int>(ipdb,m));
    }
  }
  int model_size = boost::numeric_cast<int>(MODEL.size());
  //std::cout << model_size <<std::endl;
  for(int iewald = 0; iewald < model_size; iewald++)
  {
    int ipdb = MODEL[iewald].first;
    int m = MODEL[iewald].second;
    Molecule& pdb_mol = PDB.Coords(ipdb);
    for(unsigned i = 0; i < pdb_mol.nAtoms(m); i++)
    {
      dvect3 orthXYZi = pdb_mol.getOrthXYZ(m, i);
      cctbx::eltbx::tiny_pse::table elementi(pdb_mol.getElement(m,i));
      for(unsigned j = 0; j < i; j++)
      {
        if(pdb_mol.getResNum(m,i) == pdb_mol.getResNum(m,j))
        {
          dvect3 orthXYZj = pdb_mol.getOrthXYZ(m, j);
          cctbx::eltbx::tiny_pse::table elementj(pdb_mol.getElement(m,j));
          for(unsigned r = 0; r < NREFL; r++)
            if (selected[r])
            {
              floatType dphi_arg =(Q1*(orthXYZi - orthXYZj))*miller[r];
              dphi_arg *= scitbx::constants::two_pi;
              totvar_move[r] += 2*elementi.atomic_number()*elementj.atomic_number()*std::cos(dphi_arg);
            }
        }
      }
    }
    /*for(unsigned i = 0; i < pdb_mol.nAtoms(m); i++)
    {
      cctbx::eltbx::tiny_pse::table element(pdb_mol.getElement(m,i));
      for(unsigned r = 0; r < NREFL; r++)
        if(selected[r])
          totvar_move[r] += scitbx::fn::pow2(element.atomic_number());
    }*/
  }

  /*std::string fname =  "move_var.txt";
  std::ofstream myFTFfile(fname.c_str(), std::ios::out );*/
  for (unsigned r = 0; r < NREFL; r++)
    if(selected[r])
    {
      totvar_move[r] /= model_size;
      //myFTFfile << miller[r][0] <<" " << miller[r][1] << " " << miller[r][2] << " " << totvar_move[r] <<std::endl;
    }
  //myFTFfile.close();

  /*fname = "tot_move_var.txt";
  std::ofstream myMovefile(fname.c_str(), std::ios::out );

  floatType tot_sph(0);
  for(int iewald = 0; iewald < model_size; iewald++)
  {
    int ipdb = MODEL[iewald].first;
    int m = MODEL[iewald].second;
    Molecule& pdb_mol = PDB.Coords(ipdb);
    for(unsigned i = 0; i < pdb_mol.nAtoms(m); i++)
    {
      cctbx::eltbx::tiny_pse::table element(pdb_mol.getElement(m,i));
      tot_sph += scitbx::fn::pow2(element.atomic_number());
    }
  }
  tot_sph /= model_size;
  for(unsigned r = 0; r < NREFL; r++)
    if(selected[r])
    {
      totvar_move[r] += tot_sph;
      //myMovefile << miller[r][0] << " " << miller[r][1] << " " << miller[r][2]  << " " << totvar_move[r] << std::endl;
    }*/

  //myMovefile.close();
}

cmplx1D DataMR::getInterpE(std::string MODLID,mrgyre1D& GYRE,MapEnsemble& ensemble,floatType searchBfac)
{
  dmat33 ROT,PRtr,Rperturb;
  dvect3 TRA,CMori,CMrot,CMperturb;
  cmplx1D interpE;
  interpE.resize(SpaceGroup::NSYMP*NREFL*GYRE.size());
  for (int g = 0; g < GYRE.size(); g++)
  {
    std::string modlid = GYRE[g].modlid(MODLID);
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    PRtr = ens_modlid->principalRtr();
    CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = euler2matrixDEG(GYRE[g].EULR);
    CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform
    //orientation (for refinement)
    Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,0)*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented ensemble)
    ROT = ROT*PRtr;
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr = Q1.transpose();
    floatType sqrt_scatFactor = std::sqrt(ens_modlid->ensScat(TOTAL_SCAT)/NSYMP);
    const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      floatType repsn = epsn(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      Miller1D rhkl = rotMiller(r);
      floatType terms  = sqrt_scatFactor;
                terms *= sqrt_epsn[r];
                terms *= std::exp(-searchBfac*rssqr/4.0);
                terms /= 1;//Mult
      for (unsigned isym = 0; isym < NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3 RotSymHKL = Q1tr*rhkl[isym];
          cmplxType thisE = ens_modlid->InterpE(RotSymHKL);
                    thisE *= terms;
          interpE[(r*SpaceGroup::NSYMP+isym)*GYRE.size()+g] = thisE;
         //store InterpE*fraction_scattering
        }
      }
    }
  }
  return interpE;
}

cmplx3D DataMR::getELMNxR2Serial(std::string TARGET,Output& output,int LMAX)
{
  output.logBlank(VERBOSE);
  output.logTab(1,LOGFILE,"Elmn for Data");
  if (TARGET != "CROWTHER" && TARGET != "LERF1")
  throw PhaserError(FATAL,"Unknown target " + TARGET);

  //function return
  cmplx3D elmn;

  int ZSYMM(0);

  int axis(highOrderAxis());
  if (axis==1) ZSYMM=orderX();
  if (axis==2) ZSYMM=orderY();
  if (axis==3) ZSYMM=orderZ();

  output.logTab(1,VERBOSE,"Highest symmetry of multiplicity "+itos(ZSYMM)+" around axis number "+itos(axis));

  p3Dc hkl;
  floatType intensity(0.);
  int countr = 0;
  const floatType ZERO(0.),ONE(1.),TWO(2.);

  floatType HIRES(10000);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    HIRES = std::min(static_cast<floatType>(reso(r)),HIRES);
    countr++;
  }
  output.logTab(1,VERBOSE,"Number of reflections within resolution range is "+itos(countr* SpaceGroup::NSYMP));
  countr = 0;

  const int lmax(LMAX); //change notation
//  besselx.resize(lmax+2);

  output.logTab(1,VERBOSE,"Maximum total angular momentum value is "+itos(lmax));
  output.logTab(1,VERBOSE,"Number of Symmetries is  "+itos(SpaceGroup::NSYMP));
  output.logBlank(VERBOSE);

  elmn.resize(lmax/2); // l=2 => [0],  l=lmax => [lmax/2-1]
  for (int l = 2; l <= lmax; l += 2)
  {
    int l_index = l/2-1;
    int nmax = (lmax-l+2)/2;

    elmn[l_index].resize(2*l+1);
    for (int m = 0; m <= l; m++)
    {
      elmn[l_index][l+m].resize(nmax+1);
      elmn[l_index][l-m].resize(nmax+1);

      for (int n = 1; n <= nmax; n++)
      {
        elmn[l_index][l+m][n] = 0;
        elmn[l_index][l-m][n] = 0;
      }
    }
  }

  output.startClock();
  floatType V(0),cweight(0);

  HKL_clustered HKL_list;

  bool TARGET_CROWTHER(TARGET == "CROWTHER");
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    if (reso(r) > HIRES)
    {
      floatType Esqr = fn::pow2(F[r]/sqrt_epsnSigmaN[r]);
      floatType varE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
      Miller1D  rhkl = rotMiller(r);
      /* EMfixSqr is (SigmaA*Ecalc)^2 for sum of known components,
         i.e. this accounts for fraction of scattering and coordinate error.
         sum_Esqr_known will be zero unless SOLU 3DIM is resurrected .
         variance (V) is SigmaN'/SigmaN from FRF paper:
         SigmaN'/SigmaN = 1 + (sigmaA*Ecalc)^2 - sigmaA^2
         Note that totvar_known is basically sigmaA^2 (possibly corrected by G-function)
         Comparing what is here with equation (18) of the FRF paper, the
         observed coefficients have been multiplied by SigmaN but then the
         calculated coefficients are computed in terms of E-values, which
         is where the other factor of SigmaN goes.
      */
      cweight = cent(r) ? ONE : TWO;
      floatType EMfixSqr(EM_known[r].real()*EM_known[r].real()+EM_known[r].imag()*EM_known[r].imag());
      V = PTNCS.EPSFAC[r] - totvar_known[r] + EMfixSqr + cweight*varE;
      intensity = TARGET_CROWTHER ? TWO*Esqr : TWO*cweight*(Esqr-V)/fn::pow2(V);
      if (PTNCS.use_and_present()) gfun.calcReflTerms(r);

      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          int indexh = rhkl[isym][0];
          int indexk = rhkl[isym][1];
          int indexl = rhkl[isym][2];
          dvect3 hklscaled = UnitCell::HKLtoVecS(indexh,indexk,indexl);
          if (axis==3)
          {
            hkl.x = hklscaled[0];
            hkl.y = hklscaled[1];
            hkl.z = hklscaled[2];
          }
          else if (axis==2)
          {
            hkl.y = hklscaled[0];
            hkl.z = hklscaled[1];
            hkl.x = hklscaled[2];
          }
          else
          {
            hkl.z = hklscaled[0];
            hkl.x = hklscaled[1];
            hkl.y = hklscaled[2];
          }
          HKL HKL_temp;
          HKL_temp.intensity = intensity;
          if (PTNCS.use_and_present()) HKL_temp.intensity *= gfun.refl_Gsqr[isym];
          HKL_temp.htp = toPolar(hkl);
          HKL_list.add(HKL_temp);
        }
      }
    }
  }
  //randomize the list to even up the number on each processor
  HKL_list.shuffle();

  float1D sqrt_table(3*lmax);
  for (int i=0; i<3*lmax; i++) sqrt_table[i]=sqrt(TWO*i+1);

// OpenMP reduction() clause cannot be used for parallelising the summing of
// the elmn variable as it hasn't got operator+ overloaded. So omp_elmn
// serves as an omp extra variable to that end
  cmplx1D omp_elmn_first;
  int nthreads = 1;
  int num_elmn(0);

#pragma omp parallel
  {
    int nthread = 0;
#ifdef _OPENMP
    nthreads = omp_get_num_threads();
    nthread = omp_get_thread_num();
    if (nthread != 0) output.usePhenixCallback(false);
#endif

#pragma omp single
    {
      if (nthreads > 1)
        output.logTab(1,LOGFILE,"Spreading calculation onto "+itos(nthreads)+ " threads.");
      output.logProgressBarStart(LOGFILE,"Elmn Calculation for Data",HKL_list.clustered.size()/nthreads);

      for (int m = 0; m <= lmax; m++)
      {
        if ((m/ZSYMM)*ZSYMM !=m) continue;
        for (int l = m; l <= lmax ; l += 2)
        {
          if (l%2 == 1) l++;
          if (l < 2) continue;
          int nmax = (lmax-l+2)/2;
          int l_index = l/2-1;
          for (int n = 1; n <= nmax; n++)
          {
            num_elmn++;
          } // n loop
        } // l loop
      } // m loop
      omp_elmn_first.resize(nthreads*num_elmn,0);
    }

    p3Dp htp;
    lnfactorial lnfac;

#pragma omp for
    for ( int c=0 ; c < HKL_list.clustered.size(); c++)
    {
      float1D spharmfn;
      PHASER_ASSERT(HKL_list.clustered[c].size());
      floatType x = std::cos(HKL_list.clustered[c][0].htp.theta);
      floatType sqRoot = std::sqrt(std::max(0.0,1.0-x*x));
      floatType Pmj(0),Pmm(0),Pmm1(0);
      for (int m = 0; m <= lmax; m++)
      {
        if ((m/ZSYMM)*ZSYMM !=m) continue;
        floatType Pmm_j(0),Pmm1_j(0);
        floatType Pmm(1);
        if (m>0) {
          int odd_fact=1;
          for (int i=1; i<=m; i++) {
            Pmm*= -sqRoot*odd_fact;
            odd_fact += 2;
          }
        }
        for (int l = m; l <= lmax ; l += 2)
        {
          if (l%2 == 1) l++;
          if (l < 2) continue;
          floatType normalization = std::sqrt((2*l+1)/scitbx::constants::four_pi);
          floatType scale = normalization * exp(0.5*(lnfac.get(l-scitbx::fn::absolute(m))-lnfac.get(l+scitbx::fn::absolute(m))));
          if (Pmm_j)
          {
            Pmj = (x*(2*l-3)*Pmm1_j-(l+m-2)*Pmm_j)/(l-1-m);
            Pmm = Pmm1_j;
            Pmm1 = Pmj;
            Pmj = (x*(2*l-1)*Pmm1-(l+m-1)*Pmm)/(l-m);
            Pmm_j = Pmm1;
            Pmm1_j = Pmj;
          }
          else if (l==m)
          {
            Pmj = Pmm;
          }
          else if (l==m+1)
          {
            Pmm1=x*(2*m+1)*Pmm;
            Pmj = Pmm1;
          }
          else
          {
            Pmm1=x*(2*m+1)*Pmm;
            for (int j=m+2; j<=l; j++) {
              Pmj = (x*(2*j-1)*Pmm1-(j+m-1)*Pmm)/(j-m);
              Pmm = Pmm1;
              Pmm1 = Pmj;
            }
            Pmm_j = Pmm;
            Pmm1_j = Pmm1;
          }
          spharmfn.push_back(scale*Pmj);
        } // l loop
      } // m loop

      for (int r=0; r < HKL_list.clustered[c].size(); r++)
      {
        int i(0),n_elmn(0);
        floatType pintensity = HKL_list.clustered[c][r].intensity;
        float1D besselx(lmax+2);
        floatType h = lmax*HKL_list.clustered[c][r].htp.r*HIRES; //TWOPI*htp.r*b;
        for (int u = 3; u<lmax+2; u+=2)
          besselx[u] = sqrt_table[u]*sphbessel(u,h)/h;

        cmplxType phi_phase_prev = ONE;
        cmplxType phi_phase_1 = std::exp(cmplxType(ZERO,ONE)*HKL_list.clustered[c][r].htp.phi);
        for (int m = 0; m <= lmax; m++)
        {
          cmplxType phi_phase_m = m ? phi_phase_prev * phi_phase_1 : phi_phase_prev;
          phi_phase_prev = phi_phase_m;
          if ((m/ZSYMM)*ZSYMM !=m) continue;
          for (int l = m; l <= lmax ; l += 2)
          {
            if (l%2 == 1) l++;
            if (l < 2) continue;
            cmplxType Ylm = spharmfn[i++];
                      Ylm *= phi_phase_m;
                      Ylm *= pintensity;
            int nmax = (lmax-l+2)/2;
            int l_index = l/2-1;
            for (int n = 1; n <= nmax; n++)
            {
              omp_elmn_first[nthread*num_elmn+n_elmn] += besselx[l+2*n-1]*Ylm; //note way of storing this!
              n_elmn++;
            }
          } // l loop
        } // m loop
      } // r loop for clustered seta
      if (nthread == 0) output.logProgressBarNext(LOGFILE);
    } // c loop and #pragma omp for
  }  // #pragma omp parallel
// now we're unparallelised again
  output.logProgressBarEnd(LOGFILE);
  output.logBlank(VERBOSE);

// sum up values of omp extra variable
  for (int nthread=0; nthread<nthreads; nthread++)
  {
    int n_elmn(0);
    for (int m = 0; m <= lmax; m++)
    {
      if ((m/ZSYMM)*ZSYMM !=m) continue;
      for (int l = m; l <= lmax ; l += 2)
      {
        if (l%2 == 1) l++;
        if (l < 2) continue;
        int nmax = (lmax-l+2)/2;
        int l_index = l/2-1;
        for (int n = 1; n <= nmax; n++)
        {
          elmn[l_index][l+m][n] += omp_elmn_first[nthread*num_elmn+n_elmn];
          n_elmn++;
        } // n loop
      } // l loop
    } // m loop
  }

  output.logElapsedTime(DEBUG);
  output.logBlank(DEBUG);
  for (unsigned l_index = 0; l_index < elmn.size(); l_index++)
  {
    unsigned l=2*l_index+2;
    for (int m_index = 0; m_index <= l; m_index++)
      for (unsigned n_index = 1; n_index < elmn[l_index][m_index].size(); n_index++)
      {
        floatType oddFac = (m_index & 1) ? -ONE : ONE; //bitwise AND operation - looks at last bit, very fast
        elmn[l_index][l-m_index][n_index]= oddFac*std::conj(elmn[l_index][l+m_index][n_index]);
      }
  }
  return elmn;
}

cmplx3D DataMR::getELMNxR2(std::string TARGET,Output& output,int LMAX)
{
  output.logBlank(VERBOSE);
  output.logTab(1,LOGFILE,"Elmn for Data");
  if (TARGET != "CROWTHER" && TARGET != "LERF1")
  throw PhaserError(FATAL,"Unknown target " + TARGET);

  //function return
  cmplx3D elmn;

  int ZSYMM(0);

  int axis(highOrderAxis());
  if (axis==1) ZSYMM=orderX();
  if (axis==2) ZSYMM=orderY();
  if (axis==3) ZSYMM=orderZ();

  output.logTab(1,VERBOSE,"Highest symmetry of multiplicity "+itos(ZSYMM)+" around axis number "+itos(axis));

  p3Dc hkl;
  floatType intensity(0.);
  int countr = 0;
  const floatType ZERO(0.),ONE(1.),TWO(2.);

  floatType HIRES(10000);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    HIRES = std::min(static_cast<floatType>(reso(r)),HIRES);
    countr++;
  }
  output.logTab(1,VERBOSE,"Number of reflections within resolution range is "+itos(countr* SpaceGroup::NSYMP));
  countr = 0;

  const int lmax(LMAX); //change notation
//  besselx.resize(lmax+2);

  output.logTab(1,VERBOSE,"Maximum total angular momentum value is "+itos(lmax));
  output.logTab(1,VERBOSE,"Number of Symmetries is  "+itos(SpaceGroup::NSYMP));
  output.logBlank(VERBOSE);

  elmn.resize(lmax/2); // l=2 => [0],  l=lmax => [lmax/2-1]
  for (int l = 2; l <= lmax; l += 2)
  {
    int l_index = l/2-1;
    int nmax = (lmax-l+2)/2;

    elmn[l_index].resize(2*l+1);
    for (int m = 0; m <= l; m++)
    {
      elmn[l_index][l+m].resize(nmax+1);
      elmn[l_index][l-m].resize(nmax+1);

      for (int n = 1; n <= nmax; n++)
      {
        elmn[l_index][l+m][n] = 0;
        elmn[l_index][l-m][n] = 0;
      }
    }
  }

  output.startClock();
  floatType V(0),cweight(0);

  HKL_clustered HKL_list;

  bool TARGET_CROWTHER(TARGET == "CROWTHER");
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    if (reso(r) > HIRES)
    {
      floatType Esqr = fn::pow2(F[r]/sqrt_epsnSigmaN[r]);
      floatType varE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
      Miller1D  rhkl = rotMiller(r);
      /* EMfixSqr is (SigmaA*Ecalc)^2 for sum of known components,
         i.e. this accounts for fraction of scattering and coordinate error.
         sum_Esqr_known will be zero unless SOLU 3DIM is resurrected .
         variance (V) is SigmaN'/SigmaN from FRF paper:
         SigmaN'/SigmaN = 1 + (sigmaA*Ecalc)^2 - sigmaA^2
         Note that totvar_known is basically sigmaA^2 (possibly corrected by G-function)
         Comparing what is here with equation (18) of the FRF paper, the
         observed coefficients have been multiplied by SigmaN but then the
         calculated coefficients are computed in terms of E-values, which
         is where the other factor of SigmaN goes.
      */
      cweight = cent(r) ? ONE : TWO;
      floatType EMfixSqr(EM_known[r].real()*EM_known[r].real()+EM_known[r].imag()*EM_known[r].imag());
      V = PTNCS.EPSFAC[r] - totvar_known[r] + EMfixSqr + cweight*varE;
      intensity = TARGET_CROWTHER ? TWO*Esqr : TWO*cweight*(Esqr-V)/fn::pow2(V);
      if (PTNCS.use_and_present()) gfun.calcReflTerms(r);

      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          int indexh = rhkl[isym][0];
          int indexk = rhkl[isym][1];
          int indexl = rhkl[isym][2];
          dvect3 hklscaled = UnitCell::HKLtoVecS(indexh,indexk,indexl);
          if (axis==3)
          {
            hkl.x = hklscaled[0];
            hkl.y = hklscaled[1];
            hkl.z = hklscaled[2];
          }
          else if (axis==2)
          {
            hkl.y = hklscaled[0];
            hkl.z = hklscaled[1];
            hkl.x = hklscaled[2];
          }
          else
          {
            hkl.z = hklscaled[0];
            hkl.x = hklscaled[1];
            hkl.y = hklscaled[2];
          }
          HKL HKL_temp;
          HKL_temp.intensity = intensity;
          if (PTNCS.use_and_present()) HKL_temp.intensity *= gfun.refl_Gsqr[isym];
          HKL_temp.htp = toPolar(hkl);
          HKL_list.add(HKL_temp);
        }
      }
    }
  }
  //randomize the list to even up the number on each processor
  HKL_list.shuffle();

  float1D sqrt_table(3*lmax);
  for (int i=0; i<3*lmax; i++) sqrt_table[i]=sqrt(TWO*i+1);

// OpenMP reduction() clause cannot be used for parallelising the summing of
// the elmn variable as it hasn't got operator+ overloaded. So omp_elmn
// serves as an omp extra variable to that end
  cmplx1D omp_elmn_first;
  int nthreads = 1;
  int num_elmn(0);

#pragma omp parallel
  {
    int nthread = 0;
#ifdef _OPENMP
    nthreads = omp_get_num_threads();
    nthread = omp_get_thread_num();
    if (nthread != 0) output.usePhenixCallback(false);
#endif

#pragma omp single
    {
      if (nthreads > 1)
        output.logTab(1,LOGFILE,"Spreading calculation onto "+itos(nthreads)+ " threads.");
      output.logProgressBarStart(LOGFILE,"Elmn Calculation for Data",HKL_list.clustered.size()/nthreads);

      for (int m = 0; m <= lmax; m++)
      {
        if ((m/ZSYMM)*ZSYMM !=m) continue;
        for (int l = m; l <= lmax ; l += 2)
        {
          if (l%2 == 1) l++;
          if (l < 2) continue;
          int nmax = (lmax-l+2)/2;
          int l_index = l/2-1;
          for (int n = 1; n <= nmax; n++)
          {
            num_elmn++;
          } // n loop
        } // l loop
      } // m loop
      omp_elmn_first.resize(nthreads*num_elmn,0);
    }

    p3Dp htp;
    lnfactorial lnfac;

#pragma omp for
    for ( int c=0 ; c < HKL_list.clustered.size(); c++)
    {
      float1D spharmfn;
      PHASER_ASSERT(HKL_list.clustered[c].size());
      floatType x = std::cos(HKL_list.clustered[c][0].htp.theta);
      floatType sqRoot = std::sqrt(std::max(0.0,1.0-x*x));
      floatType Pmj(0),Pmm(0),Pmm1(0);
      for (int m = 0; m <= lmax; m++)
      {
        if ((m/ZSYMM)*ZSYMM !=m) continue;
        floatType Pmm_j(0),Pmm1_j(0);
        floatType Pmm(1);
        if (m>0) {
          int odd_fact=1;
          for (int i=1; i<=m; i++) {
            Pmm*= -sqRoot*odd_fact;
            odd_fact += 2;
          }
        }
        for (int l = m; l <= lmax ; l += 2)
        {
          if (l%2 == 1) l++;
          if (l < 2) continue;
          floatType normalization = std::sqrt((2*l+1)/scitbx::constants::four_pi);
          floatType scale = normalization * exp(0.5*(lnfac.get(l-scitbx::fn::absolute(m))-lnfac.get(l+scitbx::fn::absolute(m))));
          if (Pmm_j)
          {
            Pmj = (x*(2*l-3)*Pmm1_j-(l+m-2)*Pmm_j)/(l-1-m);
            Pmm = Pmm1_j;
            Pmm1 = Pmj;
            Pmj = (x*(2*l-1)*Pmm1-(l+m-1)*Pmm)/(l-m);
            Pmm_j = Pmm1;
            Pmm1_j = Pmj;
          }
          else if (l==m)
          {
            Pmj = Pmm;
          }
          else if (l==m+1)
          {
            Pmm1=x*(2*m+1)*Pmm;
            Pmj = Pmm1;
          }
          else
          {
            Pmm1=x*(2*m+1)*Pmm;
            for (int j=m+2; j<=l; j++) {
              Pmj = (x*(2*j-1)*Pmm1-(j+m-1)*Pmm)/(j-m);
              Pmm = Pmm1;
              Pmm1 = Pmj;
            }
            Pmm_j = Pmm;
            Pmm1_j = Pmm1;
          }
          spharmfn.push_back(scale*Pmj);
        } // l loop
      } // m loop

      for (int r=0; r < HKL_list.clustered[c].size(); r++)
      {
        int i(0),n_elmn(0);
        floatType pintensity = HKL_list.clustered[c][r].intensity;
        float1D besselx(lmax+2);
        floatType h = lmax*HKL_list.clustered[c][r].htp.r*HIRES; //TWOPI*htp.r*b;
        for (int u = 3; u<lmax+2; u+=2)
          besselx[u] = sqrt_table[u]*sphbessel(u,h)/h;

        cmplxType phi_phase_prev = ONE;
        cmplxType phi_phase_1 = std::exp(cmplxType(ZERO,ONE)*HKL_list.clustered[c][r].htp.phi);
        for (int m = 0; m <= lmax; m++)
        {
          cmplxType phi_phase_m = m ? phi_phase_prev * phi_phase_1 : phi_phase_prev;
          phi_phase_prev = phi_phase_m;
          if ((m/ZSYMM)*ZSYMM !=m) continue;
          for (int l = m; l <= lmax ; l += 2)
          {
            if (l%2 == 1) l++;
            if (l < 2) continue;
            cmplxType Ylm = spharmfn[i++];
                      Ylm *= phi_phase_m;
                      Ylm *= pintensity;
            int nmax = (lmax-l+2)/2;
            int l_index = l/2-1;
            for (int n = 1; n <= nmax; n++)
            {
              omp_elmn_first[nthread*num_elmn+n_elmn] += besselx[l+2*n-1]*Ylm; //note way of storing this!
              n_elmn++;
            }
          } // l loop
        } // m loop
      } // r loop for clustered seta
      if (nthread == 0) output.logProgressBarNext(LOGFILE);
    } // c loop and #pragma omp for
  }  // #pragma omp parallel
// now we're unparallelised again
  output.logProgressBarEnd(LOGFILE);
  output.logBlank(VERBOSE);

// sum up values of omp extra variable
  for (int nthread=0; nthread<nthreads; nthread++)
  {
    int n_elmn(0);
    for (int m = 0; m <= lmax; m++)
    {
      if ((m/ZSYMM)*ZSYMM !=m) continue;
      for (int l = m; l <= lmax ; l += 2)
      {
        if (l%2 == 1) l++;
        if (l < 2) continue;
        int nmax = (lmax-l+2)/2;
        int l_index = l/2-1;
        for (int n = 1; n <= nmax; n++)
        {
          elmn[l_index][l+m][n] += omp_elmn_first[nthread*num_elmn+n_elmn];
          n_elmn++;
        } // n loop
      } // l loop
    } // m loop
  }

  output.logElapsedTime(DEBUG);
  output.logBlank(DEBUG);
  for (unsigned l_index = 0; l_index < elmn.size(); l_index++)
  {
    unsigned l=2*l_index+2;
    for (int m_index = 0; m_index <= l; m_index++)
      for (unsigned n_index = 1; n_index < elmn[l_index][m_index].size(); n_index++)
      {
        floatType oddFac = (m_index & 1) ? -ONE : ONE; //bitwise AND operation - looks at last bit, very fast
        elmn[l_index][l-m_index][n_index]= oddFac*std::conj(elmn[l_index][l+m_index][n_index]);
      }
  }
  return elmn;
}

af_cmplx DataMR::getFcalc(dmat33 ROT,dmat33 ROT2,Ensemble& MAP,af_miller P1,bool halfR,floatType searchBfac)
{
  int P1_size = P1.size();
  af_cmplx Fcalc(P1_size);

  dmat33 R = ROT*MAP.Frac2Orth();
  dmat33 Q1 = doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();

  dmat33 R2 = ROT2*MAP.Frac2Orth();
  dmat33 Q12 = doOrth2Frac(R2);
  dmat33 Q1tr2 = Q12.transpose();

  int rP1(0);
  floatType sqrt_scatFactor = std::sqrt(MAP.ensScat(TOTAL_SCAT)/NSYMP);
  const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    Miller1D rhkl = rotMiller(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    MAP.set_sqrt_ssqr(sqrtrssqr,rssqr);
    if (PTNCS.use_and_present())
    {
      gfun.calcReflTerms(r);
    }
    for (unsigned isym = 0; isym < NSYMP; isym++)
    {
      if (!duplicate(isym,rhkl))
      {
        PHASER_ASSERT(rP1 < P1_size); //paranoia
        dvect3 RotSymHKL = Q1tr*rhkl[isym];
        cmplxType thisE = MAP.InterpE(RotSymHKL);
                  thisE *= sqrt_scatFactor;
        if (searchBfac)
                  thisE *= std::exp(-searchBfac*rssqr/4.0);
        if (PTNCS.use_and_present())
        {
          floatType theta  = rhkl[isym][0]*PTNCS.TRA.VECTOR[0] +
                             rhkl[isym][1]*PTNCS.TRA.VECTOR[1] +
                             rhkl[isym][2]*PTNCS.TRA.VECTOR[2];
          if (PTNCS.NMOL > 2)
          {
            cmplxType ptncs_scat = cmplxONE; //iMOL=0
            for (floatType iMOL = 1; iMOL < PTNCS.NMOL; iMOL++)
              ptncs_scat += std::exp(TWOPII*iMOL*theta);
            thisE *= ptncs_scat;
          }
          else if (halfR) //second molecule comes from doubling at the one rotation
          {
            thisE *= cmplxONE + std::exp(TWOPII*theta);
            thisE *= gfun.refl_G[isym];
          }
          else //second molecule uses new orientation
          {
            dvect3 RotSymHKL2 = Q1tr2*rhkl[isym];
            cmplxType thisE2 = MAP.InterpE(RotSymHKL2);
                      thisE2 *= sqrt_scatFactor;
                      thisE2 *= std::exp(TWOPII*theta);
            thisE += thisE2;
          }
        }
        P1[rP1] = rhkl[isym]; //in place change
        Fcalc[rP1] = thisE;
        rP1++;
      }
    }
  }
  PHASER_ASSERT(P1_size == rP1); //paranoia
  return Fcalc;
}

af_cmplx DataMR::getFcalc(std::string MODLID,mrgyre1D& GYRE,MapEnsemble& ensemble,af_miller P1,floatType searchBfac)
{
  int P1_size = P1.size();
  af_cmplx Fcalc(P1_size);
  for (int g = 0; g < GYRE.size(); g++)
  {
    std::string modlid = GYRE[g].modlid(MODLID);
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    dmat33 ROT = euler2matrixDEG(GYRE[g].EULR);
           ROT = ROT*ens_modlid->principalRtr();
    dvect3 TRA = GYRE[g].FRAC;
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr= Q1.transpose();
    int rP1(0);
    floatType sqrt_scatFactor = std::sqrt(ens_modlid->ensScat(TOTAL_SCAT)/NSYMP);
    const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      Miller1D rhkl = rotMiller(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      for (unsigned isym = 0; isym < NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          PHASER_ASSERT(rP1 < P1_size); //paranoia
          dvect3 RotSymHKL = Q1tr*rhkl[isym];
          cmplxType thisE = ens_modlid->InterpE(RotSymHKL);
                    thisE *= sqrt_scatFactor;
          if (searchBfac)
                    thisE *= std::exp(-searchBfac*rssqr/4.0);
          floatType theta  = rhkl[isym][0]*TRA[0] +
                             rhkl[isym][1]*TRA[1] +
                             rhkl[isym][2]*TRA[2];
                    thisE *= std::exp(TWOPII*theta);
          //          thisE *= dphi(r,isym,rhkl[isym],TRA);
          P1[rP1] = rhkl[isym]; //in place change
          Fcalc[rP1] += thisE;
          rP1++;
        }
      }
    }
  }
  return Fcalc;
}

af_cmplx DataMR::getFcalcPhsTF(dmat33 ROT,Ensemble& MAP,floatType searchBfac)
{
  af_cmplx Fcalc;
  dmat33 R = ROT*MAP.Frac2Orth();
  dmat33 Q1 = doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();
  floatType sqrt_scatFactor = std::sqrt(MAP.ensScat(TOTAL_SCAT)/NSYMP);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType rssqr = ssqr(r);
    floatType sqrt_repsn = sqrt_epsn[r];
    floatType sqrtrssqr = std::sqrt(rssqr);
    MAP.set_sqrt_ssqr(sqrtrssqr,rssqr);
    Miller1D  rhkl = rotMiller(r);
    for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
    {
      if (!duplicate(isym,rhkl))
      {
        cmplxType thisE = MAP.InterpE(Q1tr*rhkl[isym]);
                  thisE *= sqrt_scatFactor;
                  thisE *= sqrt_repsn;
        if (searchBfac)
                  thisE *= std::exp(-searchBfac*rssqr/4.0);
        Fcalc.push_back(thisE);
      }
    }
  }
  return Fcalc;
}

af_cmplx DataMR::getFobsPhsTF()
{
  af_cmplx Fobs;
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    Fobs.push_back(FOM[r]*std::polar(F[r]/sqrt_epsnSigmaN[r],scitbx::deg_as_rad(PH[r])));
  return Fobs;
}

std::vector<mr_ndim> DataMR::getModls()
{
  return models_known;
}

af_cmplx DataMR::getEpart()
{
  af_cmplx Epart;
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    Epart.push_back(EM_known[r]);
  return Epart;
}

bool DataMR::haveFpart()
{
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
    if (EM_known[r] != cmplxType(0,0))
      return true;
  return false;
}

af_float DataMR::getEnat()
{
  af_float Enat;
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    Enat.push_back(E);
  }
  return Enat;
}

floatType DataMR::m_LETF1(dmat33 ROT,Ensemble& ens_modlid,bool use_rot,af_float& mI,bool halfR,floatType searchBfac)
{
  mI.clear();
  dmat33 R = ROT*ens_modlid.Frac2Orth();
  dmat33 Q1 = UnitCell::doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();
  floatType C1(0.);
  int countr(0);

  floatType scatFactor = ens_modlid.ensScat(TOTAL_SCAT)/NSYMP;
  floatType NonSearchScat = TOTAL_SCAT - ens_modlid.getEnsembleScat();
  //std::string fname =  "move_var.txt";
  //std::ofstream myFTFfile(fname.c_str(), std::ios::out );

  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    bool rcent = cent(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ens_modlid.set_sqrt_ssqr(sqrtrssqr,rssqr);
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
    if (!rcent) wtVarE *= 2.;
    floatType thisV = ens_modlid.InterpV();
    /*floatType rMoveScatFactor = totvar_move[r]/(NonSearchScat+totvar_move[r]);
              rMoveScatFactor /= NSYMP;
    if(totvar_move[r] > 0)
    {
              thisV *= rMoveScatFactor*NSYMP;
    }
    else thisV *= scatFactor*NSYMP;*/
              //thisV *= scatFactor*NSYMP;
    if (searchBfac)
              thisV *= exp(-2.0*searchBfac*ssqr(r)/4.0);
    if (PTNCS.use_and_present())
    {
      if (PTNCS.NMOL > 2 || !halfR)
        thisV *= PTNCS.NMOL*G_Vterm[r];
      else
      {
        gfun.calcReflTerms(r);
        thisV *= PTNCS.NMOL*gfun.refl_Gsqr_Vterm;
      }
    }
    //totvar_search[r] = thisV*scatFactor*NSYMP;
    totvar_search[r] = thisV*totvar_move[r]/TOTAL_SCAT;
    floatType phasedEsqr(EM_known[r].real()*EM_known[r].real()+EM_known[r].imag()*EM_known[r].imag());
    floatType eImove(0);
    const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
    floatType sum_Esqr_search(0);
    if (use_rot)
    {
      //floatType sum_Esqr_search(0);
      floatType repsn = epsn(r);
      Miller1D rhkl = rotMiller(r);
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3 RotSymHKL = Q1tr*rhkl[isym];
          cmplxType thisE = ens_modlid.InterpE(RotSymHKL);
          if (searchBfac)
                    thisE *= std::exp(-searchBfac*rssqr/4.0);
          if (PTNCS.use_and_present())
          {
            floatType theta  = rhkl[isym][0]*PTNCS.TRA.VECTOR[0] +
                               rhkl[isym][1]*PTNCS.TRA.VECTOR[1] +
                               rhkl[isym][2]*PTNCS.TRA.VECTOR[2];
            cmplxType ptncs_scat = cmplxONE;
            for (floatType iMOL = 1; iMOL < PTNCS.NMOL; iMOL++)
            {
              ptncs_scat += std::exp(TWOPII*iMOL*theta);
            }
            if (halfR && PTNCS.NMOL == 2)
              ptncs_scat *= gfun.refl_G[isym];
            thisE *= ptncs_scat;
          }
          floatType thisEsqr = thisE.real()*thisE.real() + thisE.imag()*thisE.imag();
                    //thisEsqr *= repsn;
                    thisEsqr /= repsn;
          /*if(totvar_move[r] > 0)
          {
            thisEsqr *= rMoveScatFactor;
          }
          else*/
                    thisEsqr *= scatFactor;
          sum_Esqr_search += thisEsqr;
        }
      }
      /*totvar_search[r] = sum_Esqr_search;
      totvar_known[r] /= (NonSearchScat/TOTAL_SCAT + totvar_search[r]/thisV);
      totvar_search[r] /= (NonSearchScat/TOTAL_SCAT + totvar_search[r]/thisV);

      while(PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE < 0) {
        totvar_search[r] *= 0.99;
        totvar_known[r] *= 0.99;
      }*/

      //myFTFfile << totvar_search[r] <<std::endl;

      eImove = phasedEsqr + sum_Esqr_search; // <Imove> knowing molecular transforms
    }
    else
    {
      eImove = phasedEsqr + totvar_search[r]; // <Imove> from scattering content
    }
    PHASER_ASSERT(eImove >= 0);
    floatType C(PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE);
    //floatType C(PTNCS.EPSFAC[r]*NonSearchScat/TOTAL_SCAT - totvar_known[r] + sum_Esqr_search*(1/thisV - 1));
    //floatType C(PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE + totvar_move[r]/TOTAL_SCAT*(PTNCS.EPSFAC[r] - thisV));
    floatType fracB = 1 + frac_known_B[r] + frac_search_B[r];
    C /= fracB;
    floatType sqrt_eImove = std::sqrt(eImove/fracB);
    PHASER_ASSERT(C > 0);
    floatType X = 2.*E*sqrt_eImove/C;
    floatType weight = rcent ? 2.0 : 1.0;
    X = X/weight;
    floatType rh = rcent ? std::tanh(X) : sim(X);
    floatType B1;
    if (sqrt_eImove > 0)
      B1 = (rh*E/sqrt_eImove - 1)/(weight*C);
    else
      B1 = (E*E/C - 1)/(weight*C); // Limit as sqrt_eImove -> 0
    mI.push_back(B1);
    //Compute LL for eImove, leaving out constants that are left out of Rice()
    floatType eLL = rcent ? logRelWoolfson(E,sqrt_eImove,C) : logRelRice(E,sqrt_eImove,C);
    C1 += eLL - B1*eImove;
    countr++;
  }
  //myFTFfile.close();

  return C1;
}

floatType DataMR::m_LETF1(std::string MODLID,mrgyre1D& GYRE,MapEnsemble& ensemble,bool use_rot,af_float& mI,floatType searchBfac)
{
  mI.clear();
  floatType C1(0.);
  int countr(0);
  floatType scatFactor = ensemble[MODLID].ensScat(TOTAL_SCAT)/NSYMP;
  dvect31D TRA(GYRE.size());
  dmat331D Q1tr(GYRE.size());
  float1D  scatFactor_g(GYRE.size());
  for (int g = 0; g < GYRE.size(); g++)
  {
    std::string modlid = GYRE[g].modlid(MODLID);
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    dmat33 ROT = euler2matrixDEG(GYRE[g].EULR);
           ROT = ROT*ens_modlid->principalRtr();
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
           Q1tr[g] = Q1.transpose();
           TRA[g] = GYRE[g].FRAC;
    scatFactor_g[g] = ens_modlid->ensScat(TOTAL_SCAT)/NSYMP;
  }
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    bool rcent = cent(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ensemble[MODLID].set_sqrt_ssqr(sqrtrssqr,rssqr); //same for all
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
    if (!rcent) wtVarE *= 2.;
    floatType thisV = ensemble[MODLID].InterpV();
              thisV *= scatFactor*NSYMP;
    if (searchBfac)
              thisV *= exp(-2.0*searchBfac*ssqr(r)/4.0);
    if (PTNCS.use_and_present())
    {
      thisV *= PTNCS.NMOL*G_Vterm[r];
    }
    totvar_search[r] = thisV;
    floatType phasedEsqr(EM_known[r].real()*EM_known[r].real()+EM_known[r].imag()*EM_known[r].imag());
    floatType eImove(0);
    const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
    if (use_rot)
    {
      floatType sum_Esqr_search(0);
      floatType repsn = epsn(r);
      Miller1D rhkl = rotMiller(r);
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          cmplxType thisE(0);
          for (int g = 0; g < GYRE.size(); g++)
          {
            dvect3 RotSymHKL = Q1tr[g]*rhkl[isym];
            std::string modlid = GYRE[g].modlid(MODLID);
            Ensemble* ens_modlid = &ensemble.find(modlid)->second;
            cmplxType thisEg = ens_modlid->InterpE(RotSymHKL);
            if (searchBfac)
                      thisEg *= std::exp(-searchBfac*rssqr/4.0);
            floatType theta  = rhkl[isym][0]*TRA[g][0] +
                               rhkl[isym][1]*TRA[g][1] +
                               rhkl[isym][2]*TRA[g][2];
                      thisEg *= std::exp(TWOPII*theta);
                   //   thisEg *= dphi(r,isym,rhkl[isym],TRA[g]);
                      thisEg *= std::sqrt(scatFactor_g[g]);
            thisE += thisEg;
          }
          floatType thisEsqr = thisE.real()*thisE.real() + thisE.imag()*thisE.imag();
                    thisEsqr *= repsn;
                    //thisEsqr *= scatFactor;
          sum_Esqr_search += thisEsqr;
        }
      }
      eImove = phasedEsqr + sum_Esqr_search; // <Imove> knowing molecular transforms
    }
    else
    {
      eImove = phasedEsqr + totvar_search[r]; // <Imove> from scattering content
    }
    PHASER_ASSERT(eImove >= 0);
    floatType C(PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE);
    floatType fracB = 1 + frac_known_B[r] + frac_search_B[r];
    C /= fracB;
    floatType sqrt_eImove = std::sqrt(eImove/fracB);
    PHASER_ASSERT(C > 0);
    floatType X = 2.*E*sqrt_eImove/C;
    floatType weight = rcent ? 2.0 : 1.0;
    X = X/weight;
    floatType rh = rcent ? std::tanh(X) : sim(X);
    floatType B1;
    if (sqrt_eImove > 0)
      B1 = (rh*E/sqrt_eImove - 1)/(weight*C);
    else
      B1 = (E*E/C - 1)/(weight*C); // Limit as sqrt_eImove -> 0
    mI.push_back(B1);
    //Compute LL for eImove, leaving out constants that are left out of Rice()
    floatType eLL = rcent ? logRelWoolfson(E,sqrt_eImove,C) : logRelRice(E,sqrt_eImove,C);
    C1 += eLL - B1*eImove;
    countr++;
  }
  return C1;
}

floatType DataMR::m_LETF2(dmat33 ROT,Ensemble& ens_modlid,bool use_rot,af_float& mI,af_float& mIsqr,bool halfR,floatType searchBfac)
{
  mI.clear();
  mIsqr.clear();
  dmat33 R = ROT*ens_modlid.Frac2Orth();
  dmat33 Q1 = UnitCell::doOrth2Frac(R);
  dmat33 Q1tr = Q1.transpose();
  floatType C2(0.);
  floatType scatFactor = ens_modlid.ensScat(TOTAL_SCAT)/NSYMP;
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    bool rcent = cent(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    ens_modlid.set_sqrt_ssqr(sqrtrssqr,rssqr);
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
    if (!rcent) wtVarE *= 2.;
    floatType thisV = ens_modlid.InterpV();
              thisV *= scatFactor*NSYMP;
    if (searchBfac)
              thisV *= exp(-2.0*searchBfac*ssqr(r)/4.0);
    if (PTNCS.use_and_present())
    {
      if (PTNCS.NMOL > 2 || !halfR)
        thisV *= PTNCS.NMOL*G_Vterm[r];
      else
      {
        gfun.calcReflTerms(r);
        thisV *= PTNCS.NMOL*gfun.refl_Gsqr_Vterm;
      }
    }
    totvar_search[r] = thisV;
    floatType phasedEsqr(EM_known[r].real()*EM_known[r].real()+EM_known[r].imag()*EM_known[r].imag());
    floatType I(E*E);
    floatType eImove(0);
    const cmplxType cmplxONE(1.,0.),TWOPII(0.,scitbx::constants::two_pi);
    if (use_rot)
    {
      floatType sum_Esqr_search(0);
      bool repsn = epsn(r);
      bool rssqr = ssqr(r);
      Miller1D rhkl =   rotMiller(r);
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3    RotSymHKL = Q1tr*rhkl[isym];
          cmplxType thisE = ens_modlid.InterpE(RotSymHKL);
          if (searchBfac)
                    thisE *= std::exp(-searchBfac*rssqr/4.0);
          if (PTNCS.use_and_present())
          {
            floatType theta  = rhkl[isym][0]*PTNCS.TRA.VECTOR[0] +
                               rhkl[isym][1]*PTNCS.TRA.VECTOR[1] +
                               rhkl[isym][2]*PTNCS.TRA.VECTOR[2];
            cmplxType ptncs_scat = cmplxONE;
            for (floatType iMOL = 1; iMOL < PTNCS.NMOL; iMOL++)
              ptncs_scat += std::exp(TWOPII*iMOL*theta);
            if (halfR && PTNCS.NMOL == 2)
              ptncs_scat *= gfun.refl_G[isym];
            thisE *= ptncs_scat;
          }
          floatType thisEsqr = thisE.real()*thisE.real() + thisE.imag()*thisE.imag();
                    thisEsqr *= repsn;
                    thisEsqr *= scatFactor;
          sum_Esqr_search += thisEsqr;
        }
      }
      eImove = phasedEsqr + sum_Esqr_search;
    }
    else
    {
      eImove = phasedEsqr + totvar_search[r];
    }
    floatType C(PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE);
    floatType fracB = 1 + frac_known_B[r] + frac_search_B[r];
    C /= fracB;
    floatType sqrt_eImove = std::sqrt(eImove/fracB);
    floatType X = 2.*E*sqrt_eImove/C;
    floatType weight = rcent ? 2.0 : 1.0;
    X = X/weight;
    floatType rh = rcent ? std::tanh(X) : sim(X);
    floatType LLp((rh*E/sqrt_eImove - 1)/(weight*C));
    floatType weightSqr = rcent ? 4.0 : 1.0;
    floatType LLdp((1/C/eImove/weightSqr)*((1-rh*rh)*I/(C) - rh*E/sqrt_eImove));
    mI.push_back(LLp-LLdp*eImove), mIsqr.push_back(LLdp/2.0);
    //Compute LL for eImove, leaving out constants that are left out of Rice()
    floatType eLL = rcent ? logRelWoolfson(E,sqrt_eImove,C) : logRelRice(E,sqrt_eImove,C);
    C2 += eLL - LLp*eImove + LLdp*eImove*eImove/2.;
  }
  return C2;
}

MapCoefs DataMR::getMapCoefs()
{
  MapCoefs MAPCOEFS;
  MAPCOEFS.resize(NREFL); //first call only
  std::pair<floatType,floatType> scale = ScaleFC();
  floatType scaleK(scale.first);
  floatType scaleB(scale.second);

  // Get parameters to put everything but anisotropy back into SigmaN
  floatType isotropicB = AnisoBeta2IsoB(sigmaN.ANISO,A(),B(),C(),cosAlpha(),cosBeta(),cosGamma());
  dmat6 isoBasAnisoBeta = IsoB2AnisoBeta(isotropicB,aStar(),bStar(),cStar(),cosAlphaStar(),cosBetaStar(),cosGammaStar());

  //sharpening factor is obtained from most negative eigenvalue of
  //anisotropic part of sigmaN.ANISO.
  dvect3 eigenBvals = getEigenBs(); //These apply to amplitudes, so multiply by 2
  floatType minB(2.*std::min(eigenBvals[0],std::min(eigenBvals[1],eigenBvals[2])));
  floatType resharpB = minB*Resharp.FRAC;
  dmat6 deltaAniso = IsoB2AnisoBeta(resharpB,
                                    aStar(),bStar(),cStar(),
                                    cosAlphaStar(),cosBetaStar(),cosGammaStar());
  for (unsigned n = 0; n < 6; n++)
    isoBasAnisoBeta[n] += deltaAniso[n]; //add in sharpening factor

  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r]) //case for not selected below
  {
    bool rcent = cent(r);
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]);
    if (!rcent) wtVarE *= 2.;
    // Map coeffs only calculated when relative phase known for all components.
    // However, can be classified as fixed or moving, depending on context.
    floatType sigmaA_Ecalc(std::abs(EM_known[r]+EM_search[r]));
    floatType phase = scitbx::rad_as_deg(std::arg(EM_known[r]+EM_search[r]));
    floatType V = PTNCS.EPSFAC[r] - totvar_known[r] - totvar_search[r] + wtVarE;
    PHASER_ASSERT(V > 0);
    floatType X(2*E*sigmaA_Ecalc/V);
    floatType fom = rcent ? std::tanh(X/2.) : sim(X);
    floatType fwt = rcent ? fom*E : 2*fom*E-sigmaA_Ecalc;
    fwt *= sqrt_epsnSigmaN[r];
    if (fwt >= 0)
    {
      MAPCOEFS.FWT[r] = fwt;
      MAPCOEFS.PHWT[r] = phase;
    }
    else
    {
      MAPCOEFS.FWT[r] = -fwt;
      MAPCOEFS.PHWT[r] = phase+180;
    }
    floatType delfwt(sqrt_epsnSigmaN[r]*(fom*E-sigmaA_Ecalc));
    if (delfwt >= 0)
    {
      MAPCOEFS.DELFWT[r] = delfwt;
      MAPCOEFS.PHDELWT[r] = phase;
    }
    else
    {
      MAPCOEFS.DELFWT[r] = -delfwt;
      MAPCOEFS.PHDELWT[r] = phase+180;
    }
    MAPCOEFS.FOM[r] = fom;
    floatType FC(sigmaA_Ecalc*scaleK*exp(-scaleB*ssqr(r)/4)*sqrt_epsnSigmaN[r]);
    MAPCOEFS.FC[r] = FC;
    MAPCOEFS.PHIC[r] = phase;
    floatType radphase(scitbx::deg_as_rad(phase));
    MAPCOEFS.HLA[r] = X*cos(radphase);
    MAPCOEFS.HLB[r] = X*sin(radphase);
    if (rcent)
    {
      MAPCOEFS.HLA[r] /= 2.;
      MAPCOEFS.HLB[r] /= 2.;
    }
    MAPCOEFS.HLC[r] = 0;
    MAPCOEFS.HLD[r] = 0;
  }
  else //not just selected, as list must be in sync with the miller list for output
  {
    MAPCOEFS.FWT[r] = MAPCOEFS.PHWT[r] = CCP4::ccp4_nan().f;
    MAPCOEFS.DELFWT[r] = MAPCOEFS.PHDELWT[r] = CCP4::ccp4_nan().f;
    MAPCOEFS.FC[r] = MAPCOEFS.PHIC[r] = MAPCOEFS.FOM[r] = CCP4::ccp4_nan().f;
    MAPCOEFS.HLA[r] = MAPCOEFS.HLB[r] = MAPCOEFS.HLC[r] = MAPCOEFS.HLD[r] = CCP4::ccp4_nan().f;
  }
  return MAPCOEFS;
}

floatType DataMR::logRelRice(floatType F1, floatType DF2, floatType V)
{
/* Calculate log(Rice), omitting the constant terms omitted from WilsonRice.
   Suitable for LLG calculations, where only relative values needed.
   Assuming F1 is constant observed F, this form works for both
   amplitude and intensity-based likelihoods.
*/
  PHASER_ASSERT(F1 >= 0. && DF2 >= 0. && V > 0.);
  floatType X(2*F1*DF2/V);
  return m_alogchI0.getalogI0(X)-(std::log(V)+(fn::pow2(F1)+fn::pow2(DF2))/V);
}

floatType DataMR::logRelWoolfson(floatType F1, floatType DF2, floatType V)
{
// Centric equivalent to logRelRice
  PHASER_ASSERT(F1 >= 0. && DF2 >= 0. && V > 0.);
  floatType X(F1*DF2/V),TWO(2);
  return m_alogch.getalogch(X)-(std::log(V)+(fn::pow2(F1)+fn::pow2(DF2))/V)/TWO;
}

void DataMR::initGfunction(bool& halfR)
{
  if (PTNCS.use_and_present())
  {
    gfun.init(*this,*this,PTNCS.NMOL);
    dmat33 ncsRmat(1,0,0,0,1,0,0,0,1);
    for (int dir = 0; dir < 3; dir++)
      ncsRmat = xyzRotMatDeg(dir,PTNCS.ROT.ANGLE[dir])*ncsRmat;
    if (halfR)
    {
      //convert to an angle about an axis
      dvect3     AXIS = scitbx::math::r3_rotation::axis_and_angle_from_matrix<floatType>(ncsRmat).axis;
      floatType  ANGL = scitbx::math::r3_rotation::axis_and_angle_from_matrix<floatType>(ncsRmat).angle_rad;
      //convert back to matrix, with half the rotation
      ncsRmat = scitbx::math::r3_rotation::axis_and_angle_as_matrix(AXIS,ANGL/2);
    }
    gfun.initArrays(PTNCS.TRA.VECTOR,ncsRmat,PTNCS.GFUNCTION_RADIUS,miller);
    G_Vterm.resize(NREFL);
    for (int r = 0; r < NREFL; r++)
    {
      gfun.calcReflTerms(r);
      G_Vterm[r] = gfun.refl_G_Vterm; //not correct, B-factor differences
    }
  }
}

floatType DataMR::lowerB() { return -getWilsonB() + tolB(); }
floatType DataMR::upperB() { return 50; }
floatType DataMR::tolB()   { return 1.0e-04; }

Hexagonal DataMR::calc_search_points(data_btf& BTF,floatType SAMPLING,bool isFirst,Output& output)
{
  Hexagonal search_points(*this,*this);
  if (BTF.VOLUME.full())
  {
    if (isFirst) output.logTab(1,LOGFILE,"Full Search for first ensemble");
    else output.logTab(1,LOGFILE,"Full Search for second and subsequent ensemble");
    search_points.setupAsuRegion(SAMPLING,isFirst);
    output.logTab(2,VERBOSE,"Orthogonal x, y and z ranges of search:");
    output.logTab(3,VERBOSE,"Min: " + dvtos(search_points.getOrthMin(),5,2));
    output.logTab(3,VERBOSE,"Max: " + dvtos(search_points.getOrthMax(),5,2));
  }
  else if (BTF.VOLUME.region())
  {
    Brick brick(BTF.FRAC,BTF.START,BTF.END,*this);
    search_points.setupBrickRegion(brick,SAMPLING);
    output.logTab(2,LOGFILE,"Search in Region");
    output.logTab(2,VERBOSE,"Orthogonal x, y and z ranges of search:");
    output.logTab(3,VERBOSE,"Min: " + dvtos(search_points.getOrthMin(),5,2));
    output.logTab(3,VERBOSE,"Max: " + dvtos(search_points.getOrthMax(),5,2));
  }
  else if (BTF.VOLUME.around())
  {
    dvect3 fracPoint = BTF.FRAC ? BTF.START : doOrth2Frac(BTF.START);
    dvect3 orthPoint = BTF.FRAC ? doFrac2Orth(BTF.START) : BTF.START;
    search_points.setupAroundPoint(orthPoint,SAMPLING,BTF.RANGE);
    output.logTab(2,LOGFILE,"Search around Point");
    output.logTab(3,LOGFILE,"Fractional Coordinates: " + dvtos(fracPoint,5,2));
    output.logTab(3,LOGFILE,"Orthogonal Coordinates: " + dvtos(orthPoint));
    output.logTab(3,LOGFILE,"Coordinate Range:       " + dtos(BTF.RANGE) + "A");
  }
  else if (BTF.VOLUME.line())
  {
    output.logTab(1,LOGFILE,"Search along Line");
    dvect3 tran_start = BTF.FRAC ? doFrac2Orth(BTF.START) : BTF.START;
    dvect3 tran_end   = BTF.FRAC ? doFrac2Orth(BTF.END) : BTF.END;
    search_points.setupAlongLine(tran_start,tran_end,SAMPLING);
    output.logTab(2,VERBOSE,"Orthogonal x, y and z ranges of search:");
    output.logTab(3,VERBOSE,"Min: " + dvtos(tran_start,5,2));
    output.logTab(3,VERBOSE,"Max: " + dvtos(tran_end,5,2));
  }
  return search_points;
}

int1D DataMR::get_refl_negative_variance()
{
  int1D selr;
  for (int r = NREFL-1; r >= 0; r--)
  if (selected[r])
  {
    int s = rbin(r);
    if (s == 0 || s == bin.numbins()-1)
      selr.push_back(r);
  }
  std::reverse(selr.begin(),selr.end());
  int r1 = selr[selr.size()-1];
  selr.pop_back();
  std::reverse(selr.begin(),selr.end());
  selr.push_back(r1);
  std::reverse(selr.begin(),selr.end());
  bool failsafe_check(true);
  if (PTNCS.use_and_present() || //use all the reflections for tncs present
      failsafe_check) //use all reflections anyway, belt and braces
  {
    //append the middle reflections to the selected list
    int1D selr_tncs;
    for (int r = NREFL-1; r >= 0; r--)
    if (selected[r])
    {
      int s = rbin(r);
      if (s > 0 || s < bin.numbins()-1)
        selr_tncs.push_back(r);
    }
    for (int sr = 0; sr < selr_tncs.size(); sr++)
      selr.push_back(selr_tncs[sr]);
  }
  return selr;
}

bool DataMR::negative_variance(MapEnsemble& ensemble,int1D& selr,int& ir)
{
  if (!models_known.size()) return false;
  dmat33 ROT,PRtr,Rperturb;
  dvect3 CMori,CMrot,CMperturb;
  std::vector<dvect3> TRA(models_known.size());
  std::vector<dmat33> Q1tr(models_known.size());
  float1D scatFactor(models_known.size(),0);
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    Ensemble* ens_modlid = &ensemble.find(modlid)->second;
    scatFactor[k] = ens_modlid->ensScat(TOTAL_SCAT);
    //recover original centre of mass, rotated centre of mass
    PRtr = ens_modlid->principalRtr();
    CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = models_known[k].getR();
    CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform
    //orientation (for refinement)
    Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[k][dir])*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented ensemble)
    ROT = ROT*PRtr;

    TRA[k] = models_known[k].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA[k]);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[k];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of ensemble)
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA[k] = UnitCell::doOrth2Frac(eOrthT);

    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
           Q1tr[k] = Q1.transpose();
  }
  for (; ir < selr.size(); ir++)
  {
    unsigned r = selr[ir];
    floatType rssqr = ssqr(r);
    floatType repsn = epsn(r);
    ensemble.set_ssqr(rssqr);
    floatType this_totvar_known = totvar_known[r];
    cmplxType this_EM_known = EM_known[r];
    floatType this_frac_known_B = frac_known_B[r];
    totvar_known[r] = 0;
    EM_known[r] =  0;
    frac_known_B[r] = 0;
    //save this interpolation, so not doing multiple times for multiple copies
    float1D Vknown(models_known.size(),0);
    {
      map_str_float Vmap;
      for (unsigned k = 0; k < models_known.size(); k++)
      {
        std::string modlid = models_known[k].getModlid();
        if (Vmap.find(modlid) != Vmap.end())
        {
          Vknown[k] = Vmap[modlid];
        }
        else
        {
          Ensemble* ens_modlid = &ensemble.find(modlid)->second;
          floatType thisV = ens_modlid->InterpV();
                    thisV *= scatFactor[k];
          Vmap[modlid] = thisV;
          Vknown[k] = thisV;
        }
      }
    }
    for (unsigned k = 0; k < models_known.size(); k++)
    {
      std::string modlid = models_known[k].getModlid();
      Ensemble* ens_modlid = &ensemble.find(modlid)->second;
      floatType frac = scatFactor[k];
                frac /= models_known[k].getMult();
                frac *= -1 + std::exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
      frac_known_B[r] += frac;
      floatType thisV = Vknown[k];
                thisV *= exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
                thisV /= models_known[k].getMult();
      if (PTNCS.use_and_present()) thisV *= G_Vterm[r]; //not correct, B-factor differences
      totvar_known[r] += thisV;

      Miller1D rhkl = rotMiller(r);
      floatType terms  = std::sqrt(scatFactor[k]);
                terms *= sqrt_epsn[r];
                terms *= std::exp(-models_known[k].getBfac()*rssqr/4.0);
                terms /= models_known[k].getMult();
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3    RotSymHKL = Q1tr[k]*rhkl[isym];
          cmplxType thisE =  ens_modlid->InterpE(RotSymHKL);
                    thisE *= terms;
                    thisE *= dphi(r,isym,rhkl[isym],TRA[k]);
          //no G function correction (gfun), models in separate (real) rotations
          EM_known[r] += thisE;
        }
      }
    }
    Rice_refl(r,true);
    totvar_known[r] = this_totvar_known;
    EM_known[r] =  this_EM_known;
    frac_known_B[r] = this_frac_known_B;
    if (RiceV <=0)
    {
      return true;
    }
  }
  return false;
}

} //phaser
