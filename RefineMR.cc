//(c) 2000-2014 Cambridge University Technical Services Ltd
//All rights reserved
#include <phaser/src/RefineMR.h>
#include <phaser/io/Errors.h>
#include <phaser/lib/jiffy.h>
#include <phaser/lib/maths.h>
#include <phaser/lib/between.h>
#include <phaser/lib/xyzRotMatDeg.h>
#include <cctbx/sgtbx/search_symmetry.h>
#include <scitbx/math/bessel.h>

namespace phaser {

RefineMR::RefineMR(bool use_bfactor_restraint_, floatType sigmaB_,
                   DataMR& data,Output& output,MapEnsemble* ensPtr) : RefineBase2(), DataMR(data)
{
  ensemble = ensPtr;
  ensemble_modlid.clear();
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
    ensemble_modlid.push_back(iter->first);
  use_bfactor_restraint = use_bfactor_restraint_;
  sigmaB = sigmaB_;
  use_rotref_restraint = true;
  sigmaRotRef = 0.25;
  sigmaSolpar = 0.002; //tight //std::atof(getenv("PHASER_TEST"));
  models_vrms = ensemble->getDRMS();
  models_cell = ensemble->getCELL();
}

floatType RefineMR::targetFn()
{
// returns the -log(likelihood) of the current setup
  floatType totalLL = Rice();
  floatType f = -(totalLL - LLwilson);
            f += restraint_term();
  return f;
}

floatType RefineMR::restraint_term()
{
  int m(0),i(0);
  floatType restraints(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (unsigned rot = 0; rot < 3; rot++) if (refinePar[m++]) i++;
    for (unsigned tra = 0; tra < 3; tra++) if (refinePar[m++]) i++;
    if (refinePar[m++])
    {
      if (use_bfactor_restraint)
      {
        restraints += fn::pow2(models_known[s].getBfac())/sigmaB/2.;
      }
      i++;
    }
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) i++; //VRMS
    if (refinePar[m++]) i++; //CELL
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);

  //localMLL += restraints;
  return restraints;
}

floatType RefineMR::gradientFn(TNT::Vector<floatType>& Gradient)
{
//#define PHASER_TIMINGS
#ifdef PHASER_TIMINGS
  static int countr(0);
  Output output;
  output.logTab(1,LOGFILE,"Gradient Start # " + itos(countr++));
  floatType start_clock = std::clock();
#endif
  //initialization
  Gradient.newsize(npars_ref);
  bool return_gradient(true);
  floatType localMLL = gradient_hessian(Gradient,return_gradient); //-(totalLL - LLwilson) + restraint
                       
  int m(0),i(0);
  int1D Bindex(models_known.size());
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (unsigned rot = 0; rot < 3; rot++) if (refinePar[m++]) i++;
    for (unsigned tra = 0; tra < 3; tra++) if (refinePar[m++]) i++;
    if (refinePar[m++])
    {
      Bindex[s] = i;
      if (use_bfactor_restraint)
      {
        Gradient[i] += models_known[s].getBfac()/sigmaB;
      }
      i++;
    }
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) i++; //VRMS
    if (refinePar[m++]) i++; //CELL
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);

  //average Gradient
  if (PTNCS.use_and_present())
  {
    float1D avBfac(Bindex.size(),0);
    int n(-1);
    for (unsigned s = 0; s < Bindex.size(); s++) 
    {
      if (s % PTNCS.NMOL == 0) n++;
      avBfac[n] += Gradient[Bindex[s]]/PTNCS.NMOL;
    }
    n = -1;
    for (unsigned s = 0; s < Bindex.size(); s++) 
    {
      if (s % PTNCS.NMOL == 0) n++;
      Gradient[Bindex[s]] = avBfac[n];
    } 
  } 

#ifdef PHASER_TIMINGS
  std::clock_t now_clock = std::clock();
  floatType time_for_sum = (now_clock-start_clock)/double(CLOCKS_PER_SEC);
  output.logTab(1,LOGFILE,"Gradient Time = " + dtos(time_for_sum) + " secs");
#endif
  return localMLL;
}

floatType RefineMR::hessianFn(TNT::Fortran_Matrix<floatType>& Hessian,bool& is_diagonal)
{
#ifdef PHASER_TIMINGS
  static int countr(0);
  Output output;
  output.logTab(1,LOGFILE,"Hessian Start # " + itos(countr++));
  floatType start_clock = std::clock();
#endif

  //initialization
  TNT::Vector<floatType> diag_hessian(npars_ref);
  bool return_gradient(false);
  floatType localMLL = gradient_hessian(diag_hessian,return_gradient); //-(totalLL - LLwilson) + restraint
  is_diagonal = true;

  Hessian.newsize(npars_ref,npars_ref);
  PHASER_ASSERT(Hessian.num_cols() == npars_ref);
  PHASER_ASSERT(Hessian.num_rows() == npars_ref);
  for (int ii = 0; ii < npars_ref; ii++) //not i, warning on visual-C, int i defined below
    for (int jj = 0; jj < npars_ref; jj++)
      Hessian(ii+1,jj+1) =  (ii==jj) ? diag_hessian[jj] : 0;

  int m(0),i(0);
  int1D Bindex(models_known.size());
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (unsigned rot = 0; rot < 3; rot++) if (refinePar[m++]) i++;
    for (unsigned tra = 0; tra < 3; tra++) if (refinePar[m++]) i++;
    if (refinePar[m++])
    {
      Bindex[s] = i;
      if (use_bfactor_restraint)
      {
        Hessian(i+1,i+1) += 1/sigmaB;
      }
      i++;
    }
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) i++; //VRMS
    if (refinePar[m++]) i++; //CELL
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);

  //average Hessian
  if (PTNCS.use_and_present())
  {
    float1D avBfac(Bindex.size(),0);
    int n(-1);
    for (unsigned s = 0; s < Bindex.size(); s++) 
    {
      if (s % PTNCS.NMOL == 0) n++;
      avBfac[n] += Hessian(Bindex[s]+1,Bindex[s]+1)/PTNCS.NMOL;
    }
    n = -1;
    for (unsigned s = 0; s < Bindex.size(); s++) 
    {
      if (s % PTNCS.NMOL == 0) n++;
      Hessian(Bindex[s]+1,Bindex[s]+1) = avBfac[n];
    } 
  }

#ifdef PHASER_TIMINGS
  std::clock_t now_clock = std::clock();
  floatType time_for_sum = (now_clock-start_clock)/double(CLOCKS_PER_SEC);
  output.logTab(1,LOGFILE,"Hessian Time = " + dtos(time_for_sum) + " secs");
#endif
  return localMLL;
}

floatType RefineMR::gradient_hessian(TNT::Vector<floatType>& Array,bool return_gradient)
{
  dvect31D dLL_by_drot(models_known.size(),dvect3(0,0,0));
  dvect31D dLL_by_dtra(models_known.size(),dvect3(0,0,0));
  float1D dLL_by_dB(models_known.size(),0);
  map_str_float dLL_by_dvrms;
  map_str_float dLL_by_dcell;
  for (int e = 0; e < ensemble_modlid.size(); e++) dLL_by_dvrms[ensemble_modlid[e]] = 0;
  for (int e = 0; e < ensemble_modlid.size(); e++) dLL_by_dcell[ensemble_modlid[e]] = 0;
  dvect31D d2LL_by_drot2(models_known.size(),dvect3(0,0,0));
  dvect31D d2LL_by_dtra2(models_known.size(),dvect3(0,0,0));
  float1D d2LL_by_dB2(models_known.size(),0);
  map_str_float d2LL_by_dvrms2;
  map_str_float d2LL_by_dcell2;
  for (int e = 0; e < ensemble_modlid.size(); e++) d2LL_by_dvrms2[ensemble_modlid[e]] = 0;
  for (int e = 0; e < ensemble_modlid.size(); e++) d2LL_by_dcell2[ensemble_modlid[e]] = 0;
  floatType totalLL(0);

  cvect31D dFC_by_drot(models_known.size()*NREFL,cvect3(0,0,0));
  cvect31D dFC_by_dtra(models_known.size()*NREFL,cvect3(0,0,0));
  cmplx1D  dFC_by_dcell(models_known.size()*NREFL,0);
  cvect31D d2FC_by_drot2(models_known.size()*NREFL,cvect3(0,0,0));
  cvect31D d2FC_by_dtra2(models_known.size()*NREFL,cvect3(0,0,0));
  cmplx1D  d2FC_by_dcell2(models_known.size()*NREFL,0);
  cmplx1D  Fcalc(models_known.size()*NREFL,cmplxType(0,0));

  if (!FIX_ROT || !FIX_TRA || !FIX_BFAC || !FIX_VRMS || !FIX_CELL)
  for (unsigned s = 0; s < models_known.size(); s++)
  {
    std::string modlid = models_known[s].getModlid();
    if (FAST_LAST && modlid != LAST_MODLID) continue;
    Ensemble* ens_modlid = &ensemble->find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    dmat33 PRtr = ens_modlid->principalRtr();
    dvect3 CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = models_known[s].getR();
    dvect3 CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform 
    //orientation (for refinement)
    dmat33 Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
    {
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[s][dir])*Rperturb;
    }
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    dvect3 CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented (*ensemble))
    ROT = ROT*PRtr;

    dvect3 TRA = models_known[s].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[s];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of (*ensemble))
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA = UnitCell::doOrth2Frac(eOrthT);

    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr = Q1.transpose();

    dmat33 Rtr = models_known[s].getR().transpose();
    dmat33 PR_Rtr_Dtr = ens_modlid->principalR() * Rtr * UnitCell::Orth2Frac().transpose();

    floatType sqrt_scatFactor = std::sqrt(ens_modlid->ensScat(TOTAL_SCAT)/NSYMP);
    bool return_hessian = !return_gradient;
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      int sr = s*NREFL + r;
      floatType repsn = epsn(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      floatType sqrt_repsn = sqrt_epsn[r];
      Miller1D rhkl = rotMiller(r);
      for (unsigned isym = 0; isym < SpaceGroup::NSYMP; isym++)
      {
        if (!duplicate(isym,rhkl))
        {
          dvect3    RotSymHKL = Q1tr*rhkl[isym];
          cvect3    dthisE(0,0,0);
          cmat33    hthisE(0,0,0,0,0,0,0,0,0);
          cmplxType interpE = ens_modlid->func_grad_hess_InterpE(RotSymHKL,dthisE,return_hessian,hthisE);
          floatType Bfac = models_known[s].getBfac();
          cmplxType scalefac = sqrt_repsn * std::exp(-Bfac*rssqr/4.0) * dphi(r,isym,rhkl[isym],TRA);
                    scalefac *= ens_modlid->Eterm_const();
                    scalefac *= sqrt_scatFactor;
                    scalefac /= models_known[s].getMult();
          cmplxType thisE  = scalefac * interpE;
                    dthisE = scalefac * dthisE;
                    hthisE = scalefac * hthisE;
                    Fcalc[sr] += thisE;
          dvect3 dQ1tr_h_by_drot_s(0,0,0);
          dmat33 dQ1tr_by_drot(0,0,0,0,0,0,0,0,0);
          for (int rot = 0; rot < 3; rot++)
          {
            dQ1tr_by_drot = PR_Rtr_Dtr;
            for (int dir = 2; dir >= 0; dir--)
            {
              dQ1tr_by_drot = xyzRotMatDeg(dir,models_perturbRot[s][dir],dir==rot).transpose()*dQ1tr_by_drot;
            }
            dQ1tr_by_drot = ens_modlid->Frac2Orth().transpose() * dQ1tr_by_drot;
            dQ1tr_h_by_drot_s = dQ1tr_by_drot*rhkl[isym];
            dFC_by_drot[sr][rot] += dthisE*dQ1tr_h_by_drot_s;
            if (return_hessian)
            {
              d2FC_by_drot2[sr][rot] += dQ1tr_h_by_drot_s*(hthisE*dQ1tr_h_by_drot_s);
            }
          }
          floatType cell = ens_modlid->getScaleCELL();
          dmat33 dQ1tr_by_dcell;
          for (int q = 0; q < Q1tr.size(); q++) dQ1tr_by_dcell[q] = Q1tr[q]/cell;
          dvect3 dQ1tr_h_by_dcell_s = dQ1tr_by_dcell*rhkl[isym];
          dFC_by_dcell[sr] += dthisE*dQ1tr_h_by_dcell_s; //scalar = vector.vector
          if (return_hessian)
            d2FC_by_dcell2[sr] += dQ1tr_h_by_dcell_s*(hthisE*dQ1tr_h_by_dcell_s);
          for (int tra = 0; tra < 3; tra++)
          {
            //dphi = rhkl.TRA + traMiller
            dvect3 dTRA_by_du(tra==0 ? 1 : 0, tra==1 ? 1 : 0, tra==2 ? 1 : 0);
            dvect3 dhkl(rhkl[isym][0],rhkl[isym][1],rhkl[isym][2]);
            cmplxType C = scitbx::constants::two_pi*cmplxType(0,1)*(dhkl*(UnitCell::Orth2Frac()*dTRA_by_du));
            dFC_by_dtra[sr][tra] += C*thisE;
            cmplxType C2 = fn::pow2(C);
            if (return_hessian)
            {
               d2FC_by_dtra2[sr][tra] += C2*thisE;
            }
          }
        }
      }
    }
  }

  int countr(0);
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    countr++;
    floatType reflLL(0);
    floatType fracB = 1 + frac_known_B[r] + frac_search_B[r];
    PHASER_ASSERT(fracB > 0);
    cmplxType EM(EM_known[r]+EM_search[r]);
    floatType EM_real = EM.real();
    floatType EM_imag = EM.imag();
    floatType phasedEsqr = EM_real*EM_real+EM_imag*EM_imag;
    floatType sumEsqr = phasedEsqr + sum_Esqr_search[r];
    floatType maxEsqr = std::max(phasedEsqr,max_Esqr_search[r]);

    //Now do likelihood calculation
    floatType Vterm = 0;
    Vterm = PTNCS.EPSFAC[r];
    Vterm -= totvar_known[r];
    Vterm -= totvar_search[r];
    Vterm += sumEsqr;
    Vterm -= maxEsqr;
    floatType V = Vterm/fracB;

    bool rcent = cent(r);
    floatType E = F[r]/sqrt_epsnSigmaN[r];
    floatType Esqr = fn::pow2(E);
    floatType wtVarE = fn::pow2(SIGF[r]/sqrt_epsnSigmaN[r]); 
    if (!rcent) wtVarE *= 2.;
    V += wtVarE; //centric correction has already been included
    PHASER_ASSERT(V > 0);

    floatType fc = std::sqrt(maxEsqr); //sqrt(FA*FA+FB*FB)
    floatType FCsqr = maxEsqr/fracB; //fc^2/fracB
    floatType FC = std::sqrt(FCsqr);
    floatType X = 2.0*E*FC/V;
    reflLL = -(std::log(V)+(Esqr + FCsqr)/V);
    if (rcent)
    {
      X /= 2.0;
      reflLL = reflLL/2.0 + m_alogch.getalogch(X);
    }
    else
    {
      reflLL += m_alogchI0.getalogI0(X); //acentric
    }
    totalLL += reflLL;

    if (!FIX_ROT || !FIX_TRA || !FIX_BFAC || !FIX_VRMS || !FIX_CELL)
    {
    dvect3 dLL_by_drot_s(0,0,0),d2LL_by_drot2_s(0,0,0);
    dvect3 dLL_by_dtra_s(0,0,0),d2LL_by_dtra2_s(0,0,0);
    floatType dLL_by_dB_s(0),d2LL_by_dB2_s(0);
    floatType dLL_by_dvrms_s(0),d2LL_by_dvrms2_s(0);
    floatType dLL_by_dcell_s(0),d2LL_by_dcell2_s(0);

//#define PHASER_TEST_MR_HESS
#ifdef PHASER_TEST_MR_HESS
if (getenv("PHASER_TEST_MATHEMATICA") != 0)
{
if (std::string(getenv("PHASER_TEST_MATHEMATICA")) == "centric") rcent = true;
if (std::string(getenv("PHASER_TEST_MATHEMATICA")) == "acentric") rcent = false;
}
 if (getenv("PHASER_TEST_SHIFT") == 0)
    { std::cout << "setenv PHASER_TEST_SHIFT 0.0001\n" ; std::exit(1); }
 if (getenv("PHASER_TEST_NREFL") == 0)
    { std::cout << "setenv PHASER_TEST_NREFL 1\n" ; std::exit(1); }
  floatType shift_a = std::atof(getenv("PHASER_TEST_SHIFT"));
  floatType FD_test_on_f = 0;
  floatType FD_test_forward = 0;
  floatType FD_test_backward = 0;
//#define PHASER_TEST_MR_HESS_BFAC
#ifdef PHASER_TEST_MR_HESS_BFAC
  floatType orig_a = models_known[0].getBfac();
  floatType &dL_by_dtest = dLL_by_dB_s;
  floatType &d2L_by_dtest2 = d2LL_by_dB2_s;
#else
  floatType *test_a = &(models_perturbRot[0][0]);  //or Trans
  floatType &dL_by_dtest = dLL_by_drot_s[0]; //or tra
  floatType &d2L_by_dtest2 = d2LL_by_drot2_s[0];
#endif
if (!return_gradient)
{
  calcKnown(*ensemble);  
  calcKnownVAR(*ensemble); //totvar_known
  calcKnownB(*ensemble); //frac_known_B
  floatType f_start = Ricer(r);
#ifdef PHASER_TEST_MR_HESS_BFAC
  models_known[0].setBfac(orig_a+shift_a);
#else
  (*test_a) += shift_a;
#endif
  calcKnown(*ensemble);
  calcKnownVAR(*ensemble); //totvar_known
  calcKnownB(*ensemble); //frac_known_B
  floatType f_forward = Ricer(r);
  FD_test_forward = (f_forward-f_start)/shift_a;
#ifdef PHASER_TEST_MR_HESS_BFAC
  models_known[0].setBfac(orig_a-shift_a);
#else
  (*test_a) -= shift_a;
  (*test_a) -= shift_a;
#endif
  calcKnown(*ensemble);
  calcKnownVAR(*ensemble); //totvar_known
  calcKnownB(*ensemble); //frac_known_B
  floatType f_backward = Ricer(r);
  FD_test_backward = (f_backward-f_start)/shift_a;
#ifdef PHASER_TEST_MR_HESS_BFAC
  models_known[0].setBfac(orig_a);
#else
  (*test_a) += shift_a;
#endif
  calcKnown(*ensemble);
  calcKnownVAR(*ensemble); //totvar_known
  calcKnownB(*ensemble); //frac_known_B
  FD_test_on_f = (f_forward-2*f_start+f_backward)/shift_a/shift_a;
}
#endif


    floatType bessTerm = rcent ? std::tanh(X) : scitbx::math::bessel::i1_over_i0(X);
    floatType bessTerm2 = rcent ? fn::pow2(1/std::cosh(X)) : fn::pow2(bessTerm);
    // NB: sech = 1/cosh
    floatType FA = EM_known[r].real() + EM_search[r].real();
    floatType FB = EM_known[r].imag() + EM_search[r].imag();
    floatType V2 = fn::pow2(V);
    floatType V3 = fn::pow3(V);
    floatType V4 = fn::pow4(V);
    floatType fc2 = fn::pow2(fc);
    floatType fc4 = fn::pow2(fc2);
    floatType fracB2 = fn::pow2(fracB);
    floatType fracB3 = fn::pow3(fracB);

    floatType dLL_by_dFA(0),dLL_by_dFB(0);
    floatType d2LL_by_dFA2(0),d2LL_by_dFB2(0),d2LL_by_dFA_dFB(0);
    floatType d2LL_by_dFC2(0);

    floatType dLL_by_dFC = rcent ? (1/V)*(E*bessTerm - FC) : (2/V)*(E*bessTerm - FC); 

    floatType dFC_by_dFA = 0;
    floatType dFC_by_dFB = 0;
    if (FC > 0)
    {
      dFC_by_dFA = FA/FC/fracB;
      dFC_by_dFB = FB/FC/fracB;
    }
    dLL_by_dFA = dLL_by_dFC*dFC_by_dFA;
    dLL_by_dFB = dLL_by_dFC*dFC_by_dFB;

    floatType denom = fc4*fracB*V2;
    if (!return_gradient && fc2 > 0 && fc4 > 0 && denom > 0)
    {
      d2LL_by_dFA2 = rcent ? 
        bessTerm2*FA*FA*Esqr/fc2/fracB/V2 - bessTerm*FA*FA*FC*E/fc4/V + bessTerm*FC*E/fc2/V - 1/fracB/V :
        (-2*fc2*(fc2-bessTerm*FC*E*fracB)*V-4*FA*FA*E*((-1 + bessTerm2)*fc2*E + bessTerm*FC*fracB*V))/denom;
      d2LL_by_dFB2 = rcent ? 
        bessTerm2*FB*FB*Esqr/fc2/fracB/V2 - bessTerm*FB*FB*FC*E/fc4/V + bessTerm*FC*E/fc2/V - 1/fracB/V :
        (-2*fc2*(fc2-bessTerm*FC*E*fracB)*V-4*FB*FB*E*((-1 + bessTerm2)*fc2*E + bessTerm*FC*fracB*V))/denom;
      d2LL_by_dFA_dFB = rcent ? 
        bessTerm2*FA*FB*Esqr/fc2/fracB/V2 - bessTerm*FA*FB*FC*E/fc4/V :
        -4*FA*FB*E*((-1 + bessTerm2)*fc2*E+bessTerm*FC*fracB*V)/denom;
      d2LL_by_dFC2 = rcent ?
              4*Esqr/V2 - 4*bessTerm2*Esqr/V3 -2/V-2*bessTerm*E/FC/V :
              bessTerm2*Esqr/V2 - 1/V;
    }

    floatType dLL_by_dV(0),dLL_by_dfracB(0);
    floatType d2LL_by_dV2(0),d2LL_by_dfracB2(0);
    floatType d2LL_by_dFA_dV(0),d2LL_by_dFA_dfracB(0);
    floatType d2LL_by_dFB_dV(0),d2LL_by_dFB_dfracB(0);
    floatType d2LL_by_dV_dfracB(0);
    if (!FIX_BFAC || !FIX_VRMS)
    {
      dLL_by_dV = rcent ?
         (FCsqr - 2*bessTerm*FC*E + Esqr - V)/2/V2 :
         (1/V/V)*((Esqr + FCsqr - V) - 2*E*FC*bessTerm); 
      dLL_by_dfracB = rcent ?
         (fc2 - bessTerm*FC*E*fracB)/2/fracB2/V :
         (1/fracB/V)*(FCsqr - E*FC*bessTerm);

      if (!return_gradient && V3 > 0 && V4 > 0)
      {
        d2LL_by_dV2 = rcent ?
           (2*bessTerm2*FCsqr*Esqr + V*(-2*FCsqr + 4*bessTerm*FC*E - 2*Esqr + V))/(2*V4) :
           (4*Esqr*FCsqr/V4 - 2*(FCsqr+Esqr)/V3 + 1/V2 + 2*FC*E*bessTerm/V3 - 4*Esqr*FCsqr*bessTerm2/V4);
        d2LL_by_dfracB2 =  rcent ?
          (-4*fc2*V + fc2*Esqr*bessTerm2 + 3*FC*E*fracB*V*bessTerm)/(4*fracB3*V2) :
          1/fracB3/V2*(fc2*(-(-1 + bessTerm2)*Esqr - 2*V) + bessTerm*FC*E*fracB*V);
        if (fc2 > 0)
        {
        d2LL_by_dFA_dV =  rcent ?
           FA/fracB/V2 - FA*Esqr*bessTerm2/fracB/V3 - FA*FC*E*bessTerm/fc2/V2  :
          -4*FA*Esqr/fracB/V3 + 2*FA/fracB/V2 + 4*FA*Esqr*bessTerm2/fracB/V3;
        d2LL_by_dFB_dV =  rcent ?
           FB/fracB/V2 - FB*Esqr*bessTerm2/fracB/V3 - FB*FC*E*bessTerm/fc2/V2  :
          -4*FB*Esqr/fracB/V3 + 2*FB/fracB/V2 + 4*FB*Esqr*bessTerm2/fracB/V3;
        d2LL_by_dFA_dfracB =  rcent ?
          FA/2/fc2/fracB2/V2*(-fc2*Esqr*bessTerm2 + V*(2*fc2 - FC*E*fracB*bessTerm)) :
          2*FA*((-1 + bessTerm2)*Esqr+V)/fracB2/V2;
        d2LL_by_dFB_dfracB =  rcent ?
          FB/2/fc2/fracB2/V2*(-fc2*Esqr*bessTerm2 + V*(2*fc2 - FC*E*fracB*bessTerm)) :
          2*FB*((-1 + bessTerm2)*Esqr+V)/fracB2/V2;
        }
        d2LL_by_dV_dfracB =  rcent ?   
           (-fc2*V + fc2*Esqr*bessTerm2 + FC*E*fracB*V*bessTerm)/(2*fracB2*V3):
           -fc2/fracB2/V3*(2*(-1+bessTerm2)*Esqr+V);
      }
    }

    for (int s = 0; s < models_known.size(); s++)
    {
      std::string modlid = models_known[s].getModlid();
      if (FAST_LAST && modlid != LAST_MODLID) continue;
      Ensemble* ens_modlid = &ensemble->find(modlid)->second;
      int sr = s*NREFL + r;
      if (!FIX_ROT)
      {
        dvect3 dFA_by_drot_s(dFC_by_drot[sr][0].real(),dFC_by_drot[sr][1].real(),dFC_by_drot[sr][2].real());
        dvect3 dFB_by_drot_s(dFC_by_drot[sr][0].imag(),dFC_by_drot[sr][1].imag(),dFC_by_drot[sr][2].imag());
               dLL_by_drot_s = dLL_by_dFA*dFA_by_drot_s + dLL_by_dFB*dFB_by_drot_s;
        dLL_by_drot[s] += dLL_by_drot_s;

        if (!return_gradient)
        {
          for (int rot = 0; rot < 3; rot++)
          {
            floatType d2FA_by_drot2_s = (d2FC_by_drot2[sr][rot].real());
            floatType d2FB_by_drot2_s = (d2FC_by_drot2[sr][rot].imag());
            d2LL_by_drot2_s[rot] = d2LL_by_dFA2*fn::pow2(dFA_by_drot_s[rot]) +
                                   d2LL_by_dFB2*fn::pow2(dFB_by_drot_s[rot]) +
                                   d2FA_by_drot2_s*dLL_by_dFA +
                                   d2FB_by_drot2_s*dLL_by_dFB +
                                   2*d2LL_by_dFA_dFB*dFA_by_drot_s[rot]*dFB_by_drot_s[rot];
          }
          d2LL_by_drot2[s] += d2LL_by_drot2_s;
        }
      }
      if (!FIX_TRA)
      {
        dvect3 dFA_by_dtra_s(dFC_by_dtra[sr][0].real(),dFC_by_dtra[sr][1].real(),dFC_by_dtra[sr][2].real());
        dvect3 dFB_by_dtra_s(dFC_by_dtra[sr][0].imag(),dFC_by_dtra[sr][1].imag(),dFC_by_dtra[sr][2].imag());
               dLL_by_dtra_s = dLL_by_dFA*dFA_by_dtra_s + dLL_by_dFB*dFB_by_dtra_s;
        dLL_by_dtra[s] += dLL_by_dtra_s;

        if (!return_gradient)
        {
          dvect3 d2FA_by_dtra2_s(d2FC_by_dtra2[sr][0].real(),d2FC_by_dtra2[sr][1].real(),d2FC_by_dtra2[sr][2].real());
          dvect3 d2FB_by_dtra2_s(d2FC_by_dtra2[sr][0].imag(),d2FC_by_dtra2[sr][1].imag(),d2FC_by_dtra2[sr][2].imag());
          for (int tra = 0; tra < 3; tra++)
            d2LL_by_dtra2_s[tra] = d2LL_by_dFA2*fn::pow2(dFA_by_dtra_s[tra]) +
                                   d2LL_by_dFB2*fn::pow2(dFB_by_dtra_s[tra]) +
                                   dLL_by_dFA*d2FA_by_dtra2_s[tra] +
                                   dLL_by_dFB*d2FB_by_dtra2_s[tra] +
                                   2*d2LL_by_dFA_dFB*dFA_by_dtra_s[tra]*dFB_by_dtra_s[tra];
          d2LL_by_dtra2[s] += d2LL_by_dtra2_s;
        }
      }

      if (!FIX_BFAC || !FIX_VRMS)
      {
        floatType rssqr = ssqr(r);
        floatType sqrtrssqr = std::sqrt(rssqr);
        ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
        floatType Fconst = -rssqr/4.0;
        floatType Bconst = -rssqr/2.0;
        floatType Vconst = ens_modlid->dvrms_delta_const();
        floatType FA_s = Fcalc[sr].real();
        floatType FB_s = Fcalc[sr].imag();
 
        floatType dFA_by_dB(0),dFB_by_dB(0),dV_by_dB(0),dfracB_by_dB(0),frac_known_B_s(0);
        floatType dFA_by_dvrms(0),dFB_by_dvrms(0),dV_by_dvrms(0);
        floatType thisV = ens_modlid->InterpV();
                  thisV *= ens_modlid->ensScat(TOTAL_SCAT);
                  thisV *= exp(-2.0*models_known[s].getBfac()*rssqr/4.0);
                  thisV /= models_known[s].getMult();
        floatType totvar_known_s = thisV;
  
        if (!FIX_BFAC)
        {
                    dFA_by_dB = Fconst*FA_s;
                    dFB_by_dB = Fconst*FB_s;
          floatType frac = ens_modlid->ensScat(TOTAL_SCAT);
                    frac *= std::exp(-2.0*models_known[s].getBfac()*rssqr/4.0);
                    frac /= models_known[s].getMult();
                    frac_known_B_s = frac;
          floatType dfracB_by_dB = Bconst*frac_known_B_s;
          floatType dtotvar_by_dB = Bconst*totvar_known_s;
  
                    dV_by_dB = -dfracB_by_dB*Vterm/fracB2 - dtotvar_by_dB/fracB;
                    dLL_by_dB_s = dLL_by_dFA*dFA_by_dB + dLL_by_dFB*dFB_by_dB + dLL_by_dV*dV_by_dB + dLL_by_dfracB*dfracB_by_dB;
          dLL_by_dB[s] += dLL_by_dB_s;

          if (!return_gradient)
          {
            floatType d2FA_by_dB2 = Fconst*dFA_by_dB;
            floatType d2FB_by_dB2 = Fconst*dFB_by_dB;
            floatType d2fracB_by_dB2 = Bconst*dfracB_by_dB;
            floatType d2V_by_dB2 = Bconst*Bconst*(fracB-2*frac_known_B_s)*(-fracB*totvar_known_s-frac_known_B_s*Vterm);
                      d2V_by_dB2 /= fracB3;
  
                      d2LL_by_dB2_s = d2FA_by_dB2*dLL_by_dFA +
                                      fn::pow2(dFA_by_dB)*d2LL_by_dFA2 +
                                      d2FB_by_dB2*dLL_by_dFB +
                                      2*dFA_by_dB*dFB_by_dB*d2LL_by_dFA_dFB +
                                      fn::pow2(dFB_by_dB)*d2LL_by_dFB2 +
                                      d2V_by_dB2*dLL_by_dV +
                                      2*dFA_by_dB*dV_by_dB*d2LL_by_dFA_dV +
                                      2*dFB_by_dB*dV_by_dB*d2LL_by_dFB_dV +
                                      fn::pow2(dV_by_dB)*d2LL_by_dV2 +
                                      d2fracB_by_dB2*dLL_by_dfracB +
                                      2*dFA_by_dB*dfracB_by_dB*d2LL_by_dFA_dfracB +
                                      2*dFB_by_dB*dfracB_by_dB*d2LL_by_dFB_dfracB +
                                      2*dV_by_dB*dfracB_by_dB*d2LL_by_dV_dfracB +
                                      fn::pow2(dfracB_by_dB)*d2LL_by_dfracB2;
            d2LL_by_dB2[s] += d2LL_by_dB2_s;
          }
        }
    
        if (!FIX_VRMS)
        {
          dFA_by_dvrms = Vconst*FA_s;
          dFB_by_dvrms = Vconst*FB_s;
          dV_by_dvrms = -2*Vconst*totvar_known_s; //EPSFAC - thisV term
          dLL_by_dvrms_s = dLL_by_dFA*dFA_by_dvrms + dLL_by_dFB*dFB_by_dvrms + dLL_by_dV*dV_by_dvrms;
          dLL_by_dvrms[modlid] += dLL_by_dvrms_s;

          if (!return_gradient)
          {
            floatType d2FA_by_dvrms2 = Vconst*dFA_by_dvrms;
            floatType d2FB_by_dvrms2 = Vconst*dFB_by_dvrms;
            floatType d2V_by_dvrms2 =  2*Vconst*dV_by_dvrms;
                   d2LL_by_dvrms2_s = d2FA_by_dvrms2*dLL_by_dFA +
                                      fn::pow2(dFA_by_dvrms)*d2LL_by_dFA2 +
                                      d2FB_by_dvrms2*dLL_by_dFB +
                                      2*dFA_by_dvrms*dFB_by_dvrms*d2LL_by_dFA_dFB +
                                      fn::pow2(dFB_by_dvrms)*d2LL_by_dFB2 +
                                      d2V_by_dvrms2*dLL_by_dV +
                                      2*dFA_by_dvrms*dV_by_dvrms*d2LL_by_dFA_dV +
                                      2*dFB_by_dvrms*dV_by_dvrms*d2LL_by_dFB_dV +
                                      fn::pow2(dV_by_dvrms)*d2LL_by_dV2;
            d2LL_by_dvrms2[modlid] += d2LL_by_dvrms2_s;
          }
        }
      } //loop over known

      if (!FIX_CELL)
      {
        floatType dFA_by_dcell_s = dFC_by_dcell[sr].real();
        floatType dFB_by_dcell_s = dFC_by_dcell[sr].imag();
        floatType dLL_by_dcell_s = dLL_by_dFA*dFA_by_dcell_s + dLL_by_dFB*dFB_by_dcell_s;
        dLL_by_dcell[modlid] += dLL_by_dcell_s;

        if (!return_gradient)
        {
          floatType d2FA_by_dcell2_s = d2FC_by_dcell2[sr].real();
          floatType d2FB_by_dcell2_s = d2FC_by_dcell2[sr].imag();
                    d2FA_by_dcell2_s = d2FA_by_dcell2_s/2;
                    d2FB_by_dcell2_s = d2FB_by_dcell2_s/2;
                 d2LL_by_dcell2_s = d2FA_by_dcell2_s*dLL_by_dFA +
                                    fn::pow2(dFA_by_dcell_s)*d2LL_by_dFA2 +
                                    d2FB_by_dcell2_s*dLL_by_dFB +
                                    2*dFA_by_dcell_s*dFB_by_dcell_s*d2LL_by_dFA_dFB +
                                    fn::pow2(dFB_by_dcell_s)*d2LL_by_dFB2;
          d2LL_by_dcell2[modlid] += d2LL_by_dcell2_s/2;
        }
      }

    }
    }

#ifdef PHASER_TEST_MR_HESS
if (!return_gradient)
{
  std::cout << "=== First Derivative test (r=" << r << ") ===\n";
  std::cout << "on function " << dL_by_dtest  << " fd forward "  << FD_test_forward << " fd backward " << FD_test_backward << "\n";
  PHASER_ASSERT(FD_test_forward);
  PHASER_ASSERT(FD_test_backward);
  std::cout << "Ratio forward "  << dL_by_dtest/FD_test_forward <<
                    " backward " << dL_by_dtest/FD_test_backward << "\n";
  std::cout << "=== Second Derivative test (r=" << r << ") ===\n";
  PHASER_ASSERT(FD_test_on_f);
  std::cout <<"on function " << d2L_by_dtest2 << " fd " << FD_test_on_f << "\nRatio " << d2L_by_dtest2/FD_test_on_f << "\n";
  if (countr >= std::atof(getenv("PHASER_TEST_NREFL"))) std::exit(1);
}
#endif
  }

  floatType f = -(totalLL - LLwilson);
            f += restraint_term();

  int m(0),i(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (unsigned rot = 0; rot < 3; rot++) 
      if (refinePar[m++]) Array[i++] = return_gradient ? -dLL_by_drot[s][rot] : -d2LL_by_drot2[s][rot];
    for (unsigned tra = 0; tra < 3; tra++) 
      if (refinePar[m++]) Array[i++] = return_gradient ? -dLL_by_dtra[s][tra] : -d2LL_by_dtra2[s][tra];
    if (refinePar[m++]) Array[i++] = return_gradient ? -dLL_by_dB[s] : -d2LL_by_dB2[s];
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) Array[i++] = return_gradient ? -dLL_by_dvrms[ensemble_modlid[e]] : -d2LL_by_dvrms2[ensemble_modlid[e]];
    if (refinePar[m++]) Array[i++] = return_gradient ? -dLL_by_dcell[ensemble_modlid[e]] : -d2LL_by_dcell2[ensemble_modlid[e]];
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);

  return f;
}

TNT::Vector<floatType> RefineMR::getLargeShifts()
{
  //relative sizes of these are important in finite diff grad calculation
  TNT::Vector<floatType> largeShifts(npars_ref);
  floatType dmin(HiRes());
  int m(0),i(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    std::string modlid = models_known[s].getModlid();
    Ensemble* ens_modlid = &ensemble->find(modlid)->second;
    for (unsigned rot = 0; rot < 3; rot++) 
      if (refinePar[m++])
      { 
        std::string modlid = models_known[s].getModlid();
        floatType maxperp(0);
        if (rot == 0)
          maxperp = std::max(ens_modlid->B(),ens_modlid->C());
        else if (rot == 1)
          maxperp = std::max(ens_modlid->C(),ens_modlid->A());
        else
          maxperp = std::max(ens_modlid->A(),ens_modlid->B());
        largeShifts[i++] = 4*dmin/scitbx::deg_as_rad(maxperp);
      }
    for (unsigned tra = 0; tra < 3; tra++) 
      if (refinePar[m++]) largeShifts[i++] = dmin/4.0;
    if (refinePar[m++]) largeShifts[i++] = 4.; // or higher
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) largeShifts[i++] = 1.0; //VRMS
    if (refinePar[m++]) largeShifts[i++] = 0.1; //CELL
  }
  //if smaller than 0.5 toxd will not refine ROT and TRA of 2nd peak to unique solution
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);
  return largeShifts;
}
 
void  RefineMR::applyShift(TNT::Vector<floatType>& newx)
{
  int i(0),m(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) models_perturbRot[s][rot] = newx[i++];
    for (unsigned tra = 0; tra < 3; tra++) 
      if (refinePar[m++]) models_perturbTrans[s][tra] = newx[i++];
    if (refinePar[m++]) models_known[s].setBfac(newx[i++]);
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) models_vrms[ensemble_modlid[e]] = newx[i++]; //VRMS
    if (refinePar[m++]) models_cell[ensemble_modlid[e]] = newx[i++]; //CELL
  }
  //LAST_ONLY implies refine VRMS of last only if FIX_VRMS is false (VRMS refined)
  FAST_LAST = (!FIX_ROT || !FIX_TRA || !FIX_BFAC || !FIX_VRMS) //refine something for models_known
               && LAST_ONLY && FIX_CELL;

  if (FAST_LAST)
  {
    //only set the vrms for the last ensemble
    if (!FIX_VRMS) ensemble->find(LAST_MODLID)->second.setup_vrms(models_vrms[LAST_MODLID]);
    fast_last_applyShift(); //sets search arrays, init for known in setUp
  }
  else
  {
  if (!FIX_VRMS) ensemble->setup_vrms(models_vrms);
  if (!FIX_CELL) ensemble->setup_cell(models_cell);
  calcKnown(*ensemble);
  if (!FIX_BFAC || !FIX_VRMS || !FIX_CELL) calcKnownVAR(*ensemble); //totvar_known
  if (!FIX_BFAC) calcKnownB((*ensemble).ensScat(TOTAL_SCAT),"",0); //frac_known_B
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);
}

void RefineMR::logCurrent(outStream where,Output& output)
{
  size_t len(0);
  for (unsigned s = 0; s < models_known.size(); s++)
    len = std::max(len,models_known[s].getModlid().size());

  where = VERBOSE; //override for MR
  output.logTab(1,where,"Perturbation of MR solutions (Degrees and Angstroms) and B-factor shift");
  if (!models_known.size()) output.logTab(1,where,"No solutions");
  else
  {
    for (unsigned s = 0; s < models_known.size(); s++)
    {
      output.logTabPrintf(1,where,
        "#%-2i %-*s %+6.4f %+6.4f %+6.4f    %+6.4f %+6.4f %+6.4f  %+6.2f\n",
        (s+1),len,models_known[s].getModlid().c_str(),
        models_perturbRot[s][0],models_perturbRot[s][1],models_perturbRot[s][2],
        models_perturbTrans[s][0],models_perturbTrans[s][1],models_perturbTrans[s][2],
        models_known[s].getBfac()-models_initBfac[s]);
    }
    if (!FIX_VRMS)
    for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
      output.logTab(1,where,"Ensemble " + iter->first + " VRMS delta factor: " + dtos(iter->second.getDeltaVRMS(),6,4,true));
    if (!FIX_CELL)
    for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
      output.logTab(1,where,"Ensemble " + iter->first + " Cell scale factor: " + dtos(iter->second.getScaleCELL(),6,4));
  }
  output.logBlank(where);
}

void RefineMR::logInitial(outStream where,Output& output)
{
  output.logTab(1,where,"Initial Parameters:");
  if (models_known.size())
  {
    mr_set tmp;
    tmp.KNOWN = models_known;
    for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
    {
      tmp.DRMS[iter->first] = iter->second.getDeltaVRMS();
      tmp.VRMS[iter->first] = iter->second.getInputVRMS();
      tmp.CELL[iter->first] = iter->second.getScaleCELL();
    }
    output.logTab(0,where,tmp.logfile(1,true,true));
  }
  output.logBlank(where);
}

void RefineMR::logFinal(outStream where,Output& output)
{
  output.logTab(1,where,"Final Parameters:");
  if (models_known.size())
  {
    mr_set tmp;
    tmp.KNOWN = models_known;
    for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
    {
      tmp.DRMS[iter->first] = iter->second.getDeltaVRMS();
      tmp.VRMS[iter->first] = iter->second.getInputVRMS();
      tmp.CELL[iter->first] = iter->second.getScaleCELL();
    }
    output.logTab(0,where,tmp.logfile(1,true,true));
  }
  output.logBlank(where);
}

std::string RefineMR::whatAmI(int& parameter)
{
  int i(0),m(0);
  for (unsigned s = 0; s < models_known.size(); s++)
  {
    std::string modlid = models_known[s].getModlid();
    std::string card = "ensemble " + modlid + " model #" + itos(s+1) + " ";
    if (refinePar[m++]) if (i++ == parameter) return card + " RotX";
    if (refinePar[m++]) if (i++ == parameter) return card + " RotY";
    if (refinePar[m++]) if (i++ == parameter) return card + " RotZ";
    if (refinePar[m++]) if (i++ == parameter) return card + " TraX";
    if (refinePar[m++]) if (i++ == parameter) return card + " TraY";
    if (refinePar[m++]) if (i++ == parameter) return card + " TraZ";
    if (refinePar[m++]) if (i++ == parameter) return card + " Bfac";
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    if (refinePar[m++]) if (i++ == parameter) return "ensemble " + ensemble_modlid[e] + " VRMS";
    if (refinePar[m++]) if (i++ == parameter) return "ensemble " + ensemble_modlid[e] + " CELL";
  }
  PHASER_ASSERT(i == npars_ref);
  return "Undefined parameter";
}

bool1D RefineMR::getRefineMask(protocolPtr protocol)
{
  FIX_ROT = protocol->getFIX(macm_rot);
  FIX_TRA = protocol->getFIX(macm_tra);
  FIX_BFAC = protocol->getFIX(macm_bfac);
  FIX_VRMS = protocol->getFIX(macm_vrms);
  FIX_CELL = protocol->getFIX(macm_cell);
  LAST_ONLY = protocol->getFIX(macm_last);
  LAST_MODLID = (models_known.size()) ? models_known[models_known.size()-1].getModlid() : "";
  //parameter refinement for rotref: hardwired, no user control
  bool REFINE_ON(true),REFINE_OFF(false);
  bool1D refineMask(0);
  cctbx::sgtbx::search_symmetry_flags flags(true,0,true);
  cctbx::sgtbx::space_group_type SgInfo(DataB::getCctbxSG());
  cctbx::sgtbx::structure_seminvariants semis(DataB::getCctbxSG());
  cctbx::sgtbx::search_symmetry ssym(flags,SgInfo,semis);
  af::tiny< bool, 3 > isshift(false,false,false);
  if (ssym.continuous_shifts_are_principal())
    isshift = ssym.continuous_shift_flags();

  int last = PTNCS.use_and_present() ? PTNCS.NMOL : 1;
  int known_minus_nmol(models_known.size()-last);
  LAST_KNOWN = std::max(0,known_minus_nmol);
  for (int s = 0; s < LAST_KNOWN; s++)
  {
    std::string modlid = models_known[s].getModlid();
    for (unsigned rot = 0; rot < 3; rot++) 
      (LAST_ONLY || FIX_ROT || models_known[s].getFixR() || (*ensemble)[modlid].is_one_atom()) ?
        refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
    for (unsigned tra = 0; tra < 3; tra++) 
      (LAST_ONLY || FIX_TRA || models_known[s].getFixT() || 
          (models_known.size()==1 && isshift[tra])) ?  // fix continuous shifts for only molecule
        refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
    (LAST_ONLY || FIX_BFAC || models_known[s].getFixB() ) ?
      refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
  }
  for (int s = LAST_KNOWN; s < models_known.size(); s++)
  {
    std::string modlid = models_known[s].getModlid();
    for (unsigned rot = 0; rot < 3; rot++) 
      (FIX_ROT || models_known[s].getFixR() || (*ensemble)[modlid].is_one_atom()) ?
        refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
    for (unsigned tra = 0; tra < 3; tra++) 
      (FIX_TRA || models_known[s].getFixT() ||
          (models_known.size()==1 && isshift[tra])) ?  // fix continuous shifts for only molecule
        refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
    (FIX_BFAC || models_known[s].getFixB()) ?
      refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
  }
  for (int e = 0; e < ensemble_modlid.size(); e++)
  {
    (FIX_VRMS || ensemble->find(ensemble_modlid[e])->second.is_one_atom() ||
      (!FIX_VRMS && LAST_ONLY && ensemble_modlid[e] != LAST_MODLID)) ? //VRMS
      refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
    (FIX_CELL || ensemble->find(ensemble_modlid[e])->second.from_pdb()) ? //CELL
      refineMask.push_back(REFINE_OFF) : refineMask.push_back(REFINE_ON);
  }
  //this is where npar_all and npar_ref are DEFINED
  return refineMask;
}

TNT::Vector<floatType> RefineMR::getRefinePars()
{
  TNT::Vector<floatType> pars(npars_ref);
  int i(0),m(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) pars[i++] = models_perturbRot[s][rot];
    for (int tra = 0; tra < 3; tra++)
      if (refinePar[m++]) pars[i++] = models_perturbTrans[s][tra];
    if (refinePar[m++]) pars[i++] = models_known[s].getBfac();
  }
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
  {
    if (refinePar[m++]) pars[i++] = iter->second.getDeltaVRMS(); //VRMS
    if (refinePar[m++]) pars[i++] = iter->second.getScaleCELL(); //CELL
  }
  PHASER_ASSERT(m == npars_all);
  PHASER_ASSERT(i == npars_ref);
  return pars;
}

std::vector<reparams> RefineMR::getRepar()
{
  std::vector<reparams> repar(npars_ref);
  int m(0),i(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) repar[i++].off();
    for (int tra = 0; tra < 3; tra++)
      if (refinePar[m++]) repar[i++].off();
    if (refinePar[m++]) repar[i++].off(); //B 
  }
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
  {
    if (refinePar[m++]) repar[i++].off(); //VRMS
    if (refinePar[m++]) repar[i++].off(); //CELL
  }
  PHASER_ASSERT(i == npars_ref);
  PHASER_ASSERT(m == npars_all);
  return repar;
}

std::vector<bounds>  RefineMR::getUpperBounds()
{
  std::vector<bounds> Upper(npars_ref);
  int m(0),i(0);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) Upper[i++].off();
    for (int tra = 0; tra < 3; tra++)
      if (refinePar[m++]) Upper[i++].off();
    // do not limit B values here, use correlated limit in getMaxDistSpecial
    if (refinePar[m++]) Upper[i++].on(upperB()); //B 
  }
  //if too large, molecule will get washed out and won't come back from limits
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
  {
    if (refinePar[m++]) Upper[i++].on(iter->second.upperDRMS()); //VRMS
    if (refinePar[m++]) Upper[i++].on(DEF_CELL_SCALE_MAX); //CELL
  }
  PHASER_ASSERT(i == npars_ref);
  PHASER_ASSERT(m == npars_all);
  return Upper;
}

std::vector<bounds>  RefineMR::getLowerBounds()
{
  std::vector<bounds> Lower(npars_ref);
  int m(0),i(0);
  //limit so that max factor from each term cannot take frac over 1
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) Lower[i++].off();
    for (int tra = 0; tra < 3; tra++)
      if (refinePar[m++]) Lower[i++].off();
   // do not limit B values here, use correlated limit in getMaxDistSpecial
    if (refinePar[m++]) Lower[i++].on(lowerB()); //B  
  }
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
  {
    if (refinePar[m++]) Lower[i++].on(iter->second.lowerDRMS()); //VRMS
    if (refinePar[m++]) Lower[i++].on(DEF_CELL_SCALE_MIN); //CELL
  }
  PHASER_ASSERT(i == npars_ref);
  PHASER_ASSERT(m == npars_all);
  return Lower;
}

void RefineMR::cleanUp(outStream where,Output& output) 
{
  dmat33 ROT,PRtr,Rperturb;
  dvect3 TRA,CMori,CMrot,CMperturb;
  dvect3 zero3(0,0,0);

  for (unsigned s = 0; s < models_known.size(); s++)
  {
    std::string modlid = models_known[s].getModlid();
    Ensemble* ens_modlid = &ensemble->find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    PRtr = ens_modlid->principalRtr();
    CMori = -PRtr*ens_modlid->principalT();
    ROT = models_known[s].getR();
    CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform orientation
    Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[s][dir])*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    models_known[s].setR(ROT);
    models_perturbRot[s] = zero3;

    CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented ensemble)
    //ROT = ROT*PRtr;
    TRA = models_known[s].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //apply xyz translation perturbation
    iOrthT += models_perturbTrans[s];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;

    //keep coordinates in same unit cell
    dvect3 deltaTRA = UnitCell::doOrth2Frac(iOrthT) - TRA;
    for (unsigned tra = 0; tra < 3; tra++)
    {
      while (deltaTRA[tra] <= -1) deltaTRA[tra]++;
      while (deltaTRA[tra] >  +1) deltaTRA[tra]--;
    }
    TRA += deltaTRA;

    models_known[s].setFracT(TRA);
    models_perturbTrans[s] = zero3;
  }
  output.logTab(1,where,"Reset rotation and translation perturbations to zero, update solutions");
  output.logBlank(where);
}

floatType RefineMR::getMaxDistSpecial(TNT::Vector<floatType>& x, TNT::Vector<floatType>& g, bool1D& bounded, TNT::Vector<floatType>& dist,floatType& start_distance)
// Return maximum multiple of gradient that can be shifted
// before hitting bounds for "special" parameters, i.e. aniso B-factors.
{
  if (FIX_BFAC && FIX_VRMS && FIX_CELL) return 0;
  TNT::Vector<floatType>  localg = g;
  //limit so that max factor from each term cannot take frac over 1
  int1D selr = get_refl_negative_variance();
  int ir(0);

  int m(0),i(0);
  float1D deltai;
  int1D   index;
  map_str_float mapEnsScat = ensemble->ensScat(TOTAL_SCAT);
  for (unsigned s = 0; s < models_known.size(); s++) 
  {
    std::string modlid = models_known[s].getModlid();
    PHASER_ASSERT(mapEnsScat.find(modlid) != mapEnsScat.end());
    for (int rot = 0; rot < 3; rot++)
      if (refinePar[m++]) i++;
    for (int tra = 0; tra < 3; tra++)
      if (refinePar[m++]) i++;
    if (refinePar[m++])
    {
      //check the ensemble
      if ((g[i] <= 0 && x[i] >= upperB()) || // x - d*g i.e. has run up agains upper limit
          (g[i] >= 0 && x[i] <= lowerB()))   // x - d*g i.e. has run up against lower limit
      {
        deltai.push_back(0);
        localg[i] = 0;
      }
      else deltai.push_back(g[i]);
      index.push_back(i);
      i++; //B  
    }
  }
  for (ensIter iter = ensemble->begin(); iter != ensemble->end(); iter++)
  {
    if (refinePar[m++])
    {
      if ((g[i] <= 0 && x[i] >= iter->second.upperDRMS()) || // x - d*g i.e. has run up agains upper limit
          (g[i] >= 0 && x[i] <= iter->second.lowerDRMS()))   // x - d*g i.e. has run up against lower limit
      {
        deltai.push_back(0);
        localg[i] = 0;
      }
      else deltai.push_back(g[i]);
      index.push_back(i);
      i++; //VRMS 
    }
    if (refinePar[m++])
    {
      if ((g[i] <= 0 && x[i] >= DEF_CELL_SCALE_MAX) || // x - d*g i.e. has run up agains upper limit
          (g[i] >= 0 && x[i] <= DEF_CELL_SCALE_MIN))   // x - d*g i.e. has run up against lower limit
      {
        deltai.push_back(0);
        localg[i] = 0;
      }
      else deltai.push_back(g[i]);
      index.push_back(i);
      i++; //CELL 
    }
  }

  //PHASER_ASSERT(max_step >= 0);
  // Could try using big but not limit
  floatType max_step(std::numeric_limits<floatType>::max());
  for (unsigned s = 0; s < deltai.size(); s++) 
  if (deltai[s] != 0) //i.e. 0 == flag for running up against bounds
  { //don't include the ones that are up against limits in test for max_step
    int i = index[s];
    max_step = std::min(max_step,dist[i]);
    //limits for dist have already been determined
    //this limit will be less than this
  }
  PHASER_ASSERT(max_step >= 0);

  //B-factor special restraint
  floatType distance(0),step(max_step);
  //limits now between distance and distance+step
  //doesn't matter if this goes outside minB and maxB, bound limit applied separately
  TNT::Vector<floatType> newx = x;
  int nloop(0);

  //store initial values
  map_str_float   init_models_vrms = models_vrms;
  map_str_float   init_models_cell = models_cell;
  float1D         init_Bfac(models_known.size(),0);
  for (unsigned s = 0; s < models_known.size(); s++) init_Bfac[s] = models_known[s].getBfac();

  floatType kssqr = -2*fn::pow2(1/fullHiRes())/4;
  floatType jssqr = -2*fn::pow2(1/LoRes())/4;
  if (step) //otherwise distance=0 and use 0 as the limit (max_step=0 above)
  for (;;)
  {
    distance += step;
    for (int i = 0; i < x.size(); i++)
    if (distance != std::numeric_limits<floatType>::max())
      newx[i] = x[i] - distance*localg[i];
    else
    {
      if (std::fabs(localg[i]) < 1 ||
          (std::fabs(localg[i]) >= 1 && distance < std::numeric_limits<floatType>::max()/localg[i]))
        newx[i] = x[i] - distance*localg[i]; //don't shift ones against the bounds
      else if (localg[i] < 0)
        newx[i] = std::numeric_limits<floatType>::max(); //don't shift ones against the bounds
      else if (localg[i] > 0)
        newx[i] = -std::numeric_limits<floatType>::max(); //don't shift ones against the bounds
      else if (localg[i] == 0)
        newx[i] = x[i]; //don't shift ones against the bounds
    }
    { //memory, applyShift
    int i(0),m(0);
    for (unsigned s = 0; s < models_known.size(); s++) 
    {
      for (int rot = 0; rot < 3; rot++)
        if (refinePar[m++]) i++;
      for (unsigned tra = 0; tra < 3; tra++) 
        if (refinePar[m++]) i++;
      if (refinePar[m++])
      {
        floatType this_models_bfac = newx[i++]; //BFAC
                  this_models_bfac = std::max(this_models_bfac,lowerB());
                  this_models_bfac = std::min(this_models_bfac,upperB());
        models_known[s].setBfac(this_models_bfac);
      }
    }
    for (int e = 0; e < ensemble_modlid.size(); e++)
    {
      if (refinePar[m++]) models_vrms[ensemble_modlid[e]] = newx[i++]; //VRMS
      if (refinePar[m++])
      {
        floatType this_models_cell = newx[i++]; //CELL
                  this_models_cell = std::max(this_models_cell,DEF_CELL_SCALE_MIN);
                  this_models_cell = std::min(this_models_cell,DEF_CELL_SCALE_MAX);
        models_cell[ensemble_modlid[e]] = this_models_cell;
      }
    }
    if (!FIX_VRMS)
      for (int e = 0; e < ensemble_modlid.size(); e++)
        if (!FAST_LAST || (FAST_LAST && ensemble_modlid[e] == LAST_MODLID))
          ensemble->find(ensemble_modlid[e])->second.setup_vrms(models_vrms[ensemble_modlid[e]]);
    if (!FIX_CELL) ensemble->setup_cell(models_cell);
    } //end applyShift
    if ( (FAST_LAST && fast_last_negative_variance(selr,ir)) ||
        (!FAST_LAST && negative_variance(*ensemble,selr,ir)))
      distance -= step;
    else break;
    if (nloop++ > 1000)
    { 
      distance = 0;
      break;
    }
    step /= 2.; 
    if (step < 1.e-8*distance)
    {
      break; // Step still big enough to make a difference
    }
  }
  //restore intial values
  for (unsigned k = 0; k < models_known.size(); k++) models_known[k].setBfac(init_Bfac[k]);
  models_vrms = init_models_vrms;
  models_cell = init_models_cell;
  if (!FIX_VRMS)
    for (int e = 0; e < ensemble_modlid.size(); e++)
      if (!FAST_LAST || (FAST_LAST &&  ensemble_modlid[e] == LAST_MODLID))
        ensemble->find(ensemble_modlid[e])->second.setup_vrms(models_vrms[ensemble_modlid[e]]);
  if (!FIX_CELL) ensemble->setup_cell(models_cell);

  for (unsigned s = 0; s < deltai.size(); s++) 
  if (deltai[s] != 0) // 0 == flag up against limits
  { 
    int i = index[s];
    distance = std::min(distance,dist[i]); //in case limit is lower
  }
  for (unsigned s = 0; s < deltai.size(); s++) 
  if (deltai[s] != 0) // 0 == flag
  {
    int i = index[s];
    dist[i] = distance;
    bounded[i] = true;
  }
  return std::max(0.0,distance);
}

void RefineMR::setUp(outStream where,Output& output)
{
  //for numerical stability
  TNT::Vector<floatType> old_x =  getRefinePars();
  applyShift(old_x); //calling search setup fast_last_applyShift internally
  if (FAST_LAST)
  { //fill up the known arrays with the fixed components and the search arrays
    if (!FIX_VRMS) ensemble->setup_vrms(models_vrms);
    fast_last_known_setUp(); 
  }
}

void RefineMR::fast_last_known_setUp()
{
  ///always required
  for (unsigned r = 0; r < NREFL; r++)
  if (selected[r])
  {
    sum_Esqr_search[r] = 0; //because there is no sum_Esqr_known
    max_Esqr_search[r] = 0;
    EM_known[r] = cmplxType(0,0);
  }

  float1D scatFactor(models_known.size(),0),sqrt_scatFactor(models_known.size(),0);
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    scatFactor[k] = (*ensemble)[modlid].ensScat(TOTAL_SCAT);
    sqrt_scatFactor[k] = std::sqrt(scatFactor[k]/NSYMP);
  }

  std::vector<int> klooplist;
  for (int kk = 0; kk < LAST_KNOWN; kk++)
    if (FIX_VRMS || (!FIX_VRMS && models_known[kk].getModlid() != LAST_MODLID))
      klooplist.push_back(kk);

  for (unsigned kk = 0; kk < klooplist.size(); kk++)
  {
    int k = klooplist[kk];
    std::string modlid = models_known[k].getModlid();
    Ensemble* ens_modlid = &ensemble->find(modlid)->second;
    //recover original centre of mass, rotated centre of mass
    dmat33 PRtr = ens_modlid->principalRtr();
    dvect3 CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = models_known[k].getR();
    dvect3 CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform 
    //orientation (for refinement)
    dmat33 Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[k][dir])*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    dvect3 CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented (*ensemble))
    ROT = ROT*PRtr;
 
    dvect3 TRA = models_known[k].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[k];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of (*ensemble))
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA = UnitCell::doOrth2Frac(eOrthT);
     
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr = Q1.transpose();

    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      floatType repsn = epsn(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      Miller1D rhkl = rotMiller(r);
      floatType terms  = sqrt_scatFactor[k];
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

  //calcKnownVAR
  if (!FIX_BFAC || !FIX_VRMS || !FIX_CELL)
  {
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      totvar_known[r] = 0;
      if (!FIX_BFAC) frac_known_B[r] = 0;
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      //save this interpolation, so not doing multiple times for multiple copies
      map_str_float Vmap;
    //stored  if (PTNCS.use_and_present()) gfun.calcReflTerms(r);
      for (unsigned kk = 0; kk < klooplist.size(); kk++)
      {
        int k = klooplist[kk];
        std::string modlid = models_known[k].getModlid();
        Ensemble* ens_modlid = &ensemble->find(modlid)->second;
        if (Vmap.find(modlid) == Vmap.end())
        {
          ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
          floatType thisV = ens_modlid->InterpV();
                    thisV *= scatFactor[k];
          Vmap[modlid] = thisV;
        }
        floatType expterm = exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
        floatType thisV = Vmap[modlid];
                  thisV *= expterm;
                  thisV /= models_known[k].getMult();
        if (PTNCS.use_and_present())
                  thisV *= G_Vterm[r]; //not correct, B-factor differences
        totvar_known[r] += thisV;
        if (!FIX_BFAC) //calcKnownB
        {
          floatType frac = scatFactor[k];
                    frac /= models_known[k].getMult();
                    frac *= -1 + expterm;
          frac_known_B[r] += frac;
        }
      }
    }
  }
}

void RefineMR::fast_last_applyShift()
{
  PHASER_ASSERT(FAST_LAST); //a subset changes
  for (unsigned r = 0; r < NREFL; r++) //initialize only required
    if (selected[r])
      EM_search[r] = cmplxType(0,0);

  float1D scatFactor(models_known.size(),0),sqrt_scatFactor(models_known.size(),0);
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    scatFactor[k] = (*ensemble)[modlid].ensScat(TOTAL_SCAT);
    sqrt_scatFactor[k] = std::sqrt(scatFactor[k]/NSYMP);
  }

  std::vector<int> klooplist;
  for (int kk = 0; kk < LAST_KNOWN; kk++)
    if (!FIX_VRMS && models_known[kk].getModlid() == LAST_MODLID)
       klooplist.push_back(kk);
  for (int kk = LAST_KNOWN; kk < models_known.size(); kk++)
    klooplist.push_back(kk);

  Ensemble* ens_modlid = &ensemble->find(LAST_MODLID)->second;
  //calcKnown
  for (unsigned kk = 0; kk < klooplist.size(); kk++)
  {
    int k = klooplist[kk];
    std::string modlid = models_known[k].getModlid();
    PHASER_ASSERT(modlid == LAST_MODLID);
    //recover original centre of mass, rotated centre of mass
    dmat33 PRtr = ens_modlid->principalRtr();
    dvect3 CMori = -PRtr*ens_modlid->principalT();
    dmat33 ROT = models_known[k].getR();
    dvect3 CMrot = ROT*CMori;
    //apply xyz rotation perturbation to molecular transform 
    //orientation (for refinement)
    dmat33 Rperturb = ens_modlid->principalR();
    for (int dir = 0; dir < 3; dir++)
      Rperturb = xyzRotMatDeg(dir,models_perturbRot[k][dir])*Rperturb;
    Rperturb = PRtr*Rperturb;
    ROT = ROT*Rperturb;
    dvect3 CMperturb = ROT*CMori;
    //correct for PR (refer rotation to reoriented (*ensemble))
    ROT = ROT*PRtr;

    dvect3 TRA = models_known[k].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[k];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of (*ensemble))
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA = UnitCell::doOrth2Frac(eOrthT);
     
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
    dmat33 Q1tr = Q1.transpose();

    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      floatType repsn = epsn(r);
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
      Miller1D rhkl = rotMiller(r);
      floatType terms  = sqrt_scatFactor[k];
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
          EM_search[r] += thisE;
        }
      }
    }
  }
  //calcKnownVAR
  if (!FIX_BFAC || !FIX_VRMS || !FIX_CELL)
  {
    for (unsigned r = 0; r < NREFL; r++)
    if (selected[r])
    {
      totvar_search[r] = 0;
      if (!FIX_BFAC) frac_search_B[r] = 0;
      floatType rssqr = ssqr(r);
      floatType sqrtrssqr = std::sqrt(rssqr);
      //save this interpolation, so not doing multiple times for multiple copies
      map_str_float Vmap;
    //stored  if (PTNCS.use_and_present()) gfun.calcReflTerms(r);
      for (unsigned kk = 0; kk < klooplist.size(); kk++)
      {
        int k = klooplist[kk];
        std::string modlid = models_known[k].getModlid();
        PHASER_ASSERT(modlid == LAST_MODLID); //ens_modlid is the right one
        if (Vmap.find(modlid) == Vmap.end())
        {
          ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
          floatType thisV = ens_modlid->InterpV();
                    thisV *= scatFactor[k];
          Vmap[modlid] = thisV;
        }
        floatType expterm = exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
        floatType thisV = Vmap[modlid];
                  thisV *= expterm;
                  thisV /= models_known[k].getMult();
        if (PTNCS.use_and_present())
                  thisV *= G_Vterm[r]; //not correct, B-factor differences
        totvar_search[r] += thisV;
        if (!FIX_BFAC) //calcKnownB
        {
          floatType frac = scatFactor[k];
                    frac /= models_known[k].getMult();
                    frac *= -1 + expterm;
          frac_search_B[r] += frac;
        }
      }
    }
  }
}

bool RefineMR::fast_last_negative_variance(int1D& selr,int& ir)
{
  PHASER_ASSERT(FAST_LAST); //a subset changes
  if (!models_known.size()) return false;
  std::vector<dvect3> TRA(models_known.size(),dvect3(0,0,0));
  std::vector<dmat33> Q1tr(models_known.size(),dmat33(1,0,0,0,1,0,0,0,1));
  float1D scatFactor(models_known.size());
  for (unsigned k = 0; k < models_known.size(); k++)
  {
    std::string modlid = models_known[k].getModlid();
    if (modlid != LAST_MODLID) continue;
    Ensemble* ens_modlid = &ensemble->find(modlid)->second;
    scatFactor[k] = ens_modlid->ensScat(TOTAL_SCAT);
    dmat33 PRtr,Rperturb;
    dvect3 CMori,CMrot,CMperturb;
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
    //correct for PR (refer rotation to reoriented (*ensemble))
    ROT = ROT*PRtr;

    TRA[k] = models_known[k].getT();
    dvect3 iOrthT = UnitCell::doFrac2Orth(TRA[k]);
    //apply xyz translation perturbation (for refinement)
    iOrthT += models_perturbTrans[k];
    //correct for rotation of centre of mass
    iOrthT += CMrot - CMperturb;
    //correct for PT (translation to origin of (*ensemble))
    dvect3 eOrthT = iOrthT - ROT*ens_modlid->principalT();
    TRA[k] = UnitCell::doOrth2Frac(eOrthT);
     
    dmat33 R = ROT*ens_modlid->Frac2Orth();
    dmat33 Q1 = UnitCell::doOrth2Frac(R);
           Q1tr[k] = Q1.transpose();
  }

  Ensemble* ens_modlid = LAST_MODLID.size() ? &ensemble->find(LAST_MODLID)->second : NULL;
  for (; ir < selr.size(); ir++)
  {
    unsigned r = selr[ir];
    floatType repsn = epsn(r);
    floatType rssqr = ssqr(r);
    floatType sqrtrssqr = std::sqrt(rssqr);
    //save this interpolation, so not doing multiple times for multiple copies
    ens_modlid->set_sqrt_ssqr(sqrtrssqr,rssqr);
    floatType thisV = ens_modlid->InterpV();
    floatType Vknown_modlid  = thisV;
    floatType this_totvar = totvar_search[r];
    cmplxType this_EM = EM_search[r];
    floatType this_frac_B = frac_search_B[r];
    totvar_search[r] = 0;
    EM_search[r] =  0;
    frac_search_B[r] = 0;
    std::vector<int> klooplist;
    if (FIX_VRMS)
      for (int kk = LAST_KNOWN; kk < models_known.size(); kk++)
        klooplist.push_back(kk);
    else
    {
      for (int kk = 0; kk < models_known.size(); kk++)
        if (models_known[kk].getModlid() == LAST_MODLID) klooplist.push_back(kk);
    }
    for (unsigned kk = 0; kk < klooplist.size(); kk++)
    {
      int k = klooplist[kk];
      std::string modlid = models_known[k].getModlid();
      PHASER_ASSERT(modlid == LAST_MODLID);
      floatType frac = scatFactor[k];
                frac /= models_known[k].getMult();
                frac *= -1 + std::exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
      frac_search_B[r] += frac;
      floatType thisV = Vknown_modlid;
                thisV *= scatFactor[k];
                thisV *= exp(-2.0*models_known[k].getBfac()*rssqr/4.0);
                thisV /= models_known[k].getMult();
      if (PTNCS.use_and_present()) thisV *= G_Vterm[r]; //not correct, B-factor differences
      totvar_search[r] += thisV;

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
          EM_search[r] += thisE;
        }
      }
    }
    Rice_refl(r,true);
    totvar_search[r] = this_totvar;
    EM_search[r] =  this_EM;
    frac_search_B[r] = this_frac_B;
    if (RiceV <=0)
    {
      return true;
    }
  }
  return false;
}

}//phaser
