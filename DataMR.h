//(c) 2000-2014 Cambridge University Technical Services Ltd
//All rights reserved
#ifndef __DataMRClass__
#define __DataMRClass__
#include <phaser/include/Phaser.h>
#include <phaser/include/mr_solution.h>
#include <phaser/include/MapCoefs.h>
#include <phaser/include/data_btf.h>
#include <phaser/src/DataB.h>
#include <phaser/src/MapEnsemble.h>
#include <phaser/src/Hexagonal.h>
#include <phaser/src/Brick.h>
#include <phaser/src/Gfunction.h>
#include <phaser/lib/maths.h>

namespace phaser {

class DataMR : public DataB
{
  public:
    //Constructor and Destructor
    DataMR(std::string,af::double6,data_refl&,Bin&,data_norm&,data_outl&,data_composition&,floatType,floatType,data_resharp&,data_tncs&);
    DataMR() { debug_flag = 0; }

    DataMR(const DataMR & init) : DataB((DataB&)init) // copy constructor necessary because of fast_cos_sin
    { setDataMR(init); }

    const DataMR& operator=(const DataMR& right)
    {
      if (&right != this)
      {
        setDataMR(right);
        (DataB&)*this = right;
      }
      return *this;
    }

    void setDataMR(const DataMR & init)
    {
       frac_known_B = init.frac_known_B;
       frac_search_B = init.frac_search_B;
       m_alogch = init.m_alogch;
       m_alogchI0 = init.m_alogchI0;

       models_known = init.models_known;
       models_perturbRot = init.models_perturbRot;
       models_perturbTrans = init.models_perturbTrans;

       EM_known = init.EM_known;
       totvar_known = init.totvar_known;
       EM_search = init.EM_search;
       sum_Esqr_search = init.sum_Esqr_search;
       max_Esqr_search = init.max_Esqr_search;
       totvar_search = init.totvar_search;
       G_Vterm = init.G_Vterm;

       debug_flag = init.debug_flag;

       gfun = init.gfun;

       initCos_Sin(); // dynamically allocated table created here
    }

  public:
    bool debug_flag;

  protected:
    //single value
    floatType RiceV;
    Gfunction gfun;

    //arrays
    Alogch   m_alogch;
    AlogchI0 m_alogchI0;

    //Data dependent on solution [a]
    std::vector<mr_ndim> models_known;
    std::vector<dvect3>  models_perturbRot;
    std::vector<dvect3>  models_perturbTrans;
    float1D              models_initBfac;

    //PtNcs related
    float1D     G_Vterm;

    //Data dependent on reflection [r]
    cmplx1D EM_known,EM_search; // calculated model sf
    float1D sum_Esqr_search,max_Esqr_search;
    float1D totvar_known,totvar_search,totvar_move;
    float1D frac_known_B,frac_search_B;

  private:
    pair_flt_flt ScaleFC();
    floatType logRelRice(floatType,floatType,floatType);
    floatType logRelWoolfson(floatType,floatType,floatType);

  public:
    void      logModls(outStream,Output&,std::string="");
    void      calcKnownB(map_str_float,std::string,floatType);
    void      initKnownMR(mr_set,MapEnsemble&);
    void      deleteLastKnown();
    void      initSearchMR();
    void      init();
    void      initGfunction(bool&);
    void      calcKnown(MapEnsemble&);
    void      calcKnownVAR(MapEnsemble&);
    void      calcKnownVAR(MapEnsemble&, map_str_pdb&);
    void      calcMoveVAR(dmat33, data_pdb&);
    void      calcMoveVAR(data_pdb&);
    void      calcSearchVAR(data_pdb&);
    void      calcSearchVAR(Ensemble&,bool,floatType);
    void      calcSearchVAR(std::string,mrgyre1D&,MapEnsemble&,floatType);
    void      calcSearchROT(dmat33,Ensemble&,floatType);
    void      calcSearchTRA1(dmat33,Ensemble&);
    void      calcSearchTRA2(dvect3,bool,floatType,cmplx1D*,cmplx1D* interpE2=NULL); //NULL pointer
    void      calcSearchTRA2(dvect3,std::string,MapEnsemble&,mrgyre1D&,cmplx1D*);

    floatType likelihoodFn() { return Rice()-LLwilson; }

    floatType Rice();
    floatType Rice_refl(int,bool);
    floatType Rfactor();
    cmplx1D   getInterpE(dmat33,Ensemble&);
    cmplx1D   getInterpE(std::string,mrgyre1D&,MapEnsemble&,floatType);
    bool      haveFpart();
    float1D   getInterpV();
    cmplx1D   getInterpE();
    void      setInterpV(float1D&);
    void      setInterpE(cmplx1D&);
    void      addInterpV(float1D&);
    void      addInterpE(cmplx1D&);

    cmplx3D    getELMNxR2(std::string,Output&,int);
    cmplx3D    getELMNxR2Serial(std::string, Output&,int);
    floatType  m_LETF1(dmat33,Ensemble&,bool,af_float&,bool,floatType);
    floatType  m_LETF1(std::string,mrgyre1D&,MapEnsemble&,bool,af_float&,floatType);
    floatType  m_LETF2(dmat33,Ensemble&,bool,af_float&,af_float&,bool,floatType);

    af_cmplx getFcalc(dmat33,dmat33,Ensemble&,af_miller,bool,floatType);
    af_cmplx getFcalc(std::string,mrgyre1D&,MapEnsemble&,af_miller,floatType);
    af_cmplx getFcalcPhsTF(dmat33,Ensemble&,floatType);
    af_cmplx getFobsPhsTF();
    af_float getEnat();
    af_cmplx getEpart();
    std::vector<mr_ndim>  getModls();
    MapCoefs              getMapCoefs();

    floatType lowerB();
    floatType upperB();
    floatType tolB();
    Hexagonal calc_search_points(data_btf&,floatType,bool,Output&);
    int1D     get_refl_negative_variance();
    bool      negative_variance(MapEnsemble&,int1D&,int&);
};

}
#endif
