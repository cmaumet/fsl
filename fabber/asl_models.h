/*   asl_model.h Kinetic curve models for ASL

      Michael Chappell - IBME & FMRIB Image Analysis Group

      Copyright (C) 2010-2011 University of Oxford */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#if !defined(asl_models_h)
#define asl_models_h

#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"

using namespace MISCMATHS;

#include "utils/tracer_plus.h"

using namespace Utilities;

namespace OXASL {

  // generic AIF model class
  class AIFModel {
  public:
    //evaluate the model
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const = 0;
    //report the number of dispersion parameters 
    virtual int NumDisp() const = 0;
    // return default priors for the parameters
    virtual ColumnVector Priors() const { return priors; }
    virtual string Name() const = 0;

  protected:
    ColumnVector priors; //list of prior means and precisions - all means first then precisions

  };

  //Specific AIF models
  class AIFModel_nodisp : public AIFModel {
    //AIFModel_nodisp() {}
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl, const ColumnVector dispparam) const;
  virtual int NumDisp() const {return 0;}
  virtual string Name() const { return "None"; }
  };

  class AIFModel_gammadisp : public AIFModel {
  public:
    AIFModel_gammadisp() { priors.ReSize(4); priors << 2 << -0.3 << 10 << 10; }
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 2;}
    virtual string Name() const { return "Gamma dispersion kernel"; }
  };

  //  double kcblood_gvf(const double ti, const double deltblood,const double taub,const double T_1b, const double s, const double p, bool casl);
  //  double kcblood_gaussdisp(const double ti, const double deltblood, const double taub, const double T_1b, const double sig1, const double sig2, bool casl);
  //  double kcblood_spatialgaussdisp(const double ti, const double deltblood, const double taub, const double T_1b, const double k, bool casl);
  //  double kcblood_gallichan(const double ti, const double deltblood, const double taub, const double T_1b, const double xdivVm, bool casl);

  // -------------
  //generic residue function model class
  class ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const = 0;
    //report the number of residue function parameters
    virtual int NumResid() const = 0;
    // return the default priors for the parameters
    virtual ColumnVector Priors() const {return residpriors;}
    virtual string Name() const = 0;

  protected:
    ColumnVector residpriors;
  };

  //specific residue function models
  class ResidModel_wellmix : public ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "Well mixed"; }
  };

  class ResidModel_simple : public ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;
    
    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "Simple"; }
  };

  class ResidModel_imperm : public ResidModel {
  public:
    ResidModel_imperm() { residpriors.ReSize(2); residpriors << 0.5 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "Impermeable"; }
  };

  class ResidModel_twocpt : public ResidModel {
  public:
    ResidModel_twocpt() { residpriors.ReSize(2); residpriors << 0.8 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "Two comparment (no backflow, no venous output)"; }
  };

  class ResidModel_spa : public ResidModel {
  public:
    ResidModel_spa() { residpriors.ReSize(6); residpriors << 0.02 << 0.03 << 2 << 1e-3 << 1e12 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 3;}
    virtual string Name() const { return "Single Pass Approximation (2 compartment, no backflow)"; }
  };

  // ------------
  //generic tissue model class
  class TissueModel {
  public:
    //evalute the model
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const = 0;
    // report the number of dipersion parameters
    virtual int NumDisp() const = 0;
    // report the number of residue function parameters (beyond the normal ones)
    virtual int NumResid() const = 0;
    // return default priors for the parameters
    virtual ColumnVector DispPriors() const { return disppriors; }
    virtual ColumnVector ResidPriors() const {return residpriors; }
    virtual string Name() const = 0;

  protected:
    ColumnVector disppriors; //list of prior means and precisions - all means first then precisions
    ColumnVector residpriors;
  };

  //specific tissue models

  class TissueModel_nodisp_simple : public TissueModel {
    virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return 0;}
    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "No dispersion | Simple"; }
  };

  class TissueModel_nodisp_wellmix : public TissueModel {
  public:
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 0;}
 virtual string Name() const { return "No dispersion | Well mixed"; }
  };

  class TissueModel_nodisp_imperm : public TissueModel {
  public:
    TissueModel_nodisp_imperm() { residpriors.ReSize(2); residpriors << 0.5 << 10; }
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 1;}
 virtual string Name() const { return "No dispersion | Impermeable"; }
  };

  class TissueModel_nodisp_2cpt : public TissueModel {
  public:
    TissueModel_nodisp_2cpt() { residpriors.ReSize(2); residpriors << 0.8 << 10; }
    virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return 0;}
    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "No dispersion | Two compartmentr (no backflow, no venous output)"; }
  };

  class TissueModel_nodisp_spa : public TissueModel {
  public:
    TissueModel_nodisp_spa() { residpriors.ReSize(6); residpriors << 0.02 << 0.03 << 2 << 1e-3 << 1e12 << 10; }
virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 3;}
 virtual string Name() const { return "No dispersion | Single Pass Approximation (2 compartment no backflow)"; }

  private:
 double Q(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const;
 double R(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const;
  };

  class TissueModel_gammadisp_wellmix : public TissueModel {
  public:
    TissueModel_gammadisp_wellmix() { disppriors.ReSize(4); disppriors << 2 << -0.3 << 10 << 10; } 
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 2;}
 virtual int NumResid() const {return 0;}
  virtual string Name() const { return "Gamma kernel dispersion | Well mixed"; }
  };


  //double kctissue_gvf(const double ti, const double delttiss,const double tau, const double T_1b, const double T_1app, const double s, const double p);
  //double kctissue_gaussdisp(const double ti, const double delttiss, const double tau, const double T_1b, const double T_1app, const double sig1, const double sig2);

  //  a general tissue model that does numerical convolution
  class TissueModel_aif_residue : public TissueModel {
  public:
    TissueModel_aif_residue(AIFModel* paifmodel, ResidModel* presidmodel)
      { aifmodel=paifmodel; residmodel = presidmodel;
	disppriors << aifmodel->Priors();
	residpriors << residmodel->Priors();
      }
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1app, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return aifmodel->NumDisp();}
    virtual int NumResid() const {return residmodel->NumResid();}
    virtual string Name() const { string name; name = "NUMERICAL CONVOLUTION - " + aifmodel->Name() + residmodel->Name(); return name; }

  protected:
    AIFModel* aifmodel;
    ResidModel* residmodel;
  };

  // useful functions
  double icgf(const double a, const double x);
  double gvf(const double t, const double s, const double p);


}

#endif
