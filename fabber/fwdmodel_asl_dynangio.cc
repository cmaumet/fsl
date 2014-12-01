/*  fwdmodel_asl_dynangio.cc -kinetic curve modelling for ASL Dynamic Angio data

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

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

#include "fwdmodel_asl_dynangio.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string DynAngioFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_dynangio.cc,v 1.2 2011/03/10 13:56:06 chappell Exp $";
}

void DynAngioFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("DynAngioFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

     // Set priors
     // abv
     prior.means(1) = 0;
     precisions(1,1) = 1e-12;

     // transit delay
     prior.means(2) = 0.1;
     precisions(2,2) = 10;
     
     // bolus duration
     prior.means(3) = seqtau;
     if (taufix) precisions(3,3) = 1e99;
     else        precisions(3,3) = 10;
     
     // T1b
     prior.means(4) = t1b;
     precisions(4,4) = 100;
     
     //dispersion parameters
     //default to 0 and known
     prior.means(disp_index()) = 0;
     prior.means(disp_index()+1) = 0;
     precisions(disp_index(),disp_index()) = 1e12;
     precisions(disp_index()+1,disp_index()+1) = 1e12;

     if (disptype=="gamma") {
     prior.means(disp_index()) = 1.1; //-1;//0.05;
     prior.means(disp_index()+1) = 0.1; //1.1; //0.7; //2;
     precisions(disp_index(),disp_index()) = 1; //10; //400;
     precisions(disp_index()+1,disp_index()+1) = 100; //1;
     }

     if (disptype=="gvf") {
       prior.means(disp_index()) = 1.1;//0.05;
       prior.means(disp_index()+1) = 0.7;
       precisions(disp_index(),disp_index()) = 1; //400;
       precisions(disp_index()+1,disp_index()+1) = 10;

     }

     if (disptype=="gauss") {
       prior.means(disp_index()) = -1.6;
       precisions(disp_index(),disp_index()) = 1;
     }
     if (disptype=="sgauss") {
       prior.means(disp_index()) = -1.4;
       precisions(disp_index(),disp_index()) = 1;
     }
     if (disptype=="gallichan") {
       prior.means(disp_index()) = 0.2;
       precisions(disp_index(),disp_index()) = 0.1;
     }

     prior.means(7) = 0;
     precisions(7,7) = 0.1;

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    posterior.means(1) = 1;
    precisions(1,1) = 1;

    // for dispersion parameters - initialise to a realtively non dispersed case
    if (disptype=="gamma" || disptype=="gvf") {
    prior.means(disp_index()) = 5;
    prior.means(disp_index()) = 1;
    }
    
    posterior.SetPrecisions(precisions);
    
}    
    
    

void DynAngioFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("DynAngioFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
     // sensible limits on transit times
  if (params(2)>timax-0.2) { paramcpy(2) = timax-0.2; }

  // parameters that are inferred - extract and give sensible names
  float abv;
  float delt;
  float tau;
  float T_1b;

  bool outresult=true;

  abv = params(1);
  delt = params(2);
  tau = params(3);
  T_1b = params(4);

  float floor;
  floor = params(7);

  if (abv < 0 || delt<0 || tau < 0 || T_1b < 0 ) {
    //delt=0;
    outresult=true;
  }

  if (outresult) {
  // Look-locker readout
  T_1b = 1/( 1/T_1b - log(cos(FA))/dti); //assume that all the blood has seen the small FAs


  // deal with different dispersion models
   ColumnVector kcblood(tis.Nrows());

    // generate the kinetic curves
    if (disptype=="none") {
      kcblood=kcblood_nodisp(tis,delt,tau,T_1b,false);
    //cout << kcblood << endl;
    }
    else if (disptype=="gamma") {
      float p;
      float s;
      s = exp( params(disp_index()) );
      //float k=exp(params(disp_index()+1));;
      //p = (k-1)/s; //possibly define that parameter as (k-1) ?
      p = params(disp_index()+1);
      if (p<0) p=0;

      kcblood=kcblood_gammadisp(tis,delt,tau,T_1b,s,p,false);
    //cout << kcblood << endl;
    }
    else if (disptype=="gvf") {
      float p;
      float s;
      //float sp = exp(params(disp_index()));
      s = exp( params(disp_index()) );
      //if (sp>10) sp=10;
      p = params(disp_index()+1); //sp/s;
      

      kcblood=kcblood_gvf(tis,delt,T_1b,s,p,false);
    //cout << kcblood << endl;
    }
    else if (disptype=="gauss") {
      float sig;
      sig = exp( params(disp_index()) );
      // we will only have diserpsion SD (leading edge)
      // assume trailing edge is related to leading edge as Hrabe did
      float sig2;      
      if (delt > 0) {
	sig2  = sig*sqrt( (delt+tau)/delt );
      }
      else sig2 = sig;

      //cout << sig << endl;
      kcblood = kcblood_gaussdisp(tis,delt,tau,T_1b,sig,sig2,false);
      //cout << kcblood.t() << endl;
    }
    else if (disptype=="sgauss") {
      float k;
      k = exp( params(disp_index()) );
      
      kcblood = kcblood_spatialgaussdisp(tis,delt,tau,T_1b,k,false);
    }
    else if (disptype=="gallichan") {
      float xV;
      xV = params(disp_index());

      kcblood = kcblood_gallichan(tis,delt,tau,T_1b,xV,false);
    }
    else {
      throw Exception("Unrecognised dispersion model ");
    }

    // Nan catching

    bool cont=true;
    int it=1;
    while (cont) {
      if (isnan(kcblood(it)) || isinf(kcblood(it)) ) { 
      
      LOG << "Warning NaN or Inf in kcblood" << endl;
      LOG << "params: " << params.t() << endl;
      LOG << "kcblood: " << kcblood.t() << endl;
      cont =false;
      kcblood=1e12;
 }
    it++;
    if (it>kcblood.Nrows()) cont=false;
  }

    // assemble the result
    result = floor + abv*kcblood;
  }
  else {
    result.ReSize(tis.Nrows());
    result = 0;
  }

  return;
}


DynAngioFwdModel::DynAngioFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      //dispersion model
      disptype=args.ReadWithDefault("disp","gamma");

      //repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.66"));

      seqtau = convertTo<double>(args.ReadWithDefault("tau","1000")); //bolus length as set by sequence (default of 1000 is effectively infinite

      // special options
      taufix = args.ReadBool("taufix");
     
      // Deal with tis
      tis.ReSize(1); //will add extra values onto end as needed
      tis(1) = atof(args.Read("ti1").c_str());
      
      while (true) //get the rest of the tis
	{
	  int N = tis.Nrows()+1;
	  string tiString = args.ReadWithDefault("ti"+stringify(N), "stop!");
	  if (tiString == "stop!") break; //we have run out of tis
	 
	  // append the new ti onto the end of the list
	  ColumnVector tmp(1);
	  tmp = convertTo<double>(tiString);
	  tis &= tmp; //vertical concatenation

	}
      timax = tis.Maximum(); //dtermine the final TI
      //determine the TI interval (assume it is even throughout)
      dti = tis(2)-tis(1);

      float fadeg = convertTo<double>(args.ReadWithDefault("fa","30"));
      FA = fadeg * M_PI/180;
	
      // add information about the parameters to the log
      LOG << "Inference using Dynamic Angio model" << endl;
      LOG << "TIs: ";
      for (int i=1; i <= tis.Nrows(); i++)
	LOG << tis(i) << " ";
      LOG << endl;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void DynAngioFwdModel::ModelUsage()
{ 
  cout << "To be added"
       << endl
    ;
}

void DynAngioFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void DynAngioFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("aBV");
  names.push_back("delt");
  names.push_back("tau");
  names.push_back("T_1b");

  names.push_back("disp1");
  names.push_back("disp2");

}



void DynAngioFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard)
{
  
}

void DynAngioFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  
  }
