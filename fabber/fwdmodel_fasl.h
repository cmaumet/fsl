/*  fwdmodel_fasl.h - multi-TI functional ASL model

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 20010 University of Oxford  */

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


#ifndef __FABBER_FWDMODEL_FASL_H
#define __FABBER_FWDMODEL_FASL_H 1

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class FASLFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const
  { return D0index()+Dbasis.Ncols() + ((stdevInvEff>0)?1:0) + ((stdevT1b>0)?1:0) + ((stdevT1>0)?1:0) +1;}
      // { return NnStart(echoTime.Ncols()+1) + (stdevT1b>0?1:0) 
    //   + (stdevInvEff>0?1:0) + (stdevDt>0?1:0) - 1; }

  virtual string ModelVersion() const;

  virtual ~FASLFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  // Constructor
  FASLFwdModel(ArgsType& args);
  // Usage info
  static void ModelUsage();


protected: // Constants
  
  float kctissue_nodisp(const float ti, const float delttiss, const float tau, const float T_1b, const float T_1app) const;
  
  int Q0index() const { return 1; }
  int M0index() const { return Q0index()+Qbasis.Ncols()+1; }
  int D0index() const {return M0index()+Mbasis.Ncols()+1; }
  //  int R0index() const { return M0index()+Mbasis.Ncols()+1; }
  //int NnStart(int te) const 
  //  { assert(te > 0 && te <= echoTime.Nrows()+1); // go one past, for getting end point
  //    return R0index()+Rbasis.Ncols()+(te-1)*Nbasis.Ncols() + 1; }
  int InvEffIndex() const { assert(stdevInvEff>0); return D0index()+Dbasis.Ncols() +1 ; }
  int T1bIndex() const { assert(stdevT1b>0); 
    return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+1; }
 int T1Index() const { assert(stdevT1>0); 
    return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+((stdevT1b>0)?1:0)+1; }
 int AIndex() const { 
   return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+((stdevT1b>0)?1:0)+((stdevT1>0)?1:0)+1; }
  //int dtIndex() const { assert(stdevDt>0);
  //  return NnStart(echoTime.Ncols()+1) + (stdevT1b>0?1:0) + (stdevInvEff>0?1:0); }

  // Slow submatrixing... but it works.
  // In fact, it's plenty fast enough.. fwd model is a small part calculation time
  ReturnMatrix QnOf(const ColumnVector& params) const
    { return params.Rows(Q0index()+1,Q0index()+Qbasis.Ncols()).Evaluate(); }
 ReturnMatrix DnOf(const ColumnVector& params) const
    { return params.Rows(D0index()+1,D0index()+Dbasis.Ncols()).Evaluate(); }
  ReturnMatrix MnOf(const ColumnVector& params) const
    { return params.Rows(M0index()+1,M0index()+Mbasis.Ncols()).Evaluate(); }
//  ReturnMatrix RnOf(const ColumnVector& params) const
//    { return params.Rows(R0index()+1,R0index()+Rbasis.Ncols()).Evaluate(); }  
//  ReturnMatrix NnOf(int te, const ColumnVector& params) const
//    { return params.Rows(NnStart(te),NnStart(te+1)-1).Evaluate(); }

  // scan parameters
  ColumnVector rho;
//ColumnVector echoTime;
ColumnVector tivec;
float tau;
  
 bool dataisdiff;

  // assumed known constants
 double fixedInvEff, fixedT1b, fixedDt, fixedT1;
 double stdevT1b, stdevInvEff, stdevDt, stdevT1; // >0 means treat these as parameters!
  
  // basis for each piece
  Matrix Qbasis; // basis in columns?
Matrix Dbasis;
  Matrix Mbasis;
//  Matrix Rbasis;
//  Matrix Nbasis;
};

#endif /* __FABBER_FWDMODEL_QUIPSS2_H */

