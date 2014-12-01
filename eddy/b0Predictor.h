/*! \file b0Predictor.h
    \brief Contains declaration of class for making silly (mean) predictions about b=0 data.

    \author Jesper Andersson
    \version 1.0b, Sep., 2012.
*/
// Declarations of class to make silly
// predictions about b0 scans.
//
// b0Predictor.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef b0Predictor_h
#define b0Predictor_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "DWIPredictionMaker.h"

namespace EDDY {

/****************************************************************//**
*
* \brief Class for making silly (mean) predictions about b=0 data.
*
* Will make predictions about b=0 scans. Since we don't expect any
* variability in these scans the prediction will simply be the mean
* of all the scans in the prediction maker. It is there mainly so
* that the same code can be used for registration of 
* diffusion-weighted and b=0 scans.
*
********************************************************************/ 
  class b0Predictor : public DWIPredictionMaker
{
public:
  b0Predictor() : _pop(true), _valid(false) {};
  ~b0Predictor() {}
  bool IsPopulated() const;                  // Returns true if all data present
  bool IsValid() const { return(_valid); }   // Returns true if ready to make predictions
  void SetNoOfScans(unsigned int n);
  void AddScan(const NEWIMAGE::volume<float>& scan, // NOT thread safe
               const DiffPara&                dp);
  void SetScan(const NEWIMAGE::volume<float>& scan, // _May_ be thread safe if used "sensibly"
               const DiffPara&                dp,
               unsigned int                   indx);
  void EvaluateModel(const NEWIMAGE::volume<float>& mask) 
  { 
    if (!IsPopulated()) throw EddyException("b0Predictor::EvaluateModel: Not ready to evaluate model");
    if (!IsValid()) _mean=mean_vol(mask); 
  }
  NEWIMAGE::volume<float> Predict(unsigned int indx) const 
  { if (!IsValid()) throw EddyException("b0Predictor::Predict: Not ready to make predictions");
    return(_mean); 
  }
  NEWIMAGE::volume<float> Predict(const DiffPara& dpars) const
  { if (!IsValid()) throw EddyException("b0Predictor::Predict: Not ready to make predictions");
    return(_mean); 
  }
  NEWIMAGE::volume4D<float> PredictAll() const
  { if (!IsValid()) throw EddyException("b0Predictor::PredictAll: Not ready to make predictions");
    NEWIMAGE::volume4D<float> ovol(_mean.xsize(),_mean.ysize(),_mean.zsize(),int(_sptrs.size()));
    for (int s=0; s<ovol.tsize(); s++) ovol[s] = _mean;
    return(ovol); 
  }

private:
  std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > > _sptrs; // Pointers to the scans
  NEWIMAGE::volume<float>                                   _mean;  // Mean image
  mutable bool                                              _pop;   // Tells if all data is present
  bool                                                      _valid; // Tells if it is ready to make predictions

  NEWIMAGE::volume<float> mean_vol(const NEWIMAGE::volume<float>& mask)
  {
    NEWIMAGE::volume<float> mvol = *_sptrs[0];
    mvol = 0.0;
    for (int k=0; k<_sptrs[0]->zsize(); k++) {
      for (int j=0; j<_sptrs[0]->ysize(); j++) {
	for (int i=0; i<_sptrs[0]->xsize(); i++) {
	  if (mask(i,j,k)) {
            float &v = mvol(i,j,k);
	    for (unsigned int s=0; s<_sptrs.size(); s++) {
              v += (*_sptrs[s])(i,j,k);
	    }
            v /= float(_sptrs.size());
	  }
	}
      }
    }
    _valid = true;
    return(mvol);
  }
};

} // End namespace EDDY

#endif // End #ifndef b0Predictor_h
