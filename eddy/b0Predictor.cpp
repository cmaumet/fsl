// Definitions of class to make silly
// predictions about b0 scans.
//
// b0Predictor.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "EddyUtils.h"
#include "b0Predictor.h"

using namespace EDDY;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class b0Predictor
//
// 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bool b0Predictor::IsPopulated() const
{
  if (_pop) return(true);
  else {
    _pop = true;
    for (unsigned int i=0; i<_sptrs.size(); i++) {
      if (!_sptrs[i]) { _pop = false; break; }
    }
  }
  return(_pop);
}

void b0Predictor::SetNoOfScans(unsigned int n)
{
  if (n == _sptrs.size()) return; // No change
  else if (n > _sptrs.size()) {   // If increasing size
    _sptrs.resize(n,boost::shared_ptr<NEWIMAGE::volume<float> >()); // New elements populated by NULL
    _pop = false;
    _valid = false;
  }
  else { // If decreasing size
    _sptrs.resize(n);
    _valid = false;
    if (_pop==false) { // _pop may potentially go from false to true
      _pop = IsPopulated();
    }
  }
  return;
}

//
// Ideally one would like to check if array is ready populated
// when adding a scan. But that would potentially make it 
// non-thread reentrant so I have choosen not to.
//
void b0Predictor::AddScan(const NEWIMAGE::volume<float>& scan,
                          const DiffPara&                dp)
{
  if (_sptrs.size() && !NEWIMAGE::samesize(*_sptrs[0],scan)) throw EddyException("b0Predictor::AddScan: Wrong image dimension");
  _sptrs.push_back(boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scan)));
  _valid = false;
}

void b0Predictor::SetScan(const NEWIMAGE::volume<float>& scan,
                          const DiffPara&                dp,
                          unsigned int                   indx)
{
  if (int(indx) > (int(_sptrs.size())-1)) throw EddyException("b0Predictor::SetScan: Invalid image index");
  if (_sptrs.size() && _sptrs[0] && !NEWIMAGE::samesize(*_sptrs[0],scan)) throw EddyException("b0Predictor::SetScan: Wrong image dimension");
  _sptrs[indx] = boost::shared_ptr<NEWIMAGE::volume<float> >(new NEWIMAGE::volume<float>(scan));
}

