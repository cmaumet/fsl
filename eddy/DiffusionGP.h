/*! \file DiffusionGP.h
    \brief Contains declaration of class for making Gaussian process based predictions about DWI data.

    \author Jesper Andersson
    \version 1.0b, Sep., 2012.
*/
// Declarations of class to make Gaussian-Process
// based predictions about diffusion data.
//
// DiffusionGP.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef DiffusionGP_h
#define DiffusionGP_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "DWIPredictionMaker.h"

namespace EDDY {

class DiffusionGP;

/****************************************************************//**
*
* \brief Class used to calculate the K-matrix and the derived 
* entities (such as the prediction vector) used for making GP
* predictions. 
*
* Class used to calculate the K-matrix and the derived 
* entities (such as the prediction vector) used for making GP
* predictions. It needs to be supplied with a "grouping" of the
* the dwi scans into shells and estimates of white matter (signal)
* variance for each shell and an estimate of noise variance.
* The class VarianceCalculator can be used to calculate those.
*
********************************************************************/ 
class GPModel
{
  friend class DiffusionGP;
public:
  /// Constructor
  GPModel(const std::vector<DiffPara>&          dpars,
	  const std::vector<double>             vars,
	  const std::vector<std::vector<int> >  grps,
	  const std::vector<double>             grpb,
	  double                                s2n);
  /// Returns the number of data points.
  unsigned int NPoints() const { return(_dpars.size()); }
  /// Returns the number of groups (of ~equal b-values)
  unsigned int NGroups() const { return(_vars.size()); }
  /// Returns the white matter variances for the different groups
  const std::vector<double>& WMVariances() const { return(_vars); }
  /// Returns the tentative noise variance
  double NoiseVariance() const { return(_s2n); }
  /// Returns the group for scan given by scanindx
  unsigned int Group(unsigned int scanindx) const { 
    if (scanindx>=_grp.size()) throw EddyException("GPModel::Group: Index out of range"); return(_grp[scanindx]); }
  /// Returns the b-value for group given by grpindx
  double GroupBValue(unsigned int grpindx) const {
    if (grpindx>=_grpb.size()) throw EddyException("GPModel::GroupBValue: Index out of range"); return(_grpb[grpindx]); }
  /// Returns the K-matrix
  const NEWMAT::Matrix& KMatrix() const { return(_K); }
  /// Returns a (smoothed) prediction for scan given by indx
  double Predict(unsigned int                indx,
		 const NEWMAT::ColumnVector& y,
		 double                      sigman2=-1.0) const;
  /// Returns a (n interpolated) prediction for observation given by dpars
  double Predict(const DiffPara&             dpars,
		 const NEWMAT::ColumnVector& y,
		 double                      sigman2=-1.0) const;
private:
  const std::vector<DiffPara>               _dpars;  /// Diffusion parameters.
  const std::vector<double>                 _vars;   /// The white matter variance of the different shells
  const std::vector<vector<int> >           _grps;   /// Indices (0-offs) to indicate what shell a scan belongs to
  const std::vector<double>                 _grpb;   /// b-values of the shells (bvals within a range of 100 = same)
  const double                              _s2n;    /// Tentative noise variance.
  std::vector<int>                          _grp;    /// Group index for each scan. Same info as in _grps, organised differently.
  NEWMAT::Matrix                            _K;      /// The resulting K-matrix (in the order given by scans)

  /// Returns the prediction vector for (observed) scan with index indx.
  NEWMAT::RowVector get_prediction_vector(unsigned int indx,
					  double        sigman2=-1.0) const;
  /// Returns the prediction vector for (unobserved) scan with diffusion parameters dpars
  NEWMAT::RowVector get_prediction_vector(const DiffPara& dpars,
					  double          sigman2=-1.0) const;
  std::vector<int>    reorganize_grps() const;
  std::vector<double> get_group_means(const NEWMAT::ColumnVector& y) const;
  void subtract_group_means(const std::vector<double>& means,
			    NEWMAT::ColumnVector&      y) const;
  void add_group_means(const std::vector<double>& means,
		       NEWMAT::ColumnVector&      y) const;
  NEWMAT::Matrix      get_iK(double sigman2) const;
  NEWMAT::Matrix      calculate_K() const;
  NEWMAT::RowVector   extract_K_row(unsigned int indx) const;
  NEWMAT::RowVector   calculate_K_row(const DiffPara&    dpar) const;  
};

/****************************************************************//**
*
* \brief Class used to make Gaussian process based predictions 
* about diffusion data.
*
* Will make predictions about observed (smoothing) and unobserved 
* (interpolation) scans using Gaussian Processes. It works by first
* creating an object, setting the # of scans one wants it to contain
* and then populate it with scans. Once that is done it is ready 
* to start making predictions. If you think of the signal (in a 
* given voxel) as points on a surface the Gaussian process will make
* predictions from the assumption that that surface should be smooth.
*
********************************************************************/ 
class DiffusionGP : public DWIPredictionMaker
{
public:
  /// Default constructor
  DiffusionGP() : _pop(true) {}
  /// Constructor that takes filenames from which to load data
  DiffusionGP(const std::string&                 scans_fname,
              const std::string&                 var_mask_fname,
              const std::vector<DiffPara>&       dpars);
  bool IsPopulated() const;                  // Returns true if all data present
  bool IsValid() const { return(_gpmod != NULL); }   // Returns true if ready to make predictions
  void SetNoOfScans(unsigned int n);
  void AddScan(const NEWIMAGE::volume<float>& scan,
               const DiffPara&                dp);
  void SetScan(const NEWIMAGE::volume<float>& scan,
               const DiffPara&                dp,
               unsigned int                   indx);
  void EvaluateModel(const NEWIMAGE::volume<float>& mask);
  NEWIMAGE::volume<float> Predict(double       sigman2,
                                  unsigned int indx) const;
  NEWIMAGE::volume<float> Predict(unsigned int indx) const { return(this->Predict(_gpmod->NoiseVariance(),indx)); }
  NEWIMAGE::volume<float> Predict(double            sigman2,
				  const DiffPara&   dpars) const;
  NEWIMAGE::volume<float> Predict(const DiffPara&   dpars) const { return(this->Predict(_gpmod->NoiseVariance(),dpars)); }
  NEWIMAGE::volume4D<float> PredictAll(double sigman2) const;
  NEWIMAGE::volume4D<float> PredictAll() const { return(this->PredictAll(_gpmod->NoiseVariance())); }

  // Here starts routines that should be used for debugging only.
  /// Returns the K-matrix of the Gaussian Process. Used for debugging.
  const NEWMAT::Matrix& GetKMatrix() const { if (_gpmod != NULL) return(_gpmod->KMatrix()); else throw EddyException("DiffusionGP::GetKMatrix: invalid predictor"); }
  /// Returns the signal variances for the different shells as assessed in white matter. Used for debugging.
  const std::vector<double>& GetWMVariances() const { if (_gpmod != NULL) return(_gpmod->WMVariances()); else throw EddyException("DiffusionGP::GetWMVariances: invalid predictor"); }
  /// Returns the noise variance as assessed outside the object. Used for debugging.
  double GetNoiseVariance() const { if (_gpmod != NULL) return(_gpmod->NoiseVariance()); else throw EddyException("DiffusionGP::GetNoiseVariance: invalid predictor"); }
  /// Returns a vector with the b-values for the different shells. Used for debugging.
  const std::vector<double>& GetbValues() const { if (_gpmod != NULL) return(_gpmod->_grpb); else throw EddyException("DiffusionGP::GetbValues: invalid predictor"); }
  /// Returns a vector of vectors of indicies. Each vector contains the indicies (zero offset) for one shell. Used for debugging.
  const std::vector<vector<int> >& Get_grps() const { if (_gpmod != NULL) return(_gpmod->_grps); else throw EddyException("DiffusionGP::Get_grps: invalid predictor"); }
  /// Returns a vector of group indicies. For example if rvec[4] is 0 it means that the fifth scan is in the first shell. Used for debugging.
  const std::vector<int>& Get_grp() const { if (_gpmod != NULL) return(_gpmod->_grp); else throw EddyException("DiffusionGP::Get_grp: invalid predictor"); }
  /// Writes various debug info to a file given by fname. Used for debugging.
  void WriteDebugInfo(const std::string& fname) const;

private:
  std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > > _sptrs;  // Pointers to the scans
  boost::shared_ptr<GPModel>                                _gpmod;  // The actual gaussian process
  std::vector<DiffPara>                                     _dpars;  // Diffusion parameters
  mutable bool                                              _pop;    // Tells if all data is present

  bool get_y(// Input
	     unsigned int           i,
	     unsigned int           j,
	     unsigned int           k,
	     // Output
	     NEWMAT::ColumnVector&  y) const;
};

/****************************************************************//**
*
* \brief Class used to calculate the variance (across g-vectors)
* within a shell. Only has static members.
*
* Contains utility routines for calculating signal variance (from
* white matter voxels) and noise variance (from outside the object).
* It contains only static declared functions so no instances of the
* class can be made. It is turned into a class only to keep the EDDY
* namespace less cluttered.
*
********************************************************************/ 
class VarianceCalculator
{
public:
  /// Calculates variance (across scans) confined to white matter.
  static std::vector<double> GetWMVariances(// Input
					    const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
					    const NEWIMAGE::volume<float>&                                     mask,
					    const std::vector<DiffPara>&                                       dpars,
					    // Output
					    std::vector<vector<int> >&                                         grps,
					    std::vector<double>&                                               grpb);

  /// Calculates noise variance based on voxels outside the object
  static double GetNoiseVariance(// Input
				 const std::vector<boost::shared_ptr<NEWIMAGE::volume<float> > >&   sptrs,
				 const NEWIMAGE::volume<float>&                                     mask,
				 const std::vector<DiffPara>&                                       dpars);
private:
  static std::vector<vector<int> > get_groups(const std::vector<DiffPara>&   dpars) { 
    std::vector<double> skrutt; return(get_groups(dpars,skrutt));
  }
  static std::vector<vector<int> > get_groups(// Input 
					      const std::vector<DiffPara>&   dpars,
					      // Output
					      std::vector<double>&           grpb);

};

// }}} End of fold.

} // End namespace EDDY

#endif // End #ifndef DiffusionGP_h
