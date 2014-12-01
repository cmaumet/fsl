/*! \file EddyHelperClasses.h
    \brief Contains declaration of classes that implements useful functionality for the eddy project.

    \author Jesper Andersson
    \version 1.0b, Sep., 2012.
*/
// Declarations of classes that implements useful
// concepts for the eddy current project.
// 
// EddyHelperClasses.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef EddyHelperClasses_h
#define EddyHelperClasses_h

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

namespace EDDY {

enum Parameters { MOVEMENT, EC, ALL };
enum ECModel { Linear, Quadratic, Cubic, Unknown };
enum ScanType { DWI, B0 , ANY };
enum FinalResampling { JAC, LSR, UNKNOWN_RESAMPLING};

/****************************************************************//**
*
* \brief This is the exception that is being thrown by routines
*  in the core code of eddy.
*
* This is the exception that is being thrown by routines
* in the core code of eddy.
*
********************************************************************/ 
class EddyException: public std::exception
{
private:
  std::string m_msg;
public:
  EddyException(const std::string& msg) throw(): m_msg(msg) { cout << what() << endl; }

  virtual const char * what() const throw() {
    return string("eddy: msg=" + m_msg).c_str();
  }

  ~EddyException() throw() {}
};

/****************************************************************//**
*
* \brief This class manages the diffusion parameters for one scan
*
* This class manages the diffusion parameters for one scan
*
********************************************************************/ 
class DiffPara
{
public:
  /// Default constructor. Sets b-vector to [1 0 0] and b-value to zero.
  DiffPara() { _bvec.ReSize(3); _bvec(1)=1; _bvec(2)=0; _bvec(3)=0; _bval = 0; }
  /// Contructs a DiffPara object from a b-vector and a b-value.
  DiffPara(const NEWMAT::ColumnVector&   bvec,
	   double                        bval) : _bvec(bvec), _bval(bval)
  {
    if (_bvec.Nrows() != 3) throw EddyException("DiffPara::DiffPara: b-vector must be three elements long");
    if (_bval < 0) throw EddyException("DiffPara::DiffPara: b-value must be non-negative");
    if (_bval) _bvec /= _bvec.NormFrobenius();
  }
  /// Prints out b-vector and b-value in formatted way  
  friend ostream& operator<<(ostream& op, const DiffPara& dp) { op << "b-vector: " << dp._bvec.t() << endl << "b-value:  " << dp._bval << endl; return(op); }
  /// Returns true if the b-value AND the direction are the same
  bool operator==(const DiffPara& rhs) const;
  /// Same as !(a==b)
  bool operator!=(const DiffPara& rhs) const { return(!(*this==rhs)); }
  /// Compares the b-values
  bool operator>(const DiffPara& rhs) const { return(this->bVal()>rhs.bVal()); }
  /// Compares the b-values
  bool operator<(const DiffPara& rhs) const { return(this->bVal()<rhs.bVal()); }
  /// Returns a normalised b-vector
  NEWMAT::ColumnVector bVec() const { return(_bvec); }
  /// Returns the b-value
  double bVal() const { return(_bval); }
private:
  NEWMAT::ColumnVector _bvec;
  double               _bval;
};

/****************************************************************//**
*
* \brief This class manages the acquisition parameters for one scan
*
* This class manages the acquisition parameters for one scan
*
********************************************************************/ 
class AcqPara
{
public:
  /// Constructs an AcqPara object from a phase-encode vector and a total readout-time (sec)
  AcqPara(const NEWMAT::ColumnVector&   pevec,
          double                        rotime);
  /// Prints out phase-encode vactor and readout-time (sec) in formatted way  
  friend ostream& operator<<(ostream& op, const AcqPara& ap) { op << "Phase-encode vector: " << ap._pevec.t() << endl << "Read-out time: " << ap._rotime; return(op); }
  /// Returns true if both PE direction and readout time are the same.
  bool operator==(const AcqPara& rh) const;
  /// Same as !(a==b)
  bool operator!=(const AcqPara& rh) const { return(!(*this == rh)); }
  /// Returns the phase-enocde vector
  NEWMAT::ColumnVector PhaseEncodeVector() const { return(_pevec); }
  /// Returns the readout-time in seconds
  double ReadOutTime() const { return(_rotime); }
private:
  NEWMAT::ColumnVector _pevec;
  double               _rotime;
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class MaskManager
//
// This class manages an And-mask.
// scan. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class MaskManager
{
public:
  MaskManager(const NEWIMAGE::volume<float>& mask) : _mask(mask) {}
  MaskManager(int xs, int ys, int zs) : _mask(xs,ys,zs) { _mask = 1.0; }
  void ResetMask() { _mask = 1.0; }
  void SetMask(const NEWIMAGE::volume<float>& mask) { if (!NEWIMAGE::samesize(_mask,mask)) throw EddyException("EDDY::MaskManager::SetMask: Wrong dimension"); else _mask = mask;}
  void UpdateMask(const NEWIMAGE::volume<float>& mask) { if (!NEWIMAGE::samesize(_mask,mask)) throw EddyException("EDDY::MaskManager::UpdateMask: Wrong dimension"); else _mask *= mask;}
  const NEWIMAGE::volume<float>& GetMask() const { return(_mask); }
private:
  NEWIMAGE::volume<float> _mask;
};

/****************************************************************//**
*
* \brief This class manages stats on slice wise differences.
*
* This class calculates and serves up information about slice-wise
*  (in observation space) statistics on the difference between 
*  an observation (scan) and the prediction. 
*
********************************************************************/ 
class DiffStats
{
public:
  DiffStats() {}
  /// Constructs a Diffstats object given a difference volume and a mask.
  DiffStats(const NEWIMAGE::volume<float>& diff, const NEWIMAGE::volume<float>& mask);
  /// Returns the mean (across all valid voxels in the volume) difference.
  double MeanDifference() const { return(mean_stat(_md)); }
  /// Returns the mean (across all valid voxels in slice sl (zero-offset)) difference.
  double MeanDifference(int sl) const { if (index_ok(sl)) return(_md[sl]); else return(0.0); }
  /// Returns a vector with the mean difference for all slices
  NEWMAT::RowVector MeanDifferenceVector() const { return(get_vector(_md)); }
  /// Returns the mean (across all valid voxels in the volume) squared difference.
  double MeanSqrDiff() const { return(mean_stat(_msd)); }
  /// Returns the mean (across all valid voxels in slice sl (zero-offset)) squared difference.
  double MeanSqrDiff(int sl) const { if (index_ok(sl)) return(_msd[sl]); else return(0.0); }
  /// Returns a vector with the mean squared difference for all slices
  NEWMAT::RowVector MeanSqrDiffVector() const { return(get_vector(_msd)); }
  /// Number of valid voxels in the whole volume (as determined by the mask passed to the constructor)
  unsigned int NVox() const { unsigned int n=0; for (int i=0; i<int(_n.size()); i++) n+=_n[i]; return(n); }
  /// Number of valid voxels in slice sl (zero offset).
  unsigned int NVox(int sl) const { if (index_ok(sl)) return(_n[sl]); else return(0); }
  /// Vector with the number of valid voxels in each slice.
  NEWMAT::RowVector NVoxVector() const { return(get_vector(_n)); }
  /// Number of slices.
  unsigned int NSlice() const { return(_n.size()); }
private:
  std::vector<double>        _md;  // Slice wise mean difference
  std::vector<double>        _msd; // Slice wise mean squared difference
  std::vector<unsigned int>  _n;   // Slice wise # of valid pixels

  bool index_ok(int sl) const 
  { if (sl<0 || sl>=int(_n.size())) throw EddyException("DiffStats::index_ok: Index out of range"); return(true); }

  double mean_stat(const std::vector<double>& stat) const 
  { double ms=0; for (int i=0; i<int(_n.size()); i++) ms+=_n[i]*stat[i]; ms/=double(NVox()); return(ms); }

  template<class T>
  NEWMAT::RowVector get_vector(const std::vector<T>& stat) const
  { NEWMAT::RowVector ov(stat.size()); for (unsigned int i=0; i<stat.size(); i++) ov(i+1) = double(stat[i]); return(ov); }
};

/****************************************************************//**
*
* \brief This class manages a set (one for each scan) of DiffStats objects.
*
* This class manages a vector of DiffStats objects (one for each scan).
* This means that it can look across scans (for a given slice) and 
* build up statistics of the statistics from the DiffStats objects.
* It can for example calculate the mean and standard deviation (across) 
* subjects of the slice-wise mean differences from the DiffStat objects.
* From that it can the determine how many standard deviations away
* a given scan and slice is from the mean and hence identify out-liers.
*
********************************************************************/ 
class DiffStatsVector
{
public:
  /// Constructs an object with n (empty) slots for DiffStats objects.
  DiffStatsVector(unsigned int n) : _n(n) { _ds = new DiffStats[_n]; }
  /// Copy constructor
  DiffStatsVector(const DiffStatsVector& rhs) : _n(rhs._n) { _ds = new DiffStats[_n]; for (unsigned int i=0; i<_n; i++) _ds[i] = rhs._ds[i]; }
  ~DiffStatsVector() { delete [] _ds; }
  /// Assignment
  DiffStatsVector& operator=(const DiffStatsVector& rhs) {
    delete [] _ds; _n = rhs._n; _ds = new DiffStats[_n]; for (unsigned int i=0; i<_n; i++) _ds[i] = rhs._ds[i]; return(*this);
  }
  /// Gives read-access to the ith (zero offset) DiffStats object in the vector.
  const DiffStats& operator[](unsigned int i) const { throw_if_oor(i); return(_ds[i]); }
  /// Gives read/write-access to the ith (zero offset) DiffStats object in the vector. This is used to populate the vector.
  DiffStats& operator[](unsigned int i) { throw_if_oor(i); return(_ds[i]); }
  /// Returns the number of DiffStats objects in the vector.
  unsigned int NScan() const { return(_n); }
  /// Returns the number of slices in each of the DiffStats objects.
  unsigned int NSlice() const { return(_ds[0].NSlice()); }
  /// Writes three files with information relevent for debugging.
  void Write(const std::string& bfname) const;
  /// Returns a vector of vectors if scan indicies for outliers. One element of the vector corresponds to the output from GetOutliersInSlice().
  std::vector<std::vector<unsigned int> > GetOutliers(double nstdev=4.0, unsigned int minn=100) const;
  /// Returns a vector of scan indicies for the outliers in slice sl (zero-offset).
  std::vector<unsigned int> GetOutliersInSlice(unsigned int sl, double nstdev=4.0, unsigned int minn=100) const;
  /// Returns a vector (NScan() long) with the number of std for each scan for a slice given by sl.
  std::vector<double> GetNoOfStdevInSlice(unsigned int sl, unsigned int minn=100) const;
private:
  unsigned int _n;   // Length of vector
  DiffStats    *_ds; // Old fashioned C vector of DiffStats objects

  double slice_mean_mean_diff(unsigned int sl, unsigned int minn) const;
  double slice_stdev_mean_diff(unsigned int sl, unsigned int minn) const;
  void throw_if_oor(unsigned int i) const { if (i >= _n) throw EddyException("DiffStatsVector::throw_if_oor: Index out of range"); }
};

/****************************************************************//**
*
* \brief This class decides and keeps track of which slices in which 
* scans should be replaced by their predictions
*
********************************************************************/ 
class ReplacementManager {
public:
  ReplacementManager(unsigned int nscan, unsigned int nsl, double nstdev=4.0, unsigned int minn=100) : _nstdev(nstdev), _minn(minn), _ovv(nsl), _nsv(nsl)
  { for (unsigned int sl=0; sl<nsl; sl++) { _ovv[sl].resize(nscan,false); _nsv[sl].resize(nscan,0.0);} }
  ReplacementManager(const DiffStatsVector& dsv, double nstdev=4.0, unsigned int minn=100);
  ~ReplacementManager() {}
  unsigned int NSlice() const { return(_ovv.size()); }
  unsigned int NScan() const { unsigned int rval = (_ovv.size()) ? _ovv[0].size() : 0; return(rval); }
  void Update(const DiffStatsVector& dsv);
  std::vector<unsigned int> OutliersInScan(unsigned int scan) const;
  void WriteReport(const std::vector<unsigned int>& i2i,
		   const std::string&               bfname) const;
  // For debugging
  void DumpOutlierMap(const std::string& fname) const;
private:
  double                            _nstdev;  // # of standard dev off to qualify as outlier
  unsigned int                      _minn;    // Min # of voxels in slice to be considered reliable
  std::vector<std::vector<bool> >   _ovv;     // _ovv[slice][scan] is an outlier if set
  std::vector<std::vector<double> > _nsv;     // _nsv[slice][scan] tells how many stdev off that slice-scan is.

  void scan_oor(unsigned int scan) const { if (!_ovv.size() || scan >= _ovv[0].size()) throw EddyException("ReplacementManager::scan_oor: Scan index out of range"); }
};

/****************************************************************//**
*
* \brief Helper class that manages a set of image coordinates in a way that 
* it enables calculation/implementation of partial derivatives of
* images w.r.t. transformation parameters.
*
********************************************************************/ 
class ImageCoordinates
{
public:
  ImageCoordinates(unsigned int xn, unsigned int yn, unsigned int zn, bool init=false) 
  : _xn(xn), _yn(yn), _zn(zn)
  {
    _x = new float[_xn*_yn*_zn];
    _y = new float[_xn*_yn*_zn];
    _z = new float[_xn*_yn*_zn];
    if (init) {
      for (unsigned int k=0, indx=0; k<_zn; k++) {
	for (unsigned int j=0; j<_yn; j++) {
	  for (unsigned int i=0; i<_xn; i++) {
	    _x[indx] = float(i);
	    _y[indx] = float(j);
	    _z[indx++] = float(k);
	  }
	}
      }
    }
  }
  ImageCoordinates(const ImageCoordinates& inp) 
  : _xn(inp._xn), _yn(inp._yn), _zn(inp._zn) 
  {
    _x = new float[_xn*_yn*_zn]; std::memcpy(_x,inp._x,_xn*_yn*_zn*sizeof(float));
    _y = new float[_xn*_yn*_zn]; std::memcpy(_y,inp._y,_xn*_yn*_zn*sizeof(float));
    _z = new float[_xn*_yn*_zn]; std::memcpy(_z,inp._z,_xn*_yn*_zn*sizeof(float));   
  }
  ~ImageCoordinates() { delete[] _x; delete[] _y; delete[] _z; }
  ImageCoordinates& operator=(const ImageCoordinates& rhs) {
    if (this == &rhs) return(*this);
    delete[] _x; delete[] _y; delete[] _z;
    _xn = rhs._xn; _yn = rhs._yn; _zn = rhs._zn; 
    _x = new float[_xn*_yn*_zn]; std::memcpy(_x,rhs._x,_xn*_yn*_zn*sizeof(float));
    _y = new float[_xn*_yn*_zn]; std::memcpy(_y,rhs._y,_xn*_yn*_zn*sizeof(float));
    _z = new float[_xn*_yn*_zn]; std::memcpy(_z,rhs._z,_xn*_yn*_zn*sizeof(float));   
    return(*this);
  }
  ImageCoordinates& operator-=(const ImageCoordinates& rhs) {
    if (_xn != rhs._xn || _yn != rhs._yn || _zn != rhs._zn) throw EddyException("ImageCoordinates::operator-= size mismatch");
    for (unsigned int i=0; i<_xn*_yn*_zn; i++) { _x[i]-=rhs._x[i]; _y[i]-=rhs._y[i]; _z[i]-=rhs._z[i]; }
    return(*this);
  }
  ImageCoordinates operator-(const ImageCoordinates& rhs) const {
    return(ImageCoordinates(*this)-=rhs);
  }
  NEWIMAGE::volume<float> operator*(const NEWIMAGE::volume4D<float>& vol) {
    if (int(_xn) != vol.xsize() || int(_yn) != vol.ysize() || int(_zn) != vol.zsize() || vol.tsize() != 3) {
      throw EddyException("ImageCoordinates::operator* size mismatch");
    }
    NEWIMAGE::volume<float> ovol = vol[0];
    for (unsigned int k=0, indx=0; k<_zn; k++) {
      for (unsigned int j=0; j<_yn; j++) {
	for (unsigned int i=0; i<_xn; i++) {
          ovol(i,j,k) = _x[indx]*vol(i,j,k,0) + _y[indx]*vol(i,j,k,1) + _z[indx]*vol(i,j,k,2);
          indx++;
	}
      }
    }
    return(ovol);
  }

  void Write(const std::string& fname) const
  {
    NEWMAT::Matrix omat(N(),3);
    for (unsigned int i=0; i<N(); i++) {
      omat(i+1,1) = x(i); omat(i+1,2) = y(i); omat(i+1,3) = z(i); 
    }
    MISCMATHS::write_ascii_matrix(fname,omat);
  }

  unsigned int N() const { return(_xn*_yn*_zn); }
  unsigned int NX() const { return(_xn); }
  unsigned int NY() const { return(_yn); }
  unsigned int NZ() const { return(_zn); }
  const float& x(unsigned int i) const { return(_x[i]); }
  const float& y(unsigned int i) const { return(_y[i]); }
  const float& z(unsigned int i) const { return(_z[i]); }
  float& x(unsigned int i) { return(_x[i]); }
  float& y(unsigned int i) { return(_y[i]); }
  float& z(unsigned int i) { return(_z[i]); }
   
private:
  unsigned int _xn;
  unsigned int _yn;
  unsigned int _zn;
  float *_x;
  float *_y;
  float *_z;
};

} // End namespace EDDY

#endif // End #ifndef EddyHelperClasses_h

////////////////////////////////////////////////
//
// Here starts Doxygen documentation
//
////////////////////////////////////////////////

/*! 
 * \fn EDDY::DiffPara::DiffPara(const NEWMAT::ColumnVector& bvec, double bval)
 * Contructs a DiffPara object from a b-vector and a b-value.
 * \param bvec ColumnVector with three elements. Will be normalised, but must have norm different from zero.
 * \param bval b-value. Must be non-negative.
 */

/*!
 * \fn EDDY::DiffPara::operator==(const DiffPara& rhs) const
 *  Will return true if calls to both EddyUtils::AreInSameShell and EddyUtils::HaveSameDirection are true.
 */ 

/*! 
 * \fn AcqPara::AcqPara(const NEWMAT::ColumnVector& pevec, double rotime)
 * Constructs an AcqPara object from phase-encode direction and total read-out time. 
 * \param pevec Normalised vector desribing the direction of the phase-encoding. At present the third element 
 * (the z-direction) must be zero (i.e. it only allows phase-encoding in the xy-plane.
 * \param rotime The time between the collection of the midpoint of the first echo and the last echo.
 */

