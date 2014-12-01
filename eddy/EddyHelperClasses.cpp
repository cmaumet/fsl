// Definitions of classes that implements useful
// concepts for the eddy current project.
// 
// EddyHelperClasses.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "EddyUtils.h"

using namespace EDDY;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class DiffPara (Diffusion Parameters)
//
// This class manages the diffusion parameters for a given
// scan. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bool DiffPara::operator==(const DiffPara& rhs) const
{ 
  return(EddyUtils::AreInSameShell(*this,rhs) && EddyUtils::HaveSameDirection(*this,rhs)); 
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class AcqPara (Acquisition Parameters)
//
// This class manages the acquisition parameters for a given
// scan. 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

AcqPara::AcqPara(const NEWMAT::ColumnVector&   pevec,
                 double                        rotime)
: _pevec(pevec), _rotime(rotime)
{
  if (pevec.Nrows() != 3) throw EddyException("AcqPara::AcqPara: Wrong dimension pe-vector");
  if (rotime < 0.01 || rotime > 0.2) throw EddyException("AcqPara::AcqPara: Unrealistic read-out time");
  int cc = 0; //Component count
  for (int i=0; i<3; i++) { if (fabs(pevec(i+1)) > 1e-6) cc++; }
  if (!cc) throw EddyException("AcqPara::AcqPara: Zero Phase-encode vector");
  if (cc > 1) throw EddyException("AcqPara::AcqPara: Oblique pe-vectors not yet implemented");
}

bool AcqPara::operator==(const AcqPara& rh) const
{
  if (fabs(this->_rotime-rh._rotime) > 1e-6) return(false);
  for (int i=0; i<3; i++) {
    if (fabs(this->_pevec(i+1)-rh._pevec(i+1)) > 1e-6) return(false);
  }
  return(true);
}

DiffStats::DiffStats(const NEWIMAGE::volume<float>& diff, const NEWIMAGE::volume<float>& mask)
: _md(diff.zsize(),0.0), _msd(diff.zsize(),0.0), _n(mask.zsize(),0)
{
  for (int k=0; k<diff.zsize(); k++) {
    for (int j=0; j<diff.ysize(); j++) {
      for (int i=0; i<diff.xsize(); i++) {
	if (mask(i,j,k)) {
          _md[k] += diff(i,j,k);
          _msd[k] += diff(i,j,k)*diff(i,j,k);
	  _n[k] += 1;
	}
      }
    }
    if (_n[k]) { _md[k] /= double(_n[k]); _msd[k] /= double(_n[k]); }
  }  
}

/*!
 * Will write three files containing the mean differences, the
 * mean squared differences and the number of valid voxels. The
 * organisation of the (text) files is such that the nth column 
 * of the mth row corresponds to the nth slice for the mth scan.
 * \param bfname Base file name from which will be created 'bfname'.MeanDifference, 
 * 'bfname'.MeanSquaredDifference and 'bfname'.NoOfVoxels.
 */
void DiffStatsVector::Write(const std::string& bfname) const
{
  std::string fname = bfname + std::string(".MeanDifference");
  std::ofstream file;
  file.open(fname.c_str(),ios::out|ios::trunc);
  for (unsigned int i=0; i<_n; i++) file << _ds[i].MeanDifferenceVector() << endl;
  file.close();

  fname = bfname + std::string(".MeanSquaredDifference");
  file.open(fname.c_str(),ios::out|ios::trunc);
  for (unsigned int i=0; i<_n; i++) file << _ds[i].MeanSqrDiffVector() << endl;
  file.close();

  fname = bfname + std::string(".NoOfVoxels");
  file.open(fname.c_str(),ios::out|ios::trunc);
  for (unsigned int i=0; i<_n; i++) file << _ds[i].NVoxVector() << endl;
  file.close();
}

/*!
 * Returns a vector of vectors where each vector is equivalent to the output
 * from DiffStatsVector::GetOutliersInSlice. So for example the ith vector of the return value
 * correponds to the output from a call DiffStatsVector::GetOutliersInSlice with sl=i.
 * \param sl Slice index. It is zero offset so that e.g. sl=4 implies the 5th slice.
 * \param nstdev Indicates how many standard deviations imply an outlier.
 * \param minn Indicates the minimum number of valid voxels in a slice for it to be considered to have robust/valid statistic (and hence for it to be considered an outlier).
 */
std::vector<std::vector<unsigned int> > DiffStatsVector::GetOutliers(double        nstdev,     // # of std to qualify as outlier
								    unsigned int  minn) const  // Min # of voxels to be considered
{
  std::vector<std::vector<unsigned int> > rvv(NSlice());
  for (unsigned int sl=0; sl<NSlice(); sl++) rvv[sl] = GetOutliersInSlice(sl,nstdev,minn);
  return(rvv);
}

/*!
 * Returns a vector of (scan) indicies for the outliers for a scan given by sl. So if
 * for example scan 10 and scan 36 (zero offset) are outliers for slice 6 (zero offset)
 * it will return a vector [10 36] given the call GetOutliersInSlice(6). If there are
 * no outliers an empty (zero size) vector will be returned.
 * \param sl Slice index. It is zero offset so that e.g. sl=4 implies the 5th slice.
 * \param nstdev Indicates how many standard deviations imply an outlier.
 * \param minn Indicates the minimum number of valid voxels in a slice for it to be considered to have robust/valid statistic (and hence for it to be considered an outlier).
 */
std::vector<unsigned int> DiffStatsVector::GetOutliersInSlice(unsigned int  sl,
							      double        nstdev,      // # of std to qualify as outlier
							      unsigned int  minn) const  // Min # of voxels to be considered
{
  std::vector<unsigned int> rv;
  double mv = slice_mean_mean_diff(sl,minn);
  double stdev = slice_stdev_mean_diff(sl,minn);
  // printf("ds.GetOutliersInSlice: Slice %d has a mean diff of %f and an stdev of%f\n",sl,mv,stdev);
  for (unsigned int s=0; s<NScan(); s++) {
    if (_ds[s].NVox(int(sl)) >= minn && // If enough voxels to give reliable estimate
	std::fabs(_ds[s].MeanDifference(int(sl))-mv) > nstdev*stdev) { // If outlier
      rv.push_back(s);
    }
  }
  return(rv);
}

std::vector<double> DiffStatsVector::GetNoOfStdevInSlice(unsigned int  sl,
							 unsigned int  minn) const  // Min # of voxels to be considered
{
  std::vector<double> rv(NScan(),0.0);
  double mv = slice_mean_mean_diff(sl,minn);
  double stdev = slice_stdev_mean_diff(sl,minn);
  // printf("ds.GetOutliersInSlice: Slice %d has a mean diff of %f and an stdev of%f\n",sl,mv,stdev);
  for (unsigned int s=0; s<NScan(); s++) {
    if (_ds[s].NVox(int(sl)) >= minn) { // If enough voxels to give reliable estimate
      rv[s] = (_ds[s].MeanDifference(int(sl))-mv)/stdev;
    }
  }
  return(rv);
}

double DiffStatsVector::slice_mean_mean_diff(unsigned int sl,
                                             unsigned int minn) const // Min # of voxels to be considered
{
  double mv = 0.0;
  unsigned int nsl = 0;
  for (unsigned int s=0; s<NScan(); s++) {
    if (_ds[s].NVox(int(sl)) >= minn) { mv += _ds[s].MeanDifference(int(sl)); nsl++; }
  }
  mv = (nsl) ? mv/double(nsl) : 0.0;
  return(mv);
}

double DiffStatsVector::slice_stdev_mean_diff(unsigned int sl,
					      unsigned int minn) const // Min # of voxels to be considered
{
  double stdev = 0.0; 
  double mv = slice_mean_mean_diff(sl,minn);
  unsigned int nsl = 0;
  for (unsigned int s=0; s<NScan(); s++) {
    if (_ds[s].NVox(int(sl)) >= minn) { stdev += MISCMATHS::Sqr(mv - _ds[s].MeanDifference(int(sl))); nsl++; }
  }
  stdev = (nsl>1) ? std::sqrt(stdev/double(nsl-1)) : 0.0;
  return(stdev);
}

ReplacementManager::ReplacementManager(const DiffStatsVector& dsv, double nstdev, unsigned int minn)
: _nstdev(nstdev), _minn(minn), _ovv(dsv.NSlice()), _nsv(dsv.NSlice())
{
  for (unsigned int sl=0; sl<dsv.NSlice(); sl++) {
    _ovv[sl].resize(dsv.NScan(),false);
    std::vector<unsigned int> oisl = dsv.GetOutliersInSlice(sl,_nstdev,_minn);
    for (unsigned int i=0; i<oisl.size(); i++) _ovv[sl][oisl[i]] = true;
    _nsv[sl] = dsv.GetNoOfStdevInSlice(sl,_minn);
  } 
}

void ReplacementManager::Update(const DiffStatsVector& dsv)
{
  if (dsv.NSlice() != _ovv.size() || dsv.NScan() != _ovv[0].size()) {
    throw EddyException("ReplacementManager::Update: Mismatched DiffStatsVector object");
  }
  for (unsigned int sl=0; sl<dsv.NSlice(); sl++) {
    std::vector<unsigned int> oisl = dsv.GetOutliersInSlice(sl,_nstdev,_minn);
    // printf("rm.Update(): slice %d has %d outliers\n",sl,oisl.size());
    for (unsigned int i=0; i<oisl.size(); i++) _ovv[sl][oisl[i]] = true;
    _nsv[sl] = dsv.GetNoOfStdevInSlice(sl,_minn);
  }   
}

std::vector<unsigned int> ReplacementManager::OutliersInScan(unsigned int scan) const
{
  scan_oor(scan);
  std::vector<unsigned int> ol;
  for (unsigned int sl=0; sl<NSlice(); sl++) if (_ovv[sl][scan]) ol.push_back(sl);
  return(ol);
}

void ReplacementManager::WriteReport(const std::vector<unsigned int>& i2i,
				     const std::string&               fname) const
{
  std::ofstream fout;
  fout.open(fname.c_str(),ios::out|ios::trunc);
  if (fout.fail()) {
    cout << "Failed to open outlier report file " << fname << endl;
    return;
  }
  unsigned int nsl = _ovv.size();
  unsigned int nscan = _ovv[0].size();
  for (unsigned int sl=0; sl<nsl; sl++) {
    for (unsigned int s=0; s<nscan; s++) {
      if (_ovv[sl][s]) {
	fout << "Slice " << sl << " in scan " << i2i[s] << " is an outlier " << _nsv[sl][s] << " standard deviations off" << endl;
      }
    }
  }
  fout.close();
  return;
}

void ReplacementManager::DumpOutlierMap(const std::string& bfname) const
{
  std::string fname = bfname + ".OutlierMap";
  std::ofstream fout;
  fout.open(fname.c_str(),ios::out|ios::trunc);
  fout << "One row per scan, one column per slice. Outlier: 1, Non-outlier: 0" << endl;
  for (unsigned scan=0; scan<NScan(); scan++) {
    for (unsigned int slice=0; slice<NSlice(); slice++) {
      fout << _ovv[slice][scan] << " ";
    }
    fout << endl;
  } 
  fout.close();

  fname = bfname + ".NoOfStdevMap";
  fout.open(fname.c_str(),ios::out|ios::trunc);
  fout << "One row per scan, one column per slice." << endl;
  for (unsigned scan=0; scan<NScan(); scan++) {
    for (unsigned int slice=0; slice<NSlice(); slice++) {
      fout << _nsv[slice][scan] << " ";
    }
    fout << endl;
  } 
  fout.close();
}
