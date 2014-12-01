/*! \file ECScanClasses.h
    \brief Declarations of classes that implements a scan or a collection of scans within the EC project.

    \author Jesper Andersson
    \version 1.0b, Sep., 2012.
*/
// Declarations of classes that implements a scan
// or a collection of scans within the EC project.
// 
// EddyScanClasses.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//
#pragma clang diagnostic ignored "-Wunknown-pragmas" // Ignore the OpenMP pragmas

#ifndef ECScanClasses_h
#define ECScanClasses_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "topup/topup_file_io.h"
#include "EddyHelperClasses.h"
#include "ECModels.h"

namespace EDDY {

enum HzModel { ModelMovements, ModelOffset };

/****************************************************************//**
*
* \brief This class manages one diffusion weighted or b=0 scan
*  for the eddy project.
*
********************************************************************/ 
class ECScan
{
public:
  /// Constructs an object from volume, acquisition parameters, diffusion parameters, EC-Model and session number.
  ECScan(const NEWIMAGE::volume<float>    ima,
         const AcqPara&                   acqp,
         const DiffPara&                  diffp,
         boost::shared_ptr<ScanECModel>   ecp,
         unsigned int                     sess) : _ima(ima), _acqp(acqp), _diffp(diffp), _ecp(ecp->Clone()), _fwhm(0.0), _sess(sess){}
  /// Copy constructor.
  ECScan(const ECScan& inp) : _ima(inp._ima), _sima(inp._sima), _acqp(inp._acqp), _diffp(inp._diffp), _ecp(inp._ecp->Clone()), _fwhm(inp._fwhm), _sess(inp._sess) {}
  virtual ~ECScan() {}
  /// Assignment.
  ECScan& operator=(const ECScan& rhs) {
    if (this == &rhs) return(*this);
    _ima=rhs._ima; _sima=rhs._sima; _acqp=rhs._acqp; _diffp=rhs._diffp; _ecp=rhs._ecp->Clone(); _fwhm=rhs._fwhm; _sess=rhs._sess;
    return(*this);
  }
  /// Session number (starting with session zero) of this scan.
  unsigned int Session() const { return(_sess); }
  /// Will return true if an offset (DC part) is included in the EC model.
  bool HasFieldOffset() const { return(_ecp->HasFieldOffset()); }
  /// Returns the field offset (in Hz).
  double GetFieldOffset() const { return(_ecp->GetFieldOffset()); }
  /// Sets the field offset. The value of ofst should be in Hz.
  void SetFieldOffset(double ofst) { _ecp->SetFieldOffset(ofst); }
  /// Returns a vector with the parameters for the linear part (if any) of the EC model. Will be three or zero elements.
  virtual NEWMAT::RowVector GetLinearParameters() const { return(_ecp->GetLinearParameters()); }
  /// Returns the original (unsmoothed and untransformed) volume.
  const NEWIMAGE::volume<float>& GetOriginalIma() const { return(_ima); }
  /// Returns the original (unsmoothed) volume after correction for rigid-body movement.
  NEWIMAGE::volume<float> GetMotionCorrectedOriginalIma(NEWIMAGE::volume<float>& omask) const { return(motion_correct(GetOriginalIma(),&omask)); }
  /// Returns the original (unsmoothed) volume after correction for rigid-body movement.
  NEWIMAGE::volume<float> GetMotionCorrectedOriginalIma() const { return(motion_correct(GetOriginalIma(),NULL)); }
  /// Returns the original (unsmoothed) volume after correction for rigid-body movement and EC distortions.
  NEWIMAGE::volume<float> GetUnwarpedOriginalIma(// Input
						 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
						 // Output
						 NEWIMAGE::volume<float>&                            omask) const;
  /// Returns the original (unsmoothed) volume after correction for rigid-body movement and EC distortions.
  NEWIMAGE::volume<float> GetUnwarpedOriginalIma(// Input
						 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc) const;
  /// Returns the smoothed volume in the original space (no corrections)
  const NEWIMAGE::volume<float>& GetIma() const { return( (_fwhm) ? _sima : _ima); } 
  /// Returns the smoothed and corrected (for movement and EC distortions) volume.
  NEWIMAGE::volume<float> GetUnwarpedIma(// Input
					 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
					 // Output
					 NEWIMAGE::volume<float>&                            omask) const;
  /// Returns the smoothed and corrected (for movement and EC distortions) volume.
  NEWIMAGE::volume<float> GetUnwarpedIma(// Input
					 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc) const;
  /// Returns the acquisition parameters for the scan
  AcqPara GetAcqPara() const { return(_acqp); }
  /// Returns the diffusion parameters for the scan
  DiffPara GetDiffPara() const { return(_diffp); }
  /// Returns a vector of displacements (in mm) from a field offset of 1Hz
  NEWMAT::ColumnVector GetHz2mmVector() const;
  /// Number of parameters for the EC model (including rigid body movement).
  unsigned int NParam(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->NParam(whichp)); }
  /// Returns all the current parameters for the EC model.
  NEWMAT::ColumnVector GetParams(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetParams(whichp)); }
  /// Number of parameters that are being optimised in the current ScanECModel.
  unsigned int NDerivs(EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->NDerivs(whichp)); }
  /// Get the parameter (of those being optimised) given by indx
  double GetDerivParam(unsigned int indx, EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetDerivParam(indx,whichp)); }
  /// Get the scale for the parameter (of those being optimised) given by indx
  double GetDerivScale(unsigned int indx, EDDY::Parameters whichp=EDDY::ALL) const { return(_ecp->GetDerivScale(indx,whichp)); }
  /// Set/Update the parameter (of those being optimised) given by indx
  void SetDerivParam(unsigned int indx, double p, EDDY::Parameters whichp=EDDY::ALL) { _ecp->SetDerivParam(indx,p,whichp); }
  /// Returns the FWHM of the smoothing that has been applied to the volume.
  double GetFWHM() const { return(_fwhm); }
  /// Returns matrix denoted \f$\mathbf{R}\f$ in paper.
  NEWMAT::Matrix ForwardMovementMatrix() const { return(_ecp->ForwardMovementMatrix(_ima)); }
  /// Returns matrix denoted \f$\mathbf{R}^{-1}\f$ in paper.
  NEWMAT::Matrix InverseMovementMatrix() const { return(_ecp->InverseMovementMatrix(_ima)); }
  /// Returns the EC-field relevant for the Observation->Model transform.
  NEWIMAGE::volume<float> ECField() const { return(_ecp->ECField(_ima)); }
  /// Returns the total field relevant for the Observation->Model transform
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask,
							 NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_scan_to_model_transform(susc,&omask,&jac));
  }
  /// Returns the total field relevant for the Observation->Model transform
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask) const {
    return(field_for_scan_to_model_transform(susc,&omask,NULL));
  }
  /// Returns the total field relevant for the Observation->Model transform
  NEWIMAGE::volume4D<float> FieldForScanToModelTransform(boost::shared_ptr<const NEWIMAGE::volume<float> > susc) const {
    return(field_for_scan_to_model_transform(susc,NULL,NULL));
  }
  /// Returns the total field relevant for the Model->Observation transform
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask,
							 NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_model_to_scan_transform(susc,&omask,&jac));
  }
  /// Returns the total field relevant for the Model->Observation transform
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(// Input
							 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							 // Output
							 NEWIMAGE::volume<float>&                          omask) const {
    return(field_for_model_to_scan_transform(susc,&omask,NULL));
  }
  /// Returns the total field relevant for the Model->Observation transform
  NEWIMAGE::volume4D<float> FieldForModelToScanTransform(boost::shared_ptr<const NEWIMAGE::volume<float> > susc) const {
    return(field_for_model_to_scan_transform(susc,NULL,NULL));
  }
  /// Returns the total field relevant for the Model->Observation transform and the Jacobian associated with it.
  NEWIMAGE::volume4D<float> FieldForModelToScanTransformWithJac(// Input
								boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								// Output
								NEWIMAGE::volume<float>&                          jac) const {
    return(field_for_model_to_scan_transform(susc,NULL,&jac));
  }
  /// Set/update the parameters for the EC model.
  void SetParams(const NEWMAT::ColumnVector& mpep, EDDY::Parameters whichp=EDDY::ALL) { _ecp->SetParams(mpep,whichp); }
  /// Smooth the image volume
  void SetFWHM(double fwhm);
  /// Replace selected slices. Useful for replacing outliers with predicted data.
  void ReplaceSlices(const NEWIMAGE::volume<float>&                     rep,
		     boost::shared_ptr<const NEWIMAGE::volume<float> >  susc,
		     const NEWIMAGE::volume<float>&                     inmask, 
		     const std::vector<unsigned int>&                   ol);
private:
  NEWIMAGE::volume<float>         _ima;   ///< The original (from disc) volume.
  NEWIMAGE::volume<float>         _sima;  ///< Smoothed version of _ima.
  AcqPara                         _acqp;  ///< The acquisition parameters associated with the volume.
  DiffPara                        _diffp; ///< The diffusion parameters associated with the volume.
  boost::shared_ptr<ScanECModel>  _ecp;   ///< The EC model.
  double                          _fwhm;  ///< The FWHM (mm) of the smoothing filter used for _ima->_sima.
  unsigned int                    _sess;  ///< Session number for this scan (relevant if data was acquired in multiple sessions).

  NEWIMAGE::volume<float> motion_correct(// Input
					 const NEWIMAGE::volume<float>&  inima,
					 // Output (optional)
					 NEWIMAGE::volume<float>         *omask) const;
  NEWIMAGE::volume<float> transform_to_model_space(// Input
						   const NEWIMAGE::volume<float>&                    inima,
						   boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
						   // Output
						   NEWIMAGE::volume<float>&                          omask) const;
  NEWIMAGE::volume4D<float> field_for_scan_to_model_transform(// Input
							      boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							      // Output
							      NEWIMAGE::volume<float>                           *omask,
							      NEWIMAGE::volume<float>                           *jac) const;
  NEWIMAGE::volume4D<float> field_for_model_to_scan_transform(// Input
							      boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							      // Output
							      NEWIMAGE::volume<float>                           *omask,
							      NEWIMAGE::volume<float>                           *jac) const;
};

/****************************************************************//**
*
* \brief This class manages a collection of scans (of type ECScan)
*  for the eddy project.
*
********************************************************************/ 
class ECScanManager
{
public:
  /// Constructor for the ECScanManager class
  ECScanManager(const std::string&               simafname,
                const std::string&               maskfname,
                const std::string&               acqfname,
                const std::string&               topupfname,
                const std::string&               bvecfname,
                const std::string&               bvalsfname,
                EDDY::ECModel                    ecmodel,
                const std::vector<unsigned int>& indicies,
		const std::vector<unsigned int>& session_indicies,
                unsigned int                     no_of_sessions);
  ~ECScanManager() {}
  /// Number of scans of type st
  unsigned int NScans(ScanType st=ANY) const;
  /// Number of LSR pairs of type st
  unsigned int NLSRPairs(ScanType st=ANY) const;
  /// Return number of sessions over which data was acquired
  unsigned int NoOfSessions() const { return(_nsess); }
  /// Return true if any scan has a PE component in the x-direction
  bool HasPEinX() const { return(has_pe_in_direction(1)); }
  /// Return true if any scan has a PE component in the y-direction
  bool HasPEinY() const { return(has_pe_in_direction(2)); }
  /// Return true if any scan has a PE component in the x- AND any scan in the y-direction
  bool HasPEinXandY() const { return(HasPEinX() && HasPEinY()); }
  /// Returns true if the scan indicated by indx (zero-offset) is diffusion weighted.
  bool IsDWI(unsigned int indx) const { if (indx>_fi.size()-1) throw EddyException("ECScanManager::IsDWI: index out of range"); else return(!_fi[indx].first); }
  /// Returns true if the scan indicated by indx (zero-offset) is a b=0 scan.
  bool IsB0(unsigned int indx) const { return(!IsDWI(indx)); }
  /// Returns scale factor that has been applied to data prior to registration process.
  double ScaleFactor() const { return(_sf); }
  /// Returns vector of Global (file) indicies for all dwi indicies.
  std::vector<unsigned int> GetDwi2GlobalIndexMapping() const;
  /// Returns vector of Global (file) indicies for all b0 indicies.
  std::vector<unsigned int> Getb02GlobalIndexMapping() const;
  /// Returns true if b0's have been spaced out in a way that can help with dwi movements.
  bool B0sAreInterspersed() const;
  /// Returns true if a susceptibilty induced off-resonance field has been set
  bool HasSuscHzOffResField() const { return(_has_topup); }
  /// Returns true if the EC model includes a field offset
  bool HasFieldOffset() const { return(Scan(0,DWI).HasFieldOffset()); }
  /// Returns true if data allows for LSR resampling
  bool CanDoLSRResampling() const;
  /// Returns the individual indicies for scans that constitute the i'th LSR pair
  std::pair<unsigned int,unsigned int> GetLSRPair(unsigned int i, ScanType st) const;
  /// Sets parameters for all scans
  void SetParameters(const NEWMAT::Matrix& pM, ScanType st=ANY);
  /// Sets parameters for all scans from values in file given by fname
  void SetParameters(const std::string& fname, ScanType st=ANY) { NEWMAT::Matrix pM = MISCMATHS::read_ascii_matrix(fname); SetParameters(pM,st); }
  /// Returns a read-only reference to a scan given by indx
  const ECScan& Scan(unsigned int indx, ScanType st=ANY) const;
  /// Returns a read-write reference to a scan given by indx
  ECScan& Scan(unsigned int indx, ScanType st=ANY);
  // Returns an "original" (no smoothing) "unwarped" into model space
  NEWIMAGE::volume<float> GetUnwarpedOrigScan(unsigned int              indx,
					      NEWIMAGE::volume<float>&  omask,
					      ScanType                  st=ANY) const; 
  NEWIMAGE::volume<float> GetUnwarpedOrigScan(unsigned int indx,
					      ScanType     st=ANY) const { 
    NEWIMAGE::volume<float> mask=_scans[0].GetIma(); mask=1.0; 
    return(GetUnwarpedOrigScan(indx,mask,st));
  }
  // Returns an image "unwarped" into model space
  NEWIMAGE::volume<float> GetUnwarpedScan(unsigned int              indx,
					  NEWIMAGE::volume<float>&  omask,
					  ScanType                  st=ANY) const; 
  NEWIMAGE::volume<float> GetUnwarpedScan(unsigned int indx,
					  ScanType     st=ANY) const { 
    NEWIMAGE::volume<float> mask=_scans[0].GetIma(); mask=1.0; 
    return(GetUnwarpedScan(indx,mask,st));
  }
  /// Resamples a matching pair of images using least-squares resampling.
  NEWIMAGE::volume<float> LSRResamplePair(// Input
					  unsigned int              i, 
					  unsigned int              j, 
					  ScanType                  st,
					  // Output
					  NEWIMAGE::volume<float>&  omask) const;
  /// Sets movement and EC parameters for a scan given by indx
  void SetScanParameters(unsigned int indx, const NEWMAT::ColumnVector& p, ScanType st=ANY) { Scan(indx,st).SetParams(p,ALL); }

  /// Smooths all scans to the requested FWHM
  void SetFWHM(double fwhm) {
    # pragma omp parallel for shared(fwhm)
    for (int i=0; i<int(_scans.size()); i++) {
      _scans[i].SetFWHM(fwhm);
    } 
    # pragma omp parallel for shared(fwhm)
    for (int i=0; i<int(_b0scans.size()); i++) {
      _b0scans[i].SetFWHM(fwhm);
    } 
  }

  /// Returns the current smoothness of (all) the scans
  double GetFWHM() const { if (_scans.size()) return(_scans[0].GetFWHM());  else if (_b0scans.size()) return(_b0scans[0].GetFWHM()); else return(0.0); }

  /// Returns the user defined mask (that was passed to the constructor).
  const NEWIMAGE::volume<float>& Mask() const { return(_mask); }

  /// Returns a pointer to susceptibility induced off-resonance field (in Hz). Returns NULL when no field set.
  const boost::shared_ptr<NEWIMAGE::volume<float> > GetSuscHzOffResField() const { return(_topup_field); }

  /// Returns off-resonance field pertaining to EC only. Movements are not included and its main use is for visualiation/demonstration.
  NEWIMAGE::volume<float> GetScanHzECOffResField(unsigned int indx, ScanType st=ANY) const { return(Scan(indx,st).ECField()); }
  /// Separate field offset (image FOV centre not coinciding with scanner iso-centre) from actual subject movement in PE direction.
  void SeparateFieldOffsetFromMovement(HzModel m=ModelMovements);
  /// Set specified scan as reference for location.
  void SetDWIReference(unsigned int ref=0) { set_reference(ref,DWI); }
  /// Set specified scan as reference for location.
  void Setb0Reference(unsigned int ref=0) { set_reference(ref,B0); }

  /// Writes distortion corrected images to disk.
  void WriteRegisteredImages(const std::string& fname, FinalResampling  resmethod, ScanType st=ANY)
  {
    if (resmethod==EDDY::JAC) write_jac_registered_images(fname,st);
    else if (resmethod==EDDY::LSR) write_lsr_registered_images(fname,st);
    else throw EddyException("ECScanManager::WriteRegisteredImages: Unknown resampling method");
  }
  /// Writes file with movement and EC parameters.
  void WriteParameterFile(const std::string& fname, ScanType st=ANY) const;
  /// Writes eddy-current induced fields to disk
  void WriteECFields(const std::string& fname, ScanType st=ANY) const;
  /// Writes original data with outliers replaced by predictions
  void WriteOutlierFreeData(const std::string& fname, ScanType st=ANY) const;
private:
  bool                                          _has_topup;    ///< Is true if object contains a valid susceptibility field
  boost::shared_ptr<NEWIMAGE::volume<float> >   _topup_field;  ///< Safe pointer to susceptibility field (in Hz).
  NEWIMAGE::volume<float>                       _mask;
  double                                        _sf;           ///< Scale factor applied to scans
  std::vector<pair<int,int> >                   _fi;           ///< Used to keep track of index into file
  std::vector<ECScan>                           _scans;        ///< Vector of diffusion weighted scans
  std::vector<ECScan>                           _b0scans;      ///< Vector of b=0 scans.
  unsigned int                                  _nsess;        ///< Numbes of sessions in which data was collected

  /// Extracts everything that can potentially be a field offset scaled in Hz.
  NEWMAT::ColumnVector hz_vector_with_everything() const;
  /// Creates design matrix explaining "difficult movements" in terms of other movements.
  NEWMAT::Matrix movement_design_matrix(bool debug=false) const;
  /// Creates design matrix explaining "difficult movements" as residuals after modelling linear EC.
  NEWMAT::Matrix offset_design_matrix(bool debug=false) const;
  /// Returns a vector of movement parameters for the b=0 scans inter/extrapolated onto the dwis.
  NEWMAT::Matrix get_b0_movement_vector(ScanType st=DWI) const;
  /// Sets scan indicated by ref as reference (position zero) for all scans of type st.
  void set_reference(unsigned int ref, ScanType st);
  double mean_of_first_b0(const NEWIMAGE::volume4D<float>&   vols,
                          const NEWIMAGE::volume<float>&     mask,
                          const NEWMAT::Matrix&              bvecs,
                          const NEWMAT::Matrix&              bvals) const;
  bool index_kosher(unsigned int indx, ScanType st) const
  {
    if (st==DWI) return(indx<_scans.size());
    else if (st==B0) return(indx<_b0scans.size());
    else return(indx<_fi.size());
  }
  void write_jac_registered_images(const std::string& fname, ScanType st);
  void write_lsr_registered_images(const std::string& fname, ScanType st);
  bool has_pe_in_direction(unsigned int dir, ScanType st=ANY) const;
};

} // End namespace EDDY

#endif // End #ifndef ECScanClasses_h

////////////////////////////////////////////////
//
// Here starts Doxygen documentation
//
////////////////////////////////////////////////

/*!
 * \var ECScan::_ecp
 * This is an object determining the model that is being used to model the eddy currents. 
 * It is stored as a (safe) pointer to a virtual base class ScanECModel. 
 * In any instantiation of an ECScan object it will contain an object of a class derived 
 * from ScanECModel, for example LinearScanECModel. This means that ECScan do not need to 
 * know which model is actually being used since it will only interact with _ecp through 
 * the interface specified by ScanECModel.
 */

/*!
 * \var ECScanManager::_fi
 * This is an array of pairs keeping track of the relationship between on the one hand
 * DWI and b0 indexing and on the other hand global indexing (corresponding to the indexing
 * on disk). So if we for example have a diffusion weighted volume that is the i'th volume
 * in the 4D file it was read from and that also is the j'th diffusion weighted volume then:
 * _fi[i-1].first==0 (to indicate that it is diffusion weighted) and _fi[i-1].second==j-1 to
 * indicate it is the j'th diffusion weighted volume. Correspondingly a b=0 volume that is
 * the l'th volume on disc and the m'th b=0 volume then: _fi[l-1].first==1 and 
 * _fi[l-1].second==m-1.
 */

