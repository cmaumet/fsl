// Declarations of classes that implements useful
// utility functions for the eddy current project.
// They are collections of statically declared
// functions that have been collected into classes 
// to make it explicit where they come from. There
// will never be any instances of theses classes.
// 
// EddyUtils.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef EddyUtils_h
#define EddyUtils_h

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "EddyHelperClasses.h"
#include "ECScanClasses.h"

namespace EDDY {

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class EddyUtils
//
// Helper Class used to perform various useful tasks for 
// the eddy current correction project.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class EddyUtils
{
private:
  /// b-values within this range considered equivalent
  static const int b_range = 100; 
  /// bladibla
  static NEWIMAGE::volume4D<float> get_partial_derivatives_in_scan_space(// Input
									 const NEWIMAGE::volume<float>&                      pred,
									 const EDDY::ECScan&                                 scan,
									 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
                                                                         EDDY::Parameters                                    whichp);

  static double param_update(// Input
			     const NEWIMAGE::volume<float>&                      pred,      // Prediction in model space
			     boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,      // Susceptibility induced off-resonance field
			     const NEWIMAGE::volume<float>&                      pmask,     
			     EDDY::Parameters                                    whichp,    // What parameters do we want to update
			     bool                                                cbs,       // Check (success of parameters) Before Set
			     // Input/output
			     EDDY::ECScan&                                       scan,      // Scan we want to register to pred
			     // Output
			     NEWMAT::ColumnVector                                *update);  // Vector of updates, optional output

  static EDDY::ImageCoordinates transform_coordinates_from_model_to_scan_space(// Input
									       const NEWIMAGE::volume<float>&                      pred,
									       const EDDY::ECScan&                                 scan,
									       boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
									       // Output
									       NEWIMAGE::volume<float>                             *omask,
									       NEWIMAGE::volume<float>                             *jac);
  // Returns coordinates into f transformed with 
  // displacement field d and affine matrix M
  static void transform_coordinates(// Input 
				    const NEWIMAGE::volume<float>&    f,
				    const NEWIMAGE::volume4D<float>&  d,
				    const NEWMAT::Matrix&             M,
				    // Output
				    ImageCoordinates&                 c,
				    NEWIMAGE::volume<float>           *omask);

  // Calculates X.t()*X where X is a matrix where each column is one of the volumes in vols
  static NEWMAT::Matrix make_XtX(const NEWIMAGE::volume4D<float>& vols,
				 const NEWIMAGE::volume<float>&   mask);

  // Calculates X.t()*y where X is a matrix where each column is one of the volumes in Xvols
  // and where y is the volume in Yvol.
  static NEWMAT::ColumnVector make_Xty(const NEWIMAGE::volume4D<float>& Xvols,
				       const NEWIMAGE::volume<float>&   Yvol,
				       const NEWIMAGE::volume<float>&   mask);

public:
  // This function has been temporarily moved into public space. Should probably be
  // moved back to private space at some stage.
  static NEWIMAGE::volume<float> transform_model_to_scan_space(// Input
							       const NEWIMAGE::volume<float>&                      pred,
							       const EDDY::ECScan&                                 scan,
							       boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
							       bool                                                jacmod,
							       // Output
							       NEWIMAGE::volume<float>&                            omask,
							       NEWIMAGE::volume<float>                             *jac,
							       NEWIMAGE::volume4D<float>                           *grad);


  // Some functions for comparing diffusion parameters
  /// Returns true if the difference in b-value is less than EddyUtils::b_range
  static bool AreInSameShell(const DiffPara& dp1,
                             const DiffPara& dp2) { return(fabs(dp1.bVal()-dp2.bVal())<double(b_range)); }
  static bool IsDiffusionWeighted(const DiffPara& dp) { return(dp.bVal() > double(b_range)); }
  static bool Isb0(const DiffPara& dp) { return(!IsDiffusionWeighted(dp)); }
  /// Returns true if the inner product of the b-vectors is greater than 0.999
  static bool HaveSameDirection(const DiffPara& dp1,
				const DiffPara& dp2) { return(NEWMAT::DotProduct(dp1.bVec(),dp2.bVec())>0.999); }

  // Random functions to set extrapolation and interpolation //
  template <class V>
  static void SetTrilinearInterp(V& vol) {
    if (vol.getinterpolationmethod() != NEWIMAGE::trilinear) vol.setinterpolationmethod(NEWIMAGE::trilinear);
    if (vol.getextrapolationmethod() != NEWIMAGE::mirror) vol.setextrapolationmethod(NEWIMAGE::mirror);
  }
  template <class V>
  static void SetSplineInterp(V& vol) {
    if (vol.getinterpolationmethod() != NEWIMAGE::spline) vol.setinterpolationmethod(NEWIMAGE::spline);
    if (vol.getsplineorder() != 3) vol.setsplineorder(3);
    if (vol.getextrapolationmethod() != NEWIMAGE::mirror) vol.setextrapolationmethod(NEWIMAGE::mirror);
  }

  // Check if a pair of ECScans can potentially be used in an LSR reconstruction
  static bool AreMatchingPair(const ECScan& s1, const ECScan& s2);
  
  // Get indicies for non-zero b-values
  static std::vector<unsigned int> GetIndiciesOfDWIs(const std::vector<DiffPara>& dpars);

  // Removes bvecs associated with zero b-values.
  static std::vector<DiffPara> GetDWIDiffParas(const std::vector<DiffPara>&   dpars);

  // Reads all diffusion weighted images from 4D volume
  static int read_DWI_volume4D(NEWIMAGE::volume4D<float>&     dwivols,
			       const std::string&             fname,
			       const std::vector<DiffPara>&   dpars);

  // Converts char mask (from the general_transform functions) to a float mask
  static NEWIMAGE::volume<float> ConvertMaskToFloat(const NEWIMAGE::volume<char>& charmask);

  // Calculates slice-wise statistics from the difference between observation and predicton in observation space
  static DiffStats GetSliceWiseStats(// Input
				     const NEWIMAGE::volume<float>&                    pred,      // Prediction in model space
				     boost::shared_ptr<const NEWIMAGE::volume<float> > susc,      // Susceptibility induced off-resonance field
				     const NEWIMAGE::volume<float>&                    pmask,     // "Data valid" mask in model space
				     const NEWIMAGE::volume<float>&                    bmask,     // Brain mask in model space
				     EDDY::ECScan&                                     scan);     // Scan corresponding to pred

  // Performs an update of the movement parameters for one scan
  static double MovParamUpdate(// Input
			       const NEWIMAGE::volume<float>&                      pred,
			       boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
			       const NEWIMAGE::volume<float>&                      pmask,
			       bool                                                cbs,
			       // Input/output
			       EDDY::ECScan&                                       scan) {
    return(param_update(pred,susc,pmask,MOVEMENT,cbs,scan,NULL));
  }

  // Performs an update of the EC parameters for one scan
  static double ECParamUpdate(// Input
			      const NEWIMAGE::volume<float>&                      pred,
			      boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
			      const NEWIMAGE::volume<float>&                      pmask,
			      bool                                                cbs,
			      // Input/output
			      EDDY::ECScan&                                       scan) {
    return(param_update(pred,susc,pmask,EC,cbs,scan,NULL));
  }

  // Performs an update of the EC parameters for one scan
  static double MovAndECParamUpdate(// Input
				    const NEWIMAGE::volume<float>&                      pred,
				    boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
				    const NEWIMAGE::volume<float>&                      pmask,
				    bool                                                cbs,
				    // Input/output
				    EDDY::ECScan&                                       scan) {
    return(param_update(pred,susc,pmask,ALL,cbs,scan,NULL));
  }

  // Transforms an image from model/prediction space to observation space
  static NEWIMAGE::volume<float> TransformModelToScanSpace(// Input
							   const NEWIMAGE::volume<float>&                      pred,
							   const EDDY::ECScan&                                 scan,
							   boost::shared_ptr<const NEWIMAGE::volume<float> >   susc) {
    NEWIMAGE::volume<float> mask(pred.xsize(),pred.ysize(),pred.zsize()); 
    NEWIMAGE::volume<float> jac(pred.xsize(),pred.ysize(),pred.zsize()); 
    return(transform_model_to_scan_space(pred,scan,susc,true,mask,&jac,NULL));
  }
  static NEWIMAGE::volume<float> TransformScanToModelSpace(// Input
							   const EDDY::ECScan&                               scan,
							   boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							   // Output
							   NEWIMAGE::volume<float>&                          omask);

  // The next two are alternate transformation routines that 
  // performs the transforms in several resampling steps.
  // They are intended for debugging.
  static NEWIMAGE::volume<float> DirectTransformScanToModelSpace(// Input
								 const EDDY::ECScan&                               scan,
								 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								 // Output
								 NEWIMAGE::volume<float>&                          omask);

  static NEWIMAGE::volume<float> DirectTransformModelToScanSpace(// Input
								 const NEWIMAGE::volume<float>&                    ima,
								 const EDDY::ECScan&                               scan,
								 boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								 // Output
								 NEWIMAGE::volume<float>&                          omask);

  // Alternate version of updating routine for writing out debug information.

  static double param_update_debug(// Input
				   const NEWIMAGE::volume<float>&                      pred,      // Prediction in model space
				   boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,      // Susceptibility induced off-resonance field
				   const NEWIMAGE::volume<float>&                      pmask,     // 
				   Parameters                                          whichp,    // What parameters do we want to update
				   bool                                                cbs,       // Check (success of parameters) Before Set
				   unsigned int                                        scindx,    // Scan index
				   unsigned int                                        iter,      // Iteration
				   unsigned int                                        level,     // Determines how much gets written
				   // Input/output
				   EDDY::ECScan&                                       scan,      // Scan we want to register to pred
				   // Output
				   NEWMAT::ColumnVector                                *rupdate); // Vector of updates, optional output

};

// }}} End of fold.

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class FieldUtils
//
// Helper Class used to perform various useful calculations
// on off-resonance and displacement fields.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class FieldUtils
{
public:
  // Rigid body transform of off-resonance field
  static NEWIMAGE::volume<float> RigidBodyTransformHzField(const NEWIMAGE::volume<float>& hzfield);

  // Some conversion routines off-resonance->displacements
  static NEWIMAGE::volume4D<float> Hz2VoxelDisplacements(const NEWIMAGE::volume<float>& hzfield,
                                                         const AcqPara&                 acqp);
  static NEWIMAGE::volume4D<float> Hz2MMDisplacements(const NEWIMAGE::volume<float>& hzfield,
                                                      const AcqPara&                 acqp);
  static NEWIMAGE::volume4D<float> Voxel2MMDisplacements(const NEWIMAGE::volume4D<float>& voxdisp) { 
    NEWIMAGE::volume4D<float> mmd=voxdisp; mmd[0] *= mmd.xdim(); mmd[1] *= mmd.ydim(); mmd[2] *= mmd.zdim(); return(mmd);
  }

  // Inverts 3D displacement field, ASSUMING it is really 1D (only non-zero displacements in one direction).
  static NEWIMAGE::volume4D<float> InvertDisplacementField(// Input
							   const NEWIMAGE::volume4D<float>& dfield,
							   const AcqPara&                   acqp,
							   const NEWIMAGE::volume<float>& inmask,
							   // Output
							   NEWIMAGE::volume<float>&       omask);

  // Inverts 1D displacement field. Input must be scaled to voxels (i.e. not mm).
  static NEWIMAGE::volume<float> InvertDisplacementField(// Input 
							 const NEWIMAGE::volume<float>& dfield,
							 const AcqPara&                 acqp,
							 const NEWIMAGE::volume<float>& inmask,
							 // Output
							 NEWIMAGE::volume<float>&       omask);

  // Calculates Jacobian of a displacement field
  static NEWIMAGE::volume<float> GetJacobian(const NEWIMAGE::volume4D<float>& dfield,
                                             const AcqPara&                   acqp);

  // Calculates Jacobian of a 1D displacement field
  static NEWIMAGE::volume<float> GetJacobianFrom1DField(const NEWIMAGE::volume<float>& dfield,
							unsigned int                   dir);
private:
};

} // End namespace EDDY

#endif // End #ifndef EddyUtils_h
