// Declarations of classes that implements useful
// utility functions for the eddy current project.
// They are collections of statically declared
// functions that have been collected into classes 
// to make it explicit where they come from. There
// will never be any instances of theses classes.
// 
// EddyUtils.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#include <cstdlib>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>
#include "newmat.h"
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .set_sform etc
#endif
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/warpfns.h"
#include "EddyHelperClasses.h"
#include "EddyUtils.h"

using namespace EDDY;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class EddyUtils
//
// Helper Class used to perform various useful tasks for 
// the eddy current correction project.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

std::vector<unsigned int> EddyUtils::GetIndiciesOfDWIs(const std::vector<DiffPara>& dpars)
{
  std::vector<unsigned int> indicies;
  for (unsigned int i=0; i<dpars.size(); i++) { if (EddyUtils::IsDiffusionWeighted(dpars[i])) indicies.push_back(i); }
  return(indicies);
}

std::vector<DiffPara> EddyUtils::GetDWIDiffParas(const std::vector<DiffPara>&   dpars)
{
  std::vector<unsigned int> indx = EddyUtils::GetIndiciesOfDWIs(dpars);
  std::vector<DiffPara> dwi_dpars;
  for (unsigned int i=0; i<indx.size(); i++) dwi_dpars.push_back(dpars[indx[i]]);
  return(dwi_dpars);
}

bool EddyUtils::AreMatchingPair(const ECScan& s1, const ECScan& s2)
{
  double dp = NEWMAT::DotProduct(s1.GetAcqPara().PhaseEncodeVector(),s1.GetAcqPara().PhaseEncodeVector());
  if (std::abs(dp + 1.0) > 1e-6) return(false);
  if (!EddyUtils::AreInSameShell(s1.GetDiffPara(),s2.GetDiffPara())) return(false);
  if (!HaveSameDirection(s1.GetDiffPara(),s2.GetDiffPara())) return(false);
  return(true);
}

int EddyUtils::read_DWI_volume4D(NEWIMAGE::volume4D<float>&     dwivols,
				 const std::string&             fname,
				 const std::vector<DiffPara>&   dpars)
{
  std::vector<unsigned int> indx = EddyUtils::GetIndiciesOfDWIs(dpars);
  NEWIMAGE::volume4D<float> tmp;
  read_volume4DROI(tmp,fname,0,0,0,0,-1,-1,-1,0);
  dwivols.reinitialize(tmp.xsize(),tmp.ysize(),tmp.zsize(),indx.size());
  for (unsigned int i=0; i<indx.size(); i++) {
    read_volume4DROI(tmp,fname,0,0,0,indx[i],-1,-1,-1,indx[i]);
    dwivols[i] = tmp[0];
  }
  return(1);
}

NEWIMAGE::volume<float> EddyUtils::ConvertMaskToFloat(const NEWIMAGE::volume<char>& charmask)
{
  NEWIMAGE::volume<float> floatmask(charmask.xsize(),charmask.ysize(),charmask.zsize());
  NEWIMAGE::copybasicproperties(charmask,floatmask);
  for (int k=0; k<charmask.zsize(); k++) {
    for (int j=0; j<charmask.ysize(); j++) {
      for (int i=0; i<charmask.xsize(); i++) {
	floatmask(i,j,k) = static_cast<float>(charmask(i,j,k));
      }
    }
  }
  return(floatmask);
}

DiffStats EddyUtils::GetSliceWiseStats(// Input
				       const NEWIMAGE::volume<float>&                    pred,      // Prediction in model space
				       boost::shared_ptr<const NEWIMAGE::volume<float> > susc,      // Susceptibility induced off-resonance field
				       const NEWIMAGE::volume<float>&                    pmask,     // "Data valid" mask in model space
				       const NEWIMAGE::volume<float>&                    bmask,     // Brain mask in model space
				       EDDY::ECScan&                                     scan)      // Scan we want to register to pred
{
  // Transform prediction into observation space
  NEWIMAGE::volume<float> pios = EddyUtils::TransformModelToScanSpace(pred,scan,susc);
  // Transform binary mask into observation space
  NEWIMAGE::volume<float> mask = pred; mask = 0.0;
  NEWIMAGE::volume<float> bios = EddyUtils::transform_model_to_scan_space(pmask*bmask,scan,susc,false,mask,NULL,NULL);
  bios.binarise(0.9); // Value above (arbitrary) 0.9 implies valid voxels
  mask *= bios; // Volume and prediction mask falls within FOV
  // Calculate slice-wise stats from difference image
  DiffStats stats(scan.GetOriginalIma()-pios,mask);
  return(stats);
}

double EddyUtils::param_update_debug(// Input
				     const NEWIMAGE::volume<float>&                      pred,      // Prediction in model space
				     boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,      // Susceptibility induced off-resonance field
				     const NEWIMAGE::volume<float>&                      pmask,     // "Data valid" mask in model space
				     Parameters                                          whichp,    // What parameters do we want to update
				     bool                                                cbs,       // Check (success of parameters) Before Set
				     unsigned int                                        scindx,    // Scan index
				     unsigned int                                        iter,      // Iteration
				     unsigned int                                        level,     // Determines how much gets written
				     // Input/output
				     EDDY::ECScan&                                       scan,      // Scan we want to register to pred
				     // Output
				     NEWMAT::ColumnVector                                *rupdate)  // Vector of updates, optional output
{
  // Transform prediction into observation space
  NEWIMAGE::volume<float> pios = EddyUtils::TransformModelToScanSpace(pred,scan,susc);
  // Transform binary mask into observation space
  NEWIMAGE::volume<float> mask1 = pred; mask1 = 0.0;
  NEWIMAGE::volume<float> bios = EddyUtils::transform_model_to_scan_space(pmask,scan,susc,false,mask1,NULL,NULL);
  NEWIMAGE::volume<float> bbios = bios; bbios.binarise(0.99);
  // bios.binarise(0.99); // Value above (arbitrary) 0.99 implies valid voxels
  NEWIMAGE::volume<float> mask = mask1;
  mask *= bbios; // Volume and prediction mask falls within FOV
  // Get partial derivatives w.r.t. to requested category of parameters in prediction space
  NEWIMAGE::volume4D<float> derivs = EddyUtils::get_partial_derivatives_in_scan_space(pred,scan,susc,whichp);
  // Calculate XtX where X is a matrix whos columns are the partial derivatives
  NEWMAT::Matrix XtX = EddyUtils::make_XtX(derivs,mask);
  // Calculate difference image between observed and predicted
  NEWIMAGE::volume<float> dima = pios-scan.GetIma();
  // Calculate Xty where y is the difference between observed and predicted. X as above.
  NEWMAT::ColumnVector Xty = EddyUtils::make_Xty(derivs,dima,mask);
  // Calculate mean sum of squares from difference image
  double mss = (dima*mask).sumsquares() / mask.sum();
  // Calculate update to parameters
  NEWMAT::ColumnVector update = -XtX.i()*Xty;
  // Get scan in model space before setting new parameters
  NEWIMAGE::volume<float> sims=scan.GetUnwarpedIma(susc);  // Scan In Model Space
  // Update parameters
  for (unsigned int i=0; i<scan.NDerivs(whichp); i++) {
    scan.SetDerivParam(i,scan.GetDerivParam(i,whichp)+update(i+1),whichp);
  }
  NEWIMAGE::volume<float> new_pios;
  NEWIMAGE::volume<float> new_bios;
  NEWIMAGE::volume<float> new_mask;
  if (cbs) {
    new_pios = EddyUtils::TransformModelToScanSpace(pred,scan,susc);
    // Transform binary mask into observation space
    new_mask = new_pios; new_mask = 0.0;
    new_bios = EddyUtils::transform_model_to_scan_space(pmask,scan,susc,false,new_mask,NULL,NULL);
    new_bios.binarise(0.99); // Value above (arbitrary) 0.99 implies valid voxels
    new_mask *= new_bios; // Volume and prediction mask falls within FOV
    double mss_au = ((new_pios-scan.GetIma()).sumsquares()) / new_mask.sum();
    if (mss_au > mss) { // Oh dear
      for (unsigned int i=0; i<scan.NDerivs(whichp); i++) {
	scan.SetDerivParam(i,scan.GetDerivParam(i,whichp)-update(i+1),whichp);
      }
    }
  }
  char fname[256];
  if (level>0) {
    sprintf(fname,"EDDY_DEBUG_masked_dima_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(dima*mask,fname);
    sprintf(fname,"EDDY_DEBUG_reverse_dima_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(pred-sims,fname);
  }
  if (level>1) {
    sprintf(fname,"EDDY_DEBUG_mask_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(mask,fname);
    sprintf(fname,"EDDY_DEBUG_pios_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(pios,fname);
    sprintf(fname,"EDDY_DEBUG_pred_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(pred,fname);
    sprintf(fname,"EDDY_DEBUG_dima_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(dima,fname);
    sprintf(fname,"EDDY_DEBUG_orig_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(scan.GetIma(),fname);
    if (cbs) {
      sprintf(fname,"EDDY_DEBUG_new_masked_dima_%02d_%04d",iter,scindx);
      NEWIMAGE::write_volume(new_pios-scan.GetIma()*new_mask,fname);
      sims = scan.GetUnwarpedIma(susc);
      sprintf(fname,"EDDY_DEBUG_new_reverse_dima_%02d_%04d",iter,scindx);
      NEWIMAGE::write_volume(pred-sims,fname);
      sprintf(fname,"EDDY_DEBUG_new_mask_%02d_%04d",iter,scindx);
      NEWIMAGE::write_volume(new_mask,fname);
      sprintf(fname,"EDDY_DEBUG_new_pios_%02d_%04d",iter,scindx);
      NEWIMAGE::write_volume(new_pios,fname);
      sprintf(fname,"EDDY_DEBUG_new_dima_%02d_%04d",iter,scindx);
      NEWIMAGE::write_volume(new_pios-scan.GetIma()*new_mask,fname);
    }
  }
  if (level>2) {
    sprintf(fname,"EDDY_DEBUG_mask1_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(mask1,fname);
    sprintf(fname,"EDDY_DEBUG_bios_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(bios,fname);
    sprintf(fname,"EDDY_DEBUG_bbios_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(bbios,fname);
    sprintf(fname,"EDDY_DEBUG_pmask_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume(pmask,fname);
    sprintf(fname,"EDDY_DEBUG_derivs_%02d_%04d",iter,scindx);
    NEWIMAGE::write_volume4D(derivs,fname);
    sprintf(fname,"EDDY_DEBUG_XtX_%02d_%04d",iter,scindx);
    MISCMATHS::write_ascii_matrix(fname,XtX);
    sprintf(fname,"EDDY_DEBUG_Xty_%02d_%04d",iter,scindx);
    MISCMATHS::write_ascii_matrix(fname,Xty);
    sprintf(fname,"EDDY_DEBUG_update_%02d_%04d",iter,scindx);
    MISCMATHS::write_ascii_matrix(fname,update);
  }
  if (rupdate) *rupdate = update;
  return(mss);
}

double EddyUtils::param_update(// Input
			       const NEWIMAGE::volume<float>&                      pred,      // Prediction in model space
			       boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,      // Susceptibility induced off-resonance field
			       const NEWIMAGE::volume<float>&                      pmask,     // 
			       Parameters                                          whichp,    // What parameters do we want to update
			       bool                                                cbs,       // Check (success of parameters) Before Set
			       // Input/output
			       EDDY::ECScan&                                       scan,      // Scan we want to register to pred
			       // Output
			       NEWMAT::ColumnVector                                *rupdate)  // Vector of updates, optional output
{
  // Transform prediction into observation space
  NEWIMAGE::volume<float> pios = EddyUtils::TransformModelToScanSpace(pred,scan,susc);
  // Transform binary mask into observation space
  NEWIMAGE::volume<float> mask = pred; mask = 0.0;
  NEWIMAGE::volume<float> bios = EddyUtils::transform_model_to_scan_space(pmask,scan,susc,false,mask,NULL,NULL);
  bios.binarise(0.99); // Value above (arbitrary) 0.99 implies valid voxels
  mask *= bios; // Volume and prediction mask falls within FOV
  // Get partial derivatives w.r.t. to requested category of parameters in prediction space
  NEWIMAGE::volume4D<float> derivs = EddyUtils::get_partial_derivatives_in_scan_space(pred,scan,susc,whichp);
  // Calculate XtX where X is a matrix whos columns are the partial derivatives
  NEWMAT::Matrix XtX = EddyUtils::make_XtX(derivs,mask);
  // Calculate difference image between observed and predicted
  NEWIMAGE::volume<float> dima = pios-scan.GetIma();
  // Calculate Xty where y is the difference between observed and predicted. X as above.
  NEWMAT::ColumnVector Xty = EddyUtils::make_Xty(derivs,dima,mask);
  // Calculate mean sum of squares from difference image
  double mss = (dima*mask).sumsquares() / mask.sum();
  // Calculate update to parameters
  NEWMAT::ColumnVector update = -XtX.i()*Xty;
  // Update parameters
  for (unsigned int i=0; i<scan.NDerivs(whichp); i++) {
    scan.SetDerivParam(i,scan.GetDerivParam(i,whichp)+update(i+1),whichp);
  }
  if (cbs) {
    pios = EddyUtils::TransformModelToScanSpace(pred,scan,susc);
    // Transform binary mask into observation space
    mask = 0.0;
    bios = EddyUtils::transform_model_to_scan_space(pmask,scan,susc,false,mask,NULL,NULL);
    bios.binarise(0.99); // Value above (arbitrary) 0.99 implies valid voxels
    mask *= bios; // Volume and prediction mask falls within FOV
    double mss_au = (((pios-scan.GetIma())*mask).sumsquares()) / mask.sum();
    if (mss_au > mss) { // Oh dear
      for (unsigned int i=0; i<scan.NDerivs(whichp); i++) {
	scan.SetDerivParam(i,scan.GetDerivParam(i,whichp)-update(i+1),whichp);
      }
    }
  }
  if (rupdate) *rupdate = update;
  return(mss);
}

NEWIMAGE::volume<float> EddyUtils::transform_model_to_scan_space(// Input
								 const NEWIMAGE::volume<float>&                      pred,
								 const EDDY::ECScan&                                 scan,
								 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
								 bool                                                jacmod,
								 // Output
								 NEWIMAGE::volume<float>&                            omask,
								 NEWIMAGE::volume<float>                             *jac,
								 NEWIMAGE::volume4D<float>                           *grad)
{
  // Get total field from scan
  if (jacmod && !jac) throw EddyException("EddyUtils::transform_model_to_scan_space: jacmod can only be used with valid jac");
  NEWIMAGE::volume4D<float> dfield;
  if (jacmod || jac) dfield = scan.FieldForModelToScanTransform(susc,omask,*jac);
  else dfield = scan.FieldForModelToScanTransform(susc,omask);
  // Get RB matrix and EC field
  NEWMAT::Matrix R = scan.ForwardMovementMatrix();
  // Transform prediction using RB, inverted Tot map and Jacobian
  NEWMAT::Matrix eye(4,4); eye=0; eye(1,1)=1.0; eye(2,2)=1.0; eye(3,3)=1.0; eye(4,4)=1.0;
  NEWIMAGE::volume<float> ovol = pred; ovol = 0.0;
  NEWIMAGE::volume<char> mask3(pred.xsize(),pred.ysize(),pred.zsize());
  NEWIMAGE::copybasicproperties(pred,mask3); mask3 = 0;
  if (grad) {
    std::vector<int> ddir(3); ddir[0] = 0; ddir[1] = 1; ddir[2] = 2;
    grad->reinitialize(pred.xsize(),pred.ysize(),pred.zsize(),3); 
    NEWIMAGE::copybasicproperties(pred,*grad);
    NEWIMAGE::raw_general_transform(pred,eye,dfield,ddir,ddir,&eye,&R,ovol,*grad,&mask3);
  }
  else NEWIMAGE::apply_warp(pred,eye,dfield,eye,R,ovol,mask3);
  omask *= EddyUtils::ConvertMaskToFloat(mask3); // Combine all masks
  EddyUtils::SetTrilinearInterp(omask);
  if (jacmod) ovol *= *jac;                      // Jacobian modulation if it was asked for
  return(ovol);
}

EDDY::ImageCoordinates EddyUtils::transform_coordinates_from_model_to_scan_space(// Input
										 const NEWIMAGE::volume<float>&                      pred,
										 const EDDY::ECScan&                                 scan,
										 boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
										 // Output
										 NEWIMAGE::volume<float>                             *omask,
										 NEWIMAGE::volume<float>                             *jac)
{
  // Get total field from scan
  NEWIMAGE::volume4D<float> dfield;
  if (omask && jac) dfield = scan.FieldForModelToScanTransform(susc,*omask,*jac);
  else if (omask) dfield = scan.FieldForModelToScanTransform(susc,*omask);
  else if (jac) dfield = scan.FieldForModelToScanTransformWithJac(susc,*jac);
  else dfield = scan.FieldForModelToScanTransform(susc);
  // Get RB matrix and EC field
  NEWMAT::Matrix R = scan.ForwardMovementMatrix();
  // Transform prediction using RB, inverted Tot map and Jacobian
  NEWMAT::IdentityMatrix eye(4);
  NEWIMAGE::volume<float> ovol = pred; ovol = 0.0;
  ImageCoordinates coord(pred.xsize(),pred.ysize(),pred.zsize()); 
  if (omask) { 
    NEWIMAGE::volume<float> mask2(pred.xsize(),pred.ysize(),pred.zsize());
    NEWIMAGE::copybasicproperties(pred,mask2); mask2 = 0;
    EddyUtils::transform_coordinates(pred,dfield,R,coord,&mask2);
    *omask *= mask2;
    EddyUtils::SetTrilinearInterp(*omask);
  }
  else EddyUtils::transform_coordinates(pred,dfield,R,coord,NULL);

  return(coord);
}


NEWIMAGE::volume4D<float> EddyUtils::get_partial_derivatives_in_scan_space(// Input
									   const NEWIMAGE::volume<float>&                      pred,
									   const EDDY::ECScan&                                 scan,
									   boost::shared_ptr<const NEWIMAGE::volume<float> >   susc,
									   EDDY::Parameters                                    whichp)
{
  NEWIMAGE::volume<float> basejac;
  NEWIMAGE::volume4D<float> grad;
  NEWIMAGE::volume<float> skrutt(pred.xsize(),pred.ysize(),pred.zsize());
  NEWIMAGE::volume<float> base = transform_model_to_scan_space(pred,scan,susc,true,skrutt,&basejac,&grad);
  ImageCoordinates basecoord = transform_coordinates_from_model_to_scan_space(pred,scan,susc,NULL,NULL);
  NEWIMAGE::volume4D<float> derivs(base.xsize(),base.ysize(),base.zsize(),scan.NDerivs(whichp));
  NEWIMAGE::volume<float> jac = pred;
  ECScan sc = scan;
  for (unsigned int i=0; i<sc.NDerivs(whichp); i++) {
    double p = sc.GetDerivParam(i,whichp);
    sc.SetDerivParam(i,p+sc.GetDerivScale(i,whichp),whichp);
    ImageCoordinates diff = transform_coordinates_from_model_to_scan_space(pred,sc,susc,NULL,&jac) - basecoord;
    derivs[i] = (diff*grad) / sc.GetDerivScale(i,whichp);
    derivs[i] += base * (jac-basejac) / sc.GetDerivScale(i,whichp);
    sc.SetDerivParam(i,p,whichp);
  }
  return(derivs);
}

/*
NEWIMAGE::volume<float> EddyUtils::TransformScanToModelSpace(// Input
							     const EDDY::ECScan&                               scan,
							     boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
							     // Output
							     NEWIMAGE::volume<float>&                          omask)
{
  // Get total field from scan
  NEWIMAGE::volume<float> jac;
  NEWIMAGE::volume4D<float> dfield = scan.FieldForScanToModelTransform(susc,omask,jac);
  // Transform prediction using inverse RB, dfield and Jacobian
  NEWMAT::Matrix iR = scan.InverseMovementMatrix();
  NEWIMAGE::volume<float> ovol = scan.GetIma(); ovol = 0.0;
  NEWIMAGE::volume<char> mask2(ovol.xsize(),ovol.ysize(),ovol.zsize()); 
  NEWIMAGE::copybasicproperties(scan.GetIma(),mask2); mask2 = 1;
  NEWIMAGE::general_transform(scan.GetIma(),iR,dfield,ovol,mask2);
  // Combine all masks
  omask *= EddyUtils::ConvertMaskToFloat(mask2);
  return(jac*ovol);
}
*/

NEWIMAGE::volume<float> EddyUtils::DirectTransformScanToModelSpace(// Input
								   const EDDY::ECScan&                               scan,
								   boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								   // Output
								   NEWIMAGE::volume<float>&                          omask)
{
  NEWIMAGE::volume<float> ima = scan.GetIma();
  NEWIMAGE::volume<float> eb = scan.ECField();
  NEWIMAGE::volume4D<float> dfield = FieldUtils::Hz2VoxelDisplacements(eb,scan.GetAcqPara());
  dfield = FieldUtils::Voxel2MMDisplacements(dfield);
  NEWMAT::Matrix eye(4,4); eye=0; eye(1,1)=1.0; eye(2,2)=1.0; eye(3,3)=1.0; eye(4,4)=1.0;
  NEWIMAGE::volume<float> ovol = ima; ovol = 0.0;
  NEWIMAGE::volume<char> mask(ima.xsize(),ima.ysize(),ima.zsize());
  NEWIMAGE::apply_warp(ima,eye,dfield,eye,eye,ovol,mask);

  NEWMAT::Matrix iR = scan.InverseMovementMatrix();
  NEWIMAGE::volume<float> tmp = ovol; ovol = 0;
  NEWIMAGE::affine_transform(tmp,iR,ovol,mask);
  
  return(ovol);
}

NEWIMAGE::volume<float> EddyUtils::DirectTransformModelToScanSpace(// Input
                                                                   const NEWIMAGE::volume<float>&                    ima,
								   const EDDY::ECScan&                               scan,
								   boost::shared_ptr<const NEWIMAGE::volume<float> > susc,
								   // Output
								   NEWIMAGE::volume<float>&                          omask)
{
  NEWMAT::Matrix R = scan.ForwardMovementMatrix();
  NEWIMAGE::volume<float> tmp = ima; tmp = 0;
  NEWIMAGE::volume<char> mask(tmp.xsize(),tmp.ysize(),tmp.zsize());
  NEWIMAGE::affine_transform(ima,R,tmp,mask);

  NEWIMAGE::volume<float> eb = scan.ECField();
  NEWIMAGE::volume4D<float> dfield = FieldUtils::Hz2VoxelDisplacements(eb,scan.GetAcqPara());
  NEWIMAGE::volume4D<float> idfield = FieldUtils::InvertDisplacementField(dfield,scan.GetAcqPara(),EddyUtils::ConvertMaskToFloat(mask),omask);
  idfield = FieldUtils::Voxel2MMDisplacements(idfield);
  NEWMAT::Matrix eye(4,4); eye=0; eye(1,1)=1.0; eye(2,2)=1.0; eye(3,3)=1.0; eye(4,4)=1.0;
  NEWIMAGE::volume<float> ovol = tmp; ovol = 0.0;
  NEWIMAGE::apply_warp(tmp,eye,idfield,eye,eye,ovol,mask);
  
  return(ovol);
}

NEWMAT::Matrix EddyUtils::make_XtX(const NEWIMAGE::volume4D<float>& vols,
				   const NEWIMAGE::volume<float>&   mask)
{
  NEWMAT::Matrix XtX(vols.tsize(),vols.tsize());
  XtX = 0.0;
  for (int r=1; r<=vols.tsize(); r++) {
    for (int c=r; c<=vols.tsize(); c++) {
      for (NEWIMAGE::volume<float>::fast_const_iterator rit=vols[r-1].fbegin(), ritend=vols[r-1].fend(), cit=vols[c-1].fbegin(), mit=mask.fbegin(); rit!=ritend; ++rit, ++cit, ++mit) {
	if (*mit) XtX(r,c) += (*rit)*(*cit);
      }
    }
  }
  for (int r=2; r<=vols.tsize(); r++) {
    for (int c=1; c<r; c++) XtX(r,c) = XtX(c,r);
  }
  return(XtX);
}

NEWMAT::ColumnVector EddyUtils::make_Xty(const NEWIMAGE::volume4D<float>& Xvols,
					 const NEWIMAGE::volume<float>&   Yvol,
					 const NEWIMAGE::volume<float>&   mask)
{
  NEWMAT::ColumnVector Xty(Xvols.tsize());
  Xty = 0.0;
  for (int r=1; r<=Xvols.tsize(); r++) {
    for (NEWIMAGE::volume<float>::fast_const_iterator Xit=Xvols[r-1].fbegin(), Xend=Xvols[r-1].fend(), Yit=Yvol.fbegin(), mit=mask.fbegin(); Xit!=Xend; ++Xit, ++Yit, ++mit) {
      if (*mit) Xty(r) += (*Xit)*(*Yit);
    }
  }
  return(Xty);
}

void EddyUtils::transform_coordinates(// Input 
				      const NEWIMAGE::volume<float>&    f,
				      const NEWIMAGE::volume4D<float>&  d,
				      const NEWMAT::Matrix&             M,
				      // Input/Output
				      ImageCoordinates&                 c,
                                      // Output
				      NEWIMAGE::volume<float>           *omask)
{
  NEWMAT::Matrix iA = d[0].sampling_mat();

  float A11=iA(1,1), A12=iA(1,2), A13=iA(1,3), A14=iA(1,4);
  float A21=iA(2,1), A22=iA(2,2), A23=iA(2,3), A24=iA(2,4);
  float A31=iA(3,1), A32=iA(3,2), A33=iA(3,3), A34=iA(3,4); 

  // Create a matrix mapping from mm-coordinates in volume i
  // to voxel coordinates in volume f. If the matrix M is empty
  // this is simply a mm->voxel mapping for volume f

  NEWMAT::Matrix iM = f.sampling_mat().i() * M.i();

  float M11=iM(1,1), M12=iM(1,2), M13=iM(1,3), M14=iM(1,4);
  float M21=iM(2,1), M22=iM(2,2), M23=iM(2,3), M24=iM(2,4);
  float M31=iM(3,1), M32=iM(3,2), M33=iM(3,3), M34=iM(3,4); 

  for (unsigned int z=0, index=0; z<c.NZ(); z++) {
    float xtmp1 = A13*z + A14;
    float ytmp1 = A23*z + A24;
    float ztmp1 = A33*z + A34;
    for (unsigned int y=0; y<c.NY(); y++) {
      float xtmp2 = xtmp1 + A12*y;
      float ytmp2 = ytmp1 + A22*y;
      float ztmp2 = ztmp1 + A32*y;
      for (unsigned int x=0; x<c.NX(); x++) {
	float o1 = xtmp2 + A11*x + d(x,y,z,0);
	float o2 = ytmp2 + A21*x + d(x,y,z,1);
	float o3 = ztmp2 + A31*x + d(x,y,z,2);
	if (omask) (*omask)(x,y,z) = 1;  // So far, so good
	c.x(index) = M11*o1 + M12*o2 + M13*o3 + M14;
	c.y(index) = M21*o1 + M22*o2 + M23*o3 + M24;
	c.z(index) = M31*o1 + M32*o2 + M33*o3 + M34;
	if (omask) (*omask)(x,y,z) *= (f.valid(c.x(index),c.y(index),c.z(index))) ? 1 : 0; // Kosher only if valid in both d and s
        index++;
      }
    }
  }
  return;
}
// }}} End of fold.

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// Class FieldUtils
//
// Helper Class used to perform various useful calculations
// on displacement fields.
//
// {{{ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

NEWIMAGE::volume4D<float> FieldUtils::Hz2VoxelDisplacements(const NEWIMAGE::volume<float>& hzfield,
                                                            const AcqPara&                 acqp)
{
  NEWIMAGE::volume4D<float> dfield(hzfield.xsize(),hzfield.ysize(),hzfield.zsize(),3);
  for (int i=0; i<3; i++) dfield[i] = float((acqp.PhaseEncodeVector())(i+1) * acqp.ReadOutTime()) * hzfield;
  return(dfield);
}

NEWIMAGE::volume4D<float> FieldUtils::Hz2MMDisplacements(const NEWIMAGE::volume<float>& hzfield,
                                                         const AcqPara&                 acqp)
{
  NEWIMAGE::volume4D<float> dfield(hzfield.xsize(),hzfield.ysize(),hzfield.zsize(),3);
  dfield[0] = float(hzfield.xdim()*(acqp.PhaseEncodeVector())(1) * acqp.ReadOutTime()) * hzfield;
  dfield[1] = float(hzfield.ydim()*(acqp.PhaseEncodeVector())(2) * acqp.ReadOutTime()) * hzfield;
  dfield[2] = float(hzfield.zdim()*(acqp.PhaseEncodeVector())(3) * acqp.ReadOutTime()) * hzfield;
  return(dfield);
}

/////////////////////////////////////////////////////////////////////
//
// Inverts a 1D displacementfield. The input field should be in units 
// of voxels and the output will be too.
//
/////////////////////////////////////////////////////////////////////
NEWIMAGE::volume<float> FieldUtils::InvertDisplacementField(// Input 
							    const NEWIMAGE::volume<float>& dfield,
                                                            const AcqPara&                 acqp,
							    const NEWIMAGE::volume<float>& inmask,
							    // Output
							    NEWIMAGE::volume<float>&       omask)
{
  NEWIMAGE::volume<float> fc = dfield;   // fc : field copy
  NEWIMAGE::volume<float> imc = inmask;  // imc: inmask copy
  // Make it so that we invert in first (x) direction
  unsigned int d=0;
  for (; d<3; d++) if ((acqp.PhaseEncodeVector())(d+1)) break;
  if (d==1) {
    fc.swapdimensions(2,1,3);
    imc.swapdimensions(2,1,3);
    omask.swapdimensions(2,1,3);
  }
  else if (d==2) {
    fc.swapdimensions(3,2,1);
    imc.swapdimensions(3,2,1);
    omask.swapdimensions(3,2,1);
  }
  NEWIMAGE::volume<float> idf = fc;    // idf : inverse displacement field
  // Do the inversion
  for (int k=0; k<idf.zsize(); k++) {
    for (int j=0; j<idf.ysize(); j++) {
      int oi=0;
      // Check to see if we think this column can be validly inverted
      for (int i=0; i<idf.xsize(); i++) {
	int ii=oi;
	for (; ii<idf.xsize() && fc(ii,j,k)+ii<i; ii++) ; // On purpose
	if (ii>0 && ii<idf.xsize()) { // If we are in valid range
	  idf(i,j,k) = ii - i - 1.0 + float(i+1-ii-fc(ii-1,j,k))/float(fc(ii,j,k)+1.0-fc(ii-1,j,k));
          if (imc(ii-1,j,k)) omask(i,j,k) = 1.0;
	  else omask(i,j,k) = 0.0;
	}
	else {
	  idf(i,j,k) = FLT_MAX;    // Tag for further processing
	  omask(i,j,k) = 0.0;
	}
	oi = std::max(0,ii-1);
      }
      // Process NaN's at beginning of column
      int ii=0;
      for (ii=0; ii<idf.xsize()-1 && idf(ii,j,k)==FLT_MAX; ii++) ; // On purpose
      for (; ii>0; ii--) idf(ii-1,j,k) = idf(ii,j,k); 
      // Process NaN's at end of column
      for (ii=idf.xsize()-1; ii>0 && idf(ii,j,k)==FLT_MAX; ii--) ; // On purpose
      for (; ii<idf.xsize()-1; ii++) idf(ii+1,j,k) = idf(ii,j,k); 
    }
  }
  // Swap back to original orientation
  if (d==1) {
    idf.swapdimensions(2,1,3);
    omask.swapdimensions(2,1,3);
  }
  else if (d==2) {
    idf.swapdimensions(3,2,1);
    omask.swapdimensions(3,2,1);
  }

  return(idf);
}

/////////////////////////////////////////////////////////////////////
//
// Inverts a "3D" displacementfield. The input field should be in units 
// of voxels and the output will be too. The current implementation
// expects displacements to be in one direction only (i.e. 1D).
//
/////////////////////////////////////////////////////////////////////
NEWIMAGE::volume4D<float> FieldUtils::InvertDisplacementField(// Input
							      const NEWIMAGE::volume4D<float>& dfield,
                                                              const AcqPara&                   acqp,
							      const NEWIMAGE::volume<float>& inmask,
							      // Output
							      NEWIMAGE::volume<float>&       omask)
{
  NEWIMAGE::volume4D<float> idfield = dfield;
  idfield = 0.0;
  unsigned int cnt=0;
  for (unsigned int i=0; i<3; i++) if ((acqp.PhaseEncodeVector())(i+1)) cnt++;
  if (cnt != 1) throw EddyException("FieldUtils::InvertDisplacementField: Phase encode vector must have exactly one non-zero component");
  unsigned int i=0;
  for (; i<3; i++) if ((acqp.PhaseEncodeVector())(i+1)) break;
  idfield[i] = InvertDisplacementField(dfield[i],acqp,inmask,omask);

  return(idfield);
}

/////////////////////////////////////////////////////////////////////
//
// Calculates the Jacobian determinant of a 3D displacement field.
// The field must be in units of voxels and in the present
// implementation it must also be inherently 1D.
//
/////////////////////////////////////////////////////////////////////
NEWIMAGE::volume<float> FieldUtils::GetJacobian(const NEWIMAGE::volume4D<float>& dfield,
                                                const AcqPara&                   acqp)
{
  unsigned int cnt=0;
  for (unsigned int i=0; i<3; i++) if ((acqp.PhaseEncodeVector())(i+1)) cnt++;
  if (cnt != 1) throw EddyException("FieldUtils::GetJacobian: Phase encode vector must have exactly one non-zero component");
  unsigned int i=0;
  for (; i<3; i++) if ((acqp.PhaseEncodeVector())(i+1)) break;

  NEWIMAGE::volume<float> jacfield = GetJacobianFrom1DField(dfield[i],i);

  return(jacfield);  
}

/////////////////////////////////////////////////////////////////////
//
// Calculates the Jacobian determinant of a 1D displacement field.
// The field must be in units of voxels.
//
/////////////////////////////////////////////////////////////////////
NEWIMAGE::volume<float> FieldUtils::GetJacobianFrom1DField(const NEWIMAGE::volume<float>& dfield,
                                                           unsigned int                   dir)
{
  // Calculate spline coefficients for displacement field
  std::vector<unsigned int>                        dim(3,0);
  dim[0] = dfield.xsize(); dim[1] = dfield.ysize(); dim[2] = dfield.zsize();
  std::vector<SPLINTERPOLATOR::ExtrapolationType>  ep(3,SPLINTERPOLATOR::Mirror);
  SPLINTERPOLATOR::Splinterpolator<float> spc(dfield.fbegin(),dim,ep,3,false);
  // Get Jacobian at voxel centres
  NEWIMAGE::volume<float> jacf = dfield;
  for (int k=0; k<dfield.zsize(); k++) {
    for (int j=0; j<dfield.ysize(); j++) {
      for (int i=0; i<dfield.xsize(); i++) {
        jacf(i,j,k) = 1.0 + spc.DerivXYZ(i,j,k,dir);
      }
    }
  }
  return(jacf);
}

// }}} End of fold.
