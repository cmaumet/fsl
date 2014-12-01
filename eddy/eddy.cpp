/*! \file eddy.cpp
    \brief Contains main() and some very high level functions for eddy
*/
#pragma clang diagnostic ignored "-Wunknown-pragmas" // Ignore the OpenMP pragmas

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "newmat.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "stack_dump.h"
#include "EddyHelperClasses.h"
#include "ECScanClasses.h"
#include "DiffusionGP.h"
#include "b0Predictor.h"
#include "EddyUtils.h"
#include "EddyCommandLineOptions.h"

using namespace EDDY;

/// Global function that registers a set of scans together
ReplacementManager Register(const EddyCommandLineOptions&  clo,     // Input
			    ScanType                       st,      // Input
			    ECScanManager&                 sm,      // Input/Output
			    unsigned int                   niter,   // Input
			    NEWMAT::Matrix&                msshist, // Output
			    NEWMAT::Matrix&                phist);  // Output

/// Global function that detect outlier slices and replaces them by their expectation
DiffStatsVector DetectAndReplaceOutliers(// Input
					 const EddyCommandLineOptions& clo,
					 ScanType                      st,
					 // Input/Output
					 ECScanManager&                sm,
					 ReplacementManager&           rm);

/// Global function that Loads up the prediction maker with unwarped scans
boost::shared_ptr<DWIPredictionMaker> LoadPredictionMaker(// Input
							  const EddyCommandLineOptions& clo,
							  ScanType                      st,
							  const ECScanManager&          sm,
							  // Output
							  NEWIMAGE::volume<float>&      mask);

/// Global function that Generates diagnostic information for subsequent analysis
void Diagnostics(const EddyCommandLineOptions&  clo,      // Input
		 unsigned int                   iter,     // Input
		 ScanType                       st,       // Input
		 const ECScanManager&           sm,       // Input
                 const double                   *mss_tmp, // Input
                 const DiffStatsVector&         stats,    // Input
		 NEWMAT::Matrix&                mss,      // Output
		 NEWMAT::Matrix&                phist);   // Output

/// The entry point of eddy.
int main(int argc, char *argv[])
{
  StackDump::Install(); // Gives us informative stack dump if/when program crashes

  // Parse comand line input
  EddyCommandLineOptions clo(argc,argv); // Command Line Options

  // Read all available info
  if (clo.Verbose()) cout << "Reading images" << endl;
  ECScanManager sm(clo.ImaFname(),clo.MaskFname(),clo.AcqpFname(),clo.TopupFname(),
                   clo.BVecsFname(),clo.BValsFname(),clo.FirstLevelModel(),
		   clo.Indicies(),clo.SessionIndicies(),clo.NoOfSessions()); // Scan Manager
  if (clo.Verbose() && clo.FWHM()) cout << "Smoothing images" << endl;
  sm.SetFWHM(clo.FWHM());
  if (clo.ResamplingMethod() == LSR) {
    if (!sm.CanDoLSRResampling()) throw EddyException("These data do not support least-squares resampling");
  }

  // Set initial parameters. This option is only for testing/debugging/personal use
  if (clo.InitFname() != std::string("")) {
    if (clo.RegisterDWI() && clo.Registerb0()) sm.SetParameters(clo.InitFname(),ANY);
    else if (clo.RegisterDWI()) sm.SetParameters(clo.InitFname(),DWI);
    else sm.SetParameters(clo.InitFname(),B0);
  }

  // Do the registration
  NEWMAT::Matrix dwi_mss, b0_mss, dwi_ph, b0_ph;
  ReplacementManager *dwi_rm;
  if (clo.NIter() && clo.RegisterDWI()) {
    cout << "Running Register" << endl;
    dwi_rm = new ReplacementManager(Register(clo,DWI,sm,clo.NIter(),dwi_mss,dwi_ph));
    // Write outlier information
    if (clo.ReplaceOutliers()) {
      cout << "Running sm.GetDwi2GlobalIndexMapping" << endl;
      std::vector<unsigned int> i2i = sm.GetDwi2GlobalIndexMapping();
      cout << "Running dwi_rm->WriteReport" << endl;
      dwi_rm->WriteReport(i2i,clo.OLReportFname());
      if (clo.WriteOutlierFreeData()) {
	cout << "Running sm.WriteOutlierFreeData" << endl;
        sm.WriteOutlierFreeData(clo.OLFreeDataFname());
      }
    }
  }
  if (clo.NIter() && clo.Registerb0()) {
    cout << "Running Register" << endl;
    ReplacementManager b0_rm = Register(clo,B0,sm,clo.NIter(),b0_mss,b0_ph);
  }

  /*
  // Write raw registration parameters (for debugging).
  cout << "Running sm.WriteParameterFile for raw parameters" << endl;
  sm.WriteParameterFile(clo.ParOutFname()+string(".raw"));
  cout << "Running sm.WriteRegisteredImages for raw parameters" << endl;
  sm.WriteRegisteredImages(clo.IOutFname()+string(".raw"),clo.ResamplingMethod());  
  */  

  // Separate field offset from subject movement in PE direction
  // if (clo.RegisterDWI()) { cout << "Running sm.SeparateFieldOffsetFromMovement" << endl; sm.SeparateFieldOffsetFromMovement(); }

  // Set reference for location
  if (clo.RegisterDWI()) { cout << "Running sm.SetDWIReference" << endl; sm.SetDWIReference(); }
  if (clo.Registerb0()) { cout << "Running sm.Setb0Reference" << endl; sm.Setb0Reference(); }

  // Write registration parameters
  cout << "Running sm.WriteParameterFile" << endl;
  if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteParameterFile(clo.ParOutFname());
  else if (clo.RegisterDWI()) sm.WriteParameterFile(clo.ParOutFname(),DWI);
  else sm.WriteParameterFile(clo.ParOutFname(),B0);

  // Write registered images
  cout << "Running sm.WriteRegisteredImages" << endl;
  if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod());
  else if (clo.RegisterDWI()) sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod(),DWI);
  else sm.WriteRegisteredImages(clo.IOutFname(),clo.ResamplingMethod(),B0);

  // Write EC fields
  cout << "Running sm.WriteECFields" << endl;
  if (clo.WriteFields()) {
    if (clo.RegisterDWI() && clo.Registerb0()) sm.WriteECFields(clo.ECFOutFname());
    else if (clo.RegisterDWI()) sm.WriteECFields(clo.ECFOutFname(),DWI);
    else sm.WriteECFields(clo.ECFOutFname(),B0);
  }

  if (clo.NIter() && clo.History()) { 
    if (clo.RegisterDWI()) {
      MISCMATHS::write_ascii_matrix(clo.DwiMssHistoryFname(),dwi_mss); 
      MISCMATHS::write_ascii_matrix(clo.DwiParHistoryFname(),dwi_ph);
    }
    if (clo.Registerb0()) {
      MISCMATHS::write_ascii_matrix(clo.B0MssHistoryFname(),b0_mss); 
      MISCMATHS::write_ascii_matrix(clo.B0ParHistoryFname(),b0_ph);
    }
  }
  exit(EXIT_SUCCESS);
}

/****************************************************************//**
*  								  
*  A global function that registers the scans in sm.
*  \param[in] clo Carries information about the command line options 
*  that eddy was invoked with.
*  \param[in] st Specifies if we should register the diffusion weighted 
*  images or the b=0 images.
*  \param[in,out] sm Collection of all scans. Will be updated by this call.
*  \param[in] niter Specifies how many iterations should be run
*  \param[out] msshist Returns the history of the mss. msshist(i,j) 
*  contains the mss for the jth scan on the ith iteration.
*  \param[out] phist Returns the history of the estimated parameters. 
*  phist(i,j) contains the jth parameter on the ith iteration
*  \return A ReplacementManager that details which slices in which 
*  scans were replaced by their expectations.
*
********************************************************************/ 

ReplacementManager Register(const EddyCommandLineOptions&  clo,     // Input
			    ScanType                       st,      // Input
			    ECScanManager&                 sm,      // Input/Output
			    unsigned int                   niter,   // Input
			    NEWMAT::Matrix&                msshist, // Output
			    NEWMAT::Matrix&                phist)   // Output
{
  msshist.ReSize(niter,sm.NScans(st));
  phist.ReSize(niter,sm.NScans(st)*sm.Scan(0,st).NParam());
  double *mss_tmp = new double[sm.NScans(st)]; 
  ReplacementManager rm(sm.NScans(st),static_cast<unsigned int>(sm.Scan(0,st).GetIma().zsize()),clo.OLNStdev(),clo.OLNVox());
  NEWIMAGE::volume<float> mask = sm.Scan(0,st).GetIma(); EddyUtils::SetTrilinearInterp(mask); mask = 1.0; // Mask in model space

  for (unsigned int iter=0; iter<niter; iter++) {
    // Detect outliers and replace them
    DiffStatsVector stats(sm.NScans(st));
    if (iter && st==DWI && clo.ReplaceOutliers()) stats = DetectAndReplaceOutliers(clo,st,sm,rm);

    // Load prediction maker in model space
    // printf("Starting to load and evaluate predictionmaker\n");
    // clock_t stime = clock();
    boost::shared_ptr<DWIPredictionMaker> pmp = LoadPredictionMaker(clo,st,sm,mask);
    // printf("It took %f sec to load and evaluate predictionmaker\n",double(stime-clock())/double(CLOCKS_PER_SEC));

    // Calculate the parameter updates
    if (clo.Verbose()) cout << "Calculating parameter updates" << endl;
    // int dbl = clo.DebugLevel();  // Can't check DebugLevel inside parallel section :(
    // bool vv = clo.VeryVerbose(); // Compiler crashes if I test clo.VeryVerbose() inside parallel loop :(
# pragma omp parallel for shared(mss_tmp, pmp)
    for (int s=0; s<int(sm.NScans(st)); s++) {
      // Get prediction in model space 
      NEWIMAGE::volume<float> pred = pmp->Predict(s);
      // Update parameters
      if (clo.DebugLevel()) {
	mss_tmp[s] = EddyUtils::param_update_debug(pred,sm.GetSuscHzOffResField(),mask,ALL,true,s,iter,clo.DebugLevel(),sm.Scan(s,st),NULL);
      }
      else mss_tmp[s] = EddyUtils::MovAndECParamUpdate(pred,sm.GetSuscHzOffResField(),mask,true,sm.Scan(s,st));
      if (clo.VeryVerbose()) printf("Iter: %d, scan: %d, mss = %f\n",iter,s,mss_tmp[s]);
    }

    // Print/collect some information that can be used for diagnostics
    Diagnostics(clo,iter,st,sm,mss_tmp,stats,msshist,phist);

    // Maybe use model based EC parameters
    // if () sm.SetPredictedECParam();
  }

  delete [] mss_tmp;
  return(rm);
}

/****************************************************************//**
*  								  
*  A global function that loads up a prediction maker with all scans 
*  of a given type. It will load it with unwarped scans (given the 
*  current estimates of the warps) as served up by sm.GetUnwarpedScan().
*  \param[in] clo Carries information about the command line options 
*  that eddy was invoked with.
*  \param[in] st Specifies if we should register the diffusion weighted 
*  images or the b=0 images. If it is set to DWI the function will return 
*  an EDDY::DiffusionGP prediction maker and if it is set to b0 it will 
*  return an EDDY::b0Predictor.
*  \param[in] sm Collection of all scans.
*  \param[out] mask Returns a mask that indicates the voxels where data 
*  is present for all input scans in sm.
*  \return A safe pointer to a DWIPredictionMaker that can be used to 
*  make predictions about what the scans should look like in undistorted space.
*
********************************************************************/ 
boost::shared_ptr<DWIPredictionMaker> LoadPredictionMaker(// Input
							  const EddyCommandLineOptions& clo,
							  ScanType                      st,
							  const ECScanManager&          sm,
							  // Output
							  NEWIMAGE::volume<float>&      mask)
{
  boost::shared_ptr<DWIPredictionMaker>  pmp;                                 // Prediction Maker Pointer
  if (st==DWI) pmp = boost::shared_ptr<DWIPredictionMaker>(new DiffusionGP);  // Gaussian Process
  else pmp = boost::shared_ptr<DWIPredictionMaker>(new b0Predictor);          // Silly mean predictor
  pmp->SetNoOfScans(sm.NScans(st));
  mask = sm.Scan(0,st).GetIma(); EddyUtils::SetTrilinearInterp(mask); mask = 1.0;

  if (clo.Verbose()) cout << "Loading prediction maker";
  // bool vv = clo.VeryVerbose(); // Compiler crashes if I test clo.VeryVerbose() inside parallel loop :(
  if (clo.VeryVerbose()) cout << endl << "Scan: ";
#pragma omp parallel for shared(pmp,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    if (clo.VeryVerbose()) printf(" %d",s);
    NEWIMAGE::volume<float> tmpmask = sm.Scan(s,st).GetIma(); 
    EddyUtils::SetTrilinearInterp(tmpmask); tmpmask = 1.0;
    pmp->SetScan(sm.GetUnwarpedScan(s,tmpmask,st),sm.Scan(s,st).GetDiffPara(),s);
#pragma omp critical
    {
      mask *= tmpmask;
    }
  }
  if (clo.Verbose()) cout << endl << "Evaluating prediction maker model" << endl;
  pmp->EvaluateModel(sm.Mask()*mask);

  return(pmp);
}

void Diagnostics(const EddyCommandLineOptions&  clo,      // Input
		 unsigned int                   iter,     // Input
		 ScanType                       st,       // Input
		 const ECScanManager&           sm,       // Input
                 const double                   *mss_tmp, // Input
                 const DiffStatsVector&         stats,    // Input
		 NEWMAT::Matrix&                mss,      // Output
		 NEWMAT::Matrix&                phist)    // Output
{
  if (clo.Verbose()) {
    double tss=0.0;
    for (unsigned int s=0; s<sm.NScans(st); s++) tss+=mss_tmp[s]; 
    cout << "Iter: " << iter << ", Total mss = " << tss/sm.NScans(st) << endl;
  }
  if (clo.History()) {
    for (unsigned int s=0; s<sm.NScans(st); s++) {
      mss(iter+1,s+1) = mss_tmp[s];
      phist.SubMatrix(iter+1,iter+1,s*sm.Scan(0,st).NParam()+1,(s+1)*sm.Scan(0,st).NParam()) = sm.Scan(s,st).GetParams().t();
    }
  }
  if (clo.WriteSliceStats()) {
    char istring[256];
    sprintf(istring,"EddySliceStatsIteration%02d",iter);
    stats.Write(string(istring));
  }
}

DiffStatsVector DetectAndReplaceOutliers(// Input
					 const EddyCommandLineOptions& clo,
					 ScanType                      st,
					 // Input/Output
					 ECScanManager&                sm,
					 ReplacementManager&           rm)
{
  // printf("Entering DetectAndReplaceOutliers\n");
  // Load up a prediction maker with unwarped and unsmoothed data
  DWIPredictionMaker  *pmp;            // Prediction maker used for the outlier detection
  if (st==DWI) pmp = new DiffusionGP;  // Gaussian Process
  else pmp = new b0Predictor;          // Silly mean predictor
  pmp->SetNoOfScans(sm.NScans(st));
  NEWIMAGE::volume<float> mask(sm.Scan(0,st).GetIma().xsize(),sm.Scan(0,st).GetIma().ysize(),sm.Scan(0,st).GetIma().zsize());
  NEWIMAGE::copybasicproperties(sm.Scan(0,st).GetIma(),mask); mask = 1.0;
  // printf("Starting to load predictionmaker\n");
  // clock_t stime = clock();
#pragma omp parallel for shared(pmp,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    NEWIMAGE::volume<float> tmpmask(sm.Scan(s,st).GetIma().xsize(),sm.Scan(s,st).GetIma().ysize(),sm.Scan(s,st).GetIma().zsize());
    NEWIMAGE::copybasicproperties(sm.Scan(s,st).GetIma(),tmpmask); tmpmask = 1.0;
    pmp->SetScan(sm.GetUnwarpedOrigScan(s,tmpmask,st),sm.Scan(s,st).GetDiffPara(),s);
#pragma omp critical
    { mask *= tmpmask; }
  }
  // printf("It took %f sec to load pm\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // printf("Starting to evaluate the pm model\n");
  // stime = clock();
  pmp->EvaluateModel(sm.Mask()*mask);
  // printf("It took %f sec to evaluate pm\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // Use these to generate slice-wise stats on difference between observation and prediction
  // printf("Starting to generate slice-wise stats\n");
  // stime = clock();
  DiffStatsVector stats(sm.NScans(st));
#pragma omp parallel for shared(pmp,mask,stats,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    NEWIMAGE::volume<float> pred = pmp->Predict(s);
    stats[s] = EddyUtils::GetSliceWiseStats(pred,sm.GetSuscHzOffResField(),mask,sm.Mask(),sm.Scan(s,st));
  }
  // printf("It took %f sec to generate stats\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // stats.Write("slice_wise_stats");
  // Detect outliers and update replacement manager
  // printf("Updating replacement manager\n");
  // stime = clock();
  rm.Update(stats);
  // printf("It took %f sec to update\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  // Replace outlier slices with their predictions
  // printf("Replacing outliers\n");
  // stime = clock();
#pragma omp parallel for shared(pmp,mask,st)
  for (int s=0; s<int(sm.NScans(st)); s++) {
    // printf("Checking scan %d for outliers\n",s);
    std::vector<unsigned int> ol = rm.OutliersInScan(s);
    if (ol.size()) { // If this scan has outlier slices
      // printf("Scan %d has %d outlier slices\n",s,ol.size());
      NEWIMAGE::volume<float> pred = pmp->Predict(s);
      sm.Scan(s,st).ReplaceSlices(pred,sm.GetSuscHzOffResField(),mask,ol);
    }
  }
  // printf("It took %f sec to replace outliers\n",double(clock()-stime)/double(CLOCKS_PER_SEC));
  //  rm.DumpOutlierMap("DebugReplacementManager");
  //  exit(EXIT_SUCCESS);
  delete pmp;
  return(stats);
}

/*! \mainpage
 * Here goes a description of the eddy project
 */
