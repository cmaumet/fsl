// Definitions of a class that does the parsing and
// sanity checking of commnad line options for the 
// "eddy" application.
// 
// EddyCommandLineOptions.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//
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
#include <cstdlib>
#include <string>
#include <vector>
#include "utils/options.h"
#include "topup/topup_file_io.h"
#include "EddyHelperClasses.h"
#include "EddyCommandLineOptions.h"

using namespace EDDY;

EddyCommandLineOptions::EddyCommandLineOptions(int  argc, char *argv[]) : 
  _title("eddy \nCopyright(c) 2011, University of Oxford (Jesper Andersson)"),
  _examples("eddy --monsoon"),
  _verbose(string("-v,--verbose"),false,string("switch on diagnostic messages"),false, Utilities::no_argument), 
  _help(string("-h,--help"),false,string("display this message"),false, Utilities::no_argument),
  _imain(string("--imain"),string(""),string("File containing all the images to estimate distortions for"),true,Utilities::requires_argument), 
  _mask(string("--mask"),string(""),string("Mask to indicate brain"),true,Utilities::requires_argument),
  _acqp(string("--acqp"),string(""),string("File containing acquisition parameters"),true,Utilities::requires_argument),
  _index(string("--index"),string(""),string("File containing indices for all volumes in --imain into --acqp and --topup"),true,Utilities::requires_argument),
  _session(string("--session"),string(""),string("File containing session indices for all volumes in --imain"),false,Utilities::requires_argument),
  _topup(string("--topup"),string(""),string("Base name for output files from topup"),false,Utilities::requires_argument),
  _bvecs(string("--bvecs"),string(""),string("File containing the b-vectors for all volumes in --imain"),true,Utilities::requires_argument),
  _bvals(string("--bvals"),string(""),string("File containing the b-values for all volumes in --imain"),true,Utilities::requires_argument),
  _fwhm(string("--fwhm"),0.0,string("FWHM for conditioning filter when estimating the parameters"),false,Utilities::requires_argument),
  _niter(string("--niter"),5,string("Number of iterations (default 5)"),false,Utilities::requires_argument),
  _out(string("--out"),string(""),string("Basename for output"),true,Utilities::requires_argument),
  _flm(string("--flm"),string("quadratic"),string("First level EC model (linear/quadratic/cubic)"),false,Utilities::requires_argument),
  _rep_ol(string("--repol"),false,string("Detect and replace outlier slices"),false, Utilities::no_argument), 
  _resamp(string("--resamp"),string("jac"),string("Final resampling method (jac/lsr)"),false, Utilities::requires_argument), 
  _ol_nstd(string("--ol_nstd"),4.0,string("Number of std off to qualify as outlier (default 4)"),false,Utilities::requires_argument),
  _ol_nvox(string("--ol_nvox"),250,string("Min # of voxels in a slice for inclusion in outlier detection (default 250)"),false,Utilities::requires_argument),
  _very_verbose(string("--very_verbose"),false,string("Switch on detailed diagnostic messages"),false, Utilities::no_argument), 
  _dwi_only(string("--dwi_only"),false,string("Only register the dwi images"),false, Utilities::no_argument), 
  _b0_only(string("--b0_only"),false,string("Only register the b0 images"),false, Utilities::no_argument), 
  _fields(string("--fields"),false,string("Write EC fields as images"),false, Utilities::no_argument), 
  _history(string("--history"),false,string("Write history of mss and parameter estimates"),false, Utilities::no_argument), 
  _write_slice_stats(string("--wss"),false,string("Write slice-wise stats for each iteration"),false, Utilities::no_argument), 
  _init(string("--init"),string(""),string("Text file with parameters for initialisation"),false,Utilities::requires_argument),
  _debug(string("--debug"),0,string("Level of debug print-outs."),false,Utilities::requires_argument),
  _rdwi(true), _rb0(true)
{
  // Do the parsing implemented by the Option<T> class.
  do_initial_parsing(argc,argv);

  // Do sanity checking and some initial reading of files
  // Image file exists?
  NEWIMAGE::volume4D<float> imahdr;
  NEWIMAGE::read_volume4D_hdr_only(imahdr,_imain.value());  // Throws if there is a problem
  NEWIMAGE::volume<float> maskhdr;
  NEWIMAGE::read_volume_hdr_only(maskhdr,_mask.value());    // Throws if there is a problem
  if (!samesize(imahdr[0],maskhdr,true)) throw EddyException("--imain and --mask images must have the same dimensions");
  // Index file exists and makes sense given image file
  NEWMAT::Matrix index = MISCMATHS::read_ascii_matrix(_index.value());
  if (std::min(index.Nrows(),index.Ncols()) != 1 || std::max(index.Nrows(),index.Ncols()) != imahdr.tsize()) {
    throw EddyException("--index must be an 1xN or Nx1 matrix where N is the number of volumes in --imain");
  }
  // Session file either dont exist, or exists and makes sense.
  if (!this->session_indicies_kosher(_session.value(),imahdr.tsize())) {
    throw EddyException("--session file does not make sense");
  }
  // File with acquisition parameters exists and makes sense
  NEWMAT::Matrix acqpM = MISCMATHS::read_ascii_matrix(_acqp.value());
  if (acqpM.Ncols() != 4) throw EddyException("Each row of the --acqp file must contain 4 values");
  if (!this->indicies_kosher(index,acqpM)) throw EddyException("Mismatch between --index and --acqp");
  // bvals and bvecs exists and are consistent with image file
  NEWMAT::Matrix bvecsM = MISCMATHS::read_ascii_matrix(_bvecs.value());
  if (std::min(bvecsM.Nrows(),bvecsM.Ncols()) != 3 || std::max(bvecsM.Nrows(),bvecsM.Ncols()) != imahdr.tsize()) {
    throw EddyException("--bvecs should contain a 3xN or Nx3 matrix where N is the number of volumes in --imain");    
  }
  NEWMAT::Matrix bvalsM = MISCMATHS::read_ascii_matrix(_bvals.value());
  if (std::min(bvalsM.Nrows(),bvalsM.Ncols()) != 1 || std::max(bvalsM.Nrows(),bvalsM.Ncols()) != imahdr.tsize()) {
    throw EddyException("--bvals should contain a 1xN or Nx1 matrix where N is the number of volumes in --imain");    
  }
  // Topup file exists if one has been specified
  if (_topup.set() && _topup.value() != string("")) {
    try { TOPUP::TopupFileReader  tfr(_topup.value()); }
    catch (...) { throw EddyException("Error when attempting to read --topup files"); }
  }
  // Valid final resampling method?
  if (_resamp.value() != string("jac") && _resamp.value() != string("lsr")) throw EddyException("Invalid --resamp parameter");
  // Valid first level model?
  if (_flm.value() != string("linear") && _flm.value() != string("quadratic") && _flm.value() != string("cubic")) throw EddyException("Invalid --flm parameter");
  // Reasonable FWHM
  if (_fwhm.value() < 0.0 || _fwhm.value() > 15.0) throw EddyException("--fwhm value outside valid range");
  // If very_verbose is set we should also set verbose
  if (_very_verbose.value() && !_verbose.value()) _verbose.set_T(true);
  // Make sure the "only" flags are mutually exclusive
  if (_dwi_only.value() && _b0_only.value()) throw EddyException("--dwi_only and --b0_only cannnot both be set");
  if (_dwi_only.value()) _rb0 = false;
  else if (_b0_only.value()) _rdwi = false;
  // Initialisation file exists if one has been specified
  if (_init.set() && _init.value() != string("")) {
    try { NEWMAT::Matrix  tmp = MISCMATHS::read_ascii_matrix(_init.value()); }
    catch (...) { throw EddyException("Error when attempting to read --init file"); }
  }
}

EDDY::FinalResampling EddyCommandLineOptions::ResamplingMethod() const
{
  if (_resamp.value() == std::string("jac")) return(EDDY::JAC);
  else if (_resamp.value() == std::string("lsr")) return(EDDY::LSR);
  else return(EDDY::UNKNOWN_RESAMPLING);
}

EDDY::ECModel EddyCommandLineOptions::FirstLevelModel() const
{
  if (_flm.value() == string("linear")) return(EDDY::Linear);
  else if (_flm.value() == string("quadratic")) return(EDDY::Quadratic);
  else if (_flm.value() == string("cubic")) return(EDDY::Cubic);
  return(EDDY::Unknown);
}

void EddyCommandLineOptions::do_initial_parsing(int argc, char *argv[])
{
  Utilities::OptionParser options(_title,_examples);

  try {
    options.add(_imain);
    options.add(_mask);
    options.add(_index);
    options.add(_session);
    options.add(_acqp);
    options.add(_topup);
    options.add(_bvecs);
    options.add(_bvals);
    options.add(_flm);
    options.add(_fwhm);
    options.add(_niter);
    options.add(_out);
    options.add(_very_verbose);
    options.add(_dwi_only);
    options.add(_b0_only);
    options.add(_fields);
    options.add(_history);
    options.add(_resamp);
    options.add(_write_slice_stats);
    options.add(_rep_ol);
    options.add(_ol_nstd);
    options.add(_ol_nvox);
    options.add(_init);
    options.add(_debug);
    options.add(_verbose);
    options.add(_help);

    int i=options.parse_command_line(argc, argv);
    if (i < argc) {
      for (; i<argc; i++) {
        cerr << "Unknown input: " << argv[i] << endl;
      }
      exit(EXIT_FAILURE);
    }

    if (_help.value() || !options.check_compulsory_arguments(true)) {
      options.usage();
      exit(EXIT_FAILURE);
    }
  }  
  catch(Utilities::X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
}

// Makes sure there is nothing obviously wrong with the session
// indicies, and stores them in _sessvec.
bool EddyCommandLineOptions::session_indicies_kosher(const std::string& fname, unsigned int ns)
{
  bool rval=true;

  if (!fname.length()) {
    _sessvec.resize(ns,1);
    _nsess = 1;
  }
  else {
    NEWMAT::Matrix session = MISCMATHS::read_ascii_matrix(fname);
    if (std::min(session.Nrows(),session.Ncols()) != 1 || std::max(session.Nrows(),session.Ncols()) != int(ns)) {
      throw EddyException("--session must be an 1xN or Nx1 matrix where N is the number of volumes in --imain");
    }
    if (session.Ncols() > session.Nrows()) session = session.t();
    _sessvec.resize(ns,0);
    unsigned int max_sess = static_cast<unsigned int>(session.Maximum());
    std::vector<unsigned int> scan_count(max_sess,0);
    for (unsigned int i=0; i<ns; i++) {
      _sessvec[i] = static_cast<unsigned int>(session(i+1,1));
      scan_count[_sessvec[i]-1]++;
    }
    for (unsigned int i=0; i<max_sess; i++) {
      if (!scan_count[i]) {
        rval=false;
	break;
      }
      else _nsess = max_sess;
    }
  }
  return(rval);
}

// Makes sure that the indicies makes sense w.r.t. pointing in to
// the acquistition parameter file, and stores them in _indvec.
bool EddyCommandLineOptions::indicies_kosher(NEWMAT::Matrix indx, NEWMAT::Matrix acqp)
{
  if (indx.Ncols() > indx.Nrows()) indx = indx.t();
  _indvec.resize(indx.Nrows());
  unsigned int max_indx = static_cast<unsigned int>(indx(1,1));
  for (int i=0; i<indx.Nrows(); i++) {
    _indvec[i] = static_cast<unsigned int>(indx(i+1,1));
    if (fabs(static_cast<double>(_indvec[i])-indx(i+1,1)) > 1e-6) return(false);
    max_indx = std::max(max_indx,_indvec[i]);
  }
  if (max_indx > static_cast<unsigned int>(acqp.Nrows())) return(false);
  else return(true);
}
