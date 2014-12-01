
#include <cstdlib>
#include <string>
#include <vector>
#include "utils/options.h"
#include "EddyHelperClasses.h"

#ifndef EddyCommandLineOptions_h
#define EddyCommandLineOptions_h

namespace EDDY {

class EddyCommandLineOptions {
public:
  EddyCommandLineOptions(int argc, char *argv[]);
  std::string ImaFname() const { return(_imain.value()); }
  std::string MaskFname() const { return(_mask.value()); }
  std::string AcqpFname() const { return(_acqp.value()); }
  std::string IndexFname() const { return(_index.value()); }
  std::string SessionFname() const { return(_session.value()); }
  std::string TopupFname() const { return(_topup.value()); }
  std::string BVecsFname() const { return(_bvecs.value()); }
  std::string BValsFname() const { return(_bvals.value()); }
  std::string ParOutFname() const { return(_out.value()+std::string(".eddy_parameters")); }
  std::string IOutFname() const { return(_out.value()); }
  std::string ECFOutFname() const { return(_out.value()+std::string(".eddy_fields")); }
  std::string OLReportFname() const { return(_out.value()+std::string(".eddy_outlier_report")); }
  std::string OLFreeDataFname() const { return(_out.value()+std::string(".eddy_outlier_free_data")); } 
  std::string DwiMssHistoryFname() const { return(_out.value()+std::string(".eddy_dwi_mss_history")); }
  std::string DwiParHistoryFname() const { return(_out.value()+std::string(".eddy_dwi_parameter_history")); }
  std::string B0MssHistoryFname() const { return(_out.value()+std::string(".eddy_b0_mss_history")); }
  std::string B0ParHistoryFname() const { return(_out.value()+std::string(".eddy_b0_parameter_history")); }
  std::vector<unsigned int> Indicies() const { return(_indvec); }
  std::vector<unsigned int> SessionIndicies() const { return(_sessvec); }
  unsigned int NoOfSessions() const { return(_nsess); }
  double FWHM() const { return(_fwhm.value()); }
  unsigned int NIter() const { return(static_cast<unsigned int>(_niter.value())); }
  unsigned int Index(unsigned int i) const { return(_indvec[i]); }
  EDDY::ECModel FirstLevelModel() const;
  bool ReplaceOutliers() const { return(_rep_ol.value()); }
  bool WriteOutlierFreeData() const { return(_rep_ol.value()); }
  double OLNStdev() const { return(_ol_nstd.value()); }
  unsigned int OLNVox() const { return(_ol_nvox.value()); }
  bool RegisterDWI() const { return(_rdwi); }
  bool Registerb0() const { return(_rb0); }
  bool WriteFields() const { return(_fields.value()); }
  bool History() const { return(_history.value()); }
  std::string InitFname() const { return(_init.value()); }
  bool Verbose() const { return(_verbose.value()); }
  bool VeryVerbose() const { return(_very_verbose.value()); }
  bool WriteSliceStats() const { return(_write_slice_stats.value()); }
  int DebugLevel() const { return(_debug.value()); }
  FinalResampling ResamplingMethod() const;
private:
  std::string                     _title;
  std::string                     _examples;
  Utilities::Option<bool>         _verbose;
  Utilities::Option<bool>         _help;
  Utilities::Option<string>       _imain;
  Utilities::Option<string>       _mask;
  Utilities::Option<string>       _acqp;
  Utilities::Option<string>       _index;
  Utilities::Option<string>       _session;
  Utilities::Option<string>       _topup;
  Utilities::Option<string>       _bvecs;
  Utilities::Option<string>       _bvals;
  Utilities::Option<float>        _fwhm;
  Utilities::Option<int>          _niter;
  Utilities::Option<string>       _out;
  Utilities::Option<string>       _flm;
  Utilities::Option<bool>         _rep_ol;
  Utilities::Option<string>       _resamp;
  Utilities::HiddenOption<float>  _ol_nstd;
  Utilities::HiddenOption<int>    _ol_nvox;
  Utilities::HiddenOption<bool>   _very_verbose;
  Utilities::HiddenOption<bool>   _dwi_only;
  Utilities::HiddenOption<bool>   _b0_only;
  Utilities::HiddenOption<bool>   _fields;
  Utilities::HiddenOption<bool>   _history;
  Utilities::HiddenOption<bool>   _write_slice_stats;
  Utilities::HiddenOption<string> _init;
  Utilities::HiddenOption<int>    _debug;
  bool                            _rdwi;
  bool                            _rb0;
  std::vector<unsigned int>       _indvec;
  std::vector<unsigned int>       _sessvec;
  unsigned int                    _nsess;

  void do_initial_parsing(int argc, char *argv[]);
  bool session_indicies_kosher(const std::string& fname, unsigned int ns);
  bool indicies_kosher(NEWMAT::Matrix indx, NEWMAT::Matrix acqp);
};

} // End namespace EDDY

#endif // End #ifndef EddyCommandLineOptions_h
