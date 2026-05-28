#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class PHCompositeNode;

class TFile;
class TTree;

class InModuleTrackContainer;

class TrkrHitSetContainer;
class TrkrHitSet;
class ActsGeometry;


// ===================================================================
// Per-module thread data
// ===================================================================
struct InModuleThreadData
{
  InModuleThreadData();

  struct LayerHitSet
  {
    LayerHitSet();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrHitSet* hitset;
  };

  struct RawHit
  {
    RawHit();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrDefs::hitkey hitkey;

    unsigned short pad;
    unsigned short tbin;
    unsigned short adc;

    double local_phi;
    double local_radius;
  };

  struct Blob
  {
    Blob();

    unsigned int layer;

    double pad;
    double tbin;
    double adc;

    unsigned int nhits;
    int used;

    std::vector<unsigned int> raw_hit_indices;
  };

  struct Track
  {
    Track();

    unsigned int track_id;

    unsigned int first_layer;
    unsigned int last_layer;

    unsigned int nblobs;
    unsigned int nrawhits;

    // Final fit parameters.
    // These are refit after disconnected pieces are connected.
    double pad_slope;
    double pad_intercept;

    double tbin_slope;
    double tbin_intercept;

    double chi2_pad;
    double chi2_tbin;

    int ndof_pad;
    int ndof_tbin;


    bool has_local_line_fit;
    bool has_radius_phi_line_fit;
    bool has_radius_phi_circle_fit;
    bool has_sagitta_fit;

    // Sagitta fit parameters in the rotated frame used by weighted_circle_fit.
    // These are always set when use_field_on_fit=true, even when the
    // conversion back to a circle center fails.  The display evaluates
    // sagitta_model_invR(r, S, x0_rot, invR) directly so it never falls
    // back to a straight radius-phi line.
    //
    //   local_x_rot =  cos(sagitta_theta)*local_x + sin(sagitta_theta)*(local_y - sagitta_b)
    //   local_y_rot = -sin(sagitta_theta)*local_x + cos(sagitta_theta)*(local_y - sagitta_b)
    //   phi_fit(r) = asin( sagitta_model_invR(r_rot, S, x0_rot, invR) / r )  ... see display
    double sagitta_S;
    double sagitta_x0;
    double sagitta_invR;
    double sagitta_theta;   // rotation angle: atan(slope of local_x vs local_y line fit)
    double sagitta_b;       // y-intercept of the pre-rotation line fit

    double radius_tbin_slope;
    double radius_tbin_intercept;
    double chi2_radius_tbin;
    int ndof_radius_tbin;

    double radius_phi_slope;
    double radius_phi_intercept;
    double chi2_radius_phi_line;
    int ndof_radius_phi_line;

    double circle_x0;
    double circle_y0;
    double circle_radius;
    double chi2_radius_phi_circle;
    int ndof_radius_phi_circle;

    double chi2_radius_phi_sagitta;
    int ndof_radius_phi_sagitta;

    // Blob indices used by this track.
    std::vector<unsigned int> blob_indices;

    // Raw-hit indices into InModuleThreadData::raw_hits used by this track.
    // These are translated to (hitsetkey, hitkey) pairs when filling the
    // output InModuleTrack objects.
    std::vector<unsigned int> raw_hit_indices;
  };

  // Detector/module identity
  unsigned int region;
  unsigned int sector;
  int side;
  TrkrDefs::hitsetkey module_key;

  ActsGeometry* tGeometry;

  // General configuration
  double pedestal;
  int verbosity;
  bool use_field_on_fit;

  // Noise rejection
  int noise_max_consecutive_timebins;
  int noise_keep_first_timebins;
  int noise_adc_tolerance;  

  // Blob building
  int blob_dt;
  int blob_dp;

  // Initial chain growing
  int search_dt;
  int search_dp;
  unsigned int min_track_blobs;

  // Track-piece connection.
  unsigned int connect_max_layer_gap;

  double connect_dp;
  double connect_dt;

  double connect_dpad_slope;
  double connect_dtbin_slope;

  // ADC weighting for fits
  double weight_power;
  double adc_weight_floor_frac;

  // Per-module containers
  std::vector<LayerHitSet> layer_hitsets;
  std::vector<RawHit> raw_hits;
  std::vector<Blob> blobs;
  std::vector<Track> tracks;
};

// ===================================================================
// Main Fun4All module
// ===================================================================
class InModuleTracks : public SubsysReco
{
 public:
  explicit InModuleTracks(const std::string& name = "InModuleTracks",
                          const std::string& filename = "InModuleTracks.root");

  virtual ~InModuleTracks();

  int Init(PHCompositeNode*);
  int InitRun(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

  void setMaxThreads(unsigned int n);


  void setUseFieldOnFit(bool b) { m_useFieldOnFit = b; }


  void setPedestal(double p)
  {
    m_pedestal = p;
  }

  void setBlobWindow(int dt, int dp)
  {
    m_blob_dt = dt;
    m_blob_dp = dp;
  }

  void setSearchWindow(int dt, int dp)
  {
    m_search_dt = dt;
    m_search_dp = dp;
  }

  // Reject long same-pad tails before blob/track building.
  // If a same-layer/same-pad consecutive-timebin run is longer than
  // max_consecutive_timebins, keep the first keep_first_timebins hits and
  // remove the following non-increasing tail.
  // Defaults: max_consecutive_timebins = 40, keep_first_timebins = 3.
  void setNoiseRejection(int max_consecutive_timebins = 10,
                         int keep_first_timebins = 3,
                         int adc_tolerance = 5)
  {
    m_noiseMaxConsecutiveTimebins = max_consecutive_timebins;
    m_noiseKeepFirstTimebins = keep_first_timebins;
    m_noiseAdcTolerance = adc_tolerance;
  }

  void setMinTrackBlobs(unsigned int n)
  {
    m_minTrackBlobs = n;
  }

  void setConnectMaxLayerGap(unsigned int n)
  {
    m_connectMaxLayerGap = n;
  }

  void setConnectWindow(double dt, double dp)
  {
    m_connect_dt = dt;
    m_connect_dp = dp;
  }

  void setConnectSlopeWindow(double dtbin_slope, double dpad_slope)
  {
    m_connect_dtbin_slope = dtbin_slope;
    m_connect_dpad_slope = dpad_slope;
  }

 private:
  int getNodes(PHCompositeNode*);
  void reset_tree_vars();
  int createNodes(PHCompositeNode*);

  std::string m_outputFileName;

  TFile* m_outputFile;
  TTree* m_tree;

  TrkrHitSetContainer* m_hits;
  ActsGeometry* m_tGeometry;

  InModuleTrackContainer* m_inModuleTrackContainer;
  int m_event;
  unsigned int m_maxThreads;

  bool m_useFieldOnFit;

  // General configuration
  double m_pedestal;

  // Noise rejection
  int m_noiseMaxConsecutiveTimebins;
  int m_noiseKeepFirstTimebins;
  int m_noiseAdcTolerance;

  // Blob building
  int m_blob_dt;
  int m_blob_dp;

  // Initial chain growing
  int m_search_dt;
  int m_search_dp;
  unsigned int m_minTrackBlobs;

  // Track-piece connection parameters
  unsigned int m_connectMaxLayerGap;

  double m_connect_dp;
  double m_connect_dt;

  double m_connect_dpad_slope;
  double m_connect_dtbin_slope;

  // Event number saved once per tree entry
  int m_tree_event;

  // One entry per found module-track.
  std::vector<unsigned int> m_tree_track_id;
  std::vector<unsigned int> m_tree_region;
  std::vector<unsigned int> m_tree_sector;
  std::vector<int> m_tree_side;

  std::vector<unsigned int> m_tree_nblobs;
  std::vector<unsigned int> m_tree_nrawhits;

  std::vector<unsigned int> m_tree_first_layer;
  std::vector<unsigned int> m_tree_last_layer;

  std::vector<double> m_tree_pad_slope;
  std::vector<double> m_tree_pad_intercept;

  std::vector<double> m_tree_tbin_slope;
  std::vector<double> m_tree_tbin_intercept;

  std::vector<double> m_tree_chi2_pad;
  std::vector<double> m_tree_chi2_tbin;

  std::vector<int> m_tree_ndof_pad;
  std::vector<int> m_tree_ndof_tbin;

  std::vector<int> m_tree_has_sagitta_fit;
  std::vector<double> m_tree_sagitta_S;
  std::vector<double> m_tree_sagitta_x0;
  std::vector<double> m_tree_sagitta_invR;
  std::vector<double> m_tree_sagitta_theta;
  std::vector<double> m_tree_sagitta_b;
  std::vector<double> m_tree_chi2_radius_phi_sagitta;
  std::vector<int> m_tree_ndof_radius_phi_sagitta;

  std::vector<double> m_tree_radius_tbin_slope;
  std::vector<double> m_tree_radius_tbin_intercept;
  std::vector<double> m_tree_chi2_radius_tbin;
  std::vector<int> m_tree_ndof_radius_tbin;

  std::vector<double> m_tree_radius_phi_slope;
  std::vector<double> m_tree_radius_phi_intercept;
  std::vector<double> m_tree_chi2_radius_phi_line;
  std::vector<int> m_tree_ndof_radius_phi_line;

  std::vector<double> m_tree_circle_x0;
  std::vector<double> m_tree_circle_y0;
  std::vector<double> m_tree_circle_radius;
  std::vector<double> m_tree_chi2_radius_phi_circle;
  std::vector<int> m_tree_ndof_radius_phi_circle;

  // Flat per-hit content for TTree reading.
  // Hits are identified by their TrkrHitSetContainer keys only;
  // no hit data is duplicated here.
  std::vector<unsigned int> m_tree_hit_event;
  std::vector<unsigned int> m_tree_hit_track_id;

  std::vector<unsigned int> m_tree_hit_region;
  std::vector<unsigned int> m_tree_hit_sector;
  std::vector<int> m_tree_hit_side;

  std::vector<unsigned int> m_tree_hit_layer;

  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};