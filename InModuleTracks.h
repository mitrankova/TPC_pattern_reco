#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrHitSet;
class ActsGeometry;
class TFile;
class TTree;

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

    // Blob indices used by this track.
    std::vector<unsigned int> blob_indices;

    // Raw-hit indices used by this track.
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

  // Blob building
  int blob_dt;
  int blob_dp;

  // Initial chain growing
  int search_dt;
  int search_dp;
  unsigned int min_track_blobs;

  // Track-piece connection.
  // These parameters are used to connect disconnected chains inside one module,
  // for example across dead regions or missing layers.
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

  void setMinTrackBlobs(unsigned int n)
  {
    m_minTrackBlobs = n;
  }

  // Maximum number of missing layers allowed between two disconnected
  // track pieces that can still be connected.
  void setConnectMaxLayerGap(unsigned int n)
  {
    m_connectMaxLayerGap = n;
  }

  // Matching window for connecting track pieces.
  // dt is the allowed tbin difference at the connection layer.
  // dp is the allowed pad difference at the connection layer.
  void setConnectWindow(double dt, double dp)
  {
    m_connect_dt = dt;
    m_connect_dp = dp;
  }

  // Matching window for slope compatibility between track pieces.
  // dtbin_slope is in tbin/layer.
  // dpad_slope is in pad/layer.
  void setConnectSlopeWindow(double dtbin_slope, double dpad_slope)
  {
    m_connect_dtbin_slope = dtbin_slope;
    m_connect_dpad_slope = dpad_slope;
  }

 private:
  int getNodes(PHCompositeNode*);
  void reset_tree_vars();

  std::string m_outputFileName;

  TFile* m_outputFile;
  TTree* m_tree;

  TrkrHitSetContainer* m_hits;
  ActsGeometry* m_tGeometry;

  int m_event;
  unsigned int m_maxThreads;

  // General configuration
  double m_pedestal;

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

  // Final fit parameters saved after possible track-piece connection.
  std::vector<double> m_tree_pad_slope;
  std::vector<double> m_tree_pad_intercept;

  std::vector<double> m_tree_tbin_slope;
  std::vector<double> m_tree_tbin_intercept;

  std::vector<double> m_tree_chi2_pad;
  std::vector<double> m_tree_chi2_tbin;

  std::vector<int> m_tree_ndof_pad;
  std::vector<int> m_tree_ndof_tbin;

  // Flat per-raw-hit/per-track content for easy TTree reading.
  std::vector<unsigned int> m_tree_hit_event;
  std::vector<unsigned int> m_tree_hit_track_id;

  std::vector<unsigned int> m_tree_hit_region;
  std::vector<unsigned int> m_tree_hit_sector;
  std::vector<int> m_tree_hit_side;

  std::vector<unsigned int> m_tree_hit_layer;

  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;

  std::vector<double> m_tree_hit_pad;
  std::vector<double> m_tree_hit_tbin;
  std::vector<double> m_tree_hit_adc;

  std::vector<double> m_tree_hit_local_phi;
  std::vector<double> m_tree_hit_local_radius;
};