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
// Per-module pthread input/output package
// One object = one TPC module: region 0..2, sector 0..11, side 0..1.
// The object owns all output from that pthread. ROOT is filled only
// after pthread_join in the main thread.
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
    unsigned short pad;
    unsigned short tbin;
    unsigned short adc;
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
  };

  struct Track
  {
    Track();
    unsigned int track_id;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nblobs;

    double pad_slope;
    double pad_intercept;
    double tbin_slope;
    double tbin_intercept;
    double chi2_pad;
    double chi2_tbin;
    int ndof_pad;
    int ndof_tbin;

    std::vector<unsigned int> blob_indices;
  };

  unsigned int region;
  unsigned int sector;
  int side;
  TrkrDefs::hitsetkey module_key;

  ActsGeometry* tGeometry;
  double pedestal;
  int verbosity;

  // Tunable parameters, equivalent in spirit to the notebook cuts.
  int blob_dt;
  int blob_dp;
  int search_dt;
  int search_dp;
  unsigned int min_track_blobs;
  double weight_power;
  double adc_weight_floor_frac;

  std::vector<LayerHitSet> layer_hitsets;
  std::vector<RawHit> raw_hits;
  std::vector<Blob> blobs;
  std::vector<Track> tracks;
};

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

  void setPedestal(double p) { m_pedestal = p; }
  void setBlobWindow(int dt, int dp) { m_blob_dt = dt; m_blob_dp = dp; }
  void setSearchWindow(int dt, int dp) { m_search_dt = dt; m_search_dp = dp; }
  void setMinTrackBlobs(unsigned int n) { m_minTrackBlobs = n; }

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

  double m_pedestal;
  int m_blob_dt;
  int m_blob_dp;
  int m_search_dt;
  int m_search_dp;
  unsigned int m_minTrackBlobs;

  int m_tree_event;

  // One entry per found module-track.
  std::vector<unsigned int> m_tree_track_id;
  std::vector<unsigned int> m_tree_region;
  std::vector<unsigned int> m_tree_sector;
  std::vector<int> m_tree_side;
  std::vector<unsigned int> m_tree_nblobs;
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

  // Flat per-blob/per-track content for easy TTree reading.
  std::vector<unsigned int> m_tree_hit_track_id;
  std::vector<unsigned int> m_tree_hit_region;
  std::vector<unsigned int> m_tree_hit_sector;
  std::vector<int> m_tree_hit_side;
  std::vector<unsigned int> m_tree_hit_layer;
  std::vector<double> m_tree_hit_pad;
  std::vector<double> m_tree_hit_tbin;
  std::vector<double> m_tree_hit_adc;
};
