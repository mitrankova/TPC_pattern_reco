#pragma once

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class FullTrackContainer;

class FullTrackVertexer : public SubsysReco
{
 public:
  explicit FullTrackVertexer(const std::string& name = "FullTrackVertexer",
                             const std::string& filename = "FullTrackVertex.root");
  ~FullTrackVertexer() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }

  // Track-quality cuts before using a full track in pair vertexing.
  void setMinTrackQuality(unsigned int min_layers,
                          unsigned int min_rawhits,
                          double max_chi2_phi_ndf,
                          double max_chi2_tbin_ndf)
  {
    m_min_layers = min_layers;
    m_min_rawhits = min_rawhits;
    m_max_chi2_phi_ndf = max_chi2_phi_ndf;
    m_max_chi2_tbin_ndf = max_chi2_tbin_ndf;
  }

  // Pair acceptance.  r limits are in cm.  This is intentionally loose by default
  // because early TPC-only straight-line fits may extrapolate roughly.
  void setPairCuts(double min_slope_diff_phi,
                   double min_slope_diff_tbin,
                   double min_vertex_r,
                   double max_vertex_r,
                   double max_pair_dtbin)
  {
    m_min_slope_diff_phi = min_slope_diff_phi;
    m_min_slope_diff_tbin = min_slope_diff_tbin;
    m_min_vertex_r = min_vertex_r;
    m_max_vertex_r = max_vertex_r;
    m_max_pair_dtbin = max_pair_dtbin;
  }

  void setOutlierCuts(double max_residual_rphi,
                      double max_residual_tbin)
  {
    m_max_residual_rphi = max_residual_rphi;
    m_max_residual_tbin = max_residual_tbin;
  }

  // Tight track-to-vertex DCA cuts used when assigning tracks to a found vertex.
  // DCA here is evaluated in the two fitted projections:
  //   rphi residual = |r_vertex * delta_phi(track, vertex)| in cm
  //   tbin residual = |tbin_track(r_vertex) - tbin_vertex| in timebins
  void setTrackVertexDcaCuts(double max_track_dca_rphi,
                             double max_track_dca_tbin)
  {
    m_max_track_dca_rphi = max_track_dca_rphi;
    m_max_track_dca_tbin = max_track_dca_tbin;
  }

  // Multiple laser/collision vertices in one event are separated mainly by drift time.
  // Pair candidates are sorted in tbin and split whenever the next candidate is farther
  // than max_timebin_gap from the running cluster center.
  void setTimebinClustering(double max_timebin_gap,
                            unsigned int min_pairs_per_vertex)
  {
    m_vertex_timebin_gap = max_timebin_gap;
    m_min_pairs_per_vertex = min_pairs_per_vertex;
  }

 private:
  struct TrackInfo
  {
    TrackInfo();

    unsigned int index;
    unsigned int track_id;
    int side;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nrawhits;

    double phi_slope;
    double phi_intercept;
    double tbin_slope;
    double tbin_intercept;
    double chi2_phi_ndf;
    double chi2_tbin_ndf;
  };

  struct PairVertex
  {
    PairVertex();

    unsigned int i;
    unsigned int j;
    int side;

    double r_phi;
    double phi;
    double x;
    double y;

    double r_tbin;
    double tbin;

    double residual_rphi;
    double residual_tbin;
    double weight;
  };

  int getNodes(PHCompositeNode* topNode);
  void reset_tree_vars();

  bool make_track_info(unsigned int index, TrackInfo& t) const;
  bool make_pair_vertex(const TrackInfo& a,
                        const TrackInfo& b,
                        PairVertex& pv) const;

  bool calculate_event_vertex(const std::vector<PairVertex>& pairs,
                              double& vx,
                              double& vy,
                              double& vr,
                              double& vphi,
                              double& vtbin,
                              double& vquality,
                              unsigned int& nused_pairs) const;

  double median(std::vector<double> v) const;
  double weighted_rms(const std::vector<double>& v,
                      const double center) const;

  std::string m_outputFileName;
  std::string m_inputNodeName;

  TFile* m_outputFile;
  TTree* m_tree;
  FullTrackContainer* m_fullTrackContainer;

  int m_event;

  unsigned int m_min_layers;
  unsigned int m_min_rawhits;
  double m_max_chi2_phi_ndf;
  double m_max_chi2_tbin_ndf;

  double m_min_slope_diff_phi;
  double m_min_slope_diff_tbin;
  double m_min_vertex_r;
  double m_max_vertex_r;
  double m_max_pair_dtbin;

  double m_max_residual_rphi;
  double m_max_residual_tbin;

  double m_max_track_dca_rphi;
  double m_max_track_dca_tbin;

  double m_vertex_timebin_gap;
  unsigned int m_min_pairs_per_vertex;

  int m_tree_event;
  unsigned int m_tree_ntracks_input;
  unsigned int m_tree_ntracks_used;
  unsigned int m_tree_npairs;
  unsigned int m_tree_npairs_used;

  int m_tree_vertex_side;
  double m_tree_vertex_x;
  double m_tree_vertex_y;
  double m_tree_vertex_r;
  double m_tree_vertex_phi;
  double m_tree_vertex_tbin;
  double m_tree_vertex_quality;

  std::vector<int> m_tree_all_vertex_side;
  std::vector<double> m_tree_all_vertex_x;
  std::vector<double> m_tree_all_vertex_y;
  std::vector<double> m_tree_all_vertex_r;
  std::vector<double> m_tree_all_vertex_phi;
  std::vector<double> m_tree_all_vertex_tbin;
  std::vector<double> m_tree_all_vertex_quality;
  std::vector<unsigned int> m_tree_all_vertex_npairs;
  std::vector<unsigned int> m_tree_all_vertex_ntracks_assigned;

  std::vector<int> m_tree_pair_side;
  std::vector<unsigned int> m_tree_pair_i;
  std::vector<unsigned int> m_tree_pair_j;
  std::vector<double> m_tree_pair_r_phi;
  std::vector<double> m_tree_pair_phi;
  std::vector<double> m_tree_pair_x;
  std::vector<double> m_tree_pair_y;
  std::vector<double> m_tree_pair_r_tbin;
  std::vector<double> m_tree_pair_tbin;
  std::vector<double> m_tree_pair_residual_rphi;
  std::vector<double> m_tree_pair_residual_tbin;
  std::vector<double> m_tree_pair_weight;
};
