#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include "TpcPadMap.h"

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class InModuleTrackContainer;
class FullTrackContainer;

class FullTrackConnector : public SubsysReco
{
 public:
  explicit FullTrackConnector(const std::string& name = "FullTrackConnector",
                              const std::string& filename = "FullTracks.root");
  ~FullTrackConnector() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }

  // Same idea as InModuleTracks::setConnectMaxLayerGap.
  void setConnectMaxLayerGap(unsigned int n) { m_connectMaxLayerGap = n; }

  // Connection window in local pad-phi [rad] and tbin.
  void setConnectWindow(double dphi, double dtbin)
  {
    m_connect_dphi = dphi;
    m_connect_dtbin = dtbin;
  }

  // Connection slope window in local pad-phi/radius and tbin/radius.
  void setConnectSlopeWindow(double dtbin_slope, double dphi_slope)
  {
    m_connect_dtbin_slope = dtbin_slope;
    m_connect_dphi_slope = dphi_slope;
  }

 public:
  struct Piece
  {
    Piece();

    unsigned int source_index;
    unsigned int source_track_id;
    unsigned int event;
    unsigned int region;
    unsigned int sector;
    int side;

    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nblobs;
    unsigned int nrawhits;

    // Fit of this source piece in local pad-phi/timebin vs local radius.
    double phi_slope;
    double phi_intercept;
    double tbin_slope;
    double tbin_intercept;

    

    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey> hitkeys;
  };

  struct Candidate
  {
    Candidate();

    unsigned int event;
    int side;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int first_sector;
    unsigned int last_sector;
    unsigned int first_region;
    unsigned int last_region;
    unsigned int nsegments;
    unsigned int nblobs;
    unsigned int nrawhits;

    double phi_slope;
    double phi_intercept;
    double tbin_slope;
    double tbin_intercept;
    double chi2_phi;
    double chi2_tbin;
    int ndof_phi;
    int ndof_tbin;

    

    std::vector<unsigned int> piece_indices;
    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey> hitkeys;
  };

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void reset_tree_vars();

  bool make_piece(unsigned int source_index, Piece& p) const;
  void connect_side_pieces(const std::vector<Piece>& pieces,
                           int side,
                           std::vector<Candidate>& output) const;
  bool refit_candidate(const std::vector<Piece>& pieces,
                       const std::vector<unsigned int>& piece_indices,
                       Candidate& c) const;
  bool candidates_can_connect(const Candidate& a,
                              const Piece& b,
                              double& score,
                              double& b_phi_intercept_shifted) const;

  std::string m_outputFileName;
  std::string m_inputNodeName;
  std::string m_outputNodeName;

  TFile* m_outputFile;
  TTree* m_tree;

  InModuleTrackContainer* m_inModuleTrackContainer;
  FullTrackContainer* m_fullTrackContainer;

  int m_event;

  TpcPadMap m_padMap;

  unsigned int m_connectMaxLayerGap;
  double m_connect_dphi;
  double m_connect_dtbin;
  double m_connect_dphi_slope;
  double m_connect_dtbin_slope;

  int m_tree_event;
  std::vector<unsigned int> m_tree_track_id;
  std::vector<int> m_tree_side;
  std::vector<unsigned int> m_tree_nsegments;
  std::vector<unsigned int> m_tree_nrawhits;
  std::vector<unsigned int> m_tree_nblobs;
  std::vector<unsigned int> m_tree_first_layer;
  std::vector<unsigned int> m_tree_last_layer;
  std::vector<unsigned int> m_tree_first_sector;
  std::vector<unsigned int> m_tree_last_sector;
  std::vector<unsigned int> m_tree_first_region;
  std::vector<unsigned int> m_tree_last_region;
  std::vector<double> m_tree_phi_slope;
  std::vector<double> m_tree_phi_intercept;
  std::vector<double> m_tree_tbin_slope;
  std::vector<double> m_tree_tbin_intercept;
  std::vector<double> m_tree_chi2_phi;
  std::vector<double> m_tree_chi2_tbin;
  std::vector<int> m_tree_ndof_phi;
  std::vector<int> m_tree_ndof_tbin;

  std::vector<unsigned int> m_tree_source_full_track_id;
  std::vector<unsigned int> m_tree_source_inmodule_track_id;
  std::vector<unsigned int> m_tree_source_region;
  std::vector<unsigned int> m_tree_source_sector;
  std::vector<int> m_tree_source_side;

  std::vector<unsigned int> m_tree_hit_full_track_id;
  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};
