#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

class TpcPadMap;

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

  void setInputNodeName(const std::string& n)  { m_inputNodeName  = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }

  void setConnectMaxLayerGap(unsigned int n)   { m_connectMaxLayerGap = n; }

  // Connection window in pad units and tbin units at the match layer.
  // "pad units" here means pads of the *source* (lower-layer) region.
  void setConnectWindow(double dpad, double dtbin)
  {
    m_connect_dpad  = dpad;
    m_connect_dtbin = dtbin;
  }

  // Slope windows: d(pad)/d(layer) and d(tbin)/d(layer).
  // Slopes from adjacent regions are rescaled to the same pad-pitch before
  // comparison, so these windows are in consistent units.
  void setConnectSlopeWindow(double dpad_slope, double dtbin_slope)
  {
    m_connect_dpad_slope  = dpad_slope;
    m_connect_dtbin_slope = dtbin_slope;
  }

 public:
  // ---------------------------------------------------------------
  // Piece: one InModuleTrack kept in raw hardware coordinates.
  // Fit parameterisation:   pad(layer)  = pad_slope  * layer + pad_intercept
  //                         tbin(layer) = tbin_slope * layer + tbin_intercept
  // ---------------------------------------------------------------
  struct Piece
  {
    Piece();

    unsigned int source_index;
    unsigned int source_track_id;
    unsigned int event;
    unsigned int region;
    unsigned int sector;
    int          side;

    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nblobs;
    unsigned int nrawhits;

    // Hardware-space linear fit (same units as InModuleTrack output).
    double pad_slope;
    double pad_intercept;
    double tbin_slope;
    double tbin_intercept;

    // Number of pads in this piece's region — cached for boundary rescaling.
    unsigned int npads_region;

    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey>    hitkeys;
  };

  // ---------------------------------------------------------------
  // Candidate: a growing full track built from one or more Pieces.
  // The global fit is also kept in hardware space using radius as the
  // independent variable so that it is well-defined across regions.
  //   pad_eff(radius) = pad_slope_r * radius + pad_intercept_r
  //   tbin   (radius) = tbin_slope_r * radius + tbin_intercept_r
  // where pad_eff is pad rescaled to "region-0 equivalent" units so
  // that a single line can span all three regions.
  // ---------------------------------------------------------------
  struct Candidate
  {
    Candidate();

    unsigned int event;
    int          side;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int first_sector;
    unsigned int last_sector;
    unsigned int first_region;
    unsigned int last_region;
    unsigned int nsegments;
    unsigned int nblobs;
    unsigned int nrawhits;

    // Global fit in (radius, pad_eff) and (radius, tbin).
    double pad_slope_r;       // d(pad_eff)/d(radius)
    double pad_intercept_r;
    double tbin_slope_r;      // d(tbin)/d(radius)
    double tbin_intercept_r;
    double chi2_pad;
    double chi2_tbin;
    int    ndof_pad;
    int    ndof_tbin;

    std::vector<unsigned int>        piece_indices;
    std::vector<TrkrDefs::hitsetkey> hitsetkeys;
    std::vector<TrkrDefs::hitkey>    hitkeys;
  };

 private:
  int  getNodes(PHCompositeNode*);
  int  createNodes(PHCompositeNode*);
  void reset_tree_vars();

  // Pad-count for region (thin wrapper so we don't scatter magic numbers).
  unsigned int npads_for_region(unsigned int region) const;

  // Rescale pad from region src to the equivalent pad in region dst.
  // This is a simple linear scale: same physical phi → same fraction of pads.
  double rescale_pad(double pad, unsigned int src_region,
                     unsigned int dst_region) const;

  bool make_piece(unsigned int source_index, Piece& p) const;

  // Refit a growing candidate (radius-space line through all piece endpoints).
  bool refit_candidate(const std::vector<Piece>& pieces,
                       const std::vector<unsigned int>& piece_indices,
                       Candidate& c) const;

  // Check whether piece b can be appended to candidate a.
  // Returns true and fills score.  No phi conversion — pure hardware space.
  bool candidates_can_connect(const Candidate& a,
                               const Piece&     b,
                               double&          score) const;

  void connect_side_pieces(const std::vector<Piece>& pieces,
                            int side,
                            std::vector<Candidate>& output) const;

  // ---------------------------------------------------------------
  std::string m_outputFileName;
  std::string m_inputNodeName;
  std::string m_outputNodeName;

  TFile* m_outputFile;
  TTree* m_tree;

  InModuleTrackContainer* m_inModuleTrackContainer;
  FullTrackContainer*     m_fullTrackContainer;

  int        m_event;
  TpcPadMap* m_padMap;

  unsigned int m_connectMaxLayerGap;
  double       m_connect_dpad;        // pad window at match point (src region units)
  double       m_connect_dtbin;
  double       m_connect_dpad_slope;  // d(pad_eff)/d(layer) slope window
  double       m_connect_dtbin_slope;

  // ---------------------------------------------------------------
  // TTree branches
  // ---------------------------------------------------------------
  int                        m_tree_event;
  std::vector<unsigned int>  m_tree_track_id;
  std::vector<int>           m_tree_side;
  std::vector<unsigned int>  m_tree_nsegments;
  std::vector<unsigned int>  m_tree_nblobs;
  std::vector<unsigned int>  m_tree_nrawhits;
  std::vector<unsigned int>  m_tree_first_layer;
  std::vector<unsigned int>  m_tree_last_layer;
  std::vector<unsigned int>  m_tree_first_sector;
  std::vector<unsigned int>  m_tree_last_sector;
  std::vector<unsigned int>  m_tree_first_region;
  std::vector<unsigned int>  m_tree_last_region;
  std::vector<double>        m_tree_pad_slope_r;
  std::vector<double>        m_tree_pad_intercept_r;
  std::vector<double>        m_tree_tbin_slope_r;
  std::vector<double>        m_tree_tbin_intercept_r;
  std::vector<double>        m_tree_chi2_pad;
  std::vector<double>        m_tree_chi2_tbin;
  std::vector<int>           m_tree_ndof_pad;
  std::vector<int>           m_tree_ndof_tbin;

  std::vector<unsigned int>  m_tree_source_full_track_id;
  std::vector<unsigned int>  m_tree_source_inmodule_track_id;
  std::vector<unsigned int>  m_tree_source_region;
  std::vector<unsigned int>  m_tree_source_sector;
  std::vector<int>           m_tree_source_side;

  std::vector<unsigned int>       m_tree_hit_full_track_id;
  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};