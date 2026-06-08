#include "FullTrackConnector.h"

#include "FullTrack.h"
#include "FullTrackv1.h"
#include "FullTrackContainer.h"
#include "FullTrackContainerv1.h"
#include "TpcPadMap.h"

#include "InModuleTrack.h"
#include "InModuleTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

// ===================================================================
// Internal helpers
// ===================================================================
namespace
{
  bool weighted_line_fit(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& w,
                         double& m, double& b,
                         double& chi2, int& ndof)
  {
    if (x.size() < 2 ||
        x.size() != y.size() ||
        x.size() != w.size()) return false;

    double S   = 0.0, Sx  = 0.0, Sy  = 0.0;
    double Sxx = 0.0, Sxy = 0.0;

    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = (w[i] > 0.0) ? w[i] : 1.0;
      S   += wi;
      Sx  += wi * x[i];
      Sy  += wi * y[i];
      Sxx += wi * x[i] * x[i];
      Sxy += wi * x[i] * y[i];
    }

    const double den = S * Sxx - Sx * Sx;
    if (std::fabs(den) < 1.0e-14) return false;

    m = (S * Sxy - Sx * Sy) / den;
    b = (Sy - m * Sx) / S;

    chi2 = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = (w[i] > 0.0) ? w[i] : 1.0;
      const double r  = y[i] - (m * x[i] + b);
      chi2 += wi * r * r;
    }

    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }

  struct PieceStartSort
  {
    const std::vector<FullTrackConnector::Piece>* pieces;
    explicit PieceStartSort(const std::vector<FullTrackConnector::Piece>* p)
      : pieces(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const FullTrackConnector::Piece& pa = (*pieces)[a];
      const FullTrackConnector::Piece& pb = (*pieces)[b];
      if (pa.first_layer != pb.first_layer) return pa.first_layer < pb.first_layer;
      if (pa.last_layer  != pb.last_layer)  return pa.last_layer  < pb.last_layer;
      if (pa.sector      != pb.sector)      return pa.sector      < pb.sector;
      return pa.source_track_id < pb.source_track_id;
    }
  };
} // anonymous namespace


// ===================================================================
// Struct constructors
// ===================================================================
FullTrackConnector::Piece::Piece()
  : source_index(0), source_track_id(0), event(0),
    region(0), sector(0), side(0),
    first_layer(0), last_layer(0), nblobs(0), nrawhits(0),
    pad_slope(0.0), pad_intercept(0.0),
    tbin_slope(0.0), tbin_intercept(0.0),
    npads_region(0)
{}

FullTrackConnector::Candidate::Candidate()
  : event(0), side(0),
    first_layer(0), last_layer(0),
    first_sector(0), last_sector(0),
    first_region(0), last_region(0),
    nsegments(0), nblobs(0), nrawhits(0),
    pad_slope_r(0.0), pad_intercept_r(0.0),
    tbin_slope_r(0.0), tbin_intercept_r(0.0),
    chi2_pad(0.0), chi2_tbin(0.0),
    ndof_pad(0), ndof_tbin(0)
{}


// ===================================================================
// Constructor / destructor
// ===================================================================
FullTrackConnector::FullTrackConnector(const std::string& name,
                                       const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_inputNodeName("INMODULETRACKS")
  , m_outputNodeName("FULLTRACKS")
  , m_outputFile(nullptr)
  , m_tree(nullptr)
  , m_inModuleTrackContainer(nullptr)
  , m_fullTrackContainer(nullptr)
  , m_event(0)
  , m_padMap(nullptr)
  , m_connectMaxLayerGap(1)
  , m_connect_dpad(6.0)      // ~same physical window as old 0.03 rad, typical pitch
  , m_connect_dtbin(8.0)
  , m_connect_dpad_slope(2.0)
  , m_connect_dtbin_slope(2.0)
{}

FullTrackConnector::~FullTrackConnector()
{
  if (m_outputFile) { delete m_outputFile; m_outputFile = nullptr; }
}


// ===================================================================
// Init / InitRun / End
// ===================================================================
int FullTrackConnector::Init(PHCompositeNode*)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("FullTracks", "Full tracks connected from InModuleTracks");
  m_tree->Branch("event",          &m_tree_event,          "event/I");
  m_tree->Branch("track_id",       &m_tree_track_id);
  m_tree->Branch("side",           &m_tree_side);
  m_tree->Branch("nsegments",      &m_tree_nsegments);
  m_tree->Branch("nblobs",         &m_tree_nblobs);
  m_tree->Branch("nrawhits",       &m_tree_nrawhits);
  m_tree->Branch("first_layer",    &m_tree_first_layer);
  m_tree->Branch("last_layer",     &m_tree_last_layer);
  m_tree->Branch("first_sector",   &m_tree_first_sector);
  m_tree->Branch("last_sector",    &m_tree_last_sector);
  m_tree->Branch("first_region",   &m_tree_first_region);
  m_tree->Branch("last_region",    &m_tree_last_region);
  // Global fit stored in (radius, pad_eff) and (radius, tbin) space.
  m_tree->Branch("pad_slope_r",      &m_tree_pad_slope_r);
  m_tree->Branch("pad_intercept_r",  &m_tree_pad_intercept_r);
  m_tree->Branch("tbin_slope_r",     &m_tree_tbin_slope_r);
  m_tree->Branch("tbin_intercept_r", &m_tree_tbin_intercept_r);
  m_tree->Branch("chi2_pad",         &m_tree_chi2_pad);
  m_tree->Branch("chi2_tbin",        &m_tree_chi2_tbin);
  m_tree->Branch("ndof_pad",         &m_tree_ndof_pad);
  m_tree->Branch("ndof_tbin",        &m_tree_ndof_tbin);

  m_tree->Branch("source_full_track_id",      &m_tree_source_full_track_id);
  m_tree->Branch("source_inmodule_track_id",  &m_tree_source_inmodule_track_id);
  m_tree->Branch("source_region",             &m_tree_source_region);
  m_tree->Branch("source_sector",             &m_tree_source_sector);
  m_tree->Branch("source_side",               &m_tree_source_side);

  m_tree->Branch("hit_full_track_id", &m_tree_hit_full_track_id);
  m_tree->Branch("hit_hitsetkey",     &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey",        &m_tree_hit_hitkey);

  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackConnector::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode)    != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  if (!m_padMap || !m_padMap->isValid())
  {
    std::cerr << Name() << "::InitRun - TPC_PADMAP missing or invalid" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_padMap->print_sector_edges();

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackConnector::End(PHCompositeNode*)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (m_tree) m_tree->Write();
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


// ===================================================================
// Node management
// ===================================================================
int FullTrackConnector::getNodes(PHCompositeNode* topNode)
{
  m_inModuleTrackContainer =
    findNode::getClass<InModuleTrackContainer>(topNode, m_inputNodeName);
  if (!m_inModuleTrackContainer)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_padMap = findNode::getClass<TpcPadMap>(topNode, "TPC_PADMAP");
  if (!m_padMap)
  {
    std::cerr << Name() << "::getNodes - missing TPC_PADMAP" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackConnector::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode =
    dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_fullTrackContainer =
    findNode::getClass<FullTrackContainer>(topNode, m_outputNodeName);
  if (!m_fullTrackContainer)
  {
    m_fullTrackContainer = new FullTrackContainerv1();
    PHIODataNode<PHObject>* node =
      new PHIODataNode<PHObject>(m_fullTrackContainer, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


// ===================================================================
// Helpers: pad-count rescaling
// ===================================================================
unsigned int FullTrackConnector::npads_for_region(unsigned int region) const
{
  // TpcPadMap::get_pads_per_sector(region) returns the number of pads in that region.
  // We cache the value in Piece::npads_region to avoid repeated calls.
  return static_cast<unsigned int>(m_padMap->get_pads_per_sector(region));
}

// Convert a pad number from src_region to equivalent pad in dst_region.
// Physical phi is conserved:  pad_src / N_src == pad_dst / N_dst
// so:                         pad_dst = pad_src * N_dst / N_src
double FullTrackConnector::rescale_pad(double pad,
                                       unsigned int src_region,
                                       unsigned int dst_region) const
{
  if (src_region == dst_region) return pad;
  const double n_src = static_cast<double>(npads_for_region(src_region));
  const double n_dst = static_cast<double>(npads_for_region(dst_region));
  if (n_src <= 0.0) return pad;
  return pad * n_dst / n_src;
}


// ===================================================================
// make_piece: translate one InModuleTrack → Piece (hardware coords)
// ===================================================================
bool FullTrackConnector::make_piece(unsigned int source_index, Piece& p) const
{
  const InModuleTrack* trk = m_inModuleTrackContainer->get_track(source_index);
  if (!trk || !trk->isValid()) return false;
  if (trk->get_last_layer() < trk->get_first_layer()) return false;

  p.source_index    = source_index;
  p.source_track_id = trk->get_track_id();
  p.event           = trk->get_event();
  p.region          = trk->get_region();
  p.sector          = trk->get_sector();
  p.side            = trk->get_side();
  p.first_layer     = trk->get_first_layer();
  p.last_layer      = trk->get_last_layer();
  p.nblobs          = trk->get_nblobs();
  p.nrawhits        = trk->get_nrawhits();

  // The InModuleTracks fit is already in (layer, pad) and (layer, tbin) —
  // just copy it directly.
  p.pad_slope      = 0;//trk->get_pad_slope();
  p.pad_intercept  = 0; //trk->get_pad_intercept();
  p.tbin_slope     = 0; // trk->get_tbin_slope();
  p.tbin_intercept = 0; // trk->get_tbin_intercept();

  p.npads_region = npads_for_region(p.region);

  p.hitsetkeys.clear();
  p.hitkeys.clear();
  for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
  {
    const InModuleTrack::HitIndex hi = trk->get_hit_index(ih);
    p.hitsetkeys.push_back(hi.first);
    p.hitkeys.push_back(hi.second);
  }

  return true;
}


// ===================================================================
// refit_candidate
//
// Build a single straight-line fit in (radius, pad_eff) and
// (radius, tbin) from all pieces.
//
// pad_eff normalises every pad to "region-0 equivalent" units so that
// the line is continuous across region boundaries:
//   pad_eff = pad * N_region0 / N_region_i
//
// Using radius (cm) as the independent variable removes the
// nonlinearity of layer→R from the fit entirely.
// ===================================================================
bool FullTrackConnector::refit_candidate(const std::vector<Piece>& pieces,
                                          const std::vector<unsigned int>& piece_indices,
                                          Candidate& c) const
{
  if (piece_indices.empty()) return false;

  // Reference pad count: normalise everything to region 0.
  const unsigned int npads_ref =
    static_cast<unsigned int>(m_padMap->get_pads_per_sector(0));
  if (npads_ref == 0) return false;
  const double npads_ref_d = static_cast<double>(npads_ref);

  std::vector<double> r_vec;
  std::vector<double> pad_eff_vec;
  std::vector<double> tbin_vec;
  std::vector<double> w_vec;

  r_vec.reserve(piece_indices.size() * 2);
  pad_eff_vec.reserve(piece_indices.size() * 2);
  tbin_vec.reserve(piece_indices.size() * 2);
  w_vec.reserve(piece_indices.size() * 2);

  c = Candidate();
  c.piece_indices = piece_indices;

  for (unsigned int ii = 0; ii < piece_indices.size(); ++ii)
  {
    const Piece& p = pieces[piece_indices[ii]];
    const double wt = (p.nrawhits > 0) ? static_cast<double>(p.nrawhits) : 1.0;

    const double r0 = m_padMap->get_local_radius(p.region, p.first_layer);
    const double r1 = m_padMap->get_local_radius(p.region, p.last_layer);
    if (r0 <= 0.0 || r1 <= 0.0) return false;

    const double l0 = static_cast<double>(p.first_layer);
    const double l1 = static_cast<double>(p.last_layer);

    const double pad0  = p.pad_slope * l0 + p.pad_intercept;
    const double pad1  = p.pad_slope * l1 + p.pad_intercept;
    const double tbin0 = p.tbin_slope * l0 + p.tbin_intercept;
    const double tbin1 = p.tbin_slope * l1 + p.tbin_intercept;

    // Rescale pad to region-0 equivalent units.
    const double scale =
      npads_ref_d / static_cast<double>(p.npads_region > 0 ? p.npads_region : npads_ref);
    const double peff0 = pad0 * scale;
    const double peff1 = pad1 * scale;

    r_vec.push_back(r0);
    pad_eff_vec.push_back(peff0);
    tbin_vec.push_back(tbin0);
    w_vec.push_back(wt);

    if (p.last_layer != p.first_layer)
    {
      r_vec.push_back(r1);
      pad_eff_vec.push_back(peff1);
      tbin_vec.push_back(tbin1);
      w_vec.push_back(wt);
    }

    // Accumulate bookkeeping.
    if (ii == 0)
    {
      c.event        = p.event;
      c.side         = p.side;
      c.first_layer  = p.first_layer;
      c.last_layer   = p.last_layer;
      c.first_sector = p.sector;
      c.last_sector  = p.sector;
      c.first_region = p.region;
      c.last_region  = p.region;
    }
    else
    {
      if (p.first_layer < c.first_layer)
      {
        c.first_layer  = p.first_layer;
        c.first_sector = p.sector;
        c.first_region = p.region;
      }
      if (p.last_layer > c.last_layer)
      {
        c.last_layer  = p.last_layer;
        c.last_sector = p.sector;
        c.last_region = p.region;
      }
    }

    c.nblobs    += p.nblobs;
    c.nrawhits  += p.nrawhits;
    for (unsigned int ih = 0; ih < p.hitsetkeys.size(); ++ih)
    {
      c.hitsetkeys.push_back(p.hitsetkeys[ih]);
      c.hitkeys.push_back(p.hitkeys[ih]);
    }
  }

  c.nsegments = static_cast<unsigned int>(piece_indices.size());

  if (!weighted_line_fit(r_vec, pad_eff_vec, w_vec,
                         c.pad_slope_r, c.pad_intercept_r,
                         c.chi2_pad, c.ndof_pad))  return false;
  if (!weighted_line_fit(r_vec, tbin_vec, w_vec,
                         c.tbin_slope_r, c.tbin_intercept_r,
                         c.chi2_tbin, c.ndof_tbin)) return false;

  return true;
}


// ===================================================================
// candidates_can_connect
//
// Everything stays in hardware units.  At the boundary between two
// regions the pad count differs, so we rescale b's pad prediction to
// a's pad-count convention before comparing.
//
// The match point is the radius mid-way between a.last_layer and
// b.first_layer (using the actual radius lookup, so the nonlinear
// layer→R is handled exactly).
// ===================================================================
bool FullTrackConnector::candidates_can_connect(const Candidate& a,
                                                 const Piece&     b,
                                                 double&          score) const
{
  score = std::numeric_limits<double>::max();

  // Sides must match; b must come after a; gap must be within limit.
  if (a.side != b.side) return false;
  if (a.last_layer >= b.first_layer) return false;

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap) return false;

  // Radii at the boundary layers (exact lookup → no linearity assumption).
  const double ra = m_padMap->get_local_radius(a.last_region,  a.last_layer);
  const double rb = m_padMap->get_local_radius(b.region,       b.first_layer);
  if (ra <= 0.0 || rb <= 0.0) return false;

  // Evaluate candidate a's global fit at the match radius.
  // The candidate fit is in pad_eff (region-0 units) vs radius.
  const double npads_ref_d =
    static_cast<double>(m_padMap->get_pads_per_sector(0));
  const double rmatch = 0.5 * (ra + rb);

  const double pad_eff_a = a.pad_slope_r * rmatch + a.pad_intercept_r;
  const double tbin_a    = a.tbin_slope_r * rmatch + a.tbin_intercept_r;

  // Evaluate piece b at its first layer using its hardware-space fit,
  // then rescale b's pad to the same region-0 equivalent units.
  const double lb     = static_cast<double>(b.first_layer);
  const double pad_b_raw  = b.pad_slope * lb + b.pad_intercept;
  const double tbin_b     = b.tbin_slope * lb + b.tbin_intercept;

  // Rescale b's pad to region-0 units for apples-to-apples comparison.
  const double n_b = static_cast<double>(b.npads_region > 0 ? b.npads_region
                                                             : m_padMap->get_pads_per_sector(b.region));
  const double pad_eff_b = (n_b > 0.0) ? pad_b_raw * npads_ref_d / n_b : pad_b_raw;

  // Slope of b, also rescaled to pad_eff units.
  const double pad_slope_eff_b = (n_b > 0.0) ? b.pad_slope * npads_ref_d / n_b
                                              : b.pad_slope;

  // Candidate a's slope in the same units (already in pad_eff/radius).
  // Convert a's pad_eff/radius slope to pad_eff/layer at the boundary
  // for a consistent comparison with b's pad/layer slope.
  // dr/dlayer is approximated from the two boundary radii and layer gap;
  // for a 1-layer gap this is just (rb - ra).
  // We compare slopes in (pad_eff / radius) directly — both are now in
  // region-0 pad units per cm, which is well-defined.
  const double dphi_a = a.pad_slope_r;          // d(pad_eff)/d(radius) for candidate
  // Convert b's d(pad_eff)/d(layer) to d(pad_eff)/d(radius) using the
  // local layer spacing at the boundary.
  // Approximate dr/dlayer: use radius difference over 1 layer if gap == 0,
  // else use (rb - ra) / (b.first_layer - a.last_layer).
  const double dlayer = static_cast<double>(b.first_layer) -
                        static_cast<double>(a.last_layer);
  const double drperlayer = (std::fabs(rb - ra) > 1.0e-6 && dlayer > 0.0)
                            ? (rb - ra) / dlayer : 1.0;
  const double dphi_b = (std::fabs(drperlayer) > 1.0e-9)
                        ? pad_slope_eff_b / drperlayer : 0.0;

  const double dtbin_slope_b = b.tbin_slope;     // tbin/layer for b
  // Convert a's d(tbin)/d(radius) → d(tbin)/d(layer) at boundary.
  const double dtbin_a_perlayer = a.tbin_slope_r * drperlayer;

  const double dpad   = std::fabs(pad_eff_a  - pad_eff_b);
  const double dtbin  = std::fabs(tbin_a     - tbin_b);
  const double dmphi  = std::fabs(dphi_a     - dphi_b);
  const double dmtbin = std::fabs(dtbin_a_perlayer - dtbin_slope_b);

  if (dpad   > m_connect_dpad)        return false;
  if (dtbin  > m_connect_dtbin)       return false;
  if (dmphi  > m_connect_dpad_slope)  return false;
  if (dmtbin > m_connect_dtbin_slope) return false;

  score = (dpad   / m_connect_dpad)        * (dpad   / m_connect_dpad)
        + (dtbin  / m_connect_dtbin)       * (dtbin  / m_connect_dtbin)
        + (dmphi  / m_connect_dpad_slope)  * (dmphi  / m_connect_dpad_slope)
        + (dmtbin / m_connect_dtbin_slope) * (dmtbin / m_connect_dtbin_slope)
        + 0.05 * static_cast<double>(gap);

  return true;
}


// ===================================================================
// connect_side_pieces  — greedy stitching, one TPC side at a time
// ===================================================================
void FullTrackConnector::connect_side_pieces(const std::vector<Piece>& pieces,
                                              int side,
                                              std::vector<Candidate>& output) const
{
  // Collect indices belonging to this side and sort by first_layer.
  std::vector<unsigned int> order;
  order.reserve(pieces.size());
  for (unsigned int i = 0; i < pieces.size(); ++i)
    if (pieces[i].side == side) order.push_back(i);

  if (order.empty()) return;

  std::sort(order.begin(), order.end(), PieceStartSort(&pieces));
  std::vector<int> used(pieces.size(), 0);

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int iseed = order[io];
    if (used[iseed]) continue;

    // Seed a new candidate from this piece.
    std::vector<unsigned int> current_indices;
    current_indices.push_back(iseed);
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current)) continue;

    // Greedy forward extension.
    bool merged_any = true;
    while (merged_any)
    {
      merged_any  = false;
      int    best_j     = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int jo = 0; jo < order.size(); ++jo)
      {
        const unsigned int j = order[jo];
        if (used[j]) continue;

        double score = 0.0;
        if (!candidates_can_connect(current, pieces[j], score)) continue;

        if (score < best_score) { best_score = score; best_j = static_cast<int>(j); }
      }

      if (best_j >= 0)
      {
        std::vector<unsigned int> trial_indices = current_indices;
        trial_indices.push_back(static_cast<unsigned int>(best_j));

        Candidate refit;
        if (refit_candidate(pieces, trial_indices, refit))
        {
          current = refit;
          current_indices.swap(trial_indices);
          used[best_j] = 1;
          merged_any   = true;
        }
      }
    }

    output.push_back(current);
  }
}


// ===================================================================
// reset_tree_vars
// ===================================================================
void FullTrackConnector::reset_tree_vars()
{
  m_tree_event = m_event;
  m_tree_track_id.clear();
  m_tree_side.clear();
  m_tree_nsegments.clear();
  m_tree_nblobs.clear();
  m_tree_nrawhits.clear();
  m_tree_first_layer.clear();
  m_tree_last_layer.clear();
  m_tree_first_sector.clear();
  m_tree_last_sector.clear();
  m_tree_first_region.clear();
  m_tree_last_region.clear();
  m_tree_pad_slope_r.clear();
  m_tree_pad_intercept_r.clear();
  m_tree_tbin_slope_r.clear();
  m_tree_tbin_intercept_r.clear();
  m_tree_chi2_pad.clear();
  m_tree_chi2_tbin.clear();
  m_tree_ndof_pad.clear();
  m_tree_ndof_tbin.clear();

  m_tree_source_full_track_id.clear();
  m_tree_source_inmodule_track_id.clear();
  m_tree_source_region.clear();
  m_tree_source_sector.clear();
  m_tree_source_side.clear();

  m_tree_hit_full_track_id.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
}


// ===================================================================
// process_event
// ===================================================================
int FullTrackConnector::process_event(PHCompositeNode*)
{
  reset_tree_vars();
  if (m_fullTrackContainer) m_fullTrackContainer->Reset();

  // Build Piece list from all InModuleTracks.
  std::vector<Piece> pieces;
  if (m_inModuleTrackContainer)
  {
    const unsigned int n = m_inModuleTrackContainer->size();
    pieces.reserve(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      Piece p;
      if (make_piece(i, p)) pieces.push_back(p);
    }
  }

  // Stitch each side independently.
  std::vector<Candidate> full_tracks;
  connect_side_pieces(pieces, 0, full_tracks);
  connect_side_pieces(pieces, 1, full_tracks);

  // Store results.
  for (unsigned int it = 0; it < full_tracks.size(); ++it)
  {
    const Candidate& c = full_tracks[it];
    const unsigned int full_id =
      m_fullTrackContainer ? m_fullTrackContainer->size() : it;

    FullTrackv1* out = new FullTrackv1();
    out->set_event(static_cast<unsigned int>(m_event));
    out->set_track_id(full_id);
    out->set_side(c.side);
    out->set_nsegments(c.nsegments);
    out->set_nblobs(c.nblobs);
    out->set_nrawhits(c.nrawhits);
    out->set_first_layer(c.first_layer);
    out->set_last_layer(c.last_layer);
    out->set_first_sector(c.first_sector);
    out->set_last_sector(c.last_sector);
    out->set_first_region(c.first_region);
    out->set_last_region(c.last_region);
    // Store the hardware-space global fit.
    // FullTrackv1 setters: reuse phi_slope/phi_intercept for pad_slope_r/pad_intercept_r.
    // These are stored in the ROOT output via the TTree branches above; the
    // FullTrackv1 object uses whatever setters are available.
    out->set_phi_slope(c.pad_slope_r);
    out->set_phi_intercept(c.pad_intercept_r);
    out->set_tbin_slope(c.tbin_slope_r);
    out->set_tbin_intercept(c.tbin_intercept_r);
    out->set_chi2_phi(c.chi2_pad);
    out->set_chi2_tbin(c.chi2_tbin);
    out->set_ndof_phi(c.ndof_pad);
    out->set_ndof_tbin(c.ndof_tbin);

    m_tree_track_id.push_back(full_id);
    m_tree_side.push_back(c.side);
    m_tree_nsegments.push_back(c.nsegments);
    m_tree_nblobs.push_back(c.nblobs);
    m_tree_nrawhits.push_back(c.nrawhits);
    m_tree_first_layer.push_back(c.first_layer);
    m_tree_last_layer.push_back(c.last_layer);
    m_tree_first_sector.push_back(c.first_sector);
    m_tree_last_sector.push_back(c.last_sector);
    m_tree_first_region.push_back(c.first_region);
    m_tree_last_region.push_back(c.last_region);
    m_tree_pad_slope_r.push_back(c.pad_slope_r);
    m_tree_pad_intercept_r.push_back(c.pad_intercept_r);
    m_tree_tbin_slope_r.push_back(c.tbin_slope_r);
    m_tree_tbin_intercept_r.push_back(c.tbin_intercept_r);
    m_tree_chi2_pad.push_back(c.chi2_pad);
    m_tree_chi2_tbin.push_back(c.chi2_tbin);
    m_tree_ndof_pad.push_back(c.ndof_pad);
    m_tree_ndof_tbin.push_back(c.ndof_tbin);

    for (unsigned int ip = 0; ip < c.piece_indices.size(); ++ip)
    {
      const Piece& p = pieces[c.piece_indices[ip]];
      out->add_source_track(p.source_track_id, p.region, p.sector);

      m_tree_source_full_track_id.push_back(full_id);
      m_tree_source_inmodule_track_id.push_back(p.source_track_id);
      m_tree_source_region.push_back(p.region);
      m_tree_source_sector.push_back(p.sector);
      m_tree_source_side.push_back(p.side);
    }

    for (unsigned int ih = 0; ih < c.hitsetkeys.size(); ++ih)
    {
      out->add_hit_index(c.hitsetkeys[ih], c.hitkeys[ih]);
      m_tree_hit_full_track_id.push_back(full_id);
      m_tree_hit_hitsetkey.push_back(
        static_cast<unsigned long long>(c.hitsetkeys[ih]));
      m_tree_hit_hitkey.push_back(
        static_cast<unsigned long long>(c.hitkeys[ih]));
    }

    if (m_fullTrackContainer) m_fullTrackContainer->add_track(out);
    else delete out;
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << "  input pieces=" << pieces.size()
              << "  full_tracks="  << full_tracks.size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}