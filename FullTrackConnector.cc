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
#include <map>
#include <utility>

namespace
{
  bool weighted_line_fit(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& w,
                         double& m,
                         double& b,
                         double& chi2,
                         int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size() || x.size() != w.size()) return false;

    double S = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = (w[i] > 0.0) ? w[i] : 1.0;
      S += wi;
      Sx += wi * x[i];
      Sy += wi * y[i];
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
      const double r = y[i] - (m * x[i] + b);
      chi2 += wi * r * r;
    }

    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }

  struct PieceStartSort
  {
    const std::vector<FullTrackConnector::Piece>* pieces;
    PieceStartSort(const std::vector<FullTrackConnector::Piece>* p) : pieces(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const FullTrackConnector::Piece& pa = (*pieces)[a];
      const FullTrackConnector::Piece& pb = (*pieces)[b];
      if (pa.first_layer != pb.first_layer) return pa.first_layer < pb.first_layer;
      if (pa.last_layer != pb.last_layer) return pa.last_layer < pb.last_layer;
      if (pa.sector != pb.sector) return pa.sector < pb.sector;
      return pa.source_track_id < pb.source_track_id;
    }
  };

  double unwrap_phi_to_reference(double phi, const double ref)
  {
    while (phi - ref > M_PI)  phi -= 2.0 * M_PI;
    while (phi - ref < -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double wrap_to_pi(double phi)
  {
    while (phi > M_PI)  phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }
}


FullTrackConnector::Piece::Piece() :
  source_index(0),
  source_track_id(0),
  event(0),
  region(0),
  sector(0),
  side(0),
  first_layer(0),
  last_layer(0),
  nblobs(0),
  nrawhits(0),
  phi_slope(0.0),
  phi_intercept(0.0),
  tbin_slope(0.0),
  tbin_intercept(0.0)
{
}

FullTrackConnector::Candidate::Candidate() :
  event(0),
  side(0),
  first_layer(0),
  last_layer(0),
  first_sector(0),
  last_sector(0),
  first_region(0),
  last_region(0),
  nsegments(0),
  nblobs(0),
  nrawhits(0),
  phi_slope(0.0),
  phi_intercept(0.0),
  tbin_slope(0.0),
  tbin_intercept(0.0),
  chi2_phi(0.0),
  chi2_tbin(0.0),
  ndof_phi(0),
  ndof_tbin(0)
{
}

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
  , m_connect_dphi(0.03)
  , m_connect_dtbin(8.0)
  , m_connect_dphi_slope(0.01)
  , m_connect_dtbin_slope(2.0)
{
}

FullTrackConnector::~FullTrackConnector()
{
  if (m_outputFile)
  {
    delete m_outputFile;
    m_outputFile = nullptr;
  }
}

int FullTrackConnector::Init(PHCompositeNode*)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("FullTracks", "Full tracks connected from InModuleTracks");
  m_tree->Branch("event", &m_tree_event, "event/I");
  m_tree->Branch("track_id", &m_tree_track_id);
  m_tree->Branch("side", &m_tree_side);
  m_tree->Branch("nsegments", &m_tree_nsegments);
  m_tree->Branch("nblobs", &m_tree_nblobs);
  m_tree->Branch("nrawhits", &m_tree_nrawhits);
  m_tree->Branch("first_layer", &m_tree_first_layer);
  m_tree->Branch("last_layer", &m_tree_last_layer);
  m_tree->Branch("first_sector", &m_tree_first_sector);
  m_tree->Branch("last_sector", &m_tree_last_sector);
  m_tree->Branch("first_region", &m_tree_first_region);
  m_tree->Branch("last_region", &m_tree_last_region);
  m_tree->Branch("phi_slope", &m_tree_phi_slope);
  m_tree->Branch("phi_intercept", &m_tree_phi_intercept);
  m_tree->Branch("tbin_slope", &m_tree_tbin_slope);
  m_tree->Branch("tbin_intercept", &m_tree_tbin_intercept);
  m_tree->Branch("chi2_phi", &m_tree_chi2_phi);
  m_tree->Branch("chi2_tbin", &m_tree_chi2_tbin);
  m_tree->Branch("ndof_phi", &m_tree_ndof_phi);
  m_tree->Branch("ndof_tbin", &m_tree_ndof_tbin);

  m_tree->Branch("source_full_track_id", &m_tree_source_full_track_id);
  m_tree->Branch("source_inmodule_track_id", &m_tree_source_inmodule_track_id);
  m_tree->Branch("source_region", &m_tree_source_region);
  m_tree->Branch("source_sector", &m_tree_source_sector);
  m_tree->Branch("source_side", &m_tree_source_side);

  m_tree->Branch("hit_full_track_id", &m_tree_hit_full_track_id);
  m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey", &m_tree_hit_hitkey);

  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackConnector::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  if (!m_padMap)
  {
    std::cerr << Name() << "::InitRun - missing TPC_PADMAP node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_padMap->isValid())
  {
    std::cerr << Name() << "::InitRun - TPC_PADMAP node is invalid" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //if (Verbosity() > 0)
 // {
    m_padMap->print_sector_edges();
 // }

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

int FullTrackConnector::getNodes(PHCompositeNode* topNode)
{
  m_inModuleTrackContainer = findNode::getClass<InModuleTrackContainer>(topNode, m_inputNodeName);
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

  m_fullTrackContainer = findNode::getClass<FullTrackContainer>(topNode, m_outputNodeName);
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
  m_tree_phi_slope.clear();
  m_tree_phi_intercept.clear();
  m_tree_tbin_slope.clear();
  m_tree_tbin_intercept.clear();
  m_tree_chi2_phi.clear();
  m_tree_chi2_tbin.clear();
  m_tree_ndof_phi.clear();
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

bool FullTrackConnector::make_piece(unsigned int source_index, Piece& p) const
{
  const InModuleTrack* trk = m_inModuleTrackContainer->get_track(source_index);
  if (!trk || !trk->isValid()) return false;
  if (trk->get_last_layer() < trk->get_first_layer()) return false;

  p.source_index = source_index;
  p.source_track_id = trk->get_track_id();
  p.event = trk->get_event();
  p.region = trk->get_region();
  p.sector = trk->get_sector();
  p.side = trk->get_side();
  p.first_layer = trk->get_first_layer();
  p.last_layer = trk->get_last_layer();
  p.nblobs = trk->get_nblobs();
  p.nrawhits = trk->get_nrawhits();

  const double l0 = static_cast<double>(p.first_layer);
  const double l1 = static_cast<double>(p.last_layer);
  const double r0 = m_padMap->get_local_radius(p.region, p.first_layer);
  const double r1 = m_padMap->get_local_radius(p.region, p.last_layer);
  if (r0 <= 0.0 || r1 <= 0.0) return false;

  const double pad0 = trk->get_pad_slope() * l0 + trk->get_pad_intercept();
  const double pad1 = trk->get_pad_slope() * l1 + trk->get_pad_intercept();

  // FullTrackConnector works with global phi for the full-track fit:
  //   x = radius [cm], y = global phi [rad], z = timebin.
  // The phi endpoints are unwrapped only internally so the line fit is continuous
  // across the +/- pi boundary.
  double phi0 = wrap_to_pi(m_padMap->get_global_phi(static_cast<unsigned int>(p.side),
                                                    p.region, pad0));

  double phi1 = wrap_to_pi(m_padMap->get_global_phi(static_cast<unsigned int>(p.side),
                                                    p.region, pad1));

  // unwrap endpoint so the piece line does not jump at +/- pi
  phi1 = unwrap_phi_to_reference(phi1, phi0);

  const double tbin0 = trk->get_tbin_slope() * l0 + trk->get_tbin_intercept();
  const double tbin1 = trk->get_tbin_slope() * l1 + trk->get_tbin_intercept();

  if (std::fabs(r1 - r0) > 1.0e-12)
  {
    p.phi_slope = (phi1 - phi0) / (r1 - r0);
    p.phi_intercept = phi0 - p.phi_slope * r0;

    p.tbin_slope = (tbin1 - tbin0) / (r1 - r0);
    p.tbin_intercept = tbin0 - p.tbin_slope * r0;
  }
  else
  {
    p.phi_slope = 0.0;
    p.phi_intercept = phi0;

    p.tbin_slope = 0.0;
    p.tbin_intercept = tbin0;
  }

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

bool FullTrackConnector::refit_candidate(const std::vector<Piece>& pieces,
                                         const std::vector<unsigned int>& piece_indices,
                                         Candidate& c) const
{
  if (piece_indices.empty()) return false;

  std::vector<double> x;
  std::vector<double> phi;
  std::vector<double> tbin;
  std::vector<double> w;
  x.reserve(piece_indices.size() * 2);
  phi.reserve(piece_indices.size() * 2);
  tbin.reserve(piece_indices.size() * 2);
  w.reserve(piece_indices.size() * 2);

  c = Candidate();
  c.piece_indices = piece_indices;



  for (unsigned int ii = 0; ii < piece_indices.size(); ++ii)
  {
    const Piece& p = pieces[piece_indices[ii]];
    const double wt = (p.nrawhits > 0) ? static_cast<double>(p.nrawhits) : 1.0;

    const double r0 = m_padMap->get_local_radius(p.region, p.first_layer);
    const double r1 = m_padMap->get_local_radius(p.region, p.last_layer);
    if (r0 <= 0.0 || r1 <= 0.0) return false;

    double phi0 = p.phi_slope * r0 + p.phi_intercept;
    double phi1 = p.phi_slope * r1 + p.phi_intercept;

    if (!phi.empty())
    {
      phi0 = unwrap_phi_to_reference(phi0, phi.back());
    }
    phi1 = unwrap_phi_to_reference(phi1, phi0);

    x.push_back(r0);
    phi.push_back(phi0);
    tbin.push_back(p.tbin_slope * r0 + p.tbin_intercept);
    w.push_back(wt);

    if (p.last_layer != p.first_layer)
    {
      x.push_back(r1);
      phi.push_back(phi1);
      tbin.push_back(p.tbin_slope * r1 + p.tbin_intercept);
      w.push_back(wt);
    }


    if (ii == 0)
    {
      c.event = p.event;
      c.side = p.side;
      c.first_layer = p.first_layer;
      c.last_layer = p.last_layer;
      c.first_sector = p.sector;
      c.last_sector = p.sector;
      c.first_region = p.region;
      c.last_region = p.region;
    }
    else
    {
      if (p.first_layer < c.first_layer)
      {
        c.first_layer = p.first_layer;
        c.first_sector = p.sector;
        c.first_region = p.region;
      }
      if (p.last_layer > c.last_layer)
      {
        c.last_layer = p.last_layer;
        c.last_sector = p.sector;
        c.last_region = p.region;
      }
    }
    c.nblobs += p.nblobs;
    c.nrawhits += p.nrawhits;
    for (unsigned int ih = 0; ih < p.hitsetkeys.size(); ++ih)
    {
      c.hitsetkeys.push_back(p.hitsetkeys[ih]);
      c.hitkeys.push_back(p.hitkeys[ih]);
    }
  }

  c.nsegments = static_cast<unsigned int>(piece_indices.size());

  if (!weighted_line_fit(x, phi, w, c.phi_slope, c.phi_intercept,
                         c.chi2_phi, c.ndof_phi)) return false;
  if (!weighted_line_fit(x, tbin, w, c.tbin_slope, c.tbin_intercept,
                         c.chi2_tbin, c.ndof_tbin)) return false;

  return true;
}

bool FullTrackConnector::candidates_can_connect(const Candidate& a,
                                                const Piece& b,
                                                double& score,
                                                double& b_phi_intercept_shifted) const
{
  score = std::numeric_limits<double>::max();
  b_phi_intercept_shifted = b.phi_intercept;

  if (a.side != b.side) return false;
  if (a.last_layer >= b.first_layer) return false;

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap) return false;

  const double ra = m_padMap->get_local_radius(a.last_region, a.last_layer);
  const double rb = m_padMap->get_local_radius(b.region, b.first_layer);
  if (ra <= 0.0 || rb <= 0.0) return false;

  const double rmatch = 0.5 * (ra + rb);

  const double phi_a = a.phi_slope * rmatch + a.phi_intercept;
  const double phi_b_raw = b.phi_slope * rmatch + b.phi_intercept;
  const double phi_b = unwrap_phi_to_reference(phi_b_raw, phi_a);
  b_phi_intercept_shifted = b.phi_intercept + (phi_b - phi_b_raw);

  const double tbin_a = a.tbin_slope * rmatch + a.tbin_intercept;
  const double tbin_b = b.tbin_slope * rmatch + b.tbin_intercept;

  const double dphi = std::fabs(phi_a - phi_b);
  const double dtbin = std::fabs(tbin_a - tbin_b);
  const double dmphi = std::fabs(a.phi_slope - b.phi_slope);
  const double dmtbin = std::fabs(a.tbin_slope - b.tbin_slope);

  if (dphi > m_connect_dphi) return false;
  if (dtbin > m_connect_dtbin) return false;
  if (dmphi > m_connect_dphi_slope) return false;
  if (dmtbin > m_connect_dtbin_slope) return false;

  score = (dphi / m_connect_dphi) * (dphi / m_connect_dphi)
        + (dtbin / m_connect_dtbin) * (dtbin / m_connect_dtbin)
        + (dmphi / m_connect_dphi_slope) * (dmphi / m_connect_dphi_slope)
        + (dmtbin / m_connect_dtbin_slope) * (dmtbin / m_connect_dtbin_slope)
        + 0.05 * static_cast<double>(gap);

  return true;
}

void FullTrackConnector::connect_side_pieces(const std::vector<Piece>& pieces,
                                             int side,
                                             std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side == side) order.push_back(i);
  }
  if (order.empty()) return;

  std::sort(order.begin(), order.end(), PieceStartSort(&pieces));
  std::vector<int> used(pieces.size(), 0);

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int iseed = order[io];
    if (used[iseed]) continue;

    std::vector<unsigned int> current_indices;
    current_indices.push_back(iseed);
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current)) continue;

    bool merged_any = true;
    while (merged_any)
    {
      merged_any = false;
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int jo = 0; jo < order.size(); ++jo)
      {
        const unsigned int j = order[jo];
        if (used[j]) continue;

        double score = 0.0;
        double shifted_intercept = 0.0;
        if (!candidates_can_connect(current, pieces[j], score, shifted_intercept)) continue;

        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
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
          merged_any = true;
        }
      }
    }

    output.push_back(current);
  }
}

int FullTrackConnector::process_event(PHCompositeNode*)
{
  reset_tree_vars();
  if (m_fullTrackContainer) m_fullTrackContainer->Reset();

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

  std::vector<Candidate> full_tracks;
  connect_side_pieces(pieces, 0, full_tracks);
  connect_side_pieces(pieces, 1, full_tracks);

  for (unsigned int it = 0; it < full_tracks.size(); ++it)
  {
    const Candidate& c = full_tracks[it];
    const unsigned int full_id = m_fullTrackContainer ? m_fullTrackContainer->size() : it;

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
    out->set_phi_slope(c.phi_slope);
    out->set_phi_intercept(c.phi_intercept);
    out->set_tbin_slope(c.tbin_slope);
    out->set_tbin_intercept(c.tbin_intercept);
    out->set_chi2_phi(c.chi2_phi);
    out->set_chi2_tbin(c.chi2_tbin);
    out->set_ndof_phi(c.ndof_phi);
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
    m_tree_phi_slope.push_back(c.phi_slope);
    m_tree_phi_intercept.push_back(c.phi_intercept);
    m_tree_tbin_slope.push_back(c.tbin_slope);
    m_tree_tbin_intercept.push_back(c.tbin_intercept);
    m_tree_chi2_phi.push_back(c.chi2_phi);
    m_tree_chi2_tbin.push_back(c.chi2_tbin);
    m_tree_ndof_phi.push_back(c.ndof_phi);
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
      m_tree_hit_hitsetkey.push_back(static_cast<unsigned long long>(c.hitsetkeys[ih]));
      m_tree_hit_hitkey.push_back(static_cast<unsigned long long>(c.hitkeys[ih]));
    }

    if (m_fullTrackContainer) m_fullTrackContainer->add_track(out);
    else delete out;
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " input pieces=" << pieces.size()
              << " full_tracks=" << full_tracks.size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
