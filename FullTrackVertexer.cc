#include "FullTrackVertexer.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

namespace
{
  bool finite_number(const double x)
  {
    return (x == x && std::fabs(x) < 1.0e30);
  }

  double wrap_to_pi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double unwrap_phi_to_reference(double phi, const double ref)
  {
    while (phi - ref > M_PI) phi -= 2.0 * M_PI;
    while (phi - ref < -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double sqr(const double x) { return x * x; }
}

FullTrackVertexer::TrackInfo::TrackInfo() :
  index(0),
  track_id(0),
  side(0),
  first_layer(0),
  last_layer(0),
  nrawhits(0),
  phi_slope(0.0),
  phi_intercept(0.0),
  tbin_slope(0.0),
  tbin_intercept(0.0),
  chi2_phi_ndf(0.0),
  chi2_tbin_ndf(0.0)
{
}

FullTrackVertexer::PairVertex::PairVertex() :
  i(0),
  j(0),
  side(0),
  r_phi(0.0),
  phi(0.0),
  x(0.0),
  y(0.0),
  r_tbin(0.0),
  tbin(0.0),
  residual_rphi(0.0),
  residual_tbin(0.0),
  weight(1.0)
{
}

FullTrackVertexer::FullTrackVertexer(const std::string& name,
                                     const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_inputNodeName("FULLTRACKS")
  , m_outputFile(nullptr)
  , m_tree(nullptr)
  , m_fullTrackContainer(nullptr)
  , m_event(0)
  , m_min_layers(16)
  , m_min_rawhits(25)
  , m_max_chi2_phi_ndf(100.0)
  , m_max_chi2_tbin_ndf(100.0)
  , m_min_slope_diff_phi(1.0e-6)
  , m_min_slope_diff_tbin(1.0e-6)
  , m_min_vertex_r(-30.0)
  , m_max_vertex_r(30.0)
  , m_max_pair_dtbin(150.0)
  , m_max_residual_rphi(2.0)
  , m_max_residual_tbin(80.0)
  , m_max_track_dca_rphi(1.5)
  , m_max_track_dca_tbin(50.0)
  , m_vertex_timebin_gap(20.0)
  , m_min_pairs_per_vertex(5)
  , m_tree_event(0)
  , m_tree_ntracks_input(0)
  , m_tree_ntracks_used(0)
  , m_tree_npairs(0)
  , m_tree_npairs_used(0)
  , m_tree_vertex_side(0)
  , m_tree_vertex_x(0.0)
  , m_tree_vertex_y(0.0)
  , m_tree_vertex_r(0.0)
  , m_tree_vertex_phi(0.0)
  , m_tree_vertex_tbin(0.0)
  , m_tree_vertex_quality(0.0)
{
}

FullTrackVertexer::~FullTrackVertexer()
{
  if (m_outputFile)
  {
    delete m_outputFile;
    m_outputFile = nullptr;
  }
}

int FullTrackVertexer::Init(PHCompositeNode*)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("FullTrackVertex", "TPC full-track pair vertex QA");
  m_tree->Branch("event", &m_tree_event, "event/I");
  m_tree->Branch("ntracks_input", &m_tree_ntracks_input, "ntracks_input/i");
  m_tree->Branch("ntracks_used", &m_tree_ntracks_used, "ntracks_used/i");
  m_tree->Branch("npairs", &m_tree_npairs, "npairs/i");
  m_tree->Branch("npairs_used", &m_tree_npairs_used, "npairs_used/i");
  m_tree->Branch("vertex_side", &m_tree_vertex_side, "vertex_side/I");
  m_tree->Branch("vertex_x", &m_tree_vertex_x, "vertex_x/D");
  m_tree->Branch("vertex_y", &m_tree_vertex_y, "vertex_y/D");
  m_tree->Branch("vertex_r", &m_tree_vertex_r, "vertex_r/D");
  m_tree->Branch("vertex_phi", &m_tree_vertex_phi, "vertex_phi/D");
  m_tree->Branch("vertex_tbin", &m_tree_vertex_tbin, "vertex_tbin/D");
  m_tree->Branch("vertex_quality", &m_tree_vertex_quality, "vertex_quality/D");

  // All vertices found in this event after clustering pair candidates in tbin.
  m_tree->Branch("all_vertex_side", &m_tree_all_vertex_side);
  m_tree->Branch("all_vertex_x", &m_tree_all_vertex_x);
  m_tree->Branch("all_vertex_y", &m_tree_all_vertex_y);
  m_tree->Branch("all_vertex_r", &m_tree_all_vertex_r);
  m_tree->Branch("all_vertex_phi", &m_tree_all_vertex_phi);
  m_tree->Branch("all_vertex_tbin", &m_tree_all_vertex_tbin);
  m_tree->Branch("all_vertex_quality", &m_tree_all_vertex_quality);
  m_tree->Branch("all_vertex_npairs", &m_tree_all_vertex_npairs);
  m_tree->Branch("all_vertex_ntracks_assigned", &m_tree_all_vertex_ntracks_assigned);

  m_tree->Branch("pair_side", &m_tree_pair_side);
  m_tree->Branch("pair_i", &m_tree_pair_i);
  m_tree->Branch("pair_j", &m_tree_pair_j);
  m_tree->Branch("pair_r_phi", &m_tree_pair_r_phi);
  m_tree->Branch("pair_phi", &m_tree_pair_phi);
  m_tree->Branch("pair_x", &m_tree_pair_x);
  m_tree->Branch("pair_y", &m_tree_pair_y);
  m_tree->Branch("pair_r_tbin", &m_tree_pair_r_tbin);
  m_tree->Branch("pair_tbin", &m_tree_pair_tbin);
  m_tree->Branch("pair_residual_rphi", &m_tree_pair_residual_rphi);
  m_tree->Branch("pair_residual_tbin", &m_tree_pair_residual_tbin);
  m_tree->Branch("pair_weight", &m_tree_pair_weight);

  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackVertexer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int FullTrackVertexer::End(PHCompositeNode*)
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

int FullTrackVertexer::getNodes(PHCompositeNode* topNode)
{
  m_fullTrackContainer = findNode::getClass<FullTrackContainer>(topNode, m_inputNodeName);
  if (!m_fullTrackContainer)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void FullTrackVertexer::reset_tree_vars()
{
  m_tree_event = m_event;
  m_tree_ntracks_input = 0;
  m_tree_ntracks_used = 0;
  m_tree_npairs = 0;
  m_tree_npairs_used = 0;
  m_tree_vertex_side = 0;
  m_tree_vertex_x = 0.0;
  m_tree_vertex_y = 0.0;
  m_tree_vertex_r = 0.0;
  m_tree_vertex_phi = 0.0;
  m_tree_vertex_tbin = 0.0;
  m_tree_vertex_quality = 0.0;

  m_tree_all_vertex_side.clear();
  m_tree_all_vertex_x.clear();
  m_tree_all_vertex_y.clear();
  m_tree_all_vertex_r.clear();
  m_tree_all_vertex_phi.clear();
  m_tree_all_vertex_tbin.clear();
  m_tree_all_vertex_quality.clear();
  m_tree_all_vertex_npairs.clear();
  m_tree_all_vertex_ntracks_assigned.clear();

  m_tree_pair_side.clear();
  m_tree_pair_i.clear();
  m_tree_pair_j.clear();
  m_tree_pair_r_phi.clear();
  m_tree_pair_phi.clear();
  m_tree_pair_x.clear();
  m_tree_pair_y.clear();
  m_tree_pair_r_tbin.clear();
  m_tree_pair_tbin.clear();
  m_tree_pair_residual_rphi.clear();
  m_tree_pair_residual_tbin.clear();
  m_tree_pair_weight.clear();
}

bool FullTrackVertexer::make_track_info(unsigned int index, TrackInfo& t) const
{
  const FullTrack* trk = m_fullTrackContainer->get_track(index);
  if (!trk || !trk->isValid()) return false;

  if (trk->get_last_layer() < trk->get_first_layer()) return false;
  const unsigned int nlayers = trk->get_last_layer() - trk->get_first_layer() + 1;
  if (nlayers < m_min_layers) return false;
  if (trk->get_nrawhits() < m_min_rawhits) return false;

  const int ndof_phi = trk->get_ndof_phi();
  const int ndof_tbin = trk->get_ndof_tbin();
  const double chi2_phi_ndf = (ndof_phi > 0) ? trk->get_chi2_phi() / static_cast<double>(ndof_phi) : 0.0;
  const double chi2_tbin_ndf = (ndof_tbin > 0) ? trk->get_chi2_tbin() / static_cast<double>(ndof_tbin) : 0.0;

  if (chi2_phi_ndf > m_max_chi2_phi_ndf) return false;
  if (chi2_tbin_ndf > m_max_chi2_tbin_ndf) return false;

  t.index = index;
  t.track_id = trk->get_track_id();
  t.side = trk->get_side();
  t.first_layer = trk->get_first_layer();
  t.last_layer = trk->get_last_layer();
  t.nrawhits = trk->get_nrawhits();
  t.phi_slope = trk->get_phi_slope();
  t.phi_intercept = trk->get_phi_intercept();
  t.tbin_slope = trk->get_tbin_slope();
  t.tbin_intercept = trk->get_tbin_intercept();
  t.chi2_phi_ndf = chi2_phi_ndf;
  t.chi2_tbin_ndf = chi2_tbin_ndf;

  if (!finite_number(t.phi_slope) || !finite_number(t.phi_intercept)) return false;
  if (!finite_number(t.tbin_slope) || !finite_number(t.tbin_intercept)) return false;

  return true;
}

bool FullTrackVertexer::make_pair_vertex(const TrackInfo& a,
                                         const TrackInfo& b,
                                         PairVertex& pv) const
{
  if (a.index == b.index) return false;

  // Vertex finding is independent for the two TPC sides.
  // Do not form South-North mixed pairs.
  if (a.side != b.side) return false;

  // If the tracks are almost parallel in either projection, the pair vertex is unstable.
  const double dphi_slope = a.phi_slope - b.phi_slope;
  const double dtbin_slope = a.tbin_slope - b.tbin_slope;
  if (std::fabs(dphi_slope) < m_min_slope_diff_phi) return false;
  if (std::fabs(dtbin_slope) < m_min_slope_diff_tbin) return false;

  // The full-track fit is phi(r) and tbin(r).  For each pair, find the
  // intersection separately in the r-phi and r-tbin projections.
  double b_phi_intercept = b.phi_intercept;
  b_phi_intercept = unwrap_phi_to_reference(b_phi_intercept, a.phi_intercept);

  const double r_phi = (b_phi_intercept - a.phi_intercept) / dphi_slope;
  const double phi_a = a.phi_slope * r_phi + a.phi_intercept;
  const double phi_b = b.phi_slope * r_phi + b_phi_intercept;
  const double phi = 0.5 * (phi_a + unwrap_phi_to_reference(phi_b, phi_a));

  const double r_tbin = (b.tbin_intercept - a.tbin_intercept) / dtbin_slope;
  const double tbin_a = a.tbin_slope * r_tbin + a.tbin_intercept;
  const double tbin_b = b.tbin_slope * r_tbin + b.tbin_intercept;
  const double tbin = 0.5 * (tbin_a + tbin_b);

  if (!finite_number(r_phi) || !finite_number(phi)) return false;
  if (!finite_number(r_tbin) || !finite_number(tbin)) return false;

  if (r_phi < m_min_vertex_r || r_phi > m_max_vertex_r) return false;
  if (r_tbin < m_min_vertex_r || r_tbin > m_max_vertex_r) return false;
  if (std::fabs(tbin_a - tbin_b) > m_max_pair_dtbin) return false;

  pv.i = a.index;
  pv.j = b.index;
  pv.side = a.side;
  pv.r_phi = r_phi;
  pv.phi = wrap_to_pi(phi);
  pv.x = r_phi * std::cos(pv.phi);
  pv.y = r_phi * std::sin(pv.phi);
  pv.r_tbin = r_tbin;
  pv.tbin = tbin;

  // Symmetric residuals around the 2-projection pair point.  These are useful
  // for outlier removal and QA.
  const double r_common = 0.5 * (r_phi + r_tbin);
  const double yphi_a = a.phi_slope * r_common + a.phi_intercept;
  const double yphi_b = b.phi_slope * r_common + b_phi_intercept;
  const double dtbin_a2 = a.tbin_slope * r_common + a.tbin_intercept;
  const double dtbin_b2 = b.tbin_slope * r_common + b.tbin_intercept;

  pv.residual_rphi = std::fabs(r_common * (yphi_a - unwrap_phi_to_reference(yphi_b, yphi_a)));
  pv.residual_tbin = std::fabs(dtbin_a2 - dtbin_b2);

  if (pv.residual_rphi > m_max_residual_rphi) return false;
  if (pv.residual_tbin > m_max_residual_tbin) return false;

  const double w_hits = std::sqrt(static_cast<double>(a.nrawhits) * static_cast<double>(b.nrawhits));
  const double w_chi2 = 1.0 / (1.0 + a.chi2_phi_ndf + b.chi2_phi_ndf + a.chi2_tbin_ndf + b.chi2_tbin_ndf);
  pv.weight = (w_hits > 0.0) ? w_hits * w_chi2 : w_chi2;
  if (!finite_number(pv.weight) || pv.weight <= 0.0) pv.weight = 1.0;

  return true;
}

double FullTrackVertexer::median(std::vector<double> v) const
{
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const unsigned int n = static_cast<unsigned int>(v.size());
  if (n % 2) return v[n / 2];
  return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double FullTrackVertexer::weighted_rms(const std::vector<double>& v,
                                       const double center) const
{
  if (v.empty()) return 0.0;
  double sum = 0.0;
  for (unsigned int i = 0; i < v.size(); ++i) sum += sqr(v[i] - center);
  return std::sqrt(sum / static_cast<double>(v.size()));
}

bool FullTrackVertexer::calculate_event_vertex(const std::vector<PairVertex>& pairs,
                                               double& vx,
                                               double& vy,
                                               double& vr,
                                               double& vphi,
                                               double& vtbin,
                                               double& vquality,
                                               unsigned int& nused_pairs) const
{
  vx = 0.0;
  vy = 0.0;
  vr = 0.0;
  vphi = 0.0;
  vtbin = 0.0;
  vquality = 0.0;
  nused_pairs = 0;

  if (pairs.empty()) return false;

  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> rs;
  std::vector<double> tbins;
  xs.reserve(pairs.size());
  ys.reserve(pairs.size());
  rs.reserve(pairs.size());
  tbins.reserve(pairs.size());

  for (unsigned int i = 0; i < pairs.size(); ++i)
  {
    xs.push_back(pairs[i].x);
    ys.push_back(pairs[i].y);
    rs.push_back(pairs[i].r_phi);
    tbins.push_back(pairs[i].tbin);
  }

  // Robust first pass: median pair vertex.
  const double x0 = median(xs);
  const double y0 = median(ys);
  const double tbin0 = median(tbins);
  const double sx = weighted_rms(xs, x0);
  const double sy = weighted_rms(ys, y0);
  const double st = weighted_rms(tbins, tbin0);

  double sw = 0.0;
  double swx = 0.0;
  double swy = 0.0;
  double swt = 0.0;

  for (unsigned int i = 0; i < pairs.size(); ++i)
  {
    const double dx = pairs[i].x - x0;
    const double dy = pairs[i].y - y0;
    const double dt = pairs[i].tbin - tbin0;

    if (sx > 0.0 && std::fabs(dx) > 3.0 * sx) continue;
    if (sy > 0.0 && std::fabs(dy) > 3.0 * sy) continue;
    if (st > 0.0 && std::fabs(dt) > 3.0 * st) continue;

    const double w = pairs[i].weight;
    sw += w;
    swx += w * pairs[i].x;
    swy += w * pairs[i].y;
    swt += w * pairs[i].tbin;
    ++nused_pairs;
  }

  if (sw <= 0.0 || nused_pairs == 0) return false;

  vx = swx / sw;
  vy = swy / sw;
  vtbin = swt / sw;
  vr = std::sqrt(vx * vx + vy * vy);
  vphi = wrap_to_pi(std::atan2(vy, vx));

  double chi2 = 0.0;
  for (unsigned int i = 0; i < pairs.size(); ++i)
  {
    const double dx = pairs[i].x - vx;
    const double dy = pairs[i].y - vy;
    const double dt = pairs[i].tbin - vtbin;
    chi2 += dx * dx + dy * dy + 0.01 * dt * dt;
  }
  vquality = std::sqrt(chi2 / static_cast<double>(pairs.size()));

  return true;
}

int FullTrackVertexer::process_event(PHCompositeNode*)
{
  reset_tree_vars();

  if (!m_fullTrackContainer)
  {
    if (m_tree) m_tree->Fill();
    ++m_event;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const unsigned int ntracks = m_fullTrackContainer->size();
  m_tree_ntracks_input = ntracks;

  std::vector<TrackInfo> tracks;
  tracks.reserve(ntracks);
  for (unsigned int i = 0; i < ntracks; ++i)
  {
    TrackInfo t;
    if (make_track_info(i, t)) tracks.push_back(t);
  }
  m_tree_ntracks_used = static_cast<unsigned int>(tracks.size());

  std::vector<PairVertex> pairs;
  for (unsigned int i = 0; i < tracks.size(); ++i)
  {
    for (unsigned int j = i + 1; j < tracks.size(); ++j)
    {
      PairVertex pv;
      if (!make_pair_vertex(tracks[i], tracks[j], pv)) continue;

      pairs.push_back(pv);

      m_tree_pair_side.push_back(pv.side);
      m_tree_pair_i.push_back(pv.i);
      m_tree_pair_j.push_back(pv.j);
      m_tree_pair_r_phi.push_back(pv.r_phi);
      m_tree_pair_phi.push_back(pv.phi);
      m_tree_pair_x.push_back(pv.x);
      m_tree_pair_y.push_back(pv.y);
      m_tree_pair_r_tbin.push_back(pv.r_tbin);
      m_tree_pair_tbin.push_back(pv.tbin);
      m_tree_pair_residual_rphi.push_back(pv.residual_rphi);
      m_tree_pair_residual_tbin.push_back(pv.residual_tbin);
      m_tree_pair_weight.push_back(pv.weight);
    }
  }

  m_tree_npairs = static_cast<unsigned int>(pairs.size());

  // Clear vertices first. Tracks will be assigned to the best same-side,
  // same-timebin vertex found from their pair cluster.
  for (unsigned int i = 0; i < ntracks; ++i)
  {
    FullTrack* trk = m_fullTrackContainer->get_track(i);
    if (!trk) continue;
    trk->set_vertex_valid(0);
    trk->set_vertex_x(0.0);
    trk->set_vertex_y(0.0);
    trk->set_vertex_r(0.0);
    trk->set_vertex_phi(0.0);
    trk->set_vertex_tbin(0.0);
    trk->set_vertex_npairs(0);
    trk->set_vertex_quality(0.0);
  }

  std::vector<double> best_track_dtbin(ntracks, std::numeric_limits<double>::max());

  // Do the vertex determination independently for each TPC side.
  // This avoids South/North mixed pairs and allows separate timebin clusters
  // on each side.
  std::sort(pairs.begin(), pairs.end(), [](const PairVertex& a, const PairVertex& b) {
    if (a.side != b.side) return a.side < b.side;
    return a.tbin < b.tbin;
  });

  bool any_vertex_ok = false;
  unsigned int best_vertex_index = 0;
  unsigned int best_vertex_npairs = 0;

  unsigned int ip = 0;
  while (ip < pairs.size())
  {
    const int side = pairs[ip].side;

    std::vector<PairVertex> side_pairs;
    while (ip < pairs.size() && pairs[ip].side == side)
    {
      side_pairs.push_back(pairs[ip]);
      ++ip;
    }

    std::vector<std::vector<PairVertex> > pair_clusters;
    for (unsigned int is = 0; is < side_pairs.size(); ++is)
    {
      if (pair_clusters.empty())
      {
        pair_clusters.push_back(std::vector<PairVertex>());
        pair_clusters.back().push_back(side_pairs[is]);
        continue;
      }

      // Running mean is more stable than comparing only to the previous pair.
      double mean_tbin = 0.0;
      for (unsigned int j = 0; j < pair_clusters.back().size(); ++j) mean_tbin += pair_clusters.back()[j].tbin;
      mean_tbin /= static_cast<double>(pair_clusters.back().size());

      if (std::fabs(side_pairs[is].tbin - mean_tbin) > m_vertex_timebin_gap)
      {
        pair_clusters.push_back(std::vector<PairVertex>());
      }
      pair_clusters.back().push_back(side_pairs[is]);
    }

    for (unsigned int ic = 0; ic < pair_clusters.size(); ++ic)
    {
      if (pair_clusters[ic].size() < m_min_pairs_per_vertex) continue;

      double vx = 0.0;
      double vy = 0.0;
      double vr = 0.0;
      double vphi = 0.0;
      double vtbin = 0.0;
      double vquality = 0.0;
      unsigned int nused_pairs = 0;

      const bool ok = calculate_event_vertex(pair_clusters[ic], vx, vy, vr, vphi, vtbin, vquality, nused_pairs);
      if (!ok || nused_pairs < m_min_pairs_per_vertex) continue;

      any_vertex_ok = true;

      m_tree_all_vertex_side.push_back(side);
      m_tree_all_vertex_x.push_back(vx);
      m_tree_all_vertex_y.push_back(vy);
      m_tree_all_vertex_r.push_back(vr);
      m_tree_all_vertex_phi.push_back(vphi);
      m_tree_all_vertex_tbin.push_back(vtbin);
      m_tree_all_vertex_quality.push_back(vquality);
      m_tree_all_vertex_npairs.push_back(nused_pairs);
      m_tree_all_vertex_ntracks_assigned.push_back(0);
      const unsigned int this_vertex_index = static_cast<unsigned int>(m_tree_all_vertex_npairs.size() - 1);

      if (nused_pairs > best_vertex_npairs)
      {
        best_vertex_npairs = nused_pairs;
        best_vertex_index = static_cast<unsigned int>(m_tree_all_vertex_npairs.size() - 1);
      }

      // Assign this same-side vertex only to tracks that participated in this
      // same-side tbin cluster. If a track appears in more than one cluster,
      // keep the closest one in tbin.
      std::set<unsigned int> cluster_track_indices;
      for (unsigned int jp = 0; jp < pair_clusters[ic].size(); ++jp)
      {
        cluster_track_indices.insert(pair_clusters[ic][jp].i);
        cluster_track_indices.insert(pair_clusters[ic][jp].j);
      }

      for (std::set<unsigned int>::const_iterator it = cluster_track_indices.begin();
           it != cluster_track_indices.end(); ++it)
      {
        const unsigned int idx = *it;
        if (idx >= ntracks) continue;

        FullTrack* trk = m_fullTrackContainer->get_track(idx);
        if (!trk) continue;
        if (trk->get_side() != side) continue;

        const double track_phi_at_vtx = trk->get_phi_slope() * vr + trk->get_phi_intercept();
        const double dphi = wrap_to_pi(track_phi_at_vtx - vphi);
        const double dca_rphi = std::fabs(vr * dphi);

        const double track_tbin_at_vtx = trk->get_tbin_slope() * vr + trk->get_tbin_intercept();
        const double dt = std::fabs(track_tbin_at_vtx - vtbin);

        // Tight DCA-like track-to-vertex requirement. This prevents a broad
        // time cluster from assigning unrelated tracks to the vertex.
        if (dca_rphi > m_max_track_dca_rphi) continue;
        if (dt > m_max_track_dca_tbin) continue;

        if (dt > best_track_dtbin[idx]) continue;

        if (!trk->get_vertex_valid() && this_vertex_index < m_tree_all_vertex_ntracks_assigned.size())
        {
          ++m_tree_all_vertex_ntracks_assigned[this_vertex_index];
        }

        best_track_dtbin[idx] = dt;
        trk->set_vertex_valid(1);
        trk->set_vertex_x(vx);
        trk->set_vertex_y(vy);
        trk->set_vertex_r(vr);
        trk->set_vertex_phi(vphi);
        trk->set_vertex_tbin(vtbin);
        trk->set_vertex_npairs(nused_pairs);
        trk->set_vertex_quality(vquality);
      }
    }
  }

  if (any_vertex_ok && best_vertex_index < m_tree_all_vertex_npairs.size())
  {
    // Backward-compatible scalar branches: largest/most-populated vertex among
    // all side/timebin vertices. Vector branches keep all side-separated vertices.
    m_tree_npairs_used = m_tree_all_vertex_npairs[best_vertex_index];
    m_tree_vertex_side = m_tree_all_vertex_side[best_vertex_index];
    m_tree_vertex_x = m_tree_all_vertex_x[best_vertex_index];
    m_tree_vertex_y = m_tree_all_vertex_y[best_vertex_index];
    m_tree_vertex_r = m_tree_all_vertex_r[best_vertex_index];
    m_tree_vertex_phi = m_tree_all_vertex_phi[best_vertex_index];
    m_tree_vertex_tbin = m_tree_all_vertex_tbin[best_vertex_index];
    m_tree_vertex_quality = m_tree_all_vertex_quality[best_vertex_index];
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    unsigned int nside0 = 0;
    unsigned int nside1 = 0;
    unsigned int nside_other = 0;
    for (unsigned int iv = 0; iv < m_tree_all_vertex_side.size(); ++iv)
    {
      if (m_tree_all_vertex_side[iv] == 0) ++nside0;
      else if (m_tree_all_vertex_side[iv] == 1) ++nside1;
      else ++nside_other;
    }

    std::cout << Name() << "::process_event - event " << m_event
              << " full_tracks=" << ntracks
              << " used_tracks=" << tracks.size()
              << " pairs=" << pairs.size()
              << " nvertices=" << m_tree_all_vertex_npairs.size()
              << " nvertices_side0=" << nside0
              << " nvertices_side1=" << nside1
              << " nvertices_other=" << nside_other
              << " best_side=" << m_tree_vertex_side
              << " used_pairs_best=" << m_tree_npairs_used
              << " vertex_ok=" << any_vertex_ok
              << " x_best=" << m_tree_vertex_x
              << " y_best=" << m_tree_vertex_y
              << " tbin_best=" << m_tree_vertex_tbin
              << std::endl;

    // Print every side/timebin vertex, not only the backward-compatible best one.
    for (unsigned int iv = 0; iv < m_tree_all_vertex_npairs.size(); ++iv)
    {
      std::cout << "  vertex " << iv
                << " side=" << m_tree_all_vertex_side[iv]
                << " npairs=" << m_tree_all_vertex_npairs[iv]
                << " ntracks_assigned=" << ((iv < m_tree_all_vertex_ntracks_assigned.size()) ? m_tree_all_vertex_ntracks_assigned[iv] : 0)
                << " x=" << m_tree_all_vertex_x[iv]
                << " y=" << m_tree_all_vertex_y[iv]
                << " r=" << m_tree_all_vertex_r[iv]
                << " phi=" << m_tree_all_vertex_phi[iv]
                << " tbin=" << m_tree_all_vertex_tbin[iv]
                << " quality=" << m_tree_all_vertex_quality[iv]
                << std::endl;
    }
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
