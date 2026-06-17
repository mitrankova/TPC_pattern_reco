#include "FullTrackVertexer.h"

#include "Fitter.h"
#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "FullTrackVertexv1.h"
#include "FullTrackVertexContainer.h"
#include "FullTrackVertexContainerv1.h"
#include "IdealPadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi <= -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  double unwrap_phi_near(double phi, const double reference)
  {
    while (phi - reference > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi - reference < -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  bool is_good_number(const double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }

  double sagitta_model_derivative(const double xrot,
                                  const double x0,
                                  const double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }

  double sagitta_phi_at_radius(const double radius,
                               const Fitter::SagittaFit& fit)
  {
    const double c = std::cos(fit.theta);
    const double s = std::sin(fit.theta);
    double yy = std::tan(fit.theta) * radius;

    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      const double xrot = c * radius + s * yy;
      const double yrot = -s * radius + c * yy;
      const double f = Fitter::sagittaModel(xrot, fit.S, fit.x0, fit.invR);
      const double g = yrot - f;
      const double df = sagitta_model_derivative(xrot, fit.x0, fit.invR);
      const double dg = c - df * s;
      if (std::fabs(dg) < 1.0e-12) break;
      const double step = g / dg;
      yy -= step;
      if (std::fabs(step) < 1.0e-10) break;
    }

    return fit.b + yy;
  }

  struct RadiusSort
  {
    const std::vector<FullTrackVertexer::HitPoint>* pts;
    explicit RadiusSort(const std::vector<FullTrackVertexer::HitPoint>* p) : pts(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const FullTrackVertexer::HitPoint& pa = (*pts)[a];
      const FullTrackVertexer::HitPoint& pb = (*pts)[b];
      if (pa.radius != pb.radius) return pa.radius < pb.radius;
      return pa.tbin < pb.tbin;
    }
  };

  struct VertexFit
  {
    bool ok {false};
    bool use_sagitta {false};
    double d0 {0.0};
    double timebin0 {0.0};
    double phi_slope {0.0};
    double phi_intercept {0.0};
    double tbin_slope {0.0};
    double tbin_intercept {0.0};
    Fitter::SagittaFit phi_sagitta;
  };

  struct CollisionTrack
  {
    VertexFit fit;
    unsigned int nlayers {0};
  };

  struct CollisionFit
  {
    bool ok {false};
    double radius {0.0};
    double phi {0.0};
    double timebin {0.0};
    double timebin_rms {0.0};
    unsigned int ntracks {0};
  };

  unsigned int count_unique_layers(const std::vector<FullTrackVertexer::HitPoint>& pts)
  {
    std::vector<unsigned int> layers;
    layers.reserve(pts.size());
    for (const auto& p : pts)
    {
      if (p.ok) layers.push_back(p.layer);
    }
    std::sort(layers.begin(), layers.end());
    layers.erase(std::unique(layers.begin(), layers.end()), layers.end());
    return static_cast<unsigned int>(layers.size());
  }

  VertexFit fit_vertex_points(const std::vector<FullTrackVertexer::HitPoint>& pts,
                              const bool use_sagitta,
                              const double weight_power,
                              const double weight_floor_frac)
  {
    VertexFit result;
    if (pts.size() < 2) return result;

    std::vector<unsigned int> order;
    order.reserve(pts.size());
    for (unsigned int i = 0; i < pts.size(); ++i) order.push_back(i);
    std::sort(order.begin(), order.end(), RadiusSort(&pts));

    double max_adc = 0.0;
    for (const auto& p : pts) max_adc = std::max(max_adc, static_cast<double>(p.adc));
    if (max_adc <= 0.0) max_adc = 1.0;

    std::vector<Fitter::FitPoint> radius_phi_points;
    std::vector<Fitter::FitPoint> radius_tbin_points;
    radius_phi_points.reserve(pts.size());
    radius_tbin_points.reserve(pts.size());

    bool first = true;
    double phi_reference = 0.0;
    for (unsigned int io = 0; io < order.size(); ++io)
    {
      const FullTrackVertexer::HitPoint& p = pts[order[io]];
      if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) continue;

      double phi = p.global_phi;
      if (first)
      {
        phi_reference = phi;
        first = false;
      }
      else
      {
        phi = unwrap_phi_near(phi, phi_reference);
        phi_reference = phi;
      }

      const double w = Fitter::adcWeight(static_cast<double>(p.adc), max_adc,
                                         weight_power, weight_floor_frac);
      radius_phi_points.emplace_back(p.radius, phi, w);
      radius_tbin_points.emplace_back(p.radius, static_cast<double>(p.tbin), w);
    }

    if (radius_phi_points.size() < 2 || radius_tbin_points.size() < 2) return result;

    const Fitter::LineFit phi_line = Fitter::fitLine(radius_phi_points);
    const Fitter::LineFit tbin_line = Fitter::fitLine(radius_tbin_points);
    if (!phi_line.ok || !tbin_line.ok) return result;

    result.d0 = phi_line.intercept;
    result.timebin0 = tbin_line.intercept;
    result.phi_slope = phi_line.slope;
    result.phi_intercept = phi_line.intercept;
    result.tbin_slope = tbin_line.slope;
    result.tbin_intercept = tbin_line.intercept;

    if (use_sagitta && radius_phi_points.size() >= 3)
    {
      result.phi_sagitta = Fitter::fitSagitta(radius_phi_points);
      result.use_sagitta = result.phi_sagitta.ok;
      if (result.use_sagitta) result.d0 = sagitta_phi_at_radius(0.0, result.phi_sagitta);
    }

    if (!is_good_number(result.d0) || !is_good_number(result.timebin0)) return result;
    result.d0 = wrap_phi(result.d0);
    result.ok = true;
    return result;
  }

  double phi_at_radius(const VertexFit& fit, const double radius)
  {
    return fit.use_sagitta ? sagitta_phi_at_radius(radius, fit.phi_sagitta)
                           : fit.phi_slope * radius + fit.phi_intercept;
  }

  CollisionFit fit_collision_point(const std::vector<CollisionTrack>& tracks)
  {
    CollisionFit result;
    if (tracks.size() < 2) return result;

    double s_mm = 0.0;
    double s_mt = 0.0;
    double s_tt = 0.0;
    double rhs_r = 0.0;
    double rhs_t = 0.0;

    for (const auto& trk : tracks)
    {
      const double w = static_cast<double>(trk.nlayers);
      const double m = trk.fit.tbin_slope;
      const double b = trk.fit.tbin_intercept;
      s_mm += w * m * m;
      s_mt -= w * m;
      s_tt += w;
      rhs_r -= w * m * b;
      rhs_t += w * b;
    }

    const double det = s_mm * s_tt - s_mt * s_mt;
    if (std::fabs(det) < 1.0e-12) return result;

    result.radius = (rhs_r * s_tt - s_mt * rhs_t) / det;
    result.timebin = (s_mm * rhs_t - s_mt * rhs_r) / det;
    if (!is_good_number(result.radius) || !is_good_number(result.timebin)) return result;

    double sum_w = 0.0;
    double sum_res2 = 0.0;
    double sum_sin = 0.0;
    double sum_cos = 0.0;
    for (const auto& trk : tracks)
    {
      const double w = static_cast<double>(trk.nlayers);
      const double tpred = trk.fit.tbin_slope * result.radius + trk.fit.tbin_intercept;
      const double residual = tpred - result.timebin;
      const double phi = wrap_phi(phi_at_radius(trk.fit, result.radius));
      sum_w += w;
      sum_res2 += w * residual * residual;
      sum_sin += w * std::sin(phi);
      sum_cos += w * std::cos(phi);
    }

    if (sum_w <= 0.0 || (sum_sin == 0.0 && sum_cos == 0.0)) return result;
    result.phi = wrap_phi(std::atan2(sum_sin, sum_cos));
    result.timebin_rms = std::sqrt(sum_res2 / sum_w);
    result.ntracks = static_cast<unsigned int>(tracks.size());
    result.ok = true;
    return result;
  }
}

FullTrackVertexer::FullTrackVertexer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("FULLTRACKS")
  , m_outputNodeName("FULLTRACKVERTICES")
  , m_fullTracks(nullptr)
  , m_vertices(nullptr)
  , m_hits(nullptr)
  , m_idealPadMap(nullptr)
  , m_fitWeightPower(1.0)
  , m_fitWeightFloorFrac(0.05)
  , m_useSagittaPhiFit(true)
  , m_collisionMinTrackLayers(10)
{
}

FullTrackVertexer::~FullTrackVertexer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int FullTrackVertexer::InitRun(PHCompositeNode* topNode)
{
  if (!createNodes(topNode)) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - cannot load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool FullTrackVertexer::getNodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<FullTrackContainer>(topNode, m_inputNodeName.c_str());
  if (!m_fullTracks)
  {
    const char* candidate_names[] = {
      "FULLTRACKS",
      "FullTracks",
      "FullTrackContainer",
      "FULLTRACKCONTAINER",
      "FULLTRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_fullTracks;
         ++i)
    {
      m_fullTracks = findNode::getClass<FullTrackContainer>(topNode, candidate_names[i]);
      if (m_fullTracks) m_inputNodeName = candidate_names[i];
    }
  }

  if (!m_fullTracks)
  {
    std::cerr << Name() << "::getNodes - missing FullTrackContainer node" << std::endl;
    return false;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return false;
  }

  return true;
}

bool FullTrackVertexer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_vertices = findNode::getClass<FullTrackVertexContainer>(topNode, m_outputNodeName.c_str());
  if (!m_vertices)
  {
    m_vertices = new FullTrackVertexContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_vertices, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return true;
}

int FullTrackVertexer::process_event(PHCompositeNode* topNode)
{
  if (!m_fullTracks || !m_hits || !m_vertices)
  {
    if (!getNodes(topNode) || !createNodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!m_idealPadMap || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::process_event - IdealPadMap is not loaded" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_vertices->Reset();

  std::vector<CollisionTrack> collision_tracks;

  const unsigned int ntracks = m_fullTracks ? m_fullTracks->size() : 0;
  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const FullTrack* trk = m_fullTracks->get_track(itrk);
    if (!trk) continue;

    std::vector<HitPoint> pts;
    pts.reserve(trk->size_hit_indices());
    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const FullTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);
      if (p.ok) pts.push_back(p);
    }

    const VertexFit fit = fit_vertex_points(pts, m_useSagittaPhiFit,
                                            m_fitWeightPower, m_fitWeightFloorFrac);
    if (!fit.ok) continue;

    FullTrackVertexv1* out = new FullTrackVertexv1();
    out->set_track_id(trk->get_track_id());
    out->set_d0(fit.d0);
    out->set_timebin0(fit.timebin0);
    m_vertices->add_vertex(out);

    const unsigned int nlayers = count_unique_layers(pts);
    if (nlayers >= m_collisionMinTrackLayers)
    {
      CollisionTrack accepted;
      accepted.fit = fit;
      accepted.nlayers = nlayers;
      collision_tracks.push_back(accepted);
    }
  }

  const CollisionFit collision = fit_collision_point(collision_tracks);
  m_vertices->set_collision_min_layers(m_collisionMinTrackLayers);
  m_vertices->set_collision_ntracks(collision.ntracks);
  if (collision.ok)
  {
    m_vertices->set_collision_vertex_valid(1);
    m_vertices->set_collision_radius(collision.radius);
    m_vertices->set_collision_phi(collision.phi);
    m_vertices->set_collision_timebin(collision.timebin);
    m_vertices->set_collision_timebin_rms(collision.timebin_rms);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - input tracks=" << ntracks
              << " vertices=" << m_vertices->size()
              << " collision_ntracks=" << m_vertices->get_collision_ntracks()
              << " collision_radius=" << m_vertices->get_collision_radius()
              << " collision_phi=" << m_vertices->get_collision_phi()
              << " collision_timebin=" << m_vertices->get_collision_timebin()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

FullTrackVertexer::HitPoint
FullTrackVertexer::make_hit_point(const TrkrDefs::hitsetkey hsk,
                                  const TrkrDefs::hitkey hk) const
{
  HitPoint p;

  TrkrHitSet* hitset = m_hits ? m_hits->findHitSet(hsk) : nullptr;
  if (!hitset) return p;

  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return p;

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = TrkrDefs::getLayer(hsk);
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();

  if (p.layer < 7 || p.layer > 54) return p;
  if (!m_idealPadMap) return p;

  p.radius = m_idealPadMap->get_radius(p.layer);
  p.global_phi = wrap_phi(m_idealPadMap->get_phi(static_cast<unsigned int>(TpcDefs::getSide(hsk)), p.layer, p.pad));

  if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) return p;

  p.ok = true;
  return p;
}
