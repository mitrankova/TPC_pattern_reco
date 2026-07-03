#include "TpcPolyHelixFitter.h"

#include <trackbase/TrackFitUtils.h>

#include <ffamodules/CDBInterface.h>
#include <phfield/PHField3DCartesian.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Fit/Fitter.h>
#include <Math/Functor.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  double wrap_pi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  std::pair<double, double> find_root(const double q_over_r, const double x0, const double y0)
  {
    if (std::abs(q_over_r) <= 0.0 || !std::isfinite(q_over_r)) return {0.0, 0.0};

    const double r = std::abs(1.0 / q_over_r);
    const double den = square(x0) + square(y0);
    if (den <= 0.0 || !std::isfinite(den)) return {0.0, 0.0};

    const double rad = square(x0) * square(r) * square(y0) + square(r) * square(square(y0));
    if (rad < 0.0) return {0.0, 0.0};

    const double miny = (std::sqrt(rad) + square(x0) * y0 + y0 * square(y0)) / den;
    const double miny2 = (-std::sqrt(rad) + square(x0) * y0 + y0 * square(y0)) / den;
    const double dx1 = square(r) - square(miny - y0);
    const double dx2 = square(r) - square(miny2 - y0);
    if (dx1 < 0.0 || dx2 < 0.0) return {0.0, 0.0};

    const double minx = std::sqrt(dx1) + x0;
    const double minx2 = -std::sqrt(dx2) + x0;
    const double x = (std::abs(minx) < std::abs(minx2)) ? minx : minx2;
    const double y = (std::abs(miny) < std::abs(miny2)) ? miny : miny2;
    return {x, y};
  }

  double get_phi(const double q_over_r, const double x0, const double y0,
                 const std::vector<TpcPolyHelixFitter::Point>& points)
  {
    const auto root = find_root(q_over_r, x0, y0);
    double phi = std::atan2(-1.0 * (x0 - root.first), y0 - root.second);

    if (points.size() >= 2)
    {
      const auto& pos0 = points[0];
      const auto& pos1 = points[1];
      const double phi0 = std::atan2(pos0.y - y0, pos0.x - x0);
      const double phi1 = std::atan2(pos1.y - y0, pos1.x - x0);
      double dphi = phi1 - phi0;
      if (dphi > M_PI) dphi -= 2.0 * M_PI;
      if (dphi < -M_PI) dphi += 2.0 * M_PI;
      if (dphi < 0.0) phi += M_PI;
    }

    return wrap_pi(phi);
  }

  double distance(const TpcPolyHelixFitter::Point& point, const std::array<double, 3>& traj)
  {
    return square(point.x - traj[0]) + square(point.y - traj[1]) + square(point.z - traj[2]);
  }
}

TpcPolyHelixFitter::TpcPolyHelixFitter()
  : m_cdb(CDBInterface::instance())
{
}

TpcPolyHelixFitter::~TpcPolyHelixFitter()
{
  delete m_field;
  m_field = nullptr;
}

double TpcPolyHelixFitter::State::p() const
{
  return std::sqrt(px * px + py * py + pz * pz);
}

bool TpcPolyHelixFitter::InitField(int verbosity)
{
  delete m_field;
  m_field = nullptr;
  m_verbosity = verbosity;

  if (!m_cdb) m_cdb = CDBInterface::instance();
  if (!m_cdb) return false;

  std::string url = m_cdb->getUrl("FIELDMAP_TRACKING");
  if (url.empty())
  {
    std::cerr << "TpcPolyHelixFitter::InitField - empty FIELDMAP_TRACKING url" << std::endl;
    return false;
  }

  PHField3DCartesian* field = new PHField3DCartesian(url, 1.0);
  field->Verbosity(verbosity);
  m_field = field;
  return true;
}

void TpcPolyHelixFitter::getFieldTesla(const double point_cm[4], double field_t[3]) const
{
  field_t[0] = 0.0;
  field_t[1] = 0.0;
  field_t[2] = 0.0;
  if (!m_field) return;

  const double point[4] = {
      point_cm[0] * CLHEP::cm,
      point_cm[1] * CLHEP::cm,
      point_cm[2] * CLHEP::cm,
      point_cm[3]};
  double bfield[3] = {0.0, 0.0, 0.0};
  m_field->GetFieldValue(point, bfield);

  field_t[0] = bfield[0] / CLHEP::tesla;
  field_t[1] = bfield[1] / CLHEP::tesla;
  field_t[2] = bfield[2] / CLHEP::tesla;
}

void TpcPolyHelixFitter::rk4Step(double pos[3], double dir[3], double ds, double p, double q) const
{
  auto deriv = [this, p, q](const double pos_in[3], const double dir_in[3], double out_pos[3], double out_dir[3]) {
    const double point[4] = {pos_in[0], pos_in[1], pos_in[2], 0.0};
    double field_t[3] = {0.0, 0.0, 0.0};
    getFieldTesla(point, field_t);

    const double k = q * kConv / p;

    out_pos[0] = dir_in[0];
    out_pos[1] = dir_in[1];
    out_pos[2] = dir_in[2];

    out_dir[0] = k * (dir_in[1] * field_t[2] - dir_in[2] * field_t[1]);
    out_dir[1] = k * (dir_in[2] * field_t[0] - dir_in[0] * field_t[2]);
    out_dir[2] = k * (dir_in[0] * field_t[1] - dir_in[1] * field_t[0]);
  };

  double k1p[3] = {0.0, 0.0, 0.0};
  double k1d[3] = {0.0, 0.0, 0.0};
  double k2p[3] = {0.0, 0.0, 0.0};
  double k2d[3] = {0.0, 0.0, 0.0};
  double k3p[3] = {0.0, 0.0, 0.0};
  double k3d[3] = {0.0, 0.0, 0.0};
  double k4p[3] = {0.0, 0.0, 0.0};
  double k4d[3] = {0.0, 0.0, 0.0};
  double tmp_pos[3] = {0.0, 0.0, 0.0};
  double tmp_dir[3] = {0.0, 0.0, 0.0};

  deriv(pos, dir, k1p, k1d);

  for (int i = 0; i < 3; ++i)
  {
    tmp_pos[i] = pos[i] + 0.5 * ds * k1p[i];
    tmp_dir[i] = dir[i] + 0.5 * ds * k1d[i];
  }
  deriv(tmp_pos, tmp_dir, k2p, k2d);

  for (int i = 0; i < 3; ++i)
  {
    tmp_pos[i] = pos[i] + 0.5 * ds * k2p[i];
    tmp_dir[i] = dir[i] + 0.5 * ds * k2d[i];
  }
  deriv(tmp_pos, tmp_dir, k3p, k3d);

  for (int i = 0; i < 3; ++i)
  {
    tmp_pos[i] = pos[i] + ds * k3p[i];
    tmp_dir[i] = dir[i] + ds * k3d[i];
  }
  deriv(tmp_pos, tmp_dir, k4p, k4d);

  for (int i = 0; i < 3; ++i)
  {
    pos[i] += (ds / 6.0) * (k1p[i] + 2.0 * k2p[i] + 2.0 * k3p[i] + k4p[i]);
    dir[i] += (ds / 6.0) * (k1d[i] + 2.0 * k2d[i] + 2.0 * k3d[i] + k4d[i]);
  }

  const double norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if (norm > 0.0)
  {
    dir[0] /= norm;
    dir[1] /= norm;
    dir[2] /= norm;
  }
}

std::vector<std::array<double, 3>> TpcPolyHelixFitter::propagate(const State& state, double length_cm) const
{
  std::vector<std::array<double, 3>> traj;
  const double p = state.p();
  if (p <= 0.0 || length_cm < 0.0 || m_stepSize <= 0.0) return traj;

  double pos[3] = {state.x, state.y, state.z};
  double dir[3] = {state.px / p, state.py / p, state.pz / p};

  traj.push_back({pos[0], pos[1], pos[2]});

  double s = 0.0;
  while (s < length_cm)
  {
    const double step = std::min(m_stepSize, length_cm - s);
    rk4Step(pos, dir, step, p, state.charge);
    s += step;
    traj.push_back({pos[0], pos[1], pos[2]});
  }

  return traj;
}

// -----------------------------------------------------------------------
// Propagates `state` BACKWARD (flip direction, propagate forward with
// the same stepper -- standard trick for reverse propagation) and finds
// the arc-length point where r^2 = x^2+y^2 is minimized, i.e. the PCA to
// the z-axis. Refines with a local quadratic fit around the best sample.
// -----------------------------------------------------------------------
bool TpcPolyHelixFitter::propagateToPCA(const State& state, State& pca_state, double maxLength_cm) const
{
  const double p = state.p();
  if (p <= 0.0 || !std::isfinite(p)) return false;

  State rev = state;
  rev.px = -state.px;
  rev.py = -state.py;
  rev.pz = -state.pz;

  const std::vector<std::array<double, 3>> traj = propagate(rev, maxLength_cm);
  if (traj.size() < 2) return false;

  size_t best_i = 0;
  double best_r2 = square(traj[0][0]) + square(traj[0][1]);
  for (size_t i = 1; i < traj.size(); ++i)
  {
    const double r2 = square(traj[i][0]) + square(traj[i][1]);
    if (r2 < best_r2)
    {
      best_r2 = r2;
      best_i = i;
    }
  }

  if (best_i + 1 >= traj.size() && Verbosity() > 0)
  {
    std::cerr << "TpcPolyHelixFitter::propagateToPCA - PCA search hit maxLength_cm="
              << maxLength_cm << " without a clear minimum; result may be unreliable"
              << std::endl;
  }

  double s_offset = static_cast<double>(best_i) * m_stepSize;
  if (best_i > 0 && best_i + 1 < traj.size())
  {
    const double r2m = square(traj[best_i - 1][0]) + square(traj[best_i - 1][1]);
    const double r2c = best_r2;
    const double r2p = square(traj[best_i + 1][0]) + square(traj[best_i + 1][1]);
    const double denom = (r2m - 2.0 * r2c + r2p);
    if (std::abs(denom) > 1.0e-12)
    {
      const double delta = 0.5 * (r2m - r2p) / denom;
      s_offset += std::clamp(delta, -1.0, 1.0) * m_stepSize;
    }
  }

  double pos[3] = {rev.x, rev.y, rev.z};
  double dir[3] = {rev.px / p, rev.py / p, rev.pz / p};
  double s = 0.0;
  while (s < s_offset)
  {
    const double step = std::min(m_stepSize, s_offset - s);
    rk4Step(pos, dir, step, p, rev.charge);
    s += step;
  }

  pca_state.x = pos[0];
  pca_state.y = pos[1];
  pca_state.z = pos[2];
  pca_state.px = -dir[0] * p;
  pca_state.py = -dir[1] * p;
  pca_state.pz = -dir[2] * p;
  pca_state.charge = state.charge;

  return true;
}

double TpcPolyHelixFitter::chi2ForParams(const std::vector<Point>& points, const double* pars) const
{
  State state;
  state.x = pars[0];
  state.y = pars[1];
  state.z = pars[2];
  const double phi = pars[3];
  const double tan_lambda = pars[4];
  const double q_over_pt = pars[5];
  if (!std::isfinite(q_over_pt)) return std::numeric_limits<double>::max();

  state.charge = (q_over_pt >= 0.0) ? 1.0 : -1.0;
  const double abs_q_over_pt = std::max(std::abs(q_over_pt), 1.0e-12);
  const double pt = 1.0 / abs_q_over_pt;
  state.px = pt * std::cos(phi);
  state.py = pt * std::sin(phi);
  state.pz = pt * tan_lambda;

  double max_dist = 0.0;
  for (const Point& hit : points)
  {
    const double dx = hit.x - state.x;
    const double dy = hit.y - state.y;
    const double dz = hit.z - state.z;
    max_dist = std::max(max_dist, std::sqrt(dx * dx + dy * dy + dz * dz));
  }

  const std::vector<std::array<double, 3>> traj = propagate(state, 1.5 * max_dist + 5.0);
  if (traj.size() < 2) return std::numeric_limits<double>::max();

  double chi2 = 0.0;
  for (const Point& hit : points)
  {
    double best = std::numeric_limits<double>::max();
    for (const std::array<double, 3>& pos : traj) best = std::min(best, distance(hit, pos));

    constexpr double sigma = 0.02;
    chi2 += best / square(sigma);
  }

  return chi2;
}

TpcPolyHelixFitter::State TpcPolyHelixFitter::makeSeedState(const std::vector<Point>& points, const FitResult& seed) const
{
  State state;
  if (points.empty()) return state;

  state.x = points.front().x;
  state.y = points.front().y;
  state.z = points.front().z;

  const double point[4] = {state.x, state.y, state.z, 0.0};
  double field_t[3] = {0.0, 0.0, 0.0};
  getFieldTesla(point, field_t);

  double bz_for_q_over_pt = field_t[2];
  if (std::abs(bz_for_q_over_pt) < 1.0e-6)
  {
    const double bmag = std::sqrt(square(field_t[0]) + square(field_t[1]) + square(field_t[2]));
    bz_for_q_over_pt = std::copysign((bmag > 1.0e-6) ? bmag : 1.4,
                                     (seed.curvature != 0.0) ? seed.curvature : 1.0);
  }

  double q_over_pt = seed.curvature / (kConv * bz_for_q_over_pt);
  if (!std::isfinite(q_over_pt) || std::abs(q_over_pt) < 1.0e-12)
  {
    q_over_pt = std::copysign(1.0e-12, (q_over_pt != 0.0) ? q_over_pt : seed.curvature);
  }

  const double pt = 1.0 / std::abs(q_over_pt);
  const double theta = std::clamp(seed.theta, 1.0e-3, M_PI - 1.0e-3);
  const double tan_lambda = 1.0 / std::tan(theta);
  state.px = pt * std::cos(seed.phi0);
  state.py = pt * std::sin(seed.phi0);
  state.pz = pt * tan_lambda;
  state.charge = (q_over_pt >= 0.0) ? 1.0 : -1.0;

  return state;
}

TpcPolyHelixFitter::FieldFitResult TpcPolyHelixFitter::fieldFit(const std::vector<Point>& points, const FitResult& seed) const
{
  FieldFitResult result;
  if (!m_field || points.size() < 3) return result;

  const State guess = makeSeedState(points, seed);
  const double p0 = guess.p();
  if (p0 <= 0.0 || !std::isfinite(p0)) return result;

  const double phi0 = std::atan2(guess.py, guess.px);
  const double pt0 = std::sqrt(square(guess.px) + square(guess.py));
  if (pt0 <= 0.0 || !std::isfinite(pt0)) return result;
  const double tan_lambda0 = guess.pz / pt0;
  const double q_over_pt0 = guess.charge / pt0;
  double start_pars[6] = {guess.x, guess.y, guess.z, phi0, tan_lambda0, q_over_pt0};

  ROOT::Math::Functor fcn(
      [this, &points](const double* pars) { return chi2ForParams(points, pars); }, 6);

  ROOT::Fit::Fitter fitter;
  fitter.Config().MinimizerOptions().SetMaxIterations(m_maxIterations);
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(m_maxIterations * 20);
  fitter.SetFCN(fcn, start_pars);

  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("z0");
  fitter.Config().ParSettings(3).SetName("phi");
  fitter.Config().ParSettings(4).SetName("tanLambda");
  fitter.Config().ParSettings(5).SetName("qOverPt");

  if (!fitter.FitFCN()) return result;

  const ROOT::Fit::FitResult& fit_result = fitter.Result();
  const double* pars = fit_result.GetParams();
  const double q_over_pt = pars[5];
  const double abs_q_over_pt = std::max(std::abs(q_over_pt), 1.0e-12);
  const double pt = 1.0 / abs_q_over_pt;

  result.state.x = pars[0];
  result.state.y = pars[1];
  result.state.z = pars[2];
  result.state.px = pt * std::cos(pars[3]);
  result.state.py = pt * std::sin(pars[3]);
  result.state.pz = pt * pars[4];
  result.state.charge = (q_over_pt >= 0.0) ? 1.0 : -1.0;
  result.chi2 = fit_result.Chi2();
  result.ndf = static_cast<int>(points.size()) - 6;
  result.valid = true;
  return result;
}

bool TpcPolyHelixFitter::algebraicFit(const std::vector<Point>& points, FitResult& fit_result) const
{
  fit_result = FitResult();
  if (points.size() < 3) return false;

  TrackFitUtils::position_vector_t xy_positions;
  TrackFitUtils::position_vector_t rz_positions;
  xy_positions.reserve(points.size());
  rz_positions.reserve(points.size());
  for (const Point& point : points)
  {
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z)) continue;
    xy_positions.emplace_back(point.x, point.y);
    rz_positions.emplace_back(std::sqrt(square(point.x) + square(point.y)), point.z);
  }

  if (xy_positions.size() < 3 || rz_positions.size() < 2) return false;

  const auto circle = TrackFitUtils::circle_fit_by_taubin(xy_positions);
  const double radius = std::get<0>(circle);
  const double x0 = std::get<1>(circle);
  const double y0 = std::get<2>(circle);
  if (!std::isfinite(radius) || !std::isfinite(x0) || !std::isfinite(y0) || radius <= 0.0) return false;

  double q_over_r = 1.0 / radius;
  const auto& firstpos = xy_positions.front();
  const auto& lastpos = xy_positions.back();
  const double firstphi = std::atan2(firstpos.second, firstpos.first);
  const double lastphi = std::atan2(lastpos.second, lastpos.first);
  double dphi = lastphi - firstphi;
  if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;
  if (dphi < -M_PI) dphi = 2.0 * M_PI + dphi;
  if (dphi > 0.0) q_over_r *= -1.0;

  const auto line = TrackFitUtils::line_fit(rz_positions);
  const double slope = std::get<0>(line);
  const double intercept = std::get<1>(line);
  if (!std::isfinite(slope) || !std::isfinite(intercept)) return false;

  fit_result.radius = radius;
  fit_result.x0 = x0;
  fit_result.y0 = y0;
  fit_result.curvature = q_over_r;
  fit_result.slope = slope;
  fit_result.z0 = intercept;
  fit_result.theta = std::atan2(1.0, slope);
  fit_result.phi0 = get_phi(q_over_r, x0, y0, points);
  const auto root = find_root(q_over_r, x0, y0);
  fit_result.x = root.first;
  fit_result.y = root.second;
  fit_result.d0 = std::copysign(std::sqrt(square(root.first) + square(root.second)), q_over_r);

  fit_result.chi2_xy = 0.0;
  for (const auto& xy : xy_positions)
  {
    const double residual = std::sqrt(square(xy.first - x0) + square(xy.second - y0)) - radius;
    fit_result.chi2_xy += residual * residual;
  }

  fit_result.chi2_z = 0.0;
  for (const auto& rz : rz_positions)
  {
    const double residual = rz.second - (slope * rz.first + intercept);
    fit_result.chi2_z += residual * residual;
  }

  fit_result.ndof_xy = static_cast<int>(xy_positions.size()) - 3;
  fit_result.ndof_z = static_cast<int>(rz_positions.size()) - 2;
  fit_result.ok = true;
  return true;
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points, FitResult& fit_result) const
{
  if (!algebraicFit(points, fit_result)) return false;
  if (!m_field) return true;

  const FieldFitResult field_result = fieldFit(points, fit_result);
  if (!field_result.valid) return true;

  State pca_state;
  const bool have_pca = propagateToPCA(field_result.state, pca_state);
  const State& ref = have_pca ? pca_state : field_result.state;
  if (!have_pca && Verbosity() > 0)
  {
    std::cerr << "TpcPolyHelixFitter::fit - propagateToPCA failed, "
              << "falling back to fit reference point (d0/z0 will be unreliable)"
              << std::endl;
  }

  const double p = ref.p();
  const double pt = std::sqrt(square(ref.px) + square(ref.py));
  if (p <= 0.0 || pt <= 0.0) return true;

  const double field_point[4] = {ref.x, ref.y, ref.z, 0.0};
  double field_t[3] = {0.0, 0.0, 0.0};
  getFieldTesla(field_point, field_t);

  fit_result.x = ref.x;
  fit_result.y = ref.y;
  fit_result.z0 = ref.z;
  fit_result.phi0 = std::atan2(ref.py, ref.px);
  fit_result.theta = std::acos(std::clamp(ref.pz / p, -1.0, 1.0));
  fit_result.slope = (std::sin(fit_result.theta) != 0.0) ? std::cos(fit_result.theta) / std::sin(fit_result.theta) : 0.0;
  fit_result.curvature = ref.charge * kConv * field_t[2] / pt;
  if (std::abs(fit_result.curvature) > 1.0e-12)
  {
    fit_result.radius = std::abs(1.0 / fit_result.curvature);
  }
  fit_result.d0 = std::copysign(std::sqrt(square(fit_result.x) + square(fit_result.y)), fit_result.curvature);
  fit_result.pt = pt;
  fit_result.p = p;
  fit_result.chi2_xy = field_result.chi2;
  fit_result.chi2_z = 0.0;
  fit_result.ndof_xy = field_result.ndf;
  fit_result.ndof_z = 0;
  fit_result.ok = true;
  return true;
}
