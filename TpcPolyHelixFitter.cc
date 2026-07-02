#include "TpcPolyHelixFitter.h"

#include "Fitter.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
  constexpr int kIrlsIterations = 2;
  // Tukey biweight tuning constant giving ~95% efficiency for Gaussian residuals.
  constexpr double kTukeyC = 4.685;
  // MAD-to-sigma conversion factor for a Gaussian distribution.
  constexpr double kMadToSigma = 1.4826;
  constexpr double kMinRobustScaleCm = 0.01;
  constexpr double kRejectedWeight = 1.0e-9;

  double wrap_pi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double unwrap_near(double phi, const double ref)
  {
    while (phi - ref > M_PI) phi -= 2.0 * M_PI;
    while (phi - ref < -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  bool valid_seed(const TpcPolyHelixFitter::Seed& seed)
  {
    return seed.valid && std::isfinite(seed.x) && std::isfinite(seed.y) && std::isfinite(seed.z) &&
           std::isfinite(seed.px) && std::isfinite(seed.py) && std::isfinite(seed.pz) &&
           std::hypot(seed.px, seed.py) > 1.0e-12;
  }

  // Robust re-weighting shared by the circle and z-line IRLS loops: combines each
  // point's external weight with a Tukey biweight derived from the current residual.
  void update_robust_weights(const std::vector<double>& resid, const std::vector<double>& base_w,
                             std::vector<double>& w)
  {
    const unsigned int n = static_cast<unsigned int>(resid.size());
    std::vector<double> sorted_abs(n);
    for (unsigned int i = 0; i < n; ++i) sorted_abs[i] = std::fabs(resid[i]);
    std::sort(sorted_abs.begin(), sorted_abs.end());
    const double mad = sorted_abs[n / 2];
    double scale = kMadToSigma * mad;
    // Floor at ~typical TPC hit resolution (cm) so a mostly-perfect fit doesn't
    // make the robust scale collapse to numerical noise and flag every point as
    // an outlier (which stalls IRLS and biases the result toward the outliers).
    if (scale < kMinRobustScaleCm) scale = kMinRobustScaleCm;

    // Fitter's weighted fits treat any weight <= 0 as "unspecified" and silently
    // substitute 1.0 (full weight), so rejected outliers must get a tiny positive
    // weight rather than an exact zero or they'd be reinstated at full strength.
    w.resize(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      const double u = resid[i] / (kTukeyC * scale);
      w[i] = (std::fabs(u) < 1.0) ? base_w[i] * (1.0 - u * u) * (1.0 - u * u) : kRejectedWeight;
    }
  }

  // Lets IRLS exit as soon as reweighting stops changing anything, instead of
  // always spending the full iteration budget on data with no outliers.
  bool weights_converged(const std::vector<double>& w_old, const std::vector<double>& w_new)
  {
    for (unsigned int i = 0; i < w_old.size(); ++i)
      if (std::fabs(w_new[i] - w_old[i]) > 1.0e-3) return false;
    return true;
  }

  // Rough, singularity-free direction estimate used only to pre-rotate the point
  // cloud so the sagitta model's internal line pre-fit never sees a near-vertical slope.
  double initial_direction(const std::vector<TpcPolyHelixFitter::Point>& points,
                           const TpcPolyHelixFitter::Seed& seed)
  {
    if (valid_seed(seed)) return std::atan2(seed.py, seed.px);

    double cx = 0.0, cy = 0.0;
    for (const auto& p : points)
    {
      cx += p.x;
      cy += p.y;
    }
    cx /= static_cast<double>(points.size());
    cy /= static_cast<double>(points.size());

    unsigned int ifar = 0;
    double best = -1.0;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      const double d = std::hypot(points[i].x - cx, points[i].y - cy);
      if (d > best)
      {
        best = d;
        ifar = i;
      }
    }

    unsigned int jfar = 0;
    best = -1.0;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      const double d = std::hypot(points[i].x - points[ifar].x, points[i].y - points[ifar].y);
      if (d > best)
      {
        best = d;
        jfar = i;
      }
    }

    return std::atan2(points[jfar].y - points[ifar].y, points[jfar].x - points[ifar].x);
  }

  // Circle-as-sagitta fit: pre-rotate into a frame aligned with the track so the
  // model stays well conditioned even for near-straight (large-R) tracks, fit with
  // Fitter::weightedSagittaFit, then map the local sagitta parameters back to a
  // global circle center/radius. Outliers are down-weighted via IRLS.
  bool fit_circle_sagitta(const std::vector<TpcPolyHelixFitter::Point>& points,
                          const TpcPolyHelixFitter::Seed& seed,
                          double& xc, double& yc, double& R)
  {
    const unsigned int n = static_cast<unsigned int>(points.size());
    const double psi0 = initial_direction(points, seed);
    const double cpsi = std::cos(psi0);
    const double spsi = std::sin(psi0);

    std::vector<double> xp(n), yp(n), base_w(n), w(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      xp[i] = cpsi * points[i].x + spsi * points[i].y;
      yp[i] = -spsi * points[i].x + cpsi * points[i].y;
      base_w[i] = points[i].w > 0.0 ? points[i].w : 1.0;
      w[i] = base_w[i];
    }

    bool ok = false;
    for (int iter = 0; iter < kIrlsIterations; ++iter)
    {
      double S = 0.0, x0 = 0.0, invR = 0.0, theta = 0.0, bline = 0.0, chi2 = 0.0;
      int ndof = 0;
      ok = Fitter::weightedSagittaFit(xp, yp, w, S, x0, invR, theta, bline, chi2, ndof);
      if (!ok || std::fabs(invR) < 1.0e-12) return false;

      const double Rs = 1.0 / invR;
      const double c = std::cos(theta);
      const double s = std::sin(theta);
      const double xpc = c * x0 - s * (S - Rs);
      const double ypc = bline + s * x0 + c * (S - Rs);
      const double txc = cpsi * xpc - spsi * ypc;
      const double tyc = spsi * xpc + cpsi * ypc;
      const double tR = std::fabs(Rs);
      if (!std::isfinite(txc) || !std::isfinite(tyc) || !std::isfinite(tR) || tR <= 0.0) return false;

      xc = txc;
      yc = tyc;
      R = tR;
      if (iter == kIrlsIterations - 1) break;

      std::vector<double> resid(n);
      for (unsigned int i = 0; i < n; ++i) resid[i] = std::hypot(points[i].x - xc, points[i].y - yc) - R;

      std::vector<double> w_new;
      update_robust_weights(resid, base_w, w_new);
      if (weights_converged(w, w_new)) break;

      unsigned int nactive = 0;
      for (double wi : w_new) if (wi > 1.0e-6) ++nactive;
      if (nactive < 3) break;
      w = w_new;
    }
    return ok;
  }

  bool robust_weighted_line_fit(const std::vector<double>& x, const std::vector<double>& y,
                                const std::vector<double>& base_w,
                                double& slope, double& intercept, double& chi2, int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size() || x.size() != base_w.size()) return false;

    std::vector<double> w = base_w;
    bool ok = false;
    for (int iter = 0; iter < kIrlsIterations; ++iter)
    {
      ok = Fitter::weightedLineFit(x, y, w, slope, intercept, chi2, ndof);
      if (!ok) return false;
      if (iter == kIrlsIterations - 1) break;

      std::vector<double> resid(x.size());
      for (unsigned int i = 0; i < x.size(); ++i) resid[i] = y[i] - (slope * x[i] + intercept);

      std::vector<double> w_new;
      update_robust_weights(resid, base_w, w_new);
      if (weights_converged(w, w_new)) break;

      unsigned int nactive = 0;
      for (double wi : w_new) if (wi > 1.0e-6) ++nactive;
      if (nactive < 2) break;
      w = w_new;
    }
    return ok;
  }

  double direction_sign(const std::vector<TpcPolyHelixFitter::Point>& points,
                        const TpcPolyHelixFitter::Seed& seed,
                        const double xc, const double yc, const double R)
  {
    double sign = 1.0;
    if (valid_seed(seed))
    {
      double best = 1.0e100;
      unsigned int idx = 0;
      for (unsigned int i = 0; i < points.size(); ++i)
      {
        const double d = std::hypot(points[i].x - seed.x, points[i].y - seed.y);
        if (d < best)
        {
          best = d;
          idx = i;
        }
      }
      const double rx = points[idx].x - xc;
      const double ry = points[idx].y - yc;
      const double tx_ccw = -ry / R;
      const double ty_ccw = rx / R;
      sign = (tx_ccw * seed.px + ty_ccw * seed.py >= 0.0) ? 1.0 : -1.0;
    }
    else
    {
      const auto& first = points.front();
      const auto& mid = points[points.size() / 2];
      const auto& last = points.back();
      const double cross = (mid.x - first.x) * (last.y - mid.y) - (mid.y - first.y) * (last.x - mid.x);
      sign = (cross >= 0.0) ? 1.0 : -1.0;
    }
    return sign;
  }
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points, FitResult& fit_result)
{
  Seed seed;
  return fit(points, seed, fit_result);
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points, const Seed& seed, FitResult& fit_result)
{
  fit_result = FitResult();
  if (points.size() < 3) return false;

  double xc = 0.0;
  double yc = 0.0;
  double R = 0.0;
  if (!fit_circle_sagitta(points, seed, xc, yc, R)) return false;
  const double dc = std::hypot(xc, yc);
  if (R <= 0.0 || dc <= 1.0e-12) return false;

  const double sign = direction_sign(points, seed, xc, yc, R);
  fit_result.xc = xc;
  fit_result.yc = yc;
  fit_result.radius = R;
  fit_result.curvature = sign / R;
  fit_result.d0 = sign * (dc - R);

  const double px = xc * (1.0 - R / dc);
  const double py = yc * (1.0 - R / dc);
  const double rx = px - xc;
  const double ry = py - yc;
  fit_result.phi0 = wrap_pi(std::atan2(sign * rx / R, -sign * ry / R));

  const double phi_perigee = std::atan2(ry, rx);

  std::vector<unsigned int> order(points.size());
  for (unsigned int i = 0; i < points.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&points](unsigned int a, unsigned int b)
  {
    return std::hypot(points[a].x, points[a].y) < std::hypot(points[b].x, points[b].y);
  });

  std::vector<double> svals;
  std::vector<double> zvals;
  std::vector<double> wvals;
  svals.reserve(points.size());
  zvals.reserve(points.size());
  wvals.reserve(points.size());
  double prev_angle = phi_perigee;
  for (unsigned int idx : order)
  {
    const Point& p = points[idx];
    double angle = std::atan2(p.y - yc, p.x - xc);
    angle = unwrap_near(angle, prev_angle);
    prev_angle = angle;
    double dangle = angle - phi_perigee;
    if (sign * dangle < 0.0) dangle += sign * 2.0 * M_PI;
    svals.push_back(sign * R * dangle);
    zvals.push_back(p.z);
    wvals.push_back(p.w > 0.0 ? p.w : 1.0);

    const double resid = std::hypot(p.x - xc, p.y - yc) - R;
    fit_result.chi2_xy += resid * resid;
  }

  if (!robust_weighted_line_fit(svals, zvals, wvals, fit_result.dzds, fit_result.z0, fit_result.chi2_z, fit_result.ndof_z)) return false;
  fit_result.theta = std::atan2(1.0, fit_result.dzds);
  fit_result.ndof_xy = static_cast<int>(points.size()) - 3;
  fit_result.ok = true;
  return true;
}
