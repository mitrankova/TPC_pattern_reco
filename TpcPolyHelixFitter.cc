#include "TpcPolyHelixFitter.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
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

  bool solve_3x3(double A[3][3], double b[3], double x[3])
  {
    double M[3][4] = {
      {A[0][0], A[0][1], A[0][2], b[0]},
      {A[1][0], A[1][1], A[1][2], b[1]},
      {A[2][0], A[2][1], A[2][2], b[2]}
    };
    for (int col = 0; col < 3; ++col)
    {
      int pivot = col;
      for (int row = col + 1; row < 3; ++row)
        if (std::fabs(M[row][col]) > std::fabs(M[pivot][col])) pivot = row;
      if (std::fabs(M[pivot][col]) < 1.0e-20) return false;
      if (pivot != col)
        for (int k = col; k < 4; ++k) std::swap(M[col][k], M[pivot][k]);
      const double div = M[col][col];
      for (int k = col; k < 4; ++k) M[col][k] /= div;
      for (int row = 0; row < 3; ++row)
      {
        if (row == col) continue;
        const double factor = M[row][col];
        for (int k = col; k < 4; ++k) M[row][k] -= factor * M[col][k];
      }
    }
    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
  }

  bool line_fit(const std::vector<double>& x, const std::vector<double>& y,
                double& slope, double& intercept, double& chi2, int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size()) return false;
    double S = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      S += 1.0;
      Sx += x[i];
      Sy += y[i];
      Sxx += x[i] * x[i];
      Sxy += x[i] * y[i];
    }
    const double den = S * Sxx - Sx * Sx;
    if (std::fabs(den) < 1.0e-20) return false;
    slope = (S * Sxy - Sx * Sy) / den;
    intercept = (Sy - slope * Sx) / S;
    chi2 = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double r = y[i] - (slope * x[i] + intercept);
      chi2 += r * r;
    }
    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }

  bool algebraic_circle(const std::vector<TpcPolyHelixFitter::Point>& points,
                        double& xc, double& yc, double& R)
  {
    double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double b[3] = {0.0, 0.0, 0.0};
    for (const auto& p : points)
    {
      const double row[3] = {p.x, p.y, 1.0};
      const double rhs = -(p.x * p.x + p.y * p.y);
      for (int i = 0; i < 3; ++i)
      {
        b[i] += row[i] * rhs;
        for (int j = 0; j < 3; ++j) A[i][j] += row[i] * row[j];
      }
    }

    double sol[3] = {0.0, 0.0, 0.0};
    if (!solve_3x3(A, b, sol)) return false;
    xc = -0.5 * sol[0];
    yc = -0.5 * sol[1];
    const double r2 = xc * xc + yc * yc - sol[2];
    if (r2 <= 0.0) return false;
    R = std::sqrt(r2);
    return std::isfinite(xc) && std::isfinite(yc) && std::isfinite(R) && R > 0.0;
  }

  bool valid_seed(const TpcPolyHelixFitter::Seed& seed)
  {
    return seed.valid && std::isfinite(seed.x) && std::isfinite(seed.y) && std::isfinite(seed.z) &&
           std::isfinite(seed.px) && std::isfinite(seed.py) && std::isfinite(seed.pz) &&
           std::hypot(seed.px, seed.py) > 1.0e-12;
  }

  void seed_circle(const TpcPolyHelixFitter::Seed& seed, const double ref_xc, const double ref_yc,
                   const double ref_R, double& xc, double& yc, double& R)
  {
    xc = ref_xc;
    yc = ref_yc;
    R = ref_R;
    if (!valid_seed(seed)) return;

    const double pt = std::hypot(seed.px, seed.py);
    const double nx = -seed.py / pt;
    const double ny = seed.px / pt;
    const double cx1 = seed.x + nx * ref_R;
    const double cy1 = seed.y + ny * ref_R;
    const double cx2 = seed.x - nx * ref_R;
    const double cy2 = seed.y - ny * ref_R;
    const double d1 = std::hypot(cx1 - ref_xc, cy1 - ref_yc);
    const double d2 = std::hypot(cx2 - ref_xc, cy2 - ref_yc);
    xc = (d1 <= d2) ? cx1 : cx2;
    yc = (d1 <= d2) ? cy1 : cy2;
  }

  bool refine_circle(const std::vector<TpcPolyHelixFitter::Point>& points,
                     double& xc, double& yc, double& R)
  {
    for (int iter = 0; iter < 20; ++iter)
    {
      double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double b[3] = {0.0, 0.0, 0.0};
      double chi2 = 0.0;
      for (const auto& p : points)
      {
        const double dx = p.x - xc;
        const double dy = p.y - yc;
        const double rho = std::hypot(dx, dy);
        if (rho <= 1.0e-12) return false;
        const double resid = rho - R;
        const double J[3] = {-dx / rho, -dy / rho, -1.0};
        for (int i = 0; i < 3; ++i)
        {
          b[i] += -J[i] * resid;
          for (int j = 0; j < 3; ++j) A[i][j] += J[i] * J[j];
        }
        chi2 += resid * resid;
      }

      double step[3] = {0.0, 0.0, 0.0};
      if (!solve_3x3(A, b, step)) return false;
      double scale = 1.0;
      bool accepted = false;
      for (int trial = 0; trial < 8; ++trial)
      {
        const double txc = xc + scale * step[0];
        const double tyc = yc + scale * step[1];
        const double tR = R + scale * step[2];
        if (tR > 0.0 && std::isfinite(txc) && std::isfinite(tyc) && std::isfinite(tR))
        {
          double trial_chi2 = 0.0;
          for (const auto& p : points)
          {
            const double r = std::hypot(p.x - txc, p.y - tyc) - tR;
            trial_chi2 += r * r;
          }
          if (trial_chi2 <= chi2)
          {
            xc = txc;
            yc = tyc;
            R = tR;
            accepted = true;
            break;
          }
        }
        scale *= 0.5;
      }
      if (!accepted) return true;
      if (std::hypot(step[0], step[1]) + std::fabs(step[2]) < 1.0e-7) return true;
    }
    return true;
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

  double ref_xc = 0.0;
  double ref_yc = 0.0;
  double ref_R = 0.0;
  if (!algebraic_circle(points, ref_xc, ref_yc, ref_R)) return false;

  double xc = 0.0;
  double yc = 0.0;
  double R = 0.0;
  seed_circle(seed, ref_xc, ref_yc, ref_R, xc, yc, R);
  if (!refine_circle(points, xc, yc, R)) return false;
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
  std::vector<double> svals;
  std::vector<double> zvals;
  svals.reserve(points.size());
  zvals.reserve(points.size());
  double prev_angle = phi_perigee;
  for (const Point& p : points)
  {
    double angle = std::atan2(p.y - yc, p.x - xc);
    angle = unwrap_near(angle, prev_angle);
    prev_angle = angle;
    double dangle = angle - phi_perigee;
    if (sign * dangle < 0.0) dangle += sign * 2.0 * M_PI;
    svals.push_back(sign * R * dangle);
    zvals.push_back(p.z);

    const double resid = std::hypot(p.x - xc, p.y - yc) - R;
    fit_result.chi2_xy += resid * resid;
  }

  if (!line_fit(svals, zvals, fit_result.dzds, fit_result.z0, fit_result.chi2_z, fit_result.ndof_z)) return false;
  fit_result.theta = std::atan2(1.0, fit_result.dzds);
  fit_result.ndof_xy = static_cast<int>(points.size()) - 3;
  fit_result.ok = true;
  return true;
}
