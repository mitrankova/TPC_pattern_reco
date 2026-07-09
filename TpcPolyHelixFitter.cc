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
}

bool TpcPolyHelixFitter::fit(const std::vector<Point>& points, FitResult& fit)
{
  fit = FitResult();
  if (points.size() < 3) return false;

  double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double b[3] = {0.0, 0.0, 0.0};
  for (const Point& p : points)
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

  const double xc = -0.5 * sol[0];
  const double yc = -0.5 * sol[1];
  const double r2 = xc * xc + yc * yc - sol[2];
  if (r2 <= 0.0) return false;
  const double R = std::sqrt(r2);
  const double dc = std::hypot(xc, yc);
  if (R <= 0.0 || dc <= 1.0e-12) return false;

  const Point& first = points.front();
  const Point& mid = points[points.size() / 2];
  const Point& last = points.back();
  const double cross = (mid.x - first.x) * (last.y - mid.y) - (mid.y - first.y) * (last.x - mid.x);
  const double sign = (cross >= 0.0) ? 1.0 : -1.0;

  fit.curvature = sign / R;
  fit.d0 = sign * (dc - R);

  const double px = xc * (1.0 - R / dc);
  const double py = yc * (1.0 - R / dc);
  const double rx = px - xc;
  const double ry = py - yc;
  const double tx = -sign * ry / R;
  const double ty =  sign * rx / R;
  fit.phi0 = wrap_pi(std::atan2(ty, tx));

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
    svals.push_back(std::fabs(R * dangle));
    zvals.push_back(p.z);

    const double resid = std::hypot(p.x - xc, p.y - yc) - R;
    fit.chi2_xy += resid * resid;
  }

  double dzds = 0.0;
  if (!line_fit(svals, zvals, dzds, fit.z0, fit.chi2_z, fit.ndof_z)) return false;
  fit.theta = std::atan2(1.0, dzds);
  fit.ndof_xy = static_cast<int>(points.size()) - 3;
  fit.ok = true;
  return true;
}
