#include "Fitter.h"

#include <algorithm>
#include <cmath>

namespace
{
  bool solve_3x3(double A[3][3], double b[3], double x[3])
  {
    double M[3][4] =
    {
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
      {
        for (int k = col; k < 4; ++k)
          std::swap(M[col][k], M[pivot][k]);
      }

      const double div = M[col][col];
      for (int k = col; k < 4; ++k) M[col][k] /= div;

      for (int row = 0; row < 3; ++row)
      {
        if (row == col) continue;
        const double factor = M[row][col];
        for (int k = col; k < 4; ++k)
          M[row][k] -= factor * M[col][k];
      }
    }

    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
  }
}

namespace Fitter
{
  FitPoint::FitPoint()
    : x(0.0), y(0.0), w(1.0) {}

  FitPoint::FitPoint(double x_, double y_, double w_)
    : x(x_), y(y_), w(w_) {}

  LineFit::LineFit()
    : ok(false), slope(0.0), intercept(0.0), chi2(0.0), ndof(0) {}

  SagittaFit::SagittaFit()
    : ok(false), S(0.0), x0(0.0), invR(0.0), theta(0.0), b(0.0), chi2(0.0), ndof(0) {}

  double adcWeight(double adc, double maxadc, double power, double floor_frac)
  {
    if (adc <= 0.0) return 0.0;
    if (maxadc <= 0.0) return 1.0;

    double w = std::pow(adc / maxadc, power);
    if (w < floor_frac) w = floor_frac;
    return w;
  }

  bool weightedLineFit(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<double>& w,
                       double& m,
                       double& b,
                       double& chi2,
                       int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size() || x.size() != w.size())
      return false;

    double S = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;

    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      S   += wi;
      Sx  += wi * x[i];
      Sy  += wi * y[i];
      Sxx += wi * x[i] * x[i];
      Sxy += wi * x[i] * y[i];
    }

    const double den = S * Sxx - Sx * Sx;
    if (std::fabs(den) < 1.0e-12) return false;

    m = (S * Sxy - Sx * Sy) / den;
    b = (Sy - m * Sx) / S;

    chi2 = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double r = y[i] - (m * x[i] + b);
      chi2 += wi * r * r;
    }

    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }

  LineFit fitLine(const std::vector<FitPoint>& points)
  {
    LineFit fit;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> w;
    x.reserve(points.size());
    y.reserve(points.size());
    w.reserve(points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
    {
      x.push_back(points[i].x);
      y.push_back(points[i].y);
      w.push_back(points[i].w);
    }

    fit.ok = weightedLineFit(x, y, w, fit.slope, fit.intercept, fit.chi2, fit.ndof);
    return fit;
  }

  double sagittaModel(const double xrot,
                      const double S,
                      const double x0,
                      const double invR)
  {
    const double dx = xrot - x0;
    const double invR2 = invR * invR;
    const double dx2 = dx * dx;

    return S
      - 0.5 * invR * dx2
      - 0.125 * invR * invR2 * dx2 * dx2
      - 0.0625 * invR * invR2 * invR2 * dx2 * dx2 * dx2;
  }

  bool weightedSagittaFit(const std::vector<double>& local_x,
                          const std::vector<double>& local_y,
                          const std::vector<double>& w,
                          double& S,
                          double& x0,
                          double& invR,
                          double& theta,
                          double& bline,
                          double& chi2,
                          int& ndof)
  {
    if (local_x.size() < 3 || local_x.size() != local_y.size() || local_x.size() != w.size())
      return false;

    double mline = 0.0;
    double line_chi2 = 0.0;
    int line_ndof = 0;
    if (!weightedLineFit(local_x, local_y, w, mline, bline, line_chi2, line_ndof))
      return false;

    theta = std::atan(mline);
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    std::vector<double> xrot;
    std::vector<double> yrot;
    xrot.reserve(local_x.size());
    yrot.reserve(local_x.size());

    double sw = 0.0;
    double sx = 0.0;
    double sy = 0.0;

    for (unsigned int i = 0; i < local_x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double yy = local_y[i] - bline;
      const double xr =  c * local_x[i] + s * yy;
      const double yr = -s * local_x[i] + c * yy;

      xrot.push_back(xr);
      yrot.push_back(yr);

      sw += wi;
      sx += wi * xr;
      sy += wi * yr;
    }

    if (sw <= 0.0) return false;

    x0 = sx / sw;
    S  = sy / sw;

    std::vector<double> dx2;
    dx2.reserve(xrot.size());
    for (unsigned int i = 0; i < xrot.size(); ++i)
    {
      const double dx = xrot[i] - x0;
      dx2.push_back(dx * dx);
    }

    double q = 0.0;
    double qS = 0.0;
    double qchi2 = 0.0;
    int qndof = 0;
    if (weightedLineFit(dx2, yrot, w, q, qS, qchi2, qndof))
    {
      S = qS;
      invR = -2.0 * q;
    }
    else
    {
      invR = 0.0;
    }

    if (invR >  1.0) invR =  1.0;
    if (invR < -1.0) invR = -1.0;

    chi2 = 0.0;
    for (unsigned int i = 0; i < xrot.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double r = yrot[i] - sagittaModel(xrot[i], S, x0, invR);
      chi2 += wi * r * r;
    }

    double lambda = 1.0e-6;
    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      double A[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
      double bvec[3] = {0.0, 0.0, 0.0};

      for (unsigned int i = 0; i < xrot.size(); ++i)
      {
        const double wi = w[i] > 0.0 ? w[i] : 1.0;
        const double dx = xrot[i] - x0;
        const double dx2v = dx * dx;
        const double dx3 = dx2v * dx;
        const double dx4 = dx2v * dx2v;
        const double dx5 = dx4 * dx;
        const double dx6 = dx4 * dx2v;

        const double invR2 = invR * invR;
        const double invR3 = invR2 * invR;
        const double invR4 = invR2 * invR2;
        const double invR5 = invR4 * invR;

        const double f = S
          - 0.5 * invR * dx2v
          - 0.125 * invR3 * dx4
          - 0.0625 * invR5 * dx6;

        const double resid = yrot[i] - f;

        double J[3];
        J[0] = 1.0;
        J[1] = invR * dx + 0.5 * invR3 * dx3 + 0.375 * invR5 * dx5;
        J[2] = -0.5 * dx2v - 0.375 * invR2 * dx4 - 0.3125 * invR4 * dx6;

        for (int a = 0; a < 3; ++a)
        {
          bvec[a] += wi * J[a] * resid;
          for (int b = 0; b < 3; ++b)
            A[a][b] += wi * J[a] * J[b];
        }
      }

      A[0][0] += lambda;
      A[1][1] += lambda;
      A[2][2] += lambda;

      double delta[3] = {0.0, 0.0, 0.0};
      if (!solve_3x3(A, bvec, delta)) break;

      if (delta[0] >  10.0) delta[0] =  10.0;
      if (delta[0] < -10.0) delta[0] = -10.0;
      if (delta[1] >  10.0) delta[1] =  10.0;
      if (delta[1] < -10.0) delta[1] = -10.0;
      if (delta[2] >   0.1) delta[2] =   0.1;
      if (delta[2] <  -0.1) delta[2] =  -0.1;

      const double S_new = S + delta[0];
      const double x0_new = x0 + delta[1];
      double invR_new = invR + delta[2];

      if (invR_new >  1.0) invR_new =  1.0;
      if (invR_new < -1.0) invR_new = -1.0;

      double chi2_new = 0.0;
      for (unsigned int i = 0; i < xrot.size(); ++i)
      {
        const double wi = w[i] > 0.0 ? w[i] : 1.0;
        const double r = yrot[i] - sagittaModel(xrot[i], S_new, x0_new, invR_new);
        chi2_new += wi * r * r;
      }

      if (chi2_new <= chi2)
      {
        S = S_new;
        x0 = x0_new;
        invR = invR_new;
        chi2 = chi2_new;
        lambda *= 0.3;

        if (std::fabs(delta[0]) < 1.0e-6 &&
            std::fabs(delta[1]) < 1.0e-6 &&
            std::fabs(delta[2]) < 1.0e-8)
          break;
      }
      else
      {
        lambda *= 10.0;
      }
    }

    ndof = static_cast<int>(local_x.size()) - 3;
    return true;
  }

  SagittaFit fitSagitta(const std::vector<FitPoint>& points)
  {
    SagittaFit fit;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> w;
    x.reserve(points.size());
    y.reserve(points.size());
    w.reserve(points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
    {
      x.push_back(points[i].x);
      y.push_back(points[i].y);
      w.push_back(points[i].w);
    }

    fit.ok = weightedSagittaFit(x, y, w,
                                fit.S, fit.x0, fit.invR, fit.theta, fit.b,
                                fit.chi2, fit.ndof);
    return fit;
  }
}
