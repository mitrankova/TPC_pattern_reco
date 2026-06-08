#pragma once

#include <vector>

// ===================================================================
// Small reusable fitting helper for in-module / full-track pieces.
//
// The fitter is coordinate-agnostic: x/y can be hardware coordinates
// (layer,pad), (layer,tbin), (radius,phi), or global/local coordinates.
// The caller only has to fill FitPoint(x, y, weight).
// ===================================================================
namespace Fitter
{
  enum FitMode
  {
    FIT_LINEAR  = 0,
    FIT_SAGITTA = 1
  };

  struct FitPoint
  {
    FitPoint();
    FitPoint(double x_, double y_, double w_ = 1.0);

    double x;
    double y;
    double w;
  };

  struct LineFit
  {
    LineFit();

    bool ok;
    double slope;
    double intercept;
    double chi2;
    int ndof;
  };

  struct SagittaFit
  {
    SagittaFit();

    bool ok;
    double S;
    double x0;
    double invR;
    double theta;
    double b;
    double chi2;
    int ndof;
  };

  double adcWeight(double adc, double maxadc, double power, double floor_frac);

  bool weightedLineFit(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<double>& w,
                       double& m,
                       double& b,
                       double& chi2,
                       int& ndof);

  LineFit fitLine(const std::vector<FitPoint>& points);

  double sagittaModel(double xrot, double S, double x0, double invR);

  bool weightedSagittaFit(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const std::vector<double>& w,
                          double& S,
                          double& x0,
                          double& invR,
                          double& theta,
                          double& bline,
                          double& chi2,
                          int& ndof);

  SagittaFit fitSagitta(const std::vector<FitPoint>& points);
}
