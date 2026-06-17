#pragma once

#include <vector>

class TpcPolyHelixFitter
{
 public:
  struct Point
  {
    double x {0.0};
    double y {0.0};
    double z {0.0};
  };

  struct FitResult
  {
    bool ok {false};
    double d0 {0.0};
    double z0 {0.0};
    double phi0 {0.0};
    double theta {0.0};
    double curvature {0.0};
    double chi2_xy {0.0};
    double chi2_z {0.0};
    int ndof_xy {0};
    int ndof_z {0};
  };

  static bool fit(const std::vector<Point>& points, FitResult& fit);
};
