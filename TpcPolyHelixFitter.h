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

  struct Seed
  {
    bool valid {false};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double px {0.0};
    double py {0.0};
    double pz {0.0};
  };

  struct FitResult
  {
    bool ok {false};
    double d0 {0.0};
    double z0 {0.0};
    double phi0 {0.0};
    double theta {0.0};
    double curvature {0.0};
    double xc {0.0};
    double yc {0.0};
    double radius {0.0};
    double dzds {0.0};
    double chi2_xy {0.0};
    double chi2_z {0.0};
    int ndof_xy {0};
    int ndof_z {0};
  };

  static bool fit(const std::vector<Point>& points, FitResult& fit);
  static bool fit(const std::vector<Point>& points, const Seed& seed, FitResult& fit);
};
