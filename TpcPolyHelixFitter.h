#ifndef TPCPOLYHELIXFITTER_H
#define TPCPOLYHELIXFITTER_H

#include <array>
#include <vector>

class CDBInterface;
class PHField;

// -----------------------------------------------------------------------
// TpcPolyHelixFitter
//
// Two-stage charged-particle trajectory fit:
//   1) algebraicFit(): fast Taubin circle fit (xy) + line fit (rz) seed,
//      no field map needed.
//   2) fieldFit(): refines that seed against the real (non-uniform)
//      field map via RK4 propagation + Minuit2 least squares, fitting
//      (x0, y0, z0, phi, tanLambda, qOverPt) -- q/pT so the fit is
//      well-behaved across the full physical pT range, including
//      straight (high-pT) tracks.
//
// All reported track parameters (d0, z0, phi0, theta, curvature, pt, p)
// are evaluated at the point of closest approach (PCA) to the z-axis,
// not at the arbitrary hit used to seed the fit.
//
// Units: position in cm, momentum in GeV/c, field in Tesla.
// -----------------------------------------------------------------------
class TpcPolyHelixFitter
{
 public:
  struct Point
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };

  struct State
  {
    double x = 0.0, y = 0.0, z = 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    double charge = 1.0;  // +1 or -1 (units of e)

    double p() const;
  };

  struct FitResult
  {
    bool ok = false;

    // geometric (xy circle + rz line) parameters from algebraicFit()
    double radius = 0.0;
    double x0 = 0.0;  // circle center, xy
    double y0 = 0.0;
    double curvature = 0.0;  // signed 1/R, evaluated at the PCA if field fit ran
    double slope = 0.0;      // dz/dr

    // track parameters at the point of closest approach to the z-axis
    double x = 0.0;
    double y = 0.0;
    double z0 = 0.0;
    double phi0 = 0.0;
    double theta = 0.0;
    double d0 = 0.0;

    // momentum at the PCA (only meaningful once the field-map fit runs)
    double pt = 0.0;
    double p = 0.0;

    double chi2_xy = 0.0;
    double chi2_z = 0.0;
    int ndof_xy = 0;
    int ndof_z = 0;
  };

  TpcPolyHelixFitter();
  ~TpcPolyHelixFitter();

  // Loads the FIELDMAP_TRACKING CDB payload. Must succeed before fit()
  // will run the field-map refinement stage (it will still return the
  // algebraic-only result if this hasn't been called / failed).
  bool InitField(int verbosity = 0);

  bool fit(const std::vector<Point>& points, FitResult& fit_result) const;

  int Verbosity() const { return m_verbosity; }
  void Verbosity(int v) { m_verbosity = v; }

  void setStepSize(double ds_cm) { m_stepSize = ds_cm; }
  double stepSize() const { return m_stepSize; }

  void setMaxIterations(int n) { m_maxIterations = n; }

 private:
  // c_light in convenient units: GeV/(e*Tesla*cm)
  static constexpr double kConv = 0.0029979246;

  struct FieldFitResult
  {
    State state;
    double chi2 = 0.0;
    int ndf = 0;
    bool valid = false;
  };

  bool algebraicFit(const std::vector<Point>& points, FitResult& fit_result) const;

  FieldFitResult fieldFit(const std::vector<Point>& points, const FitResult& seed) const;
  State makeSeedState(const std::vector<Point>& points, const FitResult& seed) const;
  double chi2ForParams(const std::vector<Point>& points, const double* pars) const;

  std::vector<std::array<double, 3>> propagate(const State& state, double length_cm) const;

  // Propagates `state` backward to the point of closest approach to the
  // z-axis (min x^2+y^2), refined by a local quadratic fit around the
  // best sampled point. Returns false (with the input state unusable as
  // a PCA) if no clear minimum is found within maxLength_cm.
  bool propagateToPCA(const State& state, State& pca_state, double maxLength_cm = 60.0) const;

  void rk4Step(double pos[3], double dir[3], double ds, double p, double q) const;
  void getFieldTesla(const double point_cm[4], double field_t[3]) const;

  CDBInterface* m_cdb = nullptr;
  PHField* m_field = nullptr;

  double m_stepSize = 0.3;  // cm
  int m_maxIterations = 200;
  int m_verbosity = 0;
};

#endif
