#pragma once

#include "TpcPolyTrack.h"

#include <iostream>
#include <vector>

class TpcPolyTrackv1 : public TpcPolyTrack
{
 public:
  TpcPolyTrackv1();
  ~TpcPolyTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcPolyTrackv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }
  unsigned int get_source_full_track_id() const override { return m_source_full_track_id; }
  int get_side() const override { return m_side; }
  unsigned int get_nhits() const override { return static_cast<unsigned int>(m_hit_indices.size()); }
  int get_fit_status() const override { return m_fit_status; }

  double get_d0() const override { return m_d0; }
  double get_z0() const override { return m_z0; }
  double get_phi0() const override { return m_phi0; }
  double get_theta() const override { return m_theta; }
  double get_curvature() const override { return m_curvature; }
  double get_chi2_xy() const override { return m_chi2_xy; }
  double get_chi2_z() const override { return m_chi2_z; }
  int get_ndof_xy() const override { return m_ndof_xy; }
  int get_ndof_z() const override { return m_ndof_z; }

  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_source_full_track_id(unsigned int v) override { m_source_full_track_id = v; }
  void set_side(int v) override { m_side = v; }
  void set_fit_status(int v) override { m_fit_status = v; }
  void set_d0(double v) override { m_d0 = v; }
  void set_z0(double v) override { m_z0 = v; }
  void set_phi0(double v) override { m_phi0 = v; }
  void set_theta(double v) override { m_theta = v; }
  void set_curvature(double v) override { m_curvature = v; }
  void set_chi2_xy(double v) override { m_chi2_xy = v; }
  void set_chi2_z(double v) override { m_chi2_z = v; }
  void set_ndof_xy(int v) override { m_ndof_xy = v; }
  void set_ndof_z(int v) override { m_ndof_z = v; }

  void add_hit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk,
               double x, double y, double z) override
  {
    m_hit_indices.emplace_back(hsk, hk);
    m_hit_x.push_back(x);
    m_hit_y.push_back(y);
    m_hit_z.push_back(z);
  }

  unsigned int size_hits() const override { return static_cast<unsigned int>(m_hit_indices.size()); }
  HitIndex get_hit_index(unsigned int i) const override
  {
    if (i >= m_hit_indices.size()) return {0, 0};
    return m_hit_indices[i];
  }
  double get_hit_x(unsigned int i) const override { return i < m_hit_x.size() ? m_hit_x[i] : 0.0; }
  double get_hit_y(unsigned int i) const override { return i < m_hit_y.size() ? m_hit_y[i] : 0.0; }
  double get_hit_z(unsigned int i) const override { return i < m_hit_z.size() ? m_hit_z[i] : 0.0; }
  const std::vector<HitIndex>& get_hit_indices() const override { return m_hit_indices; }

 private:
  unsigned int m_event {0};
  unsigned int m_track_id {0};
  unsigned int m_source_full_track_id {0};
  int m_side {0};
  int m_fit_status {0};

  double m_d0 {0.0};
  double m_z0 {0.0};
  double m_phi0 {0.0};
  double m_theta {0.0};
  double m_curvature {0.0};
  double m_chi2_xy {0.0};
  double m_chi2_z {0.0};
  int m_ndof_xy {0};
  int m_ndof_z {0};

  std::vector<HitIndex> m_hit_indices;
  std::vector<double> m_hit_x;
  std::vector<double> m_hit_y;
  std::vector<double> m_hit_z;

  ClassDefOverride(TpcPolyTrackv1, 2)
};
