#pragma once

#include "FullTrackVertex.h"

#include <iostream>

class FullTrackVertexv1 : public FullTrackVertex
{
 public:
  FullTrackVertexv1();
  ~FullTrackVertexv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new FullTrackVertexv1(*this); }

  unsigned int get_track_id() const override { return m_track_id; }
  double get_d0() const override { return m_d0; }
  double get_timebin0() const override { return m_timebin0; }
  int get_pca_valid() const override { return m_pca_valid; }
  double get_pca_radius() const override { return m_pca_radius; }
  double get_pca_phi() const override { return m_pca_phi; }
  double get_pca_timebin() const override { return m_pca_timebin; }

  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_d0(double v) override { m_d0 = v; }
  void set_timebin0(double v) override { m_timebin0 = v; }
  void set_pca_valid(int v) override { m_pca_valid = v; }
  void set_pca_radius(double v) override { m_pca_radius = v; }
  void set_pca_phi(double v) override { m_pca_phi = v; }
  void set_pca_timebin(double v) override { m_pca_timebin = v; }

 private:
  unsigned int m_track_id {0};
  double m_d0 {0.0};
  double m_timebin0 {0.0};
  int m_pca_valid {0};
  double m_pca_radius {0.0};
  double m_pca_phi {0.0};
  double m_pca_timebin {0.0};

  ClassDefOverride(FullTrackVertexv1, 2)
};
