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

  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_d0(double v) override { m_d0 = v; }
  void set_timebin0(double v) override { m_timebin0 = v; }

 private:
  unsigned int m_track_id {0};
  double m_d0 {0.0};
  double m_timebin0 {0.0};

  ClassDefOverride(FullTrackVertexv1, 1)
};
