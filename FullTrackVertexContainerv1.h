#pragma once

#include "FullTrackVertexContainer.h"

#include <iostream>
#include <vector>

class FullTrackVertex;

class FullTrackVertexContainerv1 : public FullTrackVertexContainer
{
 public:
  FullTrackVertexContainerv1();
  ~FullTrackVertexContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_vertices.size()); }
  void add_vertex(FullTrackVertex* vtx) override { m_vertices.push_back(vtx); }

  const FullTrackVertex* get_vertex(unsigned int i) const override
  {
    if (i >= m_vertices.size()) return nullptr;
    return m_vertices[i];
  }

  FullTrackVertex* get_vertex(unsigned int i) override
  {
    if (i >= m_vertices.size()) return nullptr;
    return m_vertices[i];
  }

  int get_collision_vertex_valid() const override { return m_collision_vertex_valid; }
  double get_collision_radius() const override { return m_collision_radius; }
  double get_collision_phi() const override { return m_collision_phi; }
  double get_collision_timebin() const override { return m_collision_timebin; }
  double get_collision_timebin_rms() const override { return m_collision_timebin_rms; }
  unsigned int get_collision_ntracks() const override { return m_collision_ntracks; }
  unsigned int get_collision_min_layers() const override { return m_collision_min_layers; }

  void set_collision_vertex_valid(int v) override { m_collision_vertex_valid = v; }
  void set_collision_radius(double v) override { m_collision_radius = v; }
  void set_collision_phi(double v) override { m_collision_phi = v; }
  void set_collision_timebin(double v) override { m_collision_timebin = v; }
  void set_collision_timebin_rms(double v) override { m_collision_timebin_rms = v; }
  void set_collision_ntracks(unsigned int v) override { m_collision_ntracks = v; }
  void set_collision_min_layers(unsigned int v) override { m_collision_min_layers = v; }

 private:
  std::vector<FullTrackVertex*> m_vertices;
  int m_collision_vertex_valid {0};
  double m_collision_radius {0.0};
  double m_collision_phi {0.0};
  double m_collision_timebin {0.0};
  double m_collision_timebin_rms {0.0};
  unsigned int m_collision_ntracks {0};
  unsigned int m_collision_min_layers {0};

  ClassDefOverride(FullTrackVertexContainerv1, 3)
};
