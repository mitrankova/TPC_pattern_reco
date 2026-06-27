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
  double get_collision_radius() const override { return get_collision_radius(0); }
  double get_collision_phi() const override { return get_collision_phi(0); }
  double get_collision_timebin() const override { return get_collision_timebin(0); }
  double get_collision_timebin_rms() const override { return get_collision_timebin_rms(0); }
  unsigned int get_collision_ntracks() const override { return get_collision_ntracks(0); }
  unsigned int get_collision_min_layers() const override { return m_collision_min_layers; }

  unsigned int get_collision_vertex_count() const override { return static_cast<unsigned int>(m_collision_radius.size()); }
  double get_collision_radius(unsigned int i) const override;
  double get_collision_phi(unsigned int i) const override;
  double get_collision_timebin(unsigned int i) const override;
  double get_collision_timebin_rms(unsigned int i) const override;
  unsigned int get_collision_ntracks(unsigned int i) const override;

  void set_collision_vertex_valid(int v) override { m_collision_vertex_valid = v; }
  void set_collision_radius(double v) override;
  void set_collision_phi(double v) override;
  void set_collision_timebin(double v) override;
  void set_collision_timebin_rms(double v) override;
  void set_collision_ntracks(unsigned int v) override;
  void set_collision_min_layers(unsigned int v) override { m_collision_min_layers = v; }
  void clear_collision_vertices() override;
  void add_collision_vertex(double radius, double phi, double timebin, double timebin_rms, unsigned int ntracks) override;

 private:
  void ensure_collision_vertex(unsigned int i);

  std::vector<FullTrackVertex*> m_vertices;
  int m_collision_vertex_valid {0};
  std::vector<double> m_collision_radius;
  std::vector<double> m_collision_phi;
  std::vector<double> m_collision_timebin;
  std::vector<double> m_collision_timebin_rms;
  std::vector<unsigned int> m_collision_ntracks;
  unsigned int m_collision_min_layers {0};

  ClassDefOverride(FullTrackVertexContainerv1, 4)
};
