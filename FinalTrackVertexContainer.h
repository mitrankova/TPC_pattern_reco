#pragma once

#include <phool/PHObject.h>

#include <iostream>

class FinalTrackVertex;

class FinalTrackVertexContainer : public PHObject
{
 public:
  FinalTrackVertexContainer() = default;
  ~FinalTrackVertexContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_vertex(FinalTrackVertex*) {}
  virtual const FinalTrackVertex* get_vertex(unsigned int) const { return nullptr; }
  virtual FinalTrackVertex* get_vertex(unsigned int) { return nullptr; }

  virtual int get_collision_vertex_valid() const { return 0; }
  virtual unsigned int get_collision_vertex_count() const { return 0; }
  virtual double get_collision_x(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_y(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_z(unsigned int = 0) const { return 0.0; }
  virtual double get_collision_z_rms(unsigned int = 0) const { return 0.0; }
  virtual unsigned int get_collision_ntracks(unsigned int = 0) const { return 0; }
  virtual unsigned int get_collision_min_clusters() const { return 0; }

  virtual void set_collision_vertex_valid(int) {}
  virtual void set_collision_min_clusters(unsigned int) {}
  virtual void clear_collision_vertices() {}
  virtual void add_collision_vertex(double, double, double, double, unsigned int) {}

 private:
  ClassDefOverride(FinalTrackVertexContainer, 0)
};
