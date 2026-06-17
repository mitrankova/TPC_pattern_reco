#pragma once

#include <phool/PHObject.h>

#include <iostream>

class FullTrackVertex;

class FullTrackVertexContainer : public PHObject
{
 public:
  FullTrackVertexContainer() = default;
  ~FullTrackVertexContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_vertex(FullTrackVertex*) {}
  virtual const FullTrackVertex* get_vertex(unsigned int) const { return nullptr; }
  virtual FullTrackVertex* get_vertex(unsigned int) { return nullptr; }

  virtual int get_collision_vertex_valid() const { return 0; }
  virtual double get_collision_radius() const { return 0.0; }
  virtual double get_collision_phi() const { return 0.0; }
  virtual double get_collision_timebin() const { return 0.0; }
  virtual double get_collision_timebin_rms() const { return 0.0; }
  virtual unsigned int get_collision_ntracks() const { return 0; }
  virtual unsigned int get_collision_min_layers() const { return 0; }

  virtual void set_collision_vertex_valid(int) {}
  virtual void set_collision_radius(double) {}
  virtual void set_collision_phi(double) {}
  virtual void set_collision_timebin(double) {}
  virtual void set_collision_timebin_rms(double) {}
  virtual void set_collision_ntracks(unsigned int) {}
  virtual void set_collision_min_layers(unsigned int) {}

 private:
  ClassDefOverride(FullTrackVertexContainer, 0)
};
