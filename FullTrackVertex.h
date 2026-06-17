#pragma once

#include <phool/PHObject.h>

#include <iostream>

class FullTrackVertex : public PHObject
{
 public:
  FullTrackVertex() = default;
  ~FullTrackVertex() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int get_track_id() const { return 0; }
  virtual double get_d0() const { return 0.0; }
  virtual double get_timebin0() const { return 0.0; }

  virtual void set_track_id(unsigned int) {}
  virtual void set_d0(double) {}
  virtual void set_timebin0(double) {}

 private:
  ClassDefOverride(FullTrackVertex, 0)
};
