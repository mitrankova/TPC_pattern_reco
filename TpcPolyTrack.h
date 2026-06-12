#pragma once

#include <trackbase/TrkrDefs.h>
#include <phool/PHObject.h>

#include <iostream>
#include <utility>
#include <vector>

class TpcPolyTrack : public PHObject
{
 public:
  TpcPolyTrack() = default;
  ~TpcPolyTrack() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_source_full_track_id() const { return 0; }
  virtual int get_side() const { return 0; }
  virtual unsigned int get_nhits() const { return 0; }
  virtual int get_fit_status() const { return 0; }

  virtual double get_d0() const { return 0.0; }
  virtual double get_z0() const { return 0.0; }
  virtual double get_phi0() const { return 0.0; }
  virtual double get_theta() const { return 0.0; }
  virtual double get_curvature() const { return 0.0; }
  virtual double get_chi2_xy() const { return 0.0; }
  virtual double get_chi2_z() const { return 0.0; }
  virtual int get_ndof_xy() const { return 0; }
  virtual int get_ndof_z() const { return 0; }

  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_source_full_track_id(unsigned int) {}
  virtual void set_side(int) {}
  virtual void set_fit_status(int) {}
  virtual void set_d0(double) {}
  virtual void set_z0(double) {}
  virtual void set_phi0(double) {}
  virtual void set_theta(double) {}
  virtual void set_curvature(double) {}
  virtual void set_chi2_xy(double) {}
  virtual void set_chi2_z(double) {}
  virtual void set_ndof_xy(int) {}
  virtual void set_ndof_z(int) {}

  virtual void add_hit(TrkrDefs::hitsetkey, TrkrDefs::hitkey,
                       double /*x*/, double /*y*/, double /*z*/) {}
  virtual unsigned int size_hits() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual double get_hit_x(unsigned int) const { return 0.0; }
  virtual double get_hit_y(unsigned int) const { return 0.0; }
  virtual double get_hit_z(unsigned int) const { return 0.0; }
  virtual const std::vector<HitIndex>& get_hit_indices() const;

 private:
  static const std::vector<HitIndex> s_empty_indices;

  ClassDefOverride(TpcPolyTrack, 0)
};
