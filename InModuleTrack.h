#pragma once

#include <trackbase/TrkrDefs.h>

#include <phool/PHObject.h>

#include <iostream>
#include <utility>
#include <vector>

// A single in-module track.
// Concrete data are stored in versioned subclasses (InModuleTrackv1, ...).
// Instead of owning InModuleTrackHit objects, tracks now store
// (hitsetkey, hitkey) index pairs that reference hits in TrkrHitSetContainer.
class InModuleTrack : public PHObject
{
 public:
  InModuleTrack() = default;
  ~InModuleTrack() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  // --- Identity ---
  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }

  virtual unsigned int get_region() const { return 0; }
  virtual unsigned int get_sector() const { return 0; }
  virtual int get_side() const { return 0; }

  // --- Topology ---
  virtual unsigned int get_nblobs() const { return 0; }
  virtual unsigned int get_nrawhits() const { return 0; }

  virtual unsigned int get_first_layer() const { return 0; }
  virtual unsigned int get_last_layer() const { return 0; }

  // --- Fit parameters ---
  virtual double get_pad_slope() const { return 0; }
  virtual double get_pad_intercept() const { return 0; }

  virtual double get_tbin_slope() const { return 0; }
  virtual double get_tbin_intercept() const { return 0; }

  virtual double get_chi2_pad() const { return 0; }
  virtual double get_chi2_tbin() const { return 0; }

  virtual int get_ndof_pad() const { return 0; }
  virtual int get_ndof_tbin() const { return 0; }

  // --- Setters ---
  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}

  virtual void set_region(unsigned int) {}
  virtual void set_sector(unsigned int) {}
  virtual void set_side(int) {}

  virtual void set_nblobs(unsigned int) {}
  virtual void set_nrawhits(unsigned int) {}

  virtual void set_first_layer(unsigned int) {}
  virtual void set_last_layer(unsigned int) {}

  virtual void set_pad_slope(double) {}
  virtual void set_pad_intercept(double) {}

  virtual void set_tbin_slope(double) {}
  virtual void set_tbin_intercept(double) {}

  virtual void set_chi2_pad(double) {}
  virtual void set_chi2_tbin(double) {}

  virtual void set_ndof_pad(int) {}
  virtual void set_ndof_tbin(int) {}

  virtual void set_radius_tbin_slope(double) {}
  virtual void set_radius_tbin_intercept(double) {}
  virtual void set_radius_phi_slope(double) {}
  virtual void set_radius_phi_intercept(double) {}

  virtual void set_circle_x0(double) {}
  virtual void set_circle_y0(double) {}
  virtual void set_circle_radius(double) {}

  // Sagitta fit in the local transverse plane.  The fit follows the
  // Python notebook convention:
  //   1) fit x/y line and rotate so that line is horizontal,
  //   2) fit y_rot = S - 1/2 invR (x_rot-x0)^2
  //                    - 1/8 invR^3 (x_rot-x0)^4
  //                    - 1/16 invR^5 (x_rot-x0)^6.
  virtual void set_has_sagitta_fit(int) {}
  virtual void set_sagitta_S(double) {}
  virtual void set_sagitta_x0(double) {}
  virtual void set_sagitta_invR(double) {}
  virtual void set_sagitta_theta(double) {}
  virtual void set_sagitta_b(double) {}

  virtual double get_radius_tbin_slope() const { return 0.0; }
  virtual double get_radius_tbin_intercept() const { return 0.0; }
  virtual double get_radius_phi_slope() const { return 0.0; }
  virtual double get_radius_phi_intercept() const { return 0.0; }

  virtual double get_circle_x0() const { return 0.0; }
  virtual double get_circle_y0() const { return 0.0; }
  virtual double get_circle_radius() const { return 0.0; }

  virtual int get_has_sagitta_fit() const { return 0; }
  virtual double get_sagitta_S() const { return 0.0; }
  virtual double get_sagitta_x0() const { return 0.0; }
  virtual double get_sagitta_invR() const { return 0.0; }
  virtual double get_sagitta_theta() const { return 0.0; }
  virtual double get_sagitta_b() const { return 0.0; }

  virtual void set_local_fit_mode(int) {}
  virtual int get_local_fit_mode() const { return -1; }

  virtual void set_chi2_radius_tbin(double) {}
  virtual void set_ndof_radius_tbin(int) {}
  virtual double get_chi2_radius_tbin() const { return 0.0; }
  virtual int get_ndof_radius_tbin() const { return 0; }

  virtual void set_chi2_radius_phi_line(double) {}
  virtual void set_ndof_radius_phi_line(int) {}
  virtual double get_chi2_radius_phi_line() const { return 0.0; }
  virtual int get_ndof_radius_phi_line() const { return 0; }

  virtual void set_chi2_radius_phi_circle(double) {}
  virtual void set_ndof_radius_phi_circle(int) {}
  virtual double get_chi2_radius_phi_circle() const { return 0.0; }
  virtual int get_ndof_radius_phi_circle() const { return 0; }

  virtual void set_chi2_radius_phi_sagitta(double) {}
  virtual void set_ndof_radius_phi_sagitta(int) {}
  virtual double get_chi2_radius_phi_sagitta() const { return 0.0; }
  virtual int get_ndof_radius_phi_sagitta() const { return 0; }

  // --- Hit index access ---
  // Hits are referenced as (hitsetkey, hitkey) pairs pointing into
  // TrkrHitSetContainer; no hit data are owned by this object.
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  virtual void add_hit_index(TrkrDefs::hitsetkey /*hsk*/, TrkrDefs::hitkey /*hk*/) {}

  virtual unsigned int size_hit_indices() const { return 0; }

  virtual HitIndex get_hit_index(unsigned int /*i*/) const
  {
    return {0, 0};
  }

  virtual const std::vector<HitIndex>& get_hit_indices() const;

 private:
  // Returned by the default get_hit_indices() implementation.
  static const std::vector<HitIndex> s_empty_indices;

  ClassDefOverride(InModuleTrack, 0)
};
