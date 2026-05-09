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
