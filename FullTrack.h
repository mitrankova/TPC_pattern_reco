#pragma once

#include <trackbase/TrkrDefs.h>
#include <phool/PHObject.h>

#include <iostream>
#include <utility>
#include <vector>

// A full TPC track built by connecting already reconstructed InModuleTrack
// pieces across modules/sectors.  Concrete payload lives in versioned
// subclasses (FullTrackv1, ...), following the same schema-evolution pattern
// as InModuleTrack/InModuleTrackv1.
class FullTrack : public PHObject
{
 public:
  FullTrack() = default;
  ~FullTrack() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  // --- Identity ---
  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual int get_side() const { return 0; }

  // --- Topology ---
  virtual unsigned int get_nsegments() const { return 0; }
  virtual unsigned int get_nblobs() const { return 0; }
  virtual unsigned int get_nrawhits() const { return 0; }
  virtual unsigned int get_first_layer() const { return 0; }
  virtual unsigned int get_last_layer() const { return 0; }
  virtual unsigned int get_first_sector() const { return 0; }
  virtual unsigned int get_last_sector() const { return 0; }
  virtual unsigned int get_first_region() const { return 0; }
  virtual unsigned int get_last_region() const { return 0; }

  // --- Full-track fit in global sector-unwrapped coordinates ---
  // global_phi is in radians.  It is built from sector + local pad using the
  // same local pad/phi calibration constants as InModuleTracks.
  virtual double get_phi_slope() const { return 0.0; }
  virtual double get_phi_intercept() const { return 0.0; }
  virtual double get_tbin_slope() const { return 0.0; }
  virtual double get_tbin_intercept() const { return 0.0; }
  virtual double get_chi2_phi() const { return 0.0; }
  virtual double get_chi2_tbin() const { return 0.0; }
  virtual int get_ndof_phi() const { return 0; }
  virtual int get_ndof_tbin() const { return 0; }

  // --- Setters ---
  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_side(int) {}
  virtual void set_nsegments(unsigned int) {}
  virtual void set_nblobs(unsigned int) {}
  virtual void set_nrawhits(unsigned int) {}
  virtual void set_first_layer(unsigned int) {}
  virtual void set_last_layer(unsigned int) {}
  virtual void set_first_sector(unsigned int) {}
  virtual void set_last_sector(unsigned int) {}
  virtual void set_first_region(unsigned int) {}
  virtual void set_last_region(unsigned int) {}
  virtual void set_phi_slope(double) {}
  virtual void set_phi_intercept(double) {}
  virtual void set_tbin_slope(double) {}
  virtual void set_tbin_intercept(double) {}
  virtual void set_chi2_phi(double) {}
  virtual void set_chi2_tbin(double) {}
  virtual void set_ndof_phi(int) {}
  virtual void set_ndof_tbin(int) {}

  // Source InModuleTrack ids used by this full track.
  virtual void add_source_track(unsigned int /*track_id*/,
                                unsigned int /*region*/,
                                unsigned int /*sector*/) {}
  virtual unsigned int size_source_tracks() const { return 0; }
  virtual unsigned int get_source_track_id(unsigned int) const { return 0; }
  virtual unsigned int get_source_region(unsigned int) const { return 0; }
  virtual unsigned int get_source_sector(unsigned int) const { return 0; }

  // Hit index references copied from source in-module tracks.
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;
  virtual void add_hit_index(TrkrDefs::hitsetkey, TrkrDefs::hitkey) {}
  virtual unsigned int size_hit_indices() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual const std::vector<HitIndex>& get_hit_indices() const;

 private:
  static const std::vector<HitIndex> s_empty_indices;

  ClassDefOverride(FullTrack, 0)
};
