#pragma once

#include "InModuleTrackHit.h"

#include <phool/PHObject.h>

#include <vector>
#include <iostream>

class InModuleTrack : public PHObject
{
 public:
  InModuleTrack();
  virtual ~InModuleTrack() {}

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const;
  PHObject* CloneMe() const { return new InModuleTrack(*this); }

  unsigned int get_event() const { return m_event; }
  unsigned int get_track_id() const { return m_track_id; }

  unsigned int get_region() const { return m_region; }
  unsigned int get_sector() const { return m_sector; }
  int get_side() const { return m_side; }

  unsigned int get_nblobs() const { return m_nblobs; }
  unsigned int get_nrawhits() const { return m_nrawhits; }

  unsigned int get_first_layer() const { return m_first_layer; }
  unsigned int get_last_layer() const { return m_last_layer; }

  double get_pad_slope() const { return m_pad_slope; }
  double get_pad_intercept() const { return m_pad_intercept; }

  double get_tbin_slope() const { return m_tbin_slope; }
  double get_tbin_intercept() const { return m_tbin_intercept; }

  double get_chi2_pad() const { return m_chi2_pad; }
  double get_chi2_tbin() const { return m_chi2_tbin; }

  int get_ndof_pad() const { return m_ndof_pad; }
  int get_ndof_tbin() const { return m_ndof_tbin; }

  void set_event(unsigned int v) { m_event = v; }
  void set_track_id(unsigned int v) { m_track_id = v; }

  void set_region(unsigned int v) { m_region = v; }
  void set_sector(unsigned int v) { m_sector = v; }
  void set_side(int v) { m_side = v; }

  void set_nblobs(unsigned int v) { m_nblobs = v; }
  void set_nrawhits(unsigned int v) { m_nrawhits = v; }

  void set_first_layer(unsigned int v) { m_first_layer = v; }
  void set_last_layer(unsigned int v) { m_last_layer = v; }

  void set_pad_slope(double v) { m_pad_slope = v; }
  void set_pad_intercept(double v) { m_pad_intercept = v; }

  void set_tbin_slope(double v) { m_tbin_slope = v; }
  void set_tbin_intercept(double v) { m_tbin_intercept = v; }

  void set_chi2_pad(double v) { m_chi2_pad = v; }
  void set_chi2_tbin(double v) { m_chi2_tbin = v; }

  void set_ndof_pad(int v) { m_ndof_pad = v; }
  void set_ndof_tbin(int v) { m_ndof_tbin = v; }

  void add_hit(const InModuleTrackHit& hit)
  {
    m_hits.push_back(hit);
  }

  unsigned int size_hits() const
  {
    return static_cast<unsigned int>(m_hits.size());
  }

  const InModuleTrackHit* get_hit(unsigned int i) const
  {
    if (i >= m_hits.size()) return 0;
    return &m_hits[i];
  }

  const std::vector<InModuleTrackHit>& get_hits() const
  {
    return m_hits;
  }

 private:
  unsigned int m_event;
  unsigned int m_track_id;

  unsigned int m_region;
  unsigned int m_sector;
  int m_side;

  unsigned int m_nblobs;
  unsigned int m_nrawhits;

  unsigned int m_first_layer;
  unsigned int m_last_layer;

  double m_pad_slope;
  double m_pad_intercept;

  double m_tbin_slope;
  double m_tbin_intercept;

  double m_chi2_pad;
  double m_chi2_tbin;

  int m_ndof_pad;
  int m_ndof_tbin;

  std::vector<InModuleTrackHit> m_hits;

  ClassDef(InModuleTrack, 1)
};