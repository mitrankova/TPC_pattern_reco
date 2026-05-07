#pragma once

#include <trackbase/TrkrDefs.h>

#include <TObject.h>

class InModuleTrackHit : public TObject
{
 public:
  InModuleTrackHit();
  virtual ~InModuleTrackHit() {}

  void Reset();

  unsigned int get_event() const { return m_event; }
  unsigned int get_track_id() const { return m_track_id; }

  unsigned int get_region() const { return m_region; }
  unsigned int get_sector() const { return m_sector; }
  int get_side() const { return m_side; }

  unsigned int get_layer() const { return m_layer; }

  TrkrDefs::hitsetkey get_hitsetkey() const { return m_hitsetkey; }
  TrkrDefs::hitkey get_hitkey() const { return m_hitkey; }

  double get_pad() const { return m_pad; }
  double get_tbin() const { return m_tbin; }
  double get_adc() const { return m_adc; }

  double get_local_phi() const { return m_local_phi; }
  double get_local_radius() const { return m_local_radius; }

  void set_event(unsigned int v) { m_event = v; }
  void set_track_id(unsigned int v) { m_track_id = v; }

  void set_region(unsigned int v) { m_region = v; }
  void set_sector(unsigned int v) { m_sector = v; }
  void set_side(int v) { m_side = v; }

  void set_layer(unsigned int v) { m_layer = v; }

  void set_hitsetkey(TrkrDefs::hitsetkey v) { m_hitsetkey = v; }
  void set_hitkey(TrkrDefs::hitkey v) { m_hitkey = v; }

  void set_pad(double v) { m_pad = v; }
  void set_tbin(double v) { m_tbin = v; }
  void set_adc(double v) { m_adc = v; }

  void set_local_phi(double v) { m_local_phi = v; }
  void set_local_radius(double v) { m_local_radius = v; }

 private:
  unsigned int m_event;
  unsigned int m_track_id;

  unsigned int m_region;
  unsigned int m_sector;
  int m_side;

  unsigned int m_layer;

  TrkrDefs::hitsetkey m_hitsetkey;
  TrkrDefs::hitkey m_hitkey;

  double m_pad;
  double m_tbin;
  double m_adc;

  double m_local_phi;
  double m_local_radius;

  ClassDef(InModuleTrackHit, 1)
};