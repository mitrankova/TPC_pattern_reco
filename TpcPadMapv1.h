#ifndef TPCPADMAPV1_H
#define TPCPADMAPV1_H

#include "TpcPadMap.h"

class TpcPadMapv1 : public TpcPadMap
{
 public:
  TpcPadMapv1();
  ~TpcPadMapv1() override {}

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;
  void Reset() override;

  bool is_loaded() const override { return m_is_loaded; }

  int load_from_cdb(int verbosity = 0) override;

  double get_local_phi(unsigned int region, double pad) const override;
  double get_local_phi(unsigned int side, unsigned int sector, unsigned int region, double pad) const override;
  double get_global_phi(unsigned int side, unsigned int region, double pad) const override;

  double get_sector_min_phi(unsigned int side, unsigned int sector, unsigned int region) const override;
  double get_sector_max_phi(unsigned int side, unsigned int sector, unsigned int region) const override;

  double get_local_radius(unsigned int region, unsigned int layer) const override;
  double get_radius(unsigned int layer) const override;

  double get_phi_bin_width(unsigned int region) const override;
  double get_phi_bin_width_for_layer(unsigned int layer) const override;
  unsigned int get_pads_per_sector(unsigned int region) const override;
  unsigned int get_total_pads(unsigned int region) const override;

  void print_sector_edges() const override;

 private:
  static const int N_SIDES = 2;
  static const int N_SECTORS = 12;
  static const int N_MODULES = 3;
  static const int N_ROWS = 16;
  static const int N_LAYERS = N_MODULES * N_ROWS;

  int layer_to_region(unsigned int layer) const;
  int layer_to_row_in_region(unsigned int region, unsigned int layer) const;

  void initialize_sector_phi_from_layer(unsigned int layer, double min_phi, double max_phi, double phi_bin_width);

  double wrap_phi(double phi) const;

  bool m_is_loaded;

  unsigned int m_nphi_bins[N_MODULES];
  unsigned int m_pads_per_sector[N_MODULES];

  double m_phi_bin_width[N_MODULES];
  double m_phi_bin_width_layer[N_LAYERS];
  double m_layer_radius[N_LAYERS];
  double m_sec_min_phi[N_SIDES][N_SECTORS][N_MODULES];
  double m_sec_max_phi[N_SIDES][N_SECTORS][N_MODULES];

  ClassDefOverride(TpcPadMapv1, 1)
};

#endif
