#ifndef TPCPADMAP_H
#define TPCPADMAP_H

#include <phool/PHObject.h>

#include <iostream>

class TpcPadMap : public PHObject
{
 public:
  TpcPadMap() {}
  ~TpcPadMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "TpcPadMap base class" << std::endl;
  }

  int isValid() const override { return 0; }
  void Reset() override {}

  virtual bool is_loaded() const { return false; }

  virtual int load_from_cdb(int /*verbosity*/ = 0) { return -1; }

  virtual double get_local_phi(unsigned int /*region*/, double /*pad*/) const { return 0.0; }

  virtual double get_local_phi(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*region*/, double /*pad*/) const { return 0.0; }

  virtual double get_global_phi(unsigned int /*side*/, unsigned int /*region*/, double /*pad*/) const { return 0.0; }

  virtual double get_sector_min_phi(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*region*/) const { return 0.0; }

  virtual double get_sector_max_phi(unsigned int /*side*/, unsigned int /*sector*/, unsigned int /*region*/) const { return 0.0; }

  virtual double get_local_radius(unsigned int /*region*/, unsigned int /*layer*/) const { return 0.0; }

  virtual double get_radius(unsigned int /*layer*/) const { return 0.0; }

  virtual double get_phi_bin_width(unsigned int /*region*/) const { return 0.0; }

  virtual double get_phi_bin_width_for_layer(unsigned int /*layer*/) const { return 0.0; }

  virtual unsigned int get_pads_per_sector(unsigned int /*region*/) const { return 0; }

  virtual unsigned int get_total_pads(unsigned int /*region*/) const { return 0; }

  virtual void print_sector_edges() const {}

 private:
  ClassDefOverride(TpcPadMap, 1)
};

#endif
