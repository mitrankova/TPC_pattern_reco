#pragma once

// Local TPC pad-coordinate helper.
//
// This helper intentionally does NOT provide a global sector phi conversion.
// It only reproduces the same local coordinates used in InModuleTracks.cc:
//   local_phi = pad * phi_bin_width[region]
//   local_radius = module_radius[region][layer_in_region]
//
// Use this for InModuleTracks and InModuleTrackDisplay local coordinates.
class TpcPadMap
{
 public:
  TpcPadMap() = default;
  ~TpcPadMap() = default;

  // Kept for compatibility with modules that already call load_from_cdb().
  // The local-coordinate constants below do not require CDB loading.
  int load_from_cdb(int verbosity = 0);
  bool is_loaded() const { return true; }

  // region: 0 inner, 1 middle, 2 outer.
  // pad may be fractional for fit/display points.
  double get_local_phi(unsigned int region, double pad) const;

  // layer is the real TPC layer number, 7..54.
  double get_local_radius(unsigned int region, unsigned int layer) const;

  // Convenience layer-only radius lookup, region inferred from layer.
  double get_radius(unsigned int layer) const;

  double get_phi_bin_width(unsigned int region) const;

 private:
  static const int N_MODULES = 3;
  static const int N_ROWS = 16;

  int layer_to_region(unsigned int layer) const;
  int layer_to_row_in_region(unsigned int region, unsigned int layer) const;
};
