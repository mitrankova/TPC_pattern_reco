#include "TpcPadMap.h"

#include <iostream>

namespace
{
  static const int N_MODULES = 3;
  static const int N_ROWS = 16;

  // Same values/convention as the current InModuleTracks.cc local coordinate helper.
  static const double PHI_BIN_WIDTH[N_MODULES] =
  {
    0.0053073,
    0.003959,
    0.00265145
  };

  static const double MODULE_RADIUS[N_MODULES][N_ROWS] =
  {
    {
      29.854978828112735, 31.869737083177956, 32.43665978627038,
      33.00171100689825,  33.56863172731403,  34.133682357783,
      34.70060474122243,  35.26565540941076,  35.83257683544541,
      36.39762877363545,  36.964549975549694, 37.52960055896088,
      38.09652180558749,  38.66157293473739,  39.228495272708216,
      39.793545257944906
    },
    {
      41.65920253621078,  42.67990048015332,  43.7005755287188,
      44.7212729094545,   45.7419615067264,   46.76264656230158,
      47.78333428983602,  48.80401878201343,  49.82471910526506,
      50.8454060012135,   51.866093793785126, 52.88677964073831,
      53.90746625152035,  54.92815969895385,  55.948864895868056,
      56.9695394315422
    },
    {
      58.910963349324035, 60.00800996331871,  61.10505851260341,
      62.202104676954924, 63.29915863086735,  64.39619682986867,
      65.49324606923312,  66.59029899562653,  67.68734047670296,
      68.78439383353172,  69.88143340055497,  70.97848786511186,
      72.07553264226554,  73.17257662017182,  74.2696338511705,
      75.36667517343196
    }
  };
}

int TpcPadMap::load_from_cdb(const int verbosity)
{
  if (verbosity > 0)
  {
    std::cout << "TpcPadMap::load_from_cdb - local-only map: no CDB load needed"
              << std::endl;
  }
  return 0;
}

int TpcPadMap::layer_to_region(const unsigned int layer) const
{
  if (layer < 7 || layer > 54) return -1;
  const int region = (static_cast<int>(layer) - 7) / N_ROWS;
  if (region < 0 || region >= N_MODULES) return -1;
  return region;
}

int TpcPadMap::layer_to_row_in_region(const unsigned int region,
                                      const unsigned int layer) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return -1;
  const int row = static_cast<int>(layer) - 7 - static_cast<int>(region) * N_ROWS;
  if (row < 0 || row >= N_ROWS) return -1;
  return row;
}

double TpcPadMap::get_phi_bin_width(const unsigned int region) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;
  return PHI_BIN_WIDTH[region];
}

double TpcPadMap::get_local_phi(const unsigned int region,
                                const double pad) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;
  return pad * PHI_BIN_WIDTH[region];
}

double TpcPadMap::get_local_radius(const unsigned int region,
                                   const unsigned int layer) const
{
  const int row = layer_to_row_in_region(region, layer);
  if (row < 0) return 0.0;
  return MODULE_RADIUS[region][row];
}

double TpcPadMap::get_radius(const unsigned int layer) const
{
  const int region = layer_to_region(layer);
  if (region < 0) return 0.0;
  return get_local_radius(static_cast<unsigned int>(region), layer);
}
