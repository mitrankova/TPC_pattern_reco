#include "TpcPadMapv1.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TObject.h>

ClassImp(TpcPadMapv1)

namespace
{
  static const int N_FEE = 26;
  static const int N_CH = 256;
}

TpcPadMapv1::TpcPadMapv1()
  : m_is_loaded(false)
{
  Reset();
}

void TpcPadMapv1::identify(std::ostream& os) const
{
  os << "TpcPadMapv1, loaded = " << m_is_loaded << std::endl;
}

int TpcPadMapv1::isValid() const
{
  return m_is_loaded ? 1 : 0;
}

void TpcPadMapv1::Reset()
{
  m_is_loaded = false;

  // Same total phibins convention as PHG4TpcDetector parameters:
  // inner/mid/outer = pads_per_sector*12 = 94*12, 128*12, 192*12.
  m_nphi_bins[0] = 94U * N_SECTORS;
  m_nphi_bins[1] = 128U * N_SECTORS;
  m_nphi_bins[2] = 192U * N_SECTORS;

  for (int region = 0; region < N_MODULES; ++region)
  {
    m_pads_per_sector[region] = m_nphi_bins[region] / N_SECTORS;
    m_phi_bin_width[region] = 0.0;
  }

  for (int layer = 0; layer < N_LAYERS; ++layer)
  {
    m_phi_bin_width_layer[layer] = 0.0;
    m_layer_radius[layer] = 0.0;
  }

  for (int side = 0; side < N_SIDES; ++side)
  {
    for (int sector = 0; sector < N_SECTORS; ++sector)
    {
      for (int region = 0; region < N_MODULES; ++region)
      {
        m_sec_min_phi[side][sector][region] = 0.0;
        m_sec_max_phi[side][sector][region] = 0.0;
      }
    }
  }
}

int TpcPadMapv1::load_from_cdb(const int verbosity)
{
  Reset();

  CDBInterface* cdb = CDBInterface::instance();
  const std::string calibdir = cdb->getUrl("TPC_FEE_CHANNEL_MAP");

  if (calibdir.empty())
  {
    std::cout << "TpcPadMapv1::load_from_cdb - no TPC_FEE_CHANNEL_MAP found" << std::endl;
    return -1;
  }

  CDBTTree* cdbttree = new CDBTTree(calibdir);
  cdbttree->LoadCalibrations();

  std::vector<double> pad_phi[N_LAYERS];
  std::vector<double> pad_R[N_LAYERS];

  for (int fee = 0; fee < N_FEE; ++fee)
  {
    for (int ch = 0; ch < N_CH; ++ch)
    {
      const unsigned int key = 256U * static_cast<unsigned int>(fee) + static_cast<unsigned int>(ch);

      const int layer = cdbttree->GetIntValue(key, "layer");
      if (layer > 6 && layer < 55)
      {
        const int v_layer = layer - 7;
        pad_phi[v_layer].push_back(cdbttree->GetDoubleValue(key, "phi"));
        pad_R[v_layer].push_back(cdbttree->GetDoubleValue(key, "R"));
      }
    }
  }

  delete cdbttree;
  cdbttree = 0;

  for (int v_layer = 0; v_layer < N_LAYERS; ++v_layer)
  {
    if (pad_phi[v_layer].empty() || pad_R[v_layer].empty())
    {
      std::cout << "TpcPadMapv1::load_from_cdb - missing pad map entries for TPC layer "
                << v_layer + 7 << std::endl;
      return -1;
    }

    double radius = 0.0;
    for (std::size_t ipad = 0; ipad < pad_R[v_layer].size(); ++ipad)
    {
      radius += pad_R[v_layer][ipad];
    }
    radius /= static_cast<double>(pad_R[v_layer].size());
    radius /= 10.0;  // CDB map R is in mm; local display uses cm.
    m_layer_radius[v_layer] = radius;

    const std::vector<double>::const_iterator min_phi_iter =
        std::min_element(pad_phi[v_layer].begin(), pad_phi[v_layer].end());
    const std::vector<double>::const_iterator max_phi_iter =
        std::max_element(pad_phi[v_layer].begin(), pad_phi[v_layer].end());

    double min_phi = *min_phi_iter;
    double max_phi = *max_phi_iter;

    // Match PHG4TpcDetector convention.
    min_phi -= M_PI / 2.0;
    max_phi -= M_PI / 2.0;

    const unsigned int region = static_cast<unsigned int>(v_layer / N_ROWS);
    const unsigned int pads_per_sector = m_pads_per_sector[region];
    if (pads_per_sector < 2U)
    {
      std::cout << "TpcPadMapv1::load_from_cdb - invalid pads_per_sector for region "
                << region << std::endl;
      return -1;
    }

    const double phi_bin_width = std::fabs(max_phi - min_phi) /
                                 static_cast<double>(pads_per_sector - 1U);
    m_phi_bin_width_layer[v_layer] = phi_bin_width;

    // PHG4TpcDetector fills the region sector edge arrays inside a loop over
    // layers. This intentionally overwrites the region-level values as layers
    // progress; after the loop, each region uses the last layer in that region.
    m_phi_bin_width[region] = phi_bin_width;
    initialize_sector_phi_from_layer(static_cast<unsigned int>(v_layer),
                                     min_phi,
                                     max_phi,
                                     phi_bin_width);
  }

  m_is_loaded = true;

  if (verbosity > 0)
  {
    std::cout << "TpcPadMapv1::load_from_cdb - loaded TPC_FEE_CHANNEL_MAP from "
              << calibdir << std::endl;
    std::cout << "  region 0 pads/sector " << m_pads_per_sector[0]
              << " phi_bin_width " << m_phi_bin_width[0] << std::endl;
    std::cout << "  region 1 pads/sector " << m_pads_per_sector[1]
              << " phi_bin_width " << m_phi_bin_width[1] << std::endl;
    std::cout << "  region 2 pads/sector " << m_pads_per_sector[2]
              << " phi_bin_width " << m_phi_bin_width[2] << std::endl;
  }

  return 0;
}

void TpcPadMapv1::initialize_sector_phi_from_layer(const unsigned int layer,
                                                 const double min_phi,
                                                 const double max_phi,
                                                 const double phi_bin_width)
{
  if (layer >= static_cast<unsigned int>(N_LAYERS)) return;
  const unsigned int region = layer / N_ROWS;
  if (region >= static_cast<unsigned int>(N_MODULES)) return;

  const double sector_phi = std::fabs(max_phi - min_phi) + phi_bin_width;

  for (int side = 0; side < N_SIDES; ++side)
  {
    for (int sector = 0; sector < N_SECTORS; ++sector)
    {
      const double sector_anchor = M_PI - 2.0 * M_PI / 12.0 * static_cast<double>(sector + 1);

      if (side == 0)
      {
        m_sec_min_phi[side][sector][region] = sector_anchor + (-max_phi - phi_bin_width / 2.0);
        m_sec_max_phi[side][sector][region] = m_sec_min_phi[side][sector][region] + sector_phi;
      }
      else
      {
        m_sec_max_phi[side][sector][region] = sector_anchor + (max_phi + phi_bin_width / 2.0);
        m_sec_min_phi[side][sector][region] = m_sec_max_phi[side][sector][region] - sector_phi;
      }
    }
  }
}

int TpcPadMapv1::layer_to_region(const unsigned int layer) const
{
  if (layer < 7U || layer > 54U) return -1;
  const int region = (static_cast<int>(layer) - 7) / N_ROWS;
  if (region < 0 || region >= N_MODULES) return -1;
  return region;
}

int TpcPadMapv1::layer_to_row_in_region(const unsigned int region,
                                      const unsigned int layer) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return -1;
  const int row = static_cast<int>(layer) - 7 - static_cast<int>(region) * N_ROWS;
  if (row < 0 || row >= N_ROWS) return -1;
  return row;
}

double TpcPadMapv1::wrap_phi(const double phi) const
{
  double out = phi;
  while (out <= -M_PI) out += 2.0 * M_PI;
  while (out > M_PI) out -= 2.0 * M_PI;
  return out;
}

double TpcPadMapv1::get_phi_bin_width(const unsigned int region) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;
  return m_phi_bin_width[region];
}

double TpcPadMapv1::get_phi_bin_width_for_layer(const unsigned int layer) const
{
  if (layer < 7U || layer > 54U) return 0.0;
  return m_phi_bin_width_layer[layer - 7U];
}

unsigned int TpcPadMapv1::get_pads_per_sector(const unsigned int region) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0U;
  return m_pads_per_sector[region];
}

unsigned int TpcPadMapv1::get_total_pads(const unsigned int region) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0U;
  return m_pads_per_sector[region] * N_SECTORS;
}

double TpcPadMapv1::get_sector_min_phi(const unsigned int side,
                                     const unsigned int sector,
                                     const unsigned int region) const
{
  if (side >= static_cast<unsigned int>(N_SIDES)) return 0.0;
  if (sector >= static_cast<unsigned int>(N_SECTORS)) return 0.0;
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;
  return m_sec_min_phi[side][sector][region];
}

double TpcPadMapv1::get_sector_max_phi(const unsigned int side,
                                     const unsigned int sector,
                                     const unsigned int region) const
{
  if (side >= static_cast<unsigned int>(N_SIDES)) return 0.0;
  if (sector >= static_cast<unsigned int>(N_SECTORS)) return 0.0;
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;
  return m_sec_max_phi[side][sector][region];
}

void TpcPadMapv1::print_sector_edges() const
{
  std::cout << "TPC sector phi edges from TpcPadMapv1" << std::endl;
  std::cout << "side sector region min_phi max_phi width phi_bin_width pads_per_sector" << std::endl;

  for (unsigned int side = 0; side < static_cast<unsigned int>(N_SIDES); ++side)
  {
    for (unsigned int sector = 0; sector < static_cast<unsigned int>(N_SECTORS); ++sector)
    {
      for (unsigned int region = 0; region < static_cast<unsigned int>(N_MODULES); ++region)
      {
        const double min_phi = m_sec_min_phi[side][sector][region];
        const double max_phi = m_sec_max_phi[side][sector][region];
        const double width = max_phi - min_phi;

        std::cout << "side " << side
                  << " sector " << sector
                  << " region " << region
                  << " min_phi " << min_phi
                  << " max_phi " << max_phi
                  << " width " << width
                  << " phi_bin_width " << m_phi_bin_width[region]
                  << " pads_per_sector " << m_pads_per_sector[region]
                  << std::endl;
      }
    }
  }
}

double TpcPadMapv1::get_local_phi(const unsigned int region,
                                const double pad) const
{
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;

  const unsigned int pads_per_sector = m_pads_per_sector[region];
  if (pads_per_sector == 0U) return 0.0;

  // pad is full-circle. Fold it to the coordinate inside one sector/module.
  const double local_pad = pad - std::floor(pad / static_cast<double>(pads_per_sector))
                               * static_cast<double>(pads_per_sector);
  return local_pad * m_phi_bin_width[region];
}

double TpcPadMapv1::get_local_phi(const unsigned int side,
                                const unsigned int sector,
                                const unsigned int region,
                                const double pad) const
{
  if (side >= static_cast<unsigned int>(N_SIDES)) return 0.0;
  if (sector >= static_cast<unsigned int>(N_SECTORS)) return 0.0;
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;

  const unsigned int pads_per_sector = m_pads_per_sector[region];
  if (pads_per_sector == 0U) return 0.0;

  // The caller passes the known sector and a full-circle pad number.
  // This keeps fractional pads from fitted lines unchanged.
  double local_pad = pad - static_cast<double>(sector * pads_per_sector);

  // Protect against callers that pass a full-circle pad but with sector not
  // synchronized to pad/pads_per_sector, and also against small extrapolations.
  if (local_pad < -1.0 || local_pad > static_cast<double>(pads_per_sector) + 1.0)
  {
    local_pad = pad - std::floor(pad / static_cast<double>(pads_per_sector))
                    * static_cast<double>(pads_per_sector);
  }

  return local_pad * m_phi_bin_width[region];
}

double TpcPadMapv1::get_global_phi(const unsigned int side,
                                 const unsigned int region,
                                 const double pad) const
{
  if (side >= static_cast<unsigned int>(N_SIDES)) return 0.0;
  if (region >= static_cast<unsigned int>(N_MODULES)) return 0.0;

  const unsigned int pads_per_sector = m_pads_per_sector[region];
  const unsigned int total_pads = pads_per_sector * N_SECTORS;
  if (pads_per_sector == 0U || pad < 0.0 || pad > static_cast<double>(total_pads)) return 0.0;

  unsigned int sector = static_cast<unsigned int>(std::floor(pad / static_cast<double>(pads_per_sector)));
  if (sector >= static_cast<unsigned int>(N_SECTORS)) sector = N_SECTORS - 1U;

  const double local_bin = pad + 0.5 - static_cast<double>(sector * pads_per_sector);
  const double phi = m_sec_max_phi[side][sector][region] - local_bin * m_phi_bin_width[region];

  return wrap_phi(phi);
}

double TpcPadMapv1::get_local_radius(const unsigned int region,
                                   const unsigned int layer) const
{
  const int row = layer_to_row_in_region(region, layer);
  if (row < 0) return 0.0;
  return m_layer_radius[static_cast<int>(region) * N_ROWS + row];
}
/*
double TpcPadMapv1::get_npads(const unsigned int region) const
{
  if (region < 0 || region > 2) return 0.0;
  return m_pads_per_sector[region];
}
*/

double TpcPadMapv1::get_radius(const unsigned int layer) const
{
  const int region = layer_to_region(layer);
  if (region < 0) return 0.0;
  return get_local_radius(static_cast<unsigned int>(region), layer);
}
