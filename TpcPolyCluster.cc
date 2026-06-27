#include "TpcPolyCluster.h"

const std::vector<TpcPolyCluster::HitIndex> TpcPolyCluster::s_empty_indices;

void TpcPolyCluster::identify(std::ostream& os) const
{
  os << "TpcPolyCluster base class" << std::endl;
}

void TpcPolyCluster::Reset()
{
}

int TpcPolyCluster::isValid() const
{
  return 0;
}

const std::vector<TpcPolyCluster::HitIndex>& TpcPolyCluster::get_hit_indices() const
{
  return s_empty_indices;
}
