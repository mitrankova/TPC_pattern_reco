#pragma once

#include "TpcPolyClusterContainer.h"

#include <iostream>
#include <vector>

class TpcPolyCluster;

class TpcPolyClusterContainerv1 : public TpcPolyClusterContainer
{
 public:
  TpcPolyClusterContainerv1();
  ~TpcPolyClusterContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_clusters.size()); }
  void add_cluster(TpcPolyCluster* cluster) override { m_clusters.push_back(cluster); }
  const TpcPolyCluster* get_cluster(unsigned int i) const override
  {
    if (i >= m_clusters.size()) return nullptr;
    return m_clusters[i];
  }
  TpcPolyCluster* get_cluster(unsigned int i) override
  {
    if (i >= m_clusters.size()) return nullptr;
    return m_clusters[i];
  }

 private:
  std::vector<TpcPolyCluster*> m_clusters;

  ClassDefOverride(TpcPolyClusterContainerv1, 1)
};
