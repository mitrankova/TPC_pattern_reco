#include "TpcPolyClusterContainerv1.h"
#include "TpcPolyCluster.h"

ClassImp(TpcPolyClusterContainerv1)

TpcPolyClusterContainerv1::TpcPolyClusterContainerv1()
{
  Reset();
}

TpcPolyClusterContainerv1::~TpcPolyClusterContainerv1()
{
  Reset();
}

void TpcPolyClusterContainerv1::identify(std::ostream& os) const
{
  os << "TpcPolyClusterContainerv1 with " << m_clusters.size() << " clusters" << std::endl;
}

void TpcPolyClusterContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_clusters.size(); ++i) delete m_clusters[i];
  m_clusters.clear();
}

int TpcPolyClusterContainerv1::isValid() const
{
  return m_clusters.empty() ? 0 : 1;
}

PHObject* TpcPolyClusterContainerv1::CloneMe() const
{
  TpcPolyClusterContainerv1* copy = new TpcPolyClusterContainerv1();
  for (unsigned int i = 0; i < m_clusters.size(); ++i)
  {
    if (m_clusters[i]) copy->m_clusters.push_back(static_cast<TpcPolyCluster*>(m_clusters[i]->CloneMe()));
  }
  return copy;
}
