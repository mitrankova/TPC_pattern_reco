#include "FullTrackVertexContainerv1.h"
#include "FullTrackVertex.h"

ClassImp(FullTrackVertexContainerv1)

FullTrackVertexContainerv1::FullTrackVertexContainerv1()
{
  Reset();
}

FullTrackVertexContainerv1::~FullTrackVertexContainerv1()
{
  Reset();
}

void FullTrackVertexContainerv1::identify(std::ostream& os) const
{
  os << "FullTrackVertexContainerv1 with "
     << m_vertices.size()
     << " full-track vertices"
     << " collision_valid=" << m_collision_vertex_valid
     << " collision_radius=" << m_collision_radius
     << " collision_phi=" << m_collision_phi
     << " collision_timebin=" << m_collision_timebin
     << " collision_timebin_rms=" << m_collision_timebin_rms
     << " collision_ntracks=" << m_collision_ntracks
     << " collision_min_layers=" << m_collision_min_layers
     << std::endl;
}

void FullTrackVertexContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_vertices.size(); ++i)
  {
    delete m_vertices[i];
  }
  m_vertices.clear();
  m_collision_vertex_valid = 0;
  m_collision_radius = 0.0;
  m_collision_phi = 0.0;
  m_collision_timebin = 0.0;
  m_collision_timebin_rms = 0.0;
  m_collision_ntracks = 0;
  m_collision_min_layers = 0;
}

int FullTrackVertexContainerv1::isValid() const
{
  return m_vertices.empty() ? 0 : 1;
}

PHObject* FullTrackVertexContainerv1::CloneMe() const
{
  FullTrackVertexContainerv1* copy = new FullTrackVertexContainerv1();
  copy->m_collision_vertex_valid = m_collision_vertex_valid;
  copy->m_collision_radius = m_collision_radius;
  copy->m_collision_phi = m_collision_phi;
  copy->m_collision_timebin = m_collision_timebin;
  copy->m_collision_timebin_rms = m_collision_timebin_rms;
  copy->m_collision_ntracks = m_collision_ntracks;
  copy->m_collision_min_layers = m_collision_min_layers;
  for (unsigned int i = 0; i < m_vertices.size(); ++i)
  {
    if (m_vertices[i]) copy->m_vertices.push_back(static_cast<FullTrackVertex*>(m_vertices[i]->CloneMe()));
  }
  return copy;
}
