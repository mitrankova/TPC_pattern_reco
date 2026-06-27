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
     << " collision_vertices=" << m_collision_radius.size()
     << " first_collision_radius=" << get_collision_radius()
     << " first_collision_phi=" << get_collision_phi()
     << " first_collision_timebin=" << get_collision_timebin()
     << " first_collision_timebin_rms=" << get_collision_timebin_rms()
     << " first_collision_ntracks=" << get_collision_ntracks()
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
  clear_collision_vertices();
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


double FullTrackVertexContainerv1::get_collision_radius(const unsigned int i) const
{
  return i < m_collision_radius.size() ? m_collision_radius[i] : 0.0;
}

double FullTrackVertexContainerv1::get_collision_phi(const unsigned int i) const
{
  return i < m_collision_phi.size() ? m_collision_phi[i] : 0.0;
}

double FullTrackVertexContainerv1::get_collision_timebin(const unsigned int i) const
{
  return i < m_collision_timebin.size() ? m_collision_timebin[i] : 0.0;
}

double FullTrackVertexContainerv1::get_collision_timebin_rms(const unsigned int i) const
{
  return i < m_collision_timebin_rms.size() ? m_collision_timebin_rms[i] : 0.0;
}

unsigned int FullTrackVertexContainerv1::get_collision_ntracks(const unsigned int i) const
{
  return i < m_collision_ntracks.size() ? m_collision_ntracks[i] : 0;
}

void FullTrackVertexContainerv1::set_collision_radius(const double v)
{
  ensure_collision_vertex(0);
  m_collision_radius[0] = v;
}

void FullTrackVertexContainerv1::set_collision_phi(const double v)
{
  ensure_collision_vertex(0);
  m_collision_phi[0] = v;
}

void FullTrackVertexContainerv1::set_collision_timebin(const double v)
{
  ensure_collision_vertex(0);
  m_collision_timebin[0] = v;
}

void FullTrackVertexContainerv1::set_collision_timebin_rms(const double v)
{
  ensure_collision_vertex(0);
  m_collision_timebin_rms[0] = v;
}

void FullTrackVertexContainerv1::set_collision_ntracks(const unsigned int v)
{
  ensure_collision_vertex(0);
  m_collision_ntracks[0] = v;
}

void FullTrackVertexContainerv1::clear_collision_vertices()
{
  m_collision_radius.clear();
  m_collision_phi.clear();
  m_collision_timebin.clear();
  m_collision_timebin_rms.clear();
  m_collision_ntracks.clear();
  m_collision_vertex_valid = 0;
}

void FullTrackVertexContainerv1::add_collision_vertex(const double radius,
                                                      const double phi,
                                                      const double timebin,
                                                      const double timebin_rms,
                                                      const unsigned int ntracks)
{
  m_collision_radius.push_back(radius);
  m_collision_phi.push_back(phi);
  m_collision_timebin.push_back(timebin);
  m_collision_timebin_rms.push_back(timebin_rms);
  m_collision_ntracks.push_back(ntracks);
  m_collision_vertex_valid = m_collision_radius.empty() ? 0 : 1;
}

void FullTrackVertexContainerv1::ensure_collision_vertex(const unsigned int i)
{
  while (m_collision_radius.size() <= i)
  {
    add_collision_vertex(0.0, 0.0, 0.0, 0.0, 0);
  }
}
