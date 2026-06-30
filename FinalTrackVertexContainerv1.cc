#include "FinalTrackVertexContainerv1.h"
#include "FinalTrackVertex.h"

ClassImp(FinalTrackVertexContainerv1)

FinalTrackVertexContainerv1::FinalTrackVertexContainerv1()
{
  Reset();
}

FinalTrackVertexContainerv1::~FinalTrackVertexContainerv1()
{
  Reset();
}

void FinalTrackVertexContainerv1::identify(std::ostream& os) const
{
  os << "FinalTrackVertexContainerv1 with " << m_vertices.size()
     << " vertices and " << m_collision_x.size()
     << " collision vertices" << std::endl;
}

void FinalTrackVertexContainerv1::Reset()
{
  for (FinalTrackVertex* vtx : m_vertices) delete vtx;
  m_vertices.clear();
  m_collision_vertex_valid = 0;
  clear_collision_vertices();
}

int FinalTrackVertexContainerv1::isValid() const
{
  return 1;
}

PHObject* FinalTrackVertexContainerv1::CloneMe() const
{
  FinalTrackVertexContainerv1* copy = new FinalTrackVertexContainerv1();
  copy->m_collision_vertex_valid = m_collision_vertex_valid;
  copy->m_collision_x = m_collision_x;
  copy->m_collision_y = m_collision_y;
  copy->m_collision_z = m_collision_z;
  copy->m_collision_z_rms = m_collision_z_rms;
  copy->m_collision_ntracks = m_collision_ntracks;
  copy->m_collision_min_clusters = m_collision_min_clusters;
  for (FinalTrackVertex* vtx : m_vertices)
  {
    if (vtx) copy->m_vertices.push_back(static_cast<FinalTrackVertex*>(vtx->CloneMe()));
  }
  return copy;
}

const FinalTrackVertex* FinalTrackVertexContainerv1::get_vertex(unsigned int i) const
{
  if (i >= m_vertices.size()) return nullptr;
  return m_vertices[i];
}

FinalTrackVertex* FinalTrackVertexContainerv1::get_vertex(unsigned int i)
{
  if (i >= m_vertices.size()) return nullptr;
  return m_vertices[i];
}

double FinalTrackVertexContainerv1::get_collision_x(unsigned int i) const
{
  if (i >= m_collision_x.size()) return 0.0;
  return m_collision_x[i];
}

double FinalTrackVertexContainerv1::get_collision_y(unsigned int i) const
{
  if (i >= m_collision_y.size()) return 0.0;
  return m_collision_y[i];
}

double FinalTrackVertexContainerv1::get_collision_z(unsigned int i) const
{
  if (i >= m_collision_z.size()) return 0.0;
  return m_collision_z[i];
}

double FinalTrackVertexContainerv1::get_collision_z_rms(unsigned int i) const
{
  if (i >= m_collision_z_rms.size()) return 0.0;
  return m_collision_z_rms[i];
}

unsigned int FinalTrackVertexContainerv1::get_collision_ntracks(unsigned int i) const
{
  if (i >= m_collision_ntracks.size()) return 0;
  return m_collision_ntracks[i];
}

void FinalTrackVertexContainerv1::clear_collision_vertices()
{
  m_collision_x.clear();
  m_collision_y.clear();
  m_collision_z.clear();
  m_collision_z_rms.clear();
  m_collision_ntracks.clear();
}

void FinalTrackVertexContainerv1::add_collision_vertex(double x, double y, double z, double z_rms, unsigned int ntracks)
{
  m_collision_x.push_back(x);
  m_collision_y.push_back(y);
  m_collision_z.push_back(z);
  m_collision_z_rms.push_back(z_rms);
  m_collision_ntracks.push_back(ntracks);
}
