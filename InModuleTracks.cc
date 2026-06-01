#include "InModuleTracks.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include "InModuleTrack.h"
#include "InModuleTrackv1.h"
#include "InModuleTrackContainer.h"
#include "InModuleTrackContainerv1.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/ActsGeometry.h>

#include <pthread.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <stdint.h>
#include <utility>
#include <vector>

// ===================================================================
// Internal helpers (anonymous namespace)
// ===================================================================
namespace
{
  struct BlobAdcSort
  {
    const std::vector<InModuleThreadData::Blob>* blobs;
    BlobAdcSort(const std::vector<InModuleThreadData::Blob>* b) : blobs(b) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      return (*blobs)[a].adc > (*blobs)[b].adc;
    }
  };

  double adc_weight(double adc, double maxadc, double power, double floor_frac)
  {
    if (adc <= 0.0) return 0.0;
    if (maxadc <= 0.0) return 1.0;

    double w = std::pow(adc / maxadc, power);
    const double floor_w = floor_frac;
    if (w < floor_w) w = floor_w;
    return w;
  }

  bool weighted_line_fit(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& w,
                         double& m,
                         double& b,
                         double& chi2,
                         int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size() || x.size() != w.size())
      return false;

    double S = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;

    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      S   += wi;
      Sx  += wi * x[i];
      Sy  += wi * y[i];
      Sxx += wi * x[i] * x[i];
      Sxy += wi * x[i] * y[i];
    }

    const double den = S * Sxx - Sx * Sx;
    if (std::fabs(den) < 1.0e-12) return false;

    m = (S * Sxy - Sx * Sy) / den;
    b = (Sy - m * Sx) / S;

    chi2 = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double r = y[i] - (m * x[i] + b);
      chi2 += wi * r * r;
    }

    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }


  double sagitta_model_fit(const double xrot,
                           const double S,
                           const double x0,
                           const double invR)
  {
    const double dx = xrot - x0;
    const double invR2 = invR * invR;
    const double dx2 = dx * dx;

    // Same truncated expansion as the Python notebook:
    // sagitta = S - 1/2 invR dx^2 - 1/8 invR^3 dx^4 - 1/16 invR^5 dx^6
    return S
      - 0.5 * invR * dx2
      - 0.125 * invR * invR2 * dx2 * dx2
      - 0.0625 * invR * invR2 * invR2 * dx2 * dx2 * dx2;
  }

  bool solve_3x3(double A[3][3], double b[3], double x[3])
  {
    double M[3][4] =
    {
      {A[0][0], A[0][1], A[0][2], b[0]},
      {A[1][0], A[1][1], A[1][2], b[1]},
      {A[2][0], A[2][1], A[2][2], b[2]}
    };

    for (int col = 0; col < 3; ++col)
    {
      int pivot = col;
      for (int row = col + 1; row < 3; ++row)
        if (std::fabs(M[row][col]) > std::fabs(M[pivot][col])) pivot = row;

      if (std::fabs(M[pivot][col]) < 1.0e-20) return false;

      if (pivot != col)
      {
        for (int k = col; k < 4; ++k)
          std::swap(M[col][k], M[pivot][k]);
      }

      const double div = M[col][col];
      for (int k = col; k < 4; ++k) M[col][k] /= div;

      for (int row = 0; row < 3; ++row)
      {
        if (row == col) continue;
        const double factor = M[row][col];
        for (int k = col; k < 4; ++k)
          M[row][k] -= factor * M[col][k];
      }
    }

    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
  }

  bool weighted_sagitta_fit(const std::vector<double>& local_x,
                            const std::vector<double>& local_y,
                            const std::vector<double>& w,
                            double& S,
                            double& x0,
                            double& invR,
                            double& theta,
                            double& bline,
                            double& chi2,
                            int& ndof)
  {
    if (local_x.size() < 3 ||
        local_x.size() != local_y.size() ||
        local_x.size() != w.size())
      return false;

    // First reproduce the notebook convention: fit a straight line in the
    // local x-y plane and rotate/shift so that this line is approximately
    // horizontal.  The sagitta is then fitted in (x_rot, y_rot).
    double mline = 0.0;
    double line_chi2 = 0.0;
    int line_ndof = 0;
    if (!weighted_line_fit(local_x, local_y, w, mline, bline, line_chi2, line_ndof))
      return false;

    theta = std::atan(mline);
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    std::vector<double> xrot;
    std::vector<double> yrot;
    xrot.reserve(local_x.size());
    yrot.reserve(local_x.size());

    double sw = 0.0;
    double sx = 0.0;
    double sy = 0.0;

    for (unsigned int i = 0; i < local_x.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double yy = local_y[i] - bline;
      const double xr =  c * local_x[i] + s * yy;
      const double yr = -s * local_x[i] + c * yy;

      xrot.push_back(xr);
      yrot.push_back(yr);

      sw += wi;
      sx += wi * xr;
      sy += wi * yr;
    }

    if (sw <= 0.0) return false;

    x0 = sx / sw;
    S  = sy / sw;

    // Initial curvature from a linear fit y_rot = S + q*(x_rot-x0)^2.
    std::vector<double> dx2;
    dx2.reserve(xrot.size());
    for (unsigned int i = 0; i < xrot.size(); ++i)
    {
      const double dx = xrot[i] - x0;
      dx2.push_back(dx * dx);
    }

    double q = 0.0;
    double qS = 0.0;
    double qchi2 = 0.0;
    int qndof = 0;
    if (weighted_line_fit(dx2, yrot, w, q, qS, qchi2, qndof))
    {
      S = qS;
      invR = -2.0 * q;
    }
    else
    {
      invR = 0.0;
    }

    // Limit pathological starting values.  TPC local distances are in cm;
    // values much above this correspond to extremely tight, unphysical arcs.
    if (invR >  1.0) invR =  1.0;
    if (invR < -1.0) invR = -1.0;

    chi2 = 0.0;
    for (unsigned int i = 0; i < xrot.size(); ++i)
    {
      const double wi = w[i] > 0.0 ? w[i] : 1.0;
      const double r = yrot[i] - sagitta_model_fit(xrot[i], S, x0, invR);
      chi2 += wi * r * r;
    }

    // Small damped Gauss-Newton refinement of (S, x0, invR).
    double lambda = 1.0e-6;
    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      double A[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
      double bvec[3] = {0.0, 0.0, 0.0};

      for (unsigned int i = 0; i < xrot.size(); ++i)
      {
        const double wi = w[i] > 0.0 ? w[i] : 1.0;
        const double dx = xrot[i] - x0;
        const double dx2v = dx * dx;
        const double dx3 = dx2v * dx;
        const double dx4 = dx2v * dx2v;
        const double dx5 = dx4 * dx;
        const double dx6 = dx4 * dx2v;

        const double invR2 = invR * invR;
        const double invR3 = invR2 * invR;
        const double invR4 = invR2 * invR2;
        const double invR5 = invR4 * invR;

        const double f = S
          - 0.5 * invR * dx2v
          - 0.125 * invR3 * dx4
          - 0.0625 * invR5 * dx6;

        const double resid = yrot[i] - f;

        double J[3];
        J[0] = 1.0;
        J[1] = invR * dx + 0.5 * invR3 * dx3 + 0.375 * invR5 * dx5;
        J[2] = -0.5 * dx2v - 0.375 * invR2 * dx4 - 0.3125 * invR4 * dx6;

        for (int a = 0; a < 3; ++a)
        {
          bvec[a] += wi * J[a] * resid;
          for (int b = 0; b < 3; ++b)
            A[a][b] += wi * J[a] * J[b];
        }
      }

      A[0][0] += lambda;
      A[1][1] += lambda;
      A[2][2] += lambda;

      double delta[3] = {0.0, 0.0, 0.0};
      if (!solve_3x3(A, bvec, delta)) break;

      // Guard against numerical jumps on nearly straight tracks.
      if (delta[0] >  10.0) delta[0] =  10.0;
      if (delta[0] < -10.0) delta[0] = -10.0;
      if (delta[1] >  10.0) delta[1] =  10.0;
      if (delta[1] < -10.0) delta[1] = -10.0;
      if (delta[2] >   0.1) delta[2] =   0.1;
      if (delta[2] <  -0.1) delta[2] =  -0.1;

      const double S_new = S + delta[0];
      const double x0_new = x0 + delta[1];
      double invR_new = invR + delta[2];

      if (invR_new >  1.0) invR_new =  1.0;
      if (invR_new < -1.0) invR_new = -1.0;

      double chi2_new = 0.0;
      for (unsigned int i = 0; i < xrot.size(); ++i)
      {
        const double wi = w[i] > 0.0 ? w[i] : 1.0;
        const double r = yrot[i] - sagitta_model_fit(xrot[i], S_new, x0_new, invR_new);
        chi2_new += wi * r * r;
      }

      if (chi2_new <= chi2)
      {
        S = S_new;
        x0 = x0_new;
        invR = invR_new;
        chi2 = chi2_new;
        lambda *= 0.3;

        if (std::fabs(delta[0]) < 1.0e-6 &&
            std::fabs(delta[1]) < 1.0e-6 &&
            std::fabs(delta[2]) < 1.0e-8)
          break;
      }
      else
      {
        lambda *= 10.0;
      }
    }

    ndof = static_cast<int>(local_x.size()) - 3;
    return true;
  }


  bool fit_track_from_blobs(const std::vector<InModuleThreadData::Blob>& blobs,
                            const std::vector<unsigned int>& idx,
                            double weight_power,
                            double floor_frac,
                            InModuleThreadData::Track& trk)
  {
    if (idx.size() < 2) return false;

    double maxadc = 0.0;
    for (unsigned int i = 0; i < idx.size(); ++i)
      if (blobs[idx[i]].adc > maxadc) maxadc = blobs[idx[i]].adc;

    std::vector<double> x, pad, tbin, w;
    x.reserve(idx.size());
    pad.reserve(idx.size());
    tbin.reserve(idx.size());
    w.reserve(idx.size());

    unsigned int first_layer = 999999, last_layer = 0;

    for (unsigned int i = 0; i < idx.size(); ++i)
    {
      const InModuleThreadData::Blob& bl = blobs[idx[i]];
      x.push_back(static_cast<double>(bl.layer));
      pad.push_back(bl.pad);
      tbin.push_back(bl.tbin);
      w.push_back(adc_weight(bl.adc, maxadc, weight_power, floor_frac));

      if (bl.layer < first_layer) first_layer = bl.layer;
      if (bl.layer > last_layer)  last_layer  = bl.layer;
    }

    double mp = 0.0, bp = 0.0, cp = 0.0;
    double mt = 0.0, bt = 0.0, ct = 0.0;
    int ndp = 0, ndt = 0;

    if (!weighted_line_fit(x, pad,  w, mp, bp, cp, ndp)) return false;
    if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

    trk.first_layer    = first_layer;
    trk.last_layer     = last_layer;
    trk.nblobs         = static_cast<unsigned int>(idx.size());
    trk.nrawhits       = 0;
    trk.pad_slope      = mp;
    trk.pad_intercept  = bp;
    trk.tbin_slope     = mt;
    trk.tbin_intercept = bt;
    trk.chi2_pad       = cp;
    trk.chi2_tbin      = ct;
    trk.ndof_pad       = ndp;
    trk.ndof_tbin      = ndt;
    trk.blob_indices   = idx;
    trk.raw_hit_indices.clear();

    return true;
  }

  void collect_raw_indices_from_blob_chain(const std::vector<InModuleThreadData::Blob>& blobs,
                                           const std::vector<unsigned int>& blob_idx,
                                           std::vector<unsigned int>& raw_idx)
  {
    raw_idx.clear();
    for (unsigned int ib = 0; ib < blob_idx.size(); ++ib)
    {
      const InModuleThreadData::Blob& bl = blobs[blob_idx[ib]];
      for (unsigned int ir = 0; ir < bl.raw_hit_indices.size(); ++ir)
        raw_idx.push_back(bl.raw_hit_indices[ir]);
    }
  }

  bool robust_fit_track_from_raw_hits(const std::vector<InModuleThreadData::RawHit>& raw_hits,
                                      const std::vector<InModuleThreadData::Blob>& blobs,
                                      const std::vector<unsigned int>& blob_idx,
                                      double weight_power,
                                      double floor_frac,
                                      bool use_field_on_fit,
                                      InModuleThreadData::Track& trk)
  {
    std::vector<unsigned int> raw_idx;
    collect_raw_indices_from_blob_chain(blobs, blob_idx, raw_idx);
    if (raw_idx.size() < 2) return false;

    double maxadc = 0.0;
    unsigned int first_layer = 999999, last_layer = 0;

    for (unsigned int i = 0; i < raw_idx.size(); ++i)
    {
      const InModuleThreadData::RawHit& rh = raw_hits[raw_idx[i]];
      if (rh.adc > maxadc) maxadc = rh.adc;
      if (rh.layer < first_layer) first_layer = rh.layer;
      if (rh.layer > last_layer)  last_layer  = rh.layer;
    }

    std::vector<double> x, pad, tbin, base_w, w;
    x.reserve(raw_idx.size());
    pad.reserve(raw_idx.size());
    tbin.reserve(raw_idx.size());
    base_w.reserve(raw_idx.size());
    w.reserve(raw_idx.size());

    std::vector<double> radius;
    std::vector<double> phi;
    std::vector<double> local_x;
    std::vector<double> local_y;

    radius.reserve(raw_idx.size());
    phi.reserve(raw_idx.size());
    local_x.reserve(raw_idx.size());
    local_y.reserve(raw_idx.size());

    for (unsigned int i = 0; i < raw_idx.size(); ++i)
    {
      const InModuleThreadData::RawHit& rh = raw_hits[raw_idx[i]];
      x.push_back(static_cast<double>(rh.layer));
      pad.push_back(static_cast<double>(rh.pad));
      tbin.push_back(static_cast<double>(rh.tbin));
      base_w.push_back(adc_weight(static_cast<double>(rh.adc), maxadc,
                                  weight_power, floor_frac));
      w.push_back(base_w.back());

      const double r  = rh.local_radius;
      const double ph = rh.local_phi;

      radius.push_back(r);
      phi.push_back(ph);

      local_x.push_back(r * std::cos(ph));
      local_y.push_back(r * std::sin(ph));
    }

    double mp = 0.0, bp = 0.0, cp = 0.0;
    double mt = 0.0, bt = 0.0, ct = 0.0;
    int ndp = 0, ndt = 0;

    const double huber_c = 2.5;
    for (unsigned int iter = 0; iter < 5; ++iter)
    {
      if (!weighted_line_fit(x, pad,  w, mp, bp, cp, ndp)) return false;
      if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

      double sw = 0.0, sp2 = 0.0, st2 = 0.0;
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        const double rp = pad[i]  - (mp * x[i] + bp);
        const double rt = tbin[i] - (mt * x[i] + bt);
        sw  += base_w[i];
        sp2 += base_w[i] * rp * rp;
        st2 += base_w[i] * rt * rt;
      }

      double scale_p = 1.0, scale_t = 1.0;
      if (sw > 0.0)
      {
        scale_p = std::sqrt(sp2 / sw);
        scale_t = std::sqrt(st2 / sw);
      }
      if (scale_p < 1.0) scale_p = 1.0;
      if (scale_t < 1.0) scale_t = 1.0;

      for (unsigned int i = 0; i < x.size(); ++i)
      {
        const double rp = pad[i]  - (mp * x[i] + bp);
        const double rt = tbin[i] - (mt * x[i] + bt);
        const double r = std::sqrt((rp / scale_p) * (rp / scale_p) +
                                   (rt / scale_t) * (rt / scale_t));
        double robust = 1.0;
        if (r > huber_c) robust = huber_c / r;
        w[i] = base_w[i] * robust;
      }
    }

    if (!weighted_line_fit(x, pad,  w, mp, bp, cp, ndp)) return false;
    if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

    trk.first_layer    = first_layer;
    trk.last_layer     = last_layer;
    trk.nblobs         = static_cast<unsigned int>(blob_idx.size());
    trk.nrawhits       = static_cast<unsigned int>(raw_idx.size());
    trk.pad_slope      = mp;
    trk.pad_intercept  = bp;
    trk.tbin_slope     = mt;
    trk.tbin_intercept = bt;
    trk.chi2_pad       = cp;
    trk.chi2_tbin      = ct;
    trk.ndof_pad       = ndp;
    trk.ndof_tbin      = ndt;
    trk.blob_indices   = blob_idx;
    trk.raw_hit_indices = raw_idx;
    double mrt = 0.0, brt = 0.0, crt = 0.0;
    int ndrt = 0;

    if (weighted_line_fit(radius, tbin, w, mrt, brt, crt, ndrt))
    {
      trk.has_local_line_fit = true;
      trk.radius_tbin_slope = mrt;
      trk.radius_tbin_intercept = brt;
      trk.chi2_radius_tbin = crt;
      trk.ndof_radius_tbin = ndrt;
    }
    if (!use_field_on_fit)
    {
      double mrp = 0.0, brp = 0.0, crp = 0.0;
      int ndrp = 0;

      if (weighted_line_fit(radius, phi, w, mrp, brp, crp, ndrp))
      {
        trk.has_radius_phi_line_fit = true;
        trk.radius_phi_slope = mrp;
        trk.radius_phi_intercept = brp;
        trk.chi2_radius_phi_line = crp;
        trk.ndof_radius_phi_line = ndrp;
      }
    }
    if (use_field_on_fit)
    {
      double sagS = 0.0;
      double sagx0 = 0.0;
      double sagInvR = 0.0;
      double sagTheta = 0.0;
      double sagB = 0.0;
      double sagChi2 = 0.0;
      int sagNdof = 0;

      // Sagitta fit requested in the (local radius, local phi) plane.
      // The fitter internally rotates this 2D plane to remove the linear trend.
      if (weighted_sagitta_fit(radius, phi, w,
                               sagS, sagx0, sagInvR, sagTheta, sagB,
                               sagChi2, sagNdof))
      {
        trk.has_sagitta_fit = true;
        trk.sagitta_S = sagS;
        trk.sagitta_x0 = sagx0;
        trk.sagitta_invR = sagInvR;
        trk.sagitta_theta = sagTheta;
        trk.sagitta_b = sagB;
        trk.chi2_radius_phi_sagitta = sagChi2;
        trk.ndof_radius_phi_sagitta = sagNdof;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------
  // Track-piece connection
  // -------------------------------------------------------------------
  struct TrackStartSort
  {
    const std::vector<InModuleThreadData::Track>* tracks;
    TrackStartSort(const std::vector<InModuleThreadData::Track>* t) : tracks(t) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const InModuleThreadData::Track& ta = (*tracks)[a];
      const InModuleThreadData::Track& tb = (*tracks)[b];
      if (ta.first_layer != tb.first_layer) return ta.first_layer < tb.first_layer;
      if (ta.last_layer  != tb.last_layer)  return ta.last_layer  < tb.last_layer;
      return ta.nrawhits > tb.nrawhits;
    }
  };

  void append_unique_blob_indices(std::vector<unsigned int>& dst,
                                  const std::vector<unsigned int>& src)
  {
    for (unsigned int i = 0; i < src.size(); ++i)
      if (std::find(dst.begin(), dst.end(), src[i]) == dst.end())
        dst.push_back(src[i]);
  }

  bool tracks_can_connect(const InModuleThreadData::Track& a,
                          const InModuleThreadData::Track& b,
                          unsigned int connect_max_layer_gap,
                          double connect_dp,
                          double connect_dt,
                          double connect_dpad_slope,
                          double connect_dtbin_slope,
                          double& score)
  {
    score = std::numeric_limits<double>::max();

    if (a.last_layer >= b.first_layer) return false;

    const unsigned int gap = b.first_layer - a.last_layer - 1;
    if (gap > connect_max_layer_gap) return false;

    const double lmatch = 0.5 * (static_cast<double>(a.last_layer) +
                                 static_cast<double>(b.first_layer));

    const double pad_a  = a.pad_slope  * lmatch + a.pad_intercept;
    const double pad_b  = b.pad_slope  * lmatch + b.pad_intercept;
    const double tbin_a = a.tbin_slope * lmatch + a.tbin_intercept;
    const double tbin_b = b.tbin_slope * lmatch + b.tbin_intercept;

    const double dp  = std::fabs(pad_a  - pad_b);
    const double dt  = std::fabs(tbin_a - tbin_b);
    const double dmp = std::fabs(a.pad_slope  - b.pad_slope);
    const double dmt = std::fabs(a.tbin_slope - b.tbin_slope);

    if (dp  > connect_dp)         return false;
    if (dt  > connect_dt)         return false;
    if (dmp > connect_dpad_slope) return false;
    if (dmt > connect_dtbin_slope) return false;

    score = (dp  / connect_dp)         * (dp  / connect_dp)
          + (dt  / connect_dt)         * (dt  / connect_dt)
          + (dmp / connect_dpad_slope) * (dmp / connect_dpad_slope)
          + (dmt / connect_dtbin_slope) * (dmt / connect_dtbin_slope)
          + 0.05 * static_cast<double>(gap);

    return true;
  }

  void connect_track_pieces_in_module(InModuleThreadData* d)
  {
    if (!d || d->tracks.size() < 2) return;

    std::vector<InModuleThreadData::Track> pieces = d->tracks;
    std::vector<InModuleThreadData::Track> output;
    std::vector<int> used(pieces.size(), 0);

    std::vector<unsigned int> order;
    order.reserve(pieces.size());
    for (unsigned int i = 0; i < pieces.size(); ++i) order.push_back(i);
    std::sort(order.begin(), order.end(), TrackStartSort(&pieces));

    for (unsigned int io = 0; io < order.size(); ++io)
    {
      const unsigned int iseed = order[io];
      if (used[iseed]) continue;

      InModuleThreadData::Track current = pieces[iseed];
      used[iseed] = 1;

      bool merged_any = true;
      while (merged_any)
      {
        merged_any = false;
        int best_j = -1;
        double best_score = std::numeric_limits<double>::max();

        for (unsigned int jo = 0; jo < order.size(); ++jo)
        {
          const unsigned int j = order[jo];
          if (used[j]) continue;

          double score = 0.0;
          if (!tracks_can_connect(current, pieces[j],
                                  d->connect_max_layer_gap,
                                  d->connect_dp,
                                  d->connect_dt,
                                  d->connect_dpad_slope,
                                  d->connect_dtbin_slope,
                                  score)) continue;

          if (score < best_score) { best_score = score; best_j = static_cast<int>(j); }
        }

        if (best_j >= 0)
        {
          append_unique_blob_indices(current.blob_indices, pieces[best_j].blob_indices);

          InModuleThreadData::Track refit;
          if (robust_fit_track_from_raw_hits(d->raw_hits, d->blobs,
                                   current.blob_indices,
                                   d->weight_power,
                                   d->adc_weight_floor_frac,
                                   d->use_field_on_fit,
                                   refit))
          {
            current  = refit;
            used[best_j] = 1;
            merged_any   = true;
          }
        }
      }

      output.push_back(current);
    }

    for (unsigned int i = 0; i < output.size(); ++i)
      output[i].track_id = i;

    if (d->verbosity > 1)
    {
      std::cout << "InModuleTracks connect pieces: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " pieces=" << pieces.size()
                << " connected_tracks=" << output.size() << std::endl;
    }

    d->tracks.swap(output);
  }

  // -------------------------------------------------------------------
  // Local hardware -> geometry conversion.
  //
  // Keep exactly the same local convention as the original
  // InModuleTracks.cc:
  //
  //   local_phi    = pad * phi_bin_width[region]
  //   local_radius = module_radius[region][layer - 7 - 16*region]
  //
  // Do not apply sector offsets here.  In-module fits must be done in the
  // local module coordinates only.
  // -------------------------------------------------------------------
  double get_local_phi(const TpcPadMap* padMap,
                       const unsigned int region,
                       const unsigned int pad)
  {
    if (!padMap) return 0.0;
    return padMap->get_local_phi(region, static_cast<double>(pad));
  }

  double get_local_radius(const TpcPadMap* padMap,
                          const unsigned int region,
                          const unsigned int layer)
  {
    if (!padMap) return 0.0;
    return padMap->get_local_radius(region, layer);
  }

  // -------------------------------------------------------------------
  // Per-module worker functions (exact originals)
  // -------------------------------------------------------------------

  struct RawHitTimeSort
  {
    const std::vector<InModuleThreadData::RawHit>* raw_hits;
    RawHitTimeSort(const std::vector<InModuleThreadData::RawHit>* h) : raw_hits(h) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const InModuleThreadData::RawHit& ha = (*raw_hits)[a];
      const InModuleThreadData::RawHit& hb = (*raw_hits)[b];
      if (ha.tbin != hb.tbin) return ha.tbin < hb.tbin;
      return ha.adc > hb.adc;
    }
  };

  void reject_long_pad_noise(InModuleThreadData* d)
  {
    if (!d) return;
    if (d->noise_max_consecutive_timebins <= 0) return;
    if (d->noise_keep_first_timebins < 0) d->noise_keep_first_timebins = 0;
    if (d->raw_hits.empty()) return;

    typedef std::pair<unsigned int, unsigned short> LayerPadKey;
    std::map<LayerPadKey, std::vector<unsigned int> > by_layer_pad;

    for (unsigned int i = 0; i < d->raw_hits.size(); ++i)
    {
      const InModuleThreadData::RawHit& rh = d->raw_hits[i];
      by_layer_pad[LayerPadKey(rh.layer, rh.pad)].push_back(i);
    }

    std::vector<int> remove(d->raw_hits.size(), 0);
    unsigned int nremoved = 0;

    for (std::map<LayerPadKey, std::vector<unsigned int> >::iterator it = by_layer_pad.begin();
         it != by_layer_pad.end(); ++it)
    {
      std::vector<unsigned int>& idx = it->second;
      if (idx.size() <= static_cast<unsigned int>(d->noise_max_consecutive_timebins)) continue;

      std::sort(idx.begin(), idx.end(), RawHitTimeSort(&d->raw_hits));

      unsigned int run_start = 0;
      while (run_start < idx.size())
      {
        unsigned int run_end = run_start;
        while (run_end + 1 < idx.size())
        {
          const unsigned short t0 = d->raw_hits[idx[run_end]].tbin;
          const unsigned short t1 = d->raw_hits[idx[run_end + 1]].tbin;
          if (static_cast<int>(t1) != static_cast<int>(t0) + 1) break;
          ++run_end;
        }

        const unsigned int run_len = run_end - run_start + 1;
        if (run_len > static_cast<unsigned int>(d->noise_max_consecutive_timebins))
        {
          const unsigned int keep_until =
            run_start + static_cast<unsigned int>(d->noise_keep_first_timebins);

          //bool in_non_increasing_tail = true;
          for (unsigned int ir = run_start; ir <= run_end; ++ir)
          {
            if (ir < keep_until) continue;

            const unsigned int cur_idx  = idx[ir];
            const unsigned int prev_idx = idx[ir - 1];

            const unsigned short cur_adc  = d->raw_hits[cur_idx].adc;
            const unsigned short prev_adc = d->raw_hits[prev_idx].adc;



            const int cur_adc_i  = static_cast<int>(cur_adc);
            const int prev_adc_i = static_cast<int>(prev_adc);
            const int tol        = d->noise_adc_tolerance;

            // remove tail if ADC is flat, decreasing, or only mildly increasing
            if (cur_adc_i <= prev_adc_i + tol)
            {
              if (!remove[cur_idx])
              {
                remove[cur_idx] = 1;
                ++nremoved;
              }
            }
            else
            {
              // Preserve a later rising signal on top of the long pad tail.
             // in_non_increasing_tail = false;
            }
          }
        }

        run_start = run_end + 1;
      }
    }

    if (nremoved == 0) return;

    std::vector<InModuleThreadData::RawHit> kept;
    kept.reserve(d->raw_hits.size() - nremoved);
    for (unsigned int i = 0; i < d->raw_hits.size(); ++i)
    {
      if (!remove[i]) kept.push_back(d->raw_hits[i]);
    }

    d->raw_hits.swap(kept);

    if (d->verbosity > 1)
    {
      std::cout << "InModuleTracks noise rejection: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " removed " << nremoved << " long same-pad tail hits"
                << std::endl;
    }
  }

  void collect_raw_hits(InModuleThreadData* d)
  {
    d->raw_hits.clear();

    for (unsigned int ihs = 0; ihs < d->layer_hitsets.size(); ++ihs)
    {
      TrkrHitSet* hitset = d->layer_hitsets[ihs].hitset;
      if (!hitset) continue;

      TrkrHitSet::ConstRange range = hitset->getHits();
      for (TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
      {
        const TrkrDefs::hitkey hitkey = hitr->first;
        const TrkrHit* hit = hitr->second;
        if (!hit) continue;

        const unsigned short pad    = TpcDefs::getPad(hitkey);
        const unsigned short tbin   = TpcDefs::getTBin(hitkey);
        const unsigned short rawAdc = hit->getAdc();
        const double fadc = static_cast<double>(rawAdc) - d->pedestal;
        if (fadc <= 0.0) continue;

        InModuleThreadData::RawHit rh;
        rh.layer        = d->layer_hitsets[ihs].layer;
        rh.hitsetkey    = d->layer_hitsets[ihs].hitsetkey;
        rh.hitkey       = hitkey;
        rh.pad          = pad;
        rh.tbin         = tbin;
        rh.adc          = static_cast<unsigned short>(fadc);
        rh.local_phi    = get_local_phi(d->padMap, d->region, pad);
        rh.local_radius = get_local_radius(d->padMap, d->region, rh.layer);
        d->raw_hits.push_back(rh);
      }
    }
  }

  void build_blobs(InModuleThreadData* d)
  {
    d->blobs.clear();
    const unsigned int n = static_cast<unsigned int>(d->raw_hits.size());
    std::vector<int> used(n, 0);

    for (unsigned int i = 0; i < n; ++i)
    {
      if (used[i]) continue;

      used[i] = 1;
      std::deque<unsigned int> q;
      q.push_back(i);

      double sw = 0.0, sp = 0.0, st = 0.0;
      unsigned int nh = 0;
      const unsigned int layer = d->raw_hits[i].layer;

      InModuleThreadData::Blob bl;

      while (!q.empty())
      {
        const unsigned int a = q.front();
        q.pop_front();

        const InModuleThreadData::RawHit& ha = d->raw_hits[a];
        bl.raw_hit_indices.push_back(a);
        const double wa = static_cast<double>(ha.adc);
        sw += wa;
        sp += wa * static_cast<double>(ha.pad);
        st += wa * static_cast<double>(ha.tbin);
        ++nh;

        for (unsigned int j = 0; j < n; ++j)
        {
          if (used[j]) continue;
          const InModuleThreadData::RawHit& hb = d->raw_hits[j];
          if (hb.layer != layer) continue;

          const int dp = std::abs(static_cast<int>(hb.pad)  - static_cast<int>(ha.pad));
          const int dt = std::abs(static_cast<int>(hb.tbin) - static_cast<int>(ha.tbin));
          if (dp <= d->blob_dp && dt <= d->blob_dt)
          {
            used[j] = 1;
            q.push_back(j);
          }
        }
      }

      if (sw <= 0.0) continue;

      bl.layer  = layer;
      bl.pad    = sp / sw;
      bl.tbin   = st / sw;
      bl.adc    = sw;
      bl.nhits  = nh;
      bl.used   = 0;
      d->blobs.push_back(bl);
    }
  }

  int find_best_blob_on_layer(const InModuleThreadData* d,
                              unsigned int target_layer,
                              double pred_pad,
                              double pred_tbin)
  {
    int best = -1;
    double best_score = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < d->blobs.size(); ++i)
    {
      const InModuleThreadData::Blob& bl = d->blobs[i];
      if (bl.used) continue;
      if (bl.layer != target_layer) continue;

      const double dp = bl.pad  - pred_pad;
      const double dt = bl.tbin - pred_tbin;
      if (std::fabs(dp) > d->search_dp) continue;
      if (std::fabs(dt) > d->search_dt) continue;

      const double score = (dp * dp) / (d->search_dp * d->search_dp + 1.0e-9)
                         + (dt * dt) / (d->search_dt * d->search_dt + 1.0e-9)
                         - 0.01 * std::log(bl.adc + 1.0);
      if (score < best_score) { best_score = score; best = static_cast<int>(i); }
    }

    return best;
  }

  void grow_one_direction(InModuleThreadData* d,
                          std::vector<unsigned int>& chain,
                          int direction)
  {
    while (true)
    {
      unsigned int edge_layer = d->blobs[chain.back()].layer;
      if (direction < 0) edge_layer = d->blobs[chain.front()].layer;

      if (direction > 0 && edge_layer >= 54) break;
      if (direction < 0 && edge_layer <= 7)  break;

      const unsigned int target_layer =
        static_cast<unsigned int>(static_cast<int>(edge_layer) + direction);

      double pred_pad  = d->blobs[chain.back()].pad;
      double pred_tbin = d->blobs[chain.back()].tbin;
      if (direction < 0)
      {
        pred_pad  = d->blobs[chain.front()].pad;
        pred_tbin = d->blobs[chain.front()].tbin;
      }

      if (chain.size() >= 2)
      {
        InModuleThreadData::Track tmp;
        if (fit_track_from_blobs(d->blobs, chain,
                                 d->weight_power, d->adc_weight_floor_frac, tmp))
        {
          pred_pad  = tmp.pad_slope  * static_cast<double>(target_layer) + tmp.pad_intercept;
          pred_tbin = tmp.tbin_slope * static_cast<double>(target_layer) + tmp.tbin_intercept;
        }
      }

      const int ibest = find_best_blob_on_layer(d, target_layer, pred_pad, pred_tbin);
      if (ibest < 0) break;

      d->blobs[ibest].used = 1;
      if (direction > 0)
        chain.push_back(static_cast<unsigned int>(ibest));
      else
        chain.insert(chain.begin(), static_cast<unsigned int>(ibest));
    }
  }

  void build_tracks_linear(InModuleThreadData* d)
  {
    d->tracks.clear();

    std::vector<unsigned int> order;
    order.reserve(d->blobs.size());
    for (unsigned int i = 0; i < d->blobs.size(); ++i) order.push_back(i);
    std::sort(order.begin(), order.end(), BlobAdcSort(&d->blobs));

    unsigned int tid = 0;
    for (unsigned int io = 0; io < order.size(); ++io)
    {
      const unsigned int seed = order[io];
      if (d->blobs[seed].used) continue;

      std::vector<unsigned int> chain;
      chain.push_back(seed);
      d->blobs[seed].used = 1;

      grow_one_direction(d, chain, +1);
      grow_one_direction(d, chain, -1);

      if (chain.size() < d->min_track_blobs)
      {
        for (unsigned int k = 0; k < chain.size(); ++k)
          d->blobs[chain[k]].used = 0;
        continue;
      }

      InModuleThreadData::Track trk;
      trk.track_id = tid;
      if (robust_fit_track_from_raw_hits(d->raw_hits, d->blobs, chain,
                                        d->weight_power,
                                        d->adc_weight_floor_frac,
                                        d->use_field_on_fit,
                                        trk))
      {
        trk.track_id = tid;
        d->tracks.push_back(trk);
        ++tid;
      }
    }
  }

  void* ProcessModule(void* arg)
  {
    InModuleThreadData* d = static_cast<InModuleThreadData*>(arg);
    if (!d) return 0;

    collect_raw_hits(d);
    reject_long_pad_noise(d);
    build_blobs(d);
    build_tracks_linear(d);
    connect_track_pieces_in_module(d);

    if (d->verbosity > 1)
    {
      std::cout << "InModuleTracks worker: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " raw_hits=" << d->raw_hits.size()
                << " blobs="    << d->blobs.size()
                << " tracks="   << d->tracks.size() << std::endl;
    }

    return 0;
  }

} // anonymous namespace

// ===================================================================
// Struct constructors
// ===================================================================
InModuleThreadData::LayerHitSet::LayerHitSet()
  : layer(0), hitsetkey(0), hitset(0) {}

InModuleThreadData::RawHit::RawHit()
  : layer(0), hitsetkey(0), hitkey(0),
    pad(0), tbin(0), adc(0),
    local_phi(0.0), local_radius(0.0) {}

InModuleThreadData::Blob::Blob()
  : layer(0), pad(0.0), tbin(0.0), adc(0.0), nhits(0), used(0) {}

InModuleThreadData::Track::Track()
  : track_id(0), first_layer(0), last_layer(0),
    nblobs(0), nrawhits(0),
    pad_slope(0.0), pad_intercept(0.0),
    tbin_slope(0.0), tbin_intercept(0.0),
    chi2_pad(0.0), chi2_tbin(0.0),
    ndof_pad(0), ndof_tbin(0),
    has_local_line_fit(false),
    has_radius_phi_line_fit(false),
    has_radius_phi_circle_fit(false),
    has_sagitta_fit(false),
    sagitta_S(0.0),
    sagitta_x0(0.0),
    sagitta_invR(0.0),
    sagitta_theta(0.0),
    sagitta_b(0.0),
    radius_tbin_slope(0.0),
    radius_tbin_intercept(0.0),
    chi2_radius_tbin(0.0),
    ndof_radius_tbin(0),
    radius_phi_slope(0.0),
    radius_phi_intercept(0.0),
    chi2_radius_phi_line(0.0),
    ndof_radius_phi_line(0),
    circle_x0(0.0),
    circle_y0(0.0),
    circle_radius(0.0),
    chi2_radius_phi_circle(0.0),
    ndof_radius_phi_circle(0),
    chi2_radius_phi_sagitta(0.0),
    ndof_radius_phi_sagitta(0) {}

InModuleThreadData::InModuleThreadData()
  : region(0), sector(0), side(0), module_key(0), tGeometry(0), padMap(0),
    pedestal(74.4), verbosity(0), use_field_on_fit(false),
    noise_max_consecutive_timebins(10), noise_keep_first_timebins(3), noise_adc_tolerance(5),
    blob_dt(2), blob_dp(2),
    search_dt(6), search_dp(6),
    min_track_blobs(4),
    connect_max_layer_gap(8),
    connect_dp(8.0), connect_dt(8.0),
    connect_dpad_slope(2.0), connect_dtbin_slope(2.0),
    weight_power(0.5), adc_weight_floor_frac(0.15) {}

// ===================================================================
// InModuleTracks
// ===================================================================
InModuleTracks::InModuleTracks(const std::string& name,
                               const std::string& filename)
  : SubsysReco(name),
    m_outputFileName(filename),
    m_outputFile(0),
    m_tree(0),
    m_hits(0),
    m_tGeometry(0),
    m_padMap(0),
    m_inModuleTrackContainer(0),
    m_event(0),
    m_maxThreads(72),
    m_useFieldOnFit(true),
    m_pedestal(0.0),
    m_noiseMaxConsecutiveTimebins(10),
    m_noiseKeepFirstTimebins(3),
    m_noiseAdcTolerance(5),
    m_blob_dt(2),
    m_blob_dp(2),
    m_search_dt(6),
    m_search_dp(6),
    m_minTrackBlobs(4),
    m_connectMaxLayerGap(8),
    m_connect_dp(8.0),
    m_connect_dt(8.0),
    m_connect_dpad_slope(2.0),
    m_connect_dtbin_slope(2.0),
    m_tree_event(0)
{
}

InModuleTracks::~InModuleTracks()
{
  if (m_outputFile)
  {
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = 0;
  }
}

void InModuleTracks::setMaxThreads(unsigned int n)
{
  m_maxThreads = (n == 0) ? 1 : n;
}

int InModuleTracks::Init(PHCompositeNode*)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("InModuleTracks", "TPC in-module straight-line pattern recognition");
  m_tree->Branch("event",       &m_tree_event, "event/I");

  m_tree->Branch("track_id",    &m_tree_track_id);
  m_tree->Branch("region",      &m_tree_region);
  m_tree->Branch("sector",      &m_tree_sector);
  m_tree->Branch("side",        &m_tree_side);
  m_tree->Branch("nblobs",      &m_tree_nblobs);
  m_tree->Branch("nrawhits",    &m_tree_nrawhits);
  m_tree->Branch("first_layer", &m_tree_first_layer);
  m_tree->Branch("last_layer",  &m_tree_last_layer);

  m_tree->Branch("pad_slope",      &m_tree_pad_slope);
  m_tree->Branch("pad_intercept",  &m_tree_pad_intercept);
  m_tree->Branch("tbin_slope",     &m_tree_tbin_slope);
  m_tree->Branch("tbin_intercept", &m_tree_tbin_intercept);
  m_tree->Branch("chi2_pad",       &m_tree_chi2_pad);
  m_tree->Branch("chi2_tbin",      &m_tree_chi2_tbin);
  m_tree->Branch("ndof_pad",       &m_tree_ndof_pad);
  m_tree->Branch("ndof_tbin",      &m_tree_ndof_tbin);

  m_tree->Branch("has_sagitta_fit",          &m_tree_has_sagitta_fit);
  m_tree->Branch("sagitta_S",                &m_tree_sagitta_S);
  m_tree->Branch("sagitta_x0",               &m_tree_sagitta_x0);
  m_tree->Branch("sagitta_invR",             &m_tree_sagitta_invR);
  m_tree->Branch("sagitta_theta",            &m_tree_sagitta_theta);
  m_tree->Branch("sagitta_b",                &m_tree_sagitta_b);
  m_tree->Branch("chi2_radius_phi_sagitta",  &m_tree_chi2_radius_phi_sagitta);
  m_tree->Branch("ndof_radius_phi_sagitta",  &m_tree_ndof_radius_phi_sagitta);

  // Per-hit branches: TrkrHitSetContainer keys only
  m_tree->Branch("hit_event",     &m_tree_hit_event);
  m_tree->Branch("hit_track_id",  &m_tree_hit_track_id);
  m_tree->Branch("hit_region",    &m_tree_hit_region);
  m_tree->Branch("hit_sector",    &m_tree_hit_sector);
  m_tree->Branch("hit_side",      &m_tree_hit_side);
  m_tree->Branch("hit_layer",     &m_tree_hit_layer);
  m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey",    &m_tree_hit_hitkey);

  std::cout << Name() << "::Init - output file " << m_outputFileName << " created" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode)    != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::End(PHCompositeNode*)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (m_tree) m_tree->Write();
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = 0;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::getNodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cerr << Name() << "::getNodes - missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_padMap = findNode::getClass<TpcPadMap>(topNode, "TPC_PADMAP");
  if (!m_padMap)
  {
    std::cerr << Name() << "::getNodes - missing TPC_PADMAP" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_padMap->isValid())
  {
    std::cerr << Name() << "::getNodes - TPC_PADMAP is invalid" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode =
    dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_inModuleTrackContainer =
    findNode::getClass<InModuleTrackContainer>(topNode, "INMODULETRACKS");

  if (!m_inModuleTrackContainer)
  {
    m_inModuleTrackContainer = new InModuleTrackContainerv1();

    PHIODataNode<PHObject>* node =
      new PHIODataNode<PHObject>(m_inModuleTrackContainer,
                                 "INMODULETRACKS", "PHObject");
    dstNode->addNode(node);

    std::cout << Name() << "::createNodes - created INMODULETRACKS node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void InModuleTracks::reset_tree_vars()
{
  m_tree_event = m_event;

  m_tree_track_id.clear();
  m_tree_region.clear();
  m_tree_sector.clear();
  m_tree_side.clear();
  m_tree_nblobs.clear();
  m_tree_nrawhits.clear();
  m_tree_first_layer.clear();
  m_tree_last_layer.clear();

  m_tree_pad_slope.clear();
  m_tree_pad_intercept.clear();
  m_tree_tbin_slope.clear();
  m_tree_tbin_intercept.clear();
  m_tree_chi2_pad.clear();
  m_tree_chi2_tbin.clear();
  m_tree_ndof_pad.clear();
  m_tree_ndof_tbin.clear();

  m_tree_has_sagitta_fit.clear();
  m_tree_sagitta_S.clear();
  m_tree_sagitta_x0.clear();
  m_tree_sagitta_invR.clear();
  m_tree_sagitta_theta.clear();
  m_tree_sagitta_b.clear();
  m_tree_chi2_radius_phi_sagitta.clear();
  m_tree_ndof_radius_phi_sagitta.clear();

  m_tree_hit_event.clear();
  m_tree_hit_track_id.clear();
  m_tree_hit_region.clear();
  m_tree_hit_sector.clear();
  m_tree_hit_side.clear();
  m_tree_hit_layer.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
}

int InModuleTracks::process_event(PHCompositeNode*)
{
  reset_tree_vars();

  if (m_inModuleTrackContainer) m_inModuleTrackContainer->Reset();

  std::vector<InModuleThreadData> tdata;
  tdata.reserve(72);

  for (unsigned int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int region = 0; region < 3; ++region)
      {
        InModuleThreadData td;
        td.region     = region;
        td.sector     = sector;
        td.side       = static_cast<int>(side);
        td.module_key = TpcDefs::genModuleHitSetKey(static_cast<uint8_t>(region),
                                                    static_cast<uint8_t>(sector),
                                                    static_cast<uint8_t>(side));
        td.tGeometry          = m_tGeometry;
        td.padMap             = m_padMap;
        td.pedestal           = m_pedestal;
        td.verbosity          = Verbosity();
        td.use_field_on_fit   = m_useFieldOnFit;
        td.noise_max_consecutive_timebins = m_noiseMaxConsecutiveTimebins;
        td.noise_keep_first_timebins      = m_noiseKeepFirstTimebins;
        td.noise_adc_tolerance = m_noiseAdcTolerance;
        td.blob_dt            = m_blob_dt;
        td.blob_dp            = m_blob_dp;
        td.search_dt          = m_search_dt;
        td.search_dp          = m_search_dp;
        td.min_track_blobs    = m_minTrackBlobs;
        td.connect_max_layer_gap = m_connectMaxLayerGap;
        td.connect_dp         = m_connect_dp;
        td.connect_dt         = m_connect_dt;
        td.connect_dpad_slope = m_connect_dpad_slope;
        td.connect_dtbin_slope = m_connect_dtbin_slope;

        for (unsigned int l = 0; l < 16; ++l)
        {
          const unsigned int layer = region * 16 + l + 7;
          const TrkrDefs::hitsetkey hitset_key = TpcDefs::genHitSetKey(layer, sector, side);
          TrkrHitSet* hitset = m_hits->findHitSet(hitset_key);
          if (!hitset) continue;

          InModuleThreadData::LayerHitSet lhs;
          lhs.layer     = layer;
          lhs.hitsetkey = hitset_key;
          lhs.hitset    = hitset;
          td.layer_hitsets.push_back(lhs);
        }

        if (!td.layer_hitsets.empty()) tdata.push_back(td);
      }
    }
  }

  std::cout << Name() << "::process_event - event " << m_event
            << " has " << tdata.size() << " non-empty modules" << std::endl;

  const unsigned int maxLive = std::max(1u,
    std::min(m_maxThreads, static_cast<unsigned int>(tdata.size())));

  for (unsigned int start = 0; start < static_cast<unsigned int>(tdata.size()); start += maxLive)
  {
    const unsigned int end   = std::min(start + maxLive,
                                        static_cast<unsigned int>(tdata.size()));
    const unsigned int nLive = end - start;

    std::vector<pthread_t> threads(nLive);
    std::vector<int>       thread_ok(nLive, 0);

    for (unsigned int i = 0; i < nLive; ++i)
    {
      const unsigned int idx = start + i;
      const int rc = pthread_create(&threads[i], 0, ProcessModule,
                                    static_cast<void*>(&tdata[idx]));
      if (rc != 0)
      {
        std::cerr << Name() << "::process_event - pthread_create failed for"
                  << " region=" << tdata[idx].region
                  << " sector=" << tdata[idx].sector
                  << " side="   << tdata[idx].side << std::endl;
      }
      else
      {
        thread_ok[i] = 1;
      }
    }

    for (unsigned int i = 0; i < nLive; ++i)
      if (thread_ok[i]) pthread_join(threads[i], 0);
  }

  // Harvest results from all modules
  for (unsigned int im = 0; im < tdata.size(); ++im)
  {
    const InModuleThreadData& td = tdata[im];

    for (unsigned int it = 0; it < td.tracks.size(); ++it)
    {
      const InModuleThreadData::Track& tr = td.tracks[it];

      const unsigned int global_track_id =
        m_inModuleTrackContainer ?
        m_inModuleTrackContainer->size() :
        static_cast<unsigned int>(m_tree_track_id.size());

      InModuleTrackv1* outTrack = new InModuleTrackv1();

      outTrack->set_event(static_cast<unsigned int>(m_event));
      outTrack->set_track_id(global_track_id);
      outTrack->set_region(td.region);
      outTrack->set_sector(td.sector);
      outTrack->set_side(td.side);
      outTrack->set_nblobs(tr.nblobs);
      outTrack->set_nrawhits(tr.nrawhits);
      outTrack->set_first_layer(tr.first_layer);
      outTrack->set_last_layer(tr.last_layer);

      outTrack->set_pad_slope(tr.pad_slope);
      outTrack->set_pad_intercept(tr.pad_intercept);
      outTrack->set_tbin_slope(tr.tbin_slope);
      outTrack->set_tbin_intercept(tr.tbin_intercept);
      outTrack->set_chi2_pad(tr.chi2_pad);
      outTrack->set_chi2_tbin(tr.chi2_tbin);
      outTrack->set_ndof_pad(tr.ndof_pad);
      outTrack->set_ndof_tbin(tr.ndof_tbin);

      outTrack->set_local_fit_mode(m_useFieldOnFit ? 1 : 0);

      outTrack->set_radius_tbin_slope(tr.radius_tbin_slope);
      outTrack->set_radius_tbin_intercept(tr.radius_tbin_intercept);
      outTrack->set_chi2_radius_tbin(tr.chi2_radius_tbin);
      outTrack->set_ndof_radius_tbin(tr.ndof_radius_tbin);

      outTrack->set_radius_phi_slope(tr.radius_phi_slope);
      outTrack->set_radius_phi_intercept(tr.radius_phi_intercept);
      outTrack->set_chi2_radius_phi_line(tr.chi2_radius_phi_line);
      outTrack->set_ndof_radius_phi_line(tr.ndof_radius_phi_line);

      outTrack->set_circle_x0(tr.circle_x0);
      outTrack->set_circle_y0(tr.circle_y0);
      outTrack->set_circle_radius(tr.circle_radius);
      outTrack->set_chi2_radius_phi_circle(tr.chi2_radius_phi_circle);
      outTrack->set_ndof_radius_phi_circle(tr.ndof_radius_phi_circle);

      outTrack->set_has_sagitta_fit(tr.has_sagitta_fit ? 1 : 0);
      outTrack->set_sagitta_S(tr.sagitta_S);
      outTrack->set_sagitta_x0(tr.sagitta_x0);
      outTrack->set_sagitta_invR(tr.sagitta_invR);
      outTrack->set_sagitta_theta(tr.sagitta_theta);
      outTrack->set_sagitta_b(tr.sagitta_b);
      outTrack->set_chi2_radius_phi_sagitta(tr.chi2_radius_phi_sagitta);
      outTrack->set_ndof_radius_phi_sagitta(tr.ndof_radius_phi_sagitta);

      // TTree track-level fill
      m_tree_track_id.push_back(global_track_id);
      m_tree_region.push_back(td.region);
      m_tree_sector.push_back(td.sector);
      m_tree_side.push_back(td.side);
      m_tree_nblobs.push_back(tr.nblobs);
      m_tree_nrawhits.push_back(tr.nrawhits);
      m_tree_first_layer.push_back(tr.first_layer);
      m_tree_last_layer.push_back(tr.last_layer);
      m_tree_pad_slope.push_back(tr.pad_slope);
      m_tree_pad_intercept.push_back(tr.pad_intercept);
      m_tree_tbin_slope.push_back(tr.tbin_slope);
      m_tree_tbin_intercept.push_back(tr.tbin_intercept);
      m_tree_chi2_pad.push_back(tr.chi2_pad);
      m_tree_chi2_tbin.push_back(tr.chi2_tbin);
      m_tree_ndof_pad.push_back(tr.ndof_pad);
      m_tree_ndof_tbin.push_back(tr.ndof_tbin);
      m_tree_has_sagitta_fit.push_back(tr.has_sagitta_fit ? 1 : 0);
      m_tree_sagitta_S.push_back(tr.sagitta_S);
      m_tree_sagitta_x0.push_back(tr.sagitta_x0);
      m_tree_sagitta_invR.push_back(tr.sagitta_invR);
      m_tree_sagitta_theta.push_back(tr.sagitta_theta);
      m_tree_sagitta_b.push_back(tr.sagitta_b);
      m_tree_chi2_radius_phi_sagitta.push_back(tr.chi2_radius_phi_sagitta);
      m_tree_ndof_radius_phi_sagitta.push_back(tr.ndof_radius_phi_sagitta);

      // Hit references: (hitsetkey, hitkey) only — no data copy
      for (unsigned int ih = 0; ih < tr.raw_hit_indices.size(); ++ih)
      {
        const InModuleThreadData::RawHit& rh = td.raw_hits[tr.raw_hit_indices[ih]];

        outTrack->add_hit_index(rh.hitsetkey, rh.hitkey);

        m_tree_hit_event.push_back(static_cast<unsigned int>(m_event));
        m_tree_hit_track_id.push_back(global_track_id);
        m_tree_hit_region.push_back(td.region);
        m_tree_hit_sector.push_back(td.sector);
        m_tree_hit_side.push_back(td.side);
        m_tree_hit_layer.push_back(rh.layer);
        m_tree_hit_hitsetkey.push_back(static_cast<unsigned long long>(rh.hitsetkey));
        m_tree_hit_hitkey.push_back(static_cast<unsigned long long>(rh.hitkey));
      }

      if (m_inModuleTrackContainer)
        m_inModuleTrackContainer->add_track(outTrack);
      else
        delete outTrack;
    }
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " tracks=" << m_tree_track_id.size()
              << " track-raw-hits=" << m_tree_hit_track_id.size() << std::endl;
  }
  if (m_inModuleTrackContainer && Verbosity() > 0)
    m_inModuleTrackContainer->identify();
  if (m_inModuleTrackContainer && Verbosity() > 1)
  {
    for (unsigned int i = 0; i < m_inModuleTrackContainer->size(); ++i)
    {
      const InModuleTrack* trk = m_inModuleTrackContainer->get_track(i);
      if (trk) trk->identify();
    }
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}