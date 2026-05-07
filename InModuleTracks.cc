#include "InModuleTracks.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <pthread.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <stdint.h>
#include <vector>

// ===================================================================
// Small C++98 helpers
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
    {
      return false;
    }

    double S = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;

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
    if (std::fabs(den) < 1.0e-12)
    {
      return false;
    }

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

  bool fit_track_from_blobs(const std::vector<InModuleThreadData::Blob>& blobs,
                            const std::vector<unsigned int>& idx,
                            double weight_power,
                            double floor_frac,
                            InModuleThreadData::Track& trk)
  {
    // This blob-center fit is used only during pattern-recognition growth
    // to predict the next layer. The final saved fit is done later on raw hits.
    if (idx.size() < 2) return false;

    double maxadc = 0.0;
    for (unsigned int i = 0; i < idx.size(); ++i)
    {
      if (blobs[idx[i]].adc > maxadc) maxadc = blobs[idx[i]].adc;
    }

    std::vector<double> x;
    std::vector<double> pad;
    std::vector<double> tbin;
    std::vector<double> w;
    x.reserve(idx.size());
    pad.reserve(idx.size());
    tbin.reserve(idx.size());
    w.reserve(idx.size());

    unsigned int first_layer = 999999;
    unsigned int last_layer = 0;

    for (unsigned int i = 0; i < idx.size(); ++i)
    {
      const InModuleThreadData::Blob& bl = blobs[idx[i]];
      x.push_back(static_cast<double>(bl.layer));
      pad.push_back(bl.pad);
      tbin.push_back(bl.tbin);
      w.push_back(adc_weight(bl.adc, maxadc, weight_power, floor_frac));

      if (bl.layer < first_layer) first_layer = bl.layer;
      if (bl.layer > last_layer) last_layer = bl.layer;
    }

    double mp = 0.0, bp = 0.0, cp = 0.0;
    double mt = 0.0, bt = 0.0, ct = 0.0;
    int ndp = 0, ndt = 0;

    if (!weighted_line_fit(x, pad, w, mp, bp, cp, ndp)) return false;
    if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

    trk.first_layer = first_layer;
    trk.last_layer = last_layer;
    trk.nblobs = static_cast<unsigned int>(idx.size());
    trk.nrawhits = 0;
    trk.pad_slope = mp;
    trk.pad_intercept = bp;
    trk.tbin_slope = mt;
    trk.tbin_intercept = bt;
    trk.chi2_pad = cp;
    trk.chi2_tbin = ct;
    trk.ndof_pad = ndp;
    trk.ndof_tbin = ndt;
    trk.blob_indices = idx;
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
      {
        raw_idx.push_back(bl.raw_hit_indices[ir]);
      }
    }
  }

  bool robust_fit_track_from_raw_hits(const std::vector<InModuleThreadData::RawHit>& raw_hits,
                                      const std::vector<InModuleThreadData::Blob>& blobs,
                                      const std::vector<unsigned int>& blob_idx,
                                      double weight_power,
                                      double floor_frac,
                                      InModuleThreadData::Track& trk)
  {
    // Final fit: use every raw ADC cell belonging to the blobs on this track.
    // ADC gives the base weight. A simple Huber-style iterative reweighting
    // down-weights raw cells with large 2D residuals in pad/tbin.
    std::vector<unsigned int> raw_idx;
    collect_raw_indices_from_blob_chain(blobs, blob_idx, raw_idx);
    if (raw_idx.size() < 2) return false;

    double maxadc = 0.0;
    unsigned int first_layer = 999999;
    unsigned int last_layer = 0;

    for (unsigned int i = 0; i < raw_idx.size(); ++i)
    {
      const InModuleThreadData::RawHit& rh = raw_hits[raw_idx[i]];
      if (rh.adc > maxadc) maxadc = rh.adc;
      if (rh.layer < first_layer) first_layer = rh.layer;
      if (rh.layer > last_layer) last_layer = rh.layer;
    }

    std::vector<double> x;
    std::vector<double> pad;
    std::vector<double> tbin;
    std::vector<double> base_w;
    std::vector<double> w;

    x.reserve(raw_idx.size());
    pad.reserve(raw_idx.size());
    tbin.reserve(raw_idx.size());
    base_w.reserve(raw_idx.size());
    w.reserve(raw_idx.size());

    for (unsigned int i = 0; i < raw_idx.size(); ++i)
    {
      const InModuleThreadData::RawHit& rh = raw_hits[raw_idx[i]];
      x.push_back(static_cast<double>(rh.layer));
      pad.push_back(static_cast<double>(rh.pad));
      tbin.push_back(static_cast<double>(rh.tbin));
      base_w.push_back(adc_weight(static_cast<double>(rh.adc), maxadc,
                                  weight_power, floor_frac));
      w.push_back(base_w.back());
    }

    double mp = 0.0, bp = 0.0, cp = 0.0;
    double mt = 0.0, bt = 0.0, ct = 0.0;
    int ndp = 0, ndt = 0;

    const double huber_c = 2.5;
    for (unsigned int iter = 0; iter < 5; ++iter)
    {
      if (!weighted_line_fit(x, pad, w, mp, bp, cp, ndp)) return false;
      if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

      double sw = 0.0;
      double sp2 = 0.0;
      double st2 = 0.0;
      for (unsigned int i = 0; i < x.size(); ++i)
      {
        const double rp = pad[i] - (mp * x[i] + bp);
        const double rt = tbin[i] - (mt * x[i] + bt);
        sw += base_w[i];
        sp2 += base_w[i] * rp * rp;
        st2 += base_w[i] * rt * rt;
      }

      double scale_p = 1.0;
      double scale_t = 1.0;
      if (sw > 0.0)
      {
        scale_p = std::sqrt(sp2 / sw);
        scale_t = std::sqrt(st2 / sw);
      }
      if (scale_p < 1.0) scale_p = 1.0;
      if (scale_t < 1.0) scale_t = 1.0;

      for (unsigned int i = 0; i < x.size(); ++i)
      {
        const double rp = pad[i] - (mp * x[i] + bp);
        const double rt = tbin[i] - (mt * x[i] + bt);
        const double r = std::sqrt((rp / scale_p) * (rp / scale_p) +
                                   (rt / scale_t) * (rt / scale_t));
        double robust = 1.0;
        if (r > huber_c) robust = huber_c / r;
        w[i] = base_w[i] * robust;
      }
    }

    if (!weighted_line_fit(x, pad, w, mp, bp, cp, ndp)) return false;
    if (!weighted_line_fit(x, tbin, w, mt, bt, ct, ndt)) return false;

    trk.first_layer = first_layer;
    trk.last_layer = last_layer;
    trk.nblobs = static_cast<unsigned int>(blob_idx.size());
    trk.nrawhits = static_cast<unsigned int>(raw_idx.size());
    trk.pad_slope = mp;
    trk.pad_intercept = bp;
    trk.tbin_slope = mt;
    trk.tbin_intercept = bt;
    trk.chi2_pad = cp;
    trk.chi2_tbin = ct;
    trk.ndof_pad = ndp;
    trk.ndof_tbin = ndt;
    trk.blob_indices = blob_idx;
    trk.raw_hit_indices = raw_idx;

    return true;
  }



  // -------------------------------------------------------------------
  // Connect disconnected track pieces inside the same TPC module.
  //
  // The first pattern-recognition pass builds short chains and stops when
  // a dead region or empty layer breaks the chain.  This second pass treats
  // those short tracks as pieces of a possible longer track.  Two pieces are
  // merged if their extrapolated pad/tbin positions agree in the empty gap
  // and their slopes are compatible.  After every accepted merge the final
  // track is refit using the raw ADC cells, not blob centers.
  // -------------------------------------------------------------------
  static const unsigned int CONNECT_MAX_LAYER_GAP = 8;  // missing layers allowed between pieces
  static const double CONNECT_DP = 8.0;                 // pad window at gap midpoint
  static const double CONNECT_DT = 8.0;                 // tbin window at gap midpoint
  static const double CONNECT_DPAD_SLOPE = 2.0;         // pad/layer slope compatibility
  static const double CONNECT_DTBIN_SLOPE = 2.0;        // tbin/layer slope compatibility

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
    {
      if (std::find(dst.begin(), dst.end(), src[i]) == dst.end())
      {
        dst.push_back(src[i]);
      }
    }
  }

  bool tracks_can_connect(const InModuleThreadData::Track& a,
                          const InModuleThreadData::Track& b,
                          double& score)
  {
    score = std::numeric_limits<double>::max();

    // Require separated pieces in layer.  This avoids merging two different
    // candidates that overlap in the same rows.
    if (a.last_layer >= b.first_layer) return false;

    const unsigned int gap = b.first_layer - a.last_layer - 1;
    if (gap > CONNECT_MAX_LAYER_GAP) return false;

    // Compare the two extrapolations in the middle of the missing region.
    const double lmatch = 0.5 * (static_cast<double>(a.last_layer) +
                                 static_cast<double>(b.first_layer));

    const double pad_a  = a.pad_slope  * lmatch + a.pad_intercept;
    const double pad_b  = b.pad_slope  * lmatch + b.pad_intercept;
    const double tbin_a = a.tbin_slope * lmatch + a.tbin_intercept;
    const double tbin_b = b.tbin_slope * lmatch + b.tbin_intercept;

    const double dp = std::fabs(pad_a - pad_b);
    const double dt = std::fabs(tbin_a - tbin_b);
    const double dmp = std::fabs(a.pad_slope  - b.pad_slope);
    const double dmt = std::fabs(a.tbin_slope - b.tbin_slope);

    if (dp  > CONNECT_DP) return false;
    if (dt  > CONNECT_DT) return false;
    if (dmp > CONNECT_DPAD_SLOPE) return false;
    if (dmt > CONNECT_DTBIN_SLOPE) return false;

    score = (dp  / CONNECT_DP) * (dp  / CONNECT_DP)
          + (dt  / CONNECT_DT) * (dt  / CONNECT_DT)
          + (dmp / CONNECT_DPAD_SLOPE) * (dmp / CONNECT_DPAD_SLOPE)
          + (dmt / CONNECT_DTBIN_SLOPE) * (dmt / CONNECT_DTBIN_SLOPE)
          + 0.05 * static_cast<double>(gap);

    return true;
  }

  void connect_track_pieces_in_module(InModuleThreadData* d)
  {
    if (!d) return;
    if (d->tracks.size() < 2) return;

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
          if (!tracks_can_connect(current, pieces[j], score)) continue;

          if (score < best_score)
          {
            best_score = score;
            best_j = static_cast<int>(j);
          }
        }

        if (best_j >= 0)
        {
          append_unique_blob_indices(current.blob_indices, pieces[best_j].blob_indices);

          InModuleThreadData::Track refit;
          if (robust_fit_track_from_raw_hits(d->raw_hits, d->blobs,
                                             current.blob_indices,
                                             d->weight_power,
                                             d->adc_weight_floor_frac,
                                             refit))
          {
            current = refit;
            used[best_j] = 1;
            merged_any = true;
          }
        }
      }

      output.push_back(current);
    }

    for (unsigned int i = 0; i < output.size(); ++i)
    {
      output[i].track_id = i;
    }

    if (d->verbosity > 1)
    {
      std::cout << "InModuleTracks connect pieces: region=" << d->region
                << " sector=" << d->sector
                << " side=" << d->side
                << " pieces=" << pieces.size()
                << " connected_tracks=" << output.size()
                << std::endl;
    }

    d->tracks.swap(output);
  }

  static const int N_MODULES = 3;
  static const int N_ROWS = 16;

  static const int Npads[N_MODULES] =
  {
    94, 128, 192
  };

  static const double phi_bin_width[N_MODULES] =
  {
    0.0053073,
    0.003959,
    0.00265145
  };

  static const double module_radius[N_MODULES][N_ROWS] =
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

  double get_local_phi(unsigned int region, unsigned int pad)
  {
    if (region >= N_MODULES) return 0.0;
    return static_cast<double>(pad) * phi_bin_width[region];
  }

  double get_local_radius(unsigned int region, unsigned int layer)
  {
    if (region >= N_MODULES) return 0.0;

    const int ilayer =
      static_cast<int>(layer) - 7 - static_cast<int>(region) * 16;

    if (ilayer < 0 || ilayer >= N_ROWS) return 0.0;

    return module_radius[region][ilayer];
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

        const unsigned short pad = TpcDefs::getPad(hitkey);
        const unsigned short tbin = TpcDefs::getTBin(hitkey);
        const unsigned short rawAdc = hit->getAdc();
        const double fadc = static_cast<double>(rawAdc) - d->pedestal;
        if (fadc <= 0.0) continue;

        InModuleThreadData::RawHit rh;
        rh.layer = d->layer_hitsets[ihs].layer;
        rh.hitsetkey = d->layer_hitsets[ihs].hitsetkey;
        rh.hitkey = hitkey;
        rh.pad = pad;
        rh.tbin = tbin;
        rh.adc = static_cast<unsigned short>(fadc);
        rh.local_phi = get_local_phi(d->region, pad);
        rh.local_radius = get_local_radius(d->region, rh.layer);
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

      double sw = 0.0;
      double sp = 0.0;
      double st = 0.0;
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

          const int dp = std::abs(static_cast<int>(hb.pad) - static_cast<int>(ha.pad));
          const int dt = std::abs(static_cast<int>(hb.tbin) - static_cast<int>(ha.tbin));
          if (dp <= d->blob_dp && dt <= d->blob_dt)
          {
            used[j] = 1;
            q.push_back(j);
          }
        }
      }

      if (sw <= 0.0) continue;

      bl.layer = layer;
      bl.pad = sp / sw;
      bl.tbin = st / sw;
      bl.adc = sw;
      bl.nhits = nh;
      bl.used = 0;
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

      const double dp = bl.pad - pred_pad;
      const double dt = bl.tbin - pred_tbin;
      if (std::fabs(dp) > d->search_dp) continue;
      if (std::fabs(dt) > d->search_dt) continue;

      const double score = (dp * dp) / (d->search_dp * d->search_dp + 1.0e-9)
                         + (dt * dt) / (d->search_dt * d->search_dt + 1.0e-9)
                         - 0.01 * std::log(bl.adc + 1.0);
      if (score < best_score)
      {
        best_score = score;
        best = static_cast<int>(i);
      }
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
      if (direction < 0 && edge_layer <= 7) break;

      const unsigned int target_layer = static_cast<unsigned int>(static_cast<int>(edge_layer) + direction);

      double pred_pad = d->blobs[chain.back()].pad;
      double pred_tbin = d->blobs[chain.back()].tbin;
      if (direction < 0)
      {
        pred_pad = d->blobs[chain.front()].pad;
        pred_tbin = d->blobs[chain.front()].tbin;
      }

      if (chain.size() >= 2)
      {
        InModuleThreadData::Track tmp;
        if (fit_track_from_blobs(d->blobs, chain, d->weight_power, d->adc_weight_floor_frac, tmp))
        {
          pred_pad = tmp.pad_slope * static_cast<double>(target_layer) + tmp.pad_intercept;
          pred_tbin = tmp.tbin_slope * static_cast<double>(target_layer) + tmp.tbin_intercept;
        }
      }

      const int ibest = find_best_blob_on_layer(d, target_layer, pred_pad, pred_tbin);
      if (ibest < 0) break;

      d->blobs[ibest].used = 1;
      if (direction > 0)
      {
        chain.push_back(static_cast<unsigned int>(ibest));
      }
      else
      {
        chain.insert(chain.begin(), static_cast<unsigned int>(ibest));
      }
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
        {
          d->blobs[chain[k]].used = 0;
        }
        continue;
      }

      InModuleThreadData::Track trk;
      trk.track_id = tid;
      if (robust_fit_track_from_raw_hits(d->raw_hits, d->blobs, chain,
                                         d->weight_power, d->adc_weight_floor_frac, trk))
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
    build_blobs(d);
    build_tracks_linear(d);
    connect_track_pieces_in_module(d);

    if (d->verbosity > 1)
    {
      std::cout << "InModuleTracks worker: region=" << d->region
                << " sector=" << d->sector
                << " side=" << d->side
                << " raw_hits=" << d->raw_hits.size()
                << " blobs=" << d->blobs.size()
                << " tracks=" << d->tracks.size()
                << std::endl;
    }

    return 0;
  }
}

// ===================================================================
// Constructors
// ===================================================================
InModuleThreadData::LayerHitSet::LayerHitSet()
  : layer(0), hitsetkey(0), hitset(0) {}

InModuleThreadData::RawHit::RawHit()
  : layer(0), hitsetkey(0), hitkey(0), pad(0), tbin(0), adc(0), local_phi(0.0), local_radius(0.0) {}

InModuleThreadData::Blob::Blob()
  : layer(0), pad(0.0), tbin(0.0), adc(0.0), nhits(0), used(0) {}

InModuleThreadData::Track::Track()
  : track_id(0), first_layer(0), last_layer(0), nblobs(0), nrawhits(0),
    pad_slope(0.0), pad_intercept(0.0), tbin_slope(0.0), tbin_intercept(0.0),
    chi2_pad(0.0), chi2_tbin(0.0), ndof_pad(0), ndof_tbin(0) {}

InModuleThreadData::InModuleThreadData()
  : region(0), sector(0), side(0), module_key(0), tGeometry(0),
    pedestal(74.4), verbosity(0), blob_dt(2), blob_dp(2),
    search_dt(6), search_dp(6), min_track_blobs(4),
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
    m_event(0),
    m_maxThreads(72),
    m_pedestal(74.4),
    m_blob_dt(2),
    m_blob_dp(2),
    m_search_dt(6),
    m_search_dp(6),
    m_minTrackBlobs(4),
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
  m_tree->Branch("event", &m_tree_event, "event/I");

  m_tree->Branch("track_id", &m_tree_track_id);
  m_tree->Branch("region", &m_tree_region);
  m_tree->Branch("sector", &m_tree_sector);
  m_tree->Branch("side", &m_tree_side);
  m_tree->Branch("nblobs", &m_tree_nblobs);
  m_tree->Branch("nrawhits", &m_tree_nrawhits);
  m_tree->Branch("first_layer", &m_tree_first_layer);
  m_tree->Branch("last_layer", &m_tree_last_layer);

  m_tree->Branch("pad_slope", &m_tree_pad_slope);
  m_tree->Branch("pad_intercept", &m_tree_pad_intercept);
  m_tree->Branch("tbin_slope", &m_tree_tbin_slope);
  m_tree->Branch("tbin_intercept", &m_tree_tbin_intercept);
  m_tree->Branch("chi2_pad", &m_tree_chi2_pad);
  m_tree->Branch("chi2_tbin", &m_tree_chi2_tbin);
  m_tree->Branch("ndof_pad", &m_tree_ndof_pad);
  m_tree->Branch("ndof_tbin", &m_tree_ndof_tbin);

  m_tree->Branch("hit_event", &m_tree_hit_event);
  m_tree->Branch("hit_track_id", &m_tree_hit_track_id);
  m_tree->Branch("hit_region", &m_tree_hit_region);
  m_tree->Branch("hit_sector", &m_tree_hit_sector);
  m_tree->Branch("hit_side", &m_tree_hit_side);
  m_tree->Branch("hit_layer", &m_tree_hit_layer);
  m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey", &m_tree_hit_hitkey);
  m_tree->Branch("hit_pad", &m_tree_hit_pad);
  m_tree->Branch("hit_tbin", &m_tree_hit_tbin);
  m_tree->Branch("hit_adc", &m_tree_hit_adc);
  m_tree->Branch("hit_local_phi", &m_tree_hit_local_phi);
  m_tree->Branch("hit_local_radius", &m_tree_hit_local_radius);

  std::cout << Name() << "::Init - output file " << m_outputFileName << " created" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
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

  return Fun4AllReturnCodes::EVENT_OK;
}

void InModuleTracks::reset_tree_vars()
{
  m_tree_event = m_event;
  m_tree_hit_event.clear();
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

  m_tree_hit_track_id.clear();
  m_tree_hit_region.clear();
  m_tree_hit_sector.clear();
  m_tree_hit_side.clear();
  m_tree_hit_layer.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
  m_tree_hit_pad.clear();
  m_tree_hit_tbin.clear();
  m_tree_hit_adc.clear();
  m_tree_hit_local_phi.clear();
m_tree_hit_local_radius.clear();
}

int InModuleTracks::process_event(PHCompositeNode*)
{
  reset_tree_vars();

  std::vector<InModuleThreadData> tdata;
  tdata.reserve(72);

  // One work package per module, not per layer.
  for (unsigned int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int region = 0; region < 3; ++region)
      {
        InModuleThreadData td;
        td.region = region;
        td.sector = sector;
        td.side = static_cast<int>(side);
        td.module_key = TpcDefs::genModuleHitSetKey(static_cast<uint8_t>(region),
                                                    static_cast<uint8_t>(sector),
                                                    static_cast<uint8_t>(side));
        td.tGeometry = m_tGeometry;
        td.pedestal = m_pedestal;
        td.verbosity = Verbosity();
        td.blob_dt = m_blob_dt;
        td.blob_dp = m_blob_dp;
        td.search_dt = m_search_dt;
        td.search_dp = m_search_dp;
        td.min_track_blobs = m_minTrackBlobs;

        for (unsigned int l = 0; l < 16; ++l)
        {
          const unsigned int layer = region * 16 + l + 7;
          const TrkrDefs::hitsetkey hitset_key = TpcDefs::genHitSetKey(layer, sector, side);
          TrkrHitSet* hitset = m_hits->findHitSet(hitset_key);
          if (!hitset) continue;

          InModuleThreadData::LayerHitSet lhs;
          lhs.layer = layer;
          lhs.hitsetkey = hitset_key;
          lhs.hitset = hitset;
          td.layer_hitsets.push_back(lhs);
        }

        if (!td.layer_hitsets.empty())
        {
          tdata.push_back(td);
        }
      }
    }
  }

  std::cout << Name() << "::process_event - event " << m_event
            << " has " << tdata.size() << " non-empty modules" << std::endl;

  const unsigned int maxLive = std::max(1u, std::min(m_maxThreads,
                         static_cast<unsigned int>(tdata.size())));

  for (unsigned int start = 0; start < static_cast<unsigned int>(tdata.size()); start += maxLive)
  {
    const unsigned int end = std::min(start + maxLive,
                                      static_cast<unsigned int>(tdata.size()));
    const unsigned int nLive = end - start;

    std::vector<pthread_t> threads(nLive);
    std::vector<int> thread_ok(nLive, 0);

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
                  << " side=" << tdata[idx].side << std::endl;
      }
      else
      {
        thread_ok[i] = 1;
      }
    }

    for (unsigned int i = 0; i < nLive; ++i)
    {
      if (thread_ok[i]) pthread_join(threads[i], 0);
    }
  }

  for (unsigned int im = 0; im < tdata.size(); ++im)
  {
    const InModuleThreadData& td = tdata[im];
    for (unsigned int it = 0; it < td.tracks.size(); ++it)
    {
      const InModuleThreadData::Track& tr = td.tracks[it];
      const unsigned int global_track_id = static_cast<unsigned int>(m_tree_track_id.size());

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

      for (unsigned int ih = 0; ih < tr.raw_hit_indices.size(); ++ih)
      {
        const InModuleThreadData::RawHit& rh = td.raw_hits[tr.raw_hit_indices[ih]];
        m_tree_hit_event.push_back(m_event);
        m_tree_hit_track_id.push_back(global_track_id);
        m_tree_hit_region.push_back(td.region);
        m_tree_hit_sector.push_back(td.sector);
        m_tree_hit_side.push_back(td.side);
        m_tree_hit_layer.push_back(rh.layer);
        m_tree_hit_hitsetkey.push_back(static_cast<unsigned long long>(rh.hitsetkey));
        m_tree_hit_hitkey.push_back(static_cast<unsigned long long>(rh.hitkey));
        m_tree_hit_pad.push_back(static_cast<double>(rh.pad));
        m_tree_hit_tbin.push_back(static_cast<double>(rh.tbin));
        m_tree_hit_adc.push_back(static_cast<double>(rh.adc));
        m_tree_hit_local_phi.push_back(rh.local_phi);
        m_tree_hit_local_radius.push_back(rh.local_radius);
      }
    }
  }

  if (m_tree) m_tree->Fill();

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " tracks=" << m_tree_track_id.size()
              << " track-raw-hits=" << m_tree_hit_track_id.size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
