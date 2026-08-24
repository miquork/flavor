// Joint four-flavor residual-response proof of principle for Run 3 AK4 PUPPI.
//
// Run with, for example:
//   root -l -b -q 'JetFlavor.C(2024,true)'
//
// The fit uses each tagged balance divided by the inclusive balance in the same
// channel.  This removes channel-wide Z/gamma calibration modes before fitting
// the flavor contrast.  Four residual responses (uds, g, c, b) are fitted
// simultaneously and centered on an explicit 100 GeV reference-mixture
// convention.  The tiny unmatched component is monitored and removed by the
// matched-component normalization.  This first implementation propagates the
// stored balance errors with their shared inclusive denominator and a clearly
// labeled diagnostic floor; tag-SF, background, and other correlated
// systematics are intentionally not yet included.

#include "tdrstyle_mod22.C"

#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace jointflavor {

constexpr int kNTags = 4;
constexpr int kNFitFlavors = 4;
constexpr int kNFlavors = 5;
constexpr int kNPar = 8;

const char *kTags[kNTags] = {"q", "g", "c", "b"};
const char *kFlavors[kNFlavors] = {"q", "g", "c", "b", "n"};
const char *kFlavorLabels[kNFitFlavors] = {"uds", "g", "c", "b"};
const int kColors[kNFitFlavors] = {kBlue + 1, kRed + 1, kGreen + 2,
                                    kMagenta + 1};

struct ResponseBin {
  int sample = -1;  // 0: Z+jet, 1: photon+jet
  double pt = 0;
  double ptLow = 0;
  double ptHigh = 0;
  double fraction[kNFlavors] = {0};
  double unmatchedFraction = 0;
  double inclusiveY = 0;
  double inclusiveEy = 0;
  double inclusiveAlpha[kNFlavors] = {0};
  double y[kNTags] = {0};
  double ey[kNTags] = {0};
  double alpha[kNTags][kNFlavors] = {{0}};
  bool valid[kNTags] = {false};
  bool inclusiveValid = false;
};

struct FitRow {
  int sample = -1;
  int tag = -1;
  double pt = 0;
  double y = 0;
  double ey = 0;
  double base = 0;
  double a[kNFitFlavors] = {0};
};

struct FlavorPoint {
  int sample = -1;
  int flavor = -1;
  double pt = 0;
  double ept = 0;
  double value = 0;
  double error = 0;
  double condition = 0;
};

struct FitResult {
  TVectorD p;
  TMatrixD cov;
  double chi2 = 0;
  double chi2Sample[2] = {0, 0};
  int ndf = 0;
  int rowsSample[2] = {0, 0};
  double covarianceScale = 1;
  double condition = 0;
  bool ok = false;
  FitResult() : p(kNPar), cov(kNPar, kNPar) {}
};

TH1D *getHistogram(TDirectory *dir, const TString &name, bool required = true) {
  TH1D *h = dir ? dynamic_cast<TH1D *>(dir->Get(name)) : nullptr;
  if (!h && required)
    std::cerr << "Missing histogram " << name << std::endl;
  return h;
}

double safeContent(TH1D *h, double pt) {
  if (!h)
    return 0;
  const int bin = h->GetXaxis()->FindFixBin(pt);
  if (bin < 1 || bin > h->GetNbinsX())
    return 0;
  const double value = h->GetBinContent(bin);
  return std::isfinite(value) ? value : 0;
}

// Minimize ||E*x-target||^2 with sum(x)=sumTarget and x>=0.  There are only
// four fitted fractions, so enumerating all active sets is both transparent and
// robust against negative fractions from a direct matrix inversion.
bool solveNonNegativeFractions(const double efficiency[kNTags][kNFlavors],
                               const double dataTagFraction[kNTags],
                               double unmatchedFraction,
                               double output[kNFitFlavors]) {
  double target[kNTags] = {0};
  for (int t = 0; t < kNTags; ++t)
    target[t] = dataTagFraction[t] -
                efficiency[t][kNFlavors - 1] * unmatchedFraction;

  const double sumTarget = std::max(0.0, 1.0 - unmatchedFraction);
  double best = std::numeric_limits<double>::max();
  bool found = false;
  double bestX[kNFitFlavors] = {0};

  for (int mask = 1; mask < (1 << kNFitFlavors); ++mask) {
    std::vector<int> active;
    for (int f = 0; f < kNFitFlavors; ++f)
      if (mask & (1 << f))
        active.push_back(f);

    const int n = static_cast<int>(active.size());
    TMatrixD kkt(n + 1, n + 1);
    TVectorD rhs(n + 1);
    kkt.Zero();
    rhs.Zero();

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double value = 0;
        for (int t = 0; t < kNTags; ++t)
          value += 2.0 * efficiency[t][active[i]] *
                   efficiency[t][active[j]];
        kkt(i, j) = value;
      }
      for (int t = 0; t < kNTags; ++t)
        rhs(i) += 2.0 * efficiency[t][active[i]] * target[t];
      kkt(i, n) = 1.0;
      kkt(n, i) = 1.0;
    }
    rhs(n) = sumTarget;

    Bool_t solved = false;
    TDecompSVD decomposition(kkt);
    TVectorD solution = decomposition.Solve(rhs, solved);
    if (!solved)
      continue;

    double candidate[kNFitFlavors] = {0};
    bool nonNegative = true;
    for (int i = 0; i < n; ++i) {
      candidate[active[i]] = solution(i);
      if (!std::isfinite(solution(i)) || solution(i) < -1e-9)
        nonNegative = false;
    }
    if (!nonNegative)
      continue;

    double objective = 0;
    for (int t = 0; t < kNTags; ++t) {
      double prediction = 0;
      for (int f = 0; f < kNFitFlavors; ++f)
        prediction += efficiency[t][f] * candidate[f];
      const double residual = prediction - target[t];
      objective += residual * residual;
    }
    if (objective < best) {
      best = objective;
      found = true;
      for (int f = 0; f < kNFitFlavors; ++f)
        bestX[f] = std::max(0.0, candidate[f]);
    }
  }

  if (!found)
    return false;
  double norm = 0;
  for (double value : bestX)
    norm += value;
  if (!(norm > 0))
    return false;
  for (int f = 0; f < kNFitFlavors; ++f)
    output[f] = bestX[f] * sumTarget / norm;
  return true;
}

bool buildResponseBin(TDirectory *mc, TDirectory *data, TDirectory *ratio,
                      const char *sampleName, int sample, int ratioBin,
                      bool useDataFractions, ResponseBin &result) {
  TH1D *reference = getHistogram(
      ratio, Form("hdm_mpfchs1_%s%s", sampleName, kTags[0]));
  if (!reference || ratioBin < 1 || ratioBin > reference->GetNbinsX())
    return false;

  result.sample = sample;
  result.pt = reference->GetBinCenter(ratioBin);
  result.ptLow = reference->GetBinLowEdge(ratioBin);
  result.ptHigh = reference->GetBinLowEdge(ratioBin + 1);

  const double minPt = (sample == 0 ? 45.0 : 40.0);
  const double maxPt = (sample == 0 ? 300.0 : 850.0);
  if (!(result.pt >= minPt && result.pt <= maxPt))
    return false;

  double inclusiveCounts[kNFlavors] = {0};
  double totalInclusive = 0;
  for (int f = 0; f < kNFlavors; ++f) {
    TH1D *h = getHistogram(
        mc, Form("counts_%si%s_a100", sampleName, kFlavors[f]));
    inclusiveCounts[f] = std::max(0.0, safeContent(h, result.pt));
    totalInclusive += inclusiveCounts[f];
  }
  if (!(totalInclusive > 0))
    return false;

  double mcFraction[kNFlavors] = {0};
  for (int f = 0; f < kNFlavors; ++f)
    mcFraction[f] = inclusiveCounts[f] / totalInclusive;

  double efficiency[kNTags][kNFlavors] = {{0}};
  for (int f = 0; f < kNFlavors; ++f) {
    double tagSum = 0;
    double tagged[kNTags] = {0};
    for (int t = 0; t < kNTags; ++t) {
      TH1D *h = getHistogram(
          mc, Form("counts_%s%s%s_a100", sampleName, kTags[t],
                   kFlavors[f]));
      tagged[t] = std::max(0.0, safeContent(h, result.pt));
      tagSum += tagged[t];
    }
    if (!(tagSum > 0))
      return false;
    for (int t = 0; t < kNTags; ++t)
      efficiency[t][f] = tagged[t] / tagSum;
  }

  double dataTagCounts[kNTags] = {0};
  double dataTagFraction[kNTags] = {0};
  double dataTotal = 0;
  for (int t = 0; t < kNTags; ++t) {
    TH1D *h =
        getHistogram(data, Form("counts_%s%s_a100", sampleName, kTags[t]));
    dataTagCounts[t] = std::max(0.0, safeContent(h, result.pt));
    dataTotal += dataTagCounts[t];
  }
  if (!(dataTotal > 0))
    return false;
  for (int t = 0; t < kNTags; ++t)
    dataTagFraction[t] = dataTagCounts[t] / dataTotal;

  for (int f = 0; f < kNFlavors; ++f)
    result.fraction[f] = mcFraction[f];
  result.unmatchedFraction = mcFraction[kNFlavors - 1];

  if (useDataFractions) {
    double fittedFraction[kNFitFlavors] = {0};
    const bool solved = solveNonNegativeFractions(
        efficiency, dataTagFraction, result.unmatchedFraction, fittedFraction);
    if (solved)
      for (int f = 0; f < kNFitFlavors; ++f)
        result.fraction[f] = fittedFraction[f];
  }

  TH1D *hInclusiveRatio =
      getHistogram(ratio, Form("hdm_mpfchs1_%si", sampleName));
  TH1D *hInclusiveMc =
      getHistogram(mc, Form("hdm_mpfchs1_%si", sampleName));
  const int inclusiveBin =
      hInclusiveRatio ? hInclusiveRatio->FindBin(result.pt) : -1;
  result.inclusiveY =
      hInclusiveRatio ? hInclusiveRatio->GetBinContent(inclusiveBin) : 0;
  result.inclusiveEy =
      hInclusiveRatio ? hInclusiveRatio->GetBinError(inclusiveBin) : 0;
  const double inclusiveMcResponse = safeContent(hInclusiveMc, result.pt);
  double inclusiveWeight[kNFlavors] = {0};
  double inclusiveWeightSum = 0;
  for (int f = 0; f < kNFlavors; ++f) {
    TH1D *component = getHistogram(
        mc, Form("hdm_mpfchs1_%si%s", sampleName, kFlavors[f]));
    double response = safeContent(component, result.pt);
    if (!(response > 0) && f == kNFlavors - 1)
      response = inclusiveMcResponse;
    inclusiveWeight[f] = result.fraction[f] * response;
    inclusiveWeightSum += inclusiveWeight[f];
  }
  result.inclusiveValid =
      result.inclusiveY > 0 && result.inclusiveEy > 0 &&
      inclusiveMcResponse > 0 && inclusiveWeightSum > 0 &&
      std::isfinite(result.inclusiveY) && std::isfinite(result.inclusiveEy);
  if (inclusiveWeightSum > 0)
    for (int f = 0; f < kNFlavors; ++f)
      result.inclusiveAlpha[f] = inclusiveWeight[f] / inclusiveWeightSum;

  for (int t = 0; t < kNTags; ++t) {
    TH1D *hRatio = getHistogram(
        ratio, Form("hdm_mpfchs1_%s%s", sampleName, kTags[t]));
    TH1D *hMcTag =
        getHistogram(mc, Form("hdm_mpfchs1_%s%s", sampleName, kTags[t]));
    const int bRatio = hRatio ? hRatio->FindBin(result.pt) : -1;
    const double y = hRatio ? hRatio->GetBinContent(bRatio) : 0;
    const double ey = hRatio ? hRatio->GetBinError(bRatio) : 0;
    const double rTag = safeContent(hMcTag, result.pt);

    double responseWeight[kNFlavors] = {0};
    double responseWeightSum = 0;
    bool valid = y > 0 && ey > 0 && rTag > 0 && result.inclusiveValid &&
                 std::isfinite(y) && std::isfinite(ey);
    for (int f = 0; f < kNFlavors; ++f) {
      TH1D *hMcComponent = getHistogram(
          mc, Form("hdm_mpfchs1_%s%s%s", sampleName, kTags[t],
                   kFlavors[f]));
      double rComponent = safeContent(hMcComponent, result.pt);
      if (!(rComponent > 0) && f == kNFlavors - 1)
        rComponent = rTag;
      responseWeight[f] =
          efficiency[t][f] * result.fraction[f] * rComponent;
      responseWeightSum += responseWeight[f];
      valid = valid && rComponent > 0 && std::isfinite(responseWeight[f]);
    }
    valid = valid && responseWeightSum > 0;
    if (responseWeightSum > 0)
      for (int f = 0; f < kNFlavors; ++f)
        result.alpha[t][f] = responseWeight[f] / responseWeightSum;
    result.y[t] = y;
    result.ey[t] = ey;
    result.valid[t] = valid;
  }
  return true;
}

std::vector<ResponseBin> loadInputs(TFile &input, bool useDataFractions) {
  TDirectory *mc = dynamic_cast<TDirectory *>(input.Get("mc/eta00-13"));
  TDirectory *data = dynamic_cast<TDirectory *>(input.Get("data/eta00-13"));
  TDirectory *ratio = dynamic_cast<TDirectory *>(input.Get("ratio/eta00-13"));
  if (!mc || !data || !ratio)
    return {};

  std::vector<ResponseBin> bins;
  const char *sampleNames[2] = {"z", "g"};
  for (int sample = 0; sample < 2; ++sample) {
    TH1D *reference = getHistogram(
        ratio, Form("hdm_mpfchs1_%s%s", sampleNames[sample], kTags[0]));
    if (!reference)
      continue;
    for (int bin = 1; bin <= reference->GetNbinsX(); ++bin) {
      ResponseBin inputBin;
      if (buildResponseBin(mc, data, ratio, sampleNames[sample], sample, bin,
                           useDataFractions, inputBin))
        bins.push_back(inputBin);
    }
  }
  return bins;
}

std::vector<FitRow> makeRows(const std::vector<ResponseBin> &bins) {
  std::vector<FitRow> rows;
  for (const ResponseBin &bin : bins) {
    for (int t = 0; t < kNTags; ++t) {
      if (!bin.valid[t])
        continue;
      FitRow row;
      row.sample = bin.sample;
      row.tag = t;
      row.pt = bin.pt;
      row.y = bin.y[t] / bin.inclusiveY;
      row.ey = row.y * std::sqrt(
                           std::pow(bin.ey[t] / bin.y[t], 2) +
                           std::pow(bin.inclusiveEy / bin.inclusiveY, 2));
      row.base = 1;
      for (int f = 0; f < kNFitFlavors; ++f)
        row.a[f] = bin.alpha[t][f] - bin.inclusiveAlpha[f];
      rows.push_back(row);
    }
  }
  return rows;
}

void fillDesign(const FitRow &row, double design[kNPar]) {
  const double shape = std::pow(row.pt / 100.0, -0.3) - 1.0;
  for (int f = 0; f < kNFitFlavors; ++f) {
    design[2 * f] = row.a[f];
    design[2 * f + 1] = row.a[f] * shape;
  }
}

double predictionPercent(const FitRow &row, const TVectorD &parameters) {
  double design[kNPar] = {0};
  fillDesign(row, design);
  double shift = 0;
  for (int p = 0; p < kNPar; ++p)
    shift += design[p] * parameters(p);
  return shift;
}

std::array<double, kNFitFlavors>
makeAnchorWeights(const std::vector<ResponseBin> &nominalBins) {
  std::array<double, kNFitFlavors> weights = {0, 0, 0, 0};
  for (int sample = 0; sample < 2; ++sample) {
    const ResponseBin *closest = nullptr;
    double distance = std::numeric_limits<double>::max();
    for (const ResponseBin &bin : nominalBins) {
      if (bin.sample != sample || !bin.inclusiveValid)
        continue;
      if (std::abs(bin.pt - 100.0) < distance) {
        closest = &bin;
        distance = std::abs(bin.pt - 100.0);
      }
    }
    if (closest)
      for (int f = 0; f < kNFitFlavors; ++f)
        weights[f] += 0.5 * closest->inclusiveAlpha[f];
  }
  double sum = 0;
  for (double weight : weights)
    sum += weight;
  if (sum > 0)
    for (double &weight : weights)
      weight /= sum;
  return weights;
}

bool observationBlock(const ResponseBin &bin, double diagnosticFloorPercent,
                      TVectorD &target, TMatrixD &covariance) {
  if (!bin.inclusiveValid)
    return false;
  for (int t = 0; t < kNTags; ++t)
    if (!bin.valid[t])
      return false;

  target.ResizeTo(kNTags);
  covariance.ResizeTo(kNTags, kNTags);
  target.Zero();
  covariance.Zero();
  double doubleRatio[kNTags] = {0};
  for (int t = 0; t < kNTags; ++t) {
    doubleRatio[t] = bin.y[t] / bin.inclusiveY;
    target(t) = 100.0 * std::log(doubleRatio[t]);
  }

  const double inclusiveRelativeVariance =
      std::pow(bin.inclusiveEy / bin.inclusiveY, 2);
  for (int t = 0; t < kNTags; ++t) {
    for (int u = 0; u < kNTags; ++u) {
      // All tagged double ratios in a native bin share the inclusive
      // denominator.  The tag/inclusive overlap covariance is unavailable in
      // the current ROOT schema; an event-level bootstrap is required later.
      covariance(t, u) = 10000.0 * inclusiveRelativeVariance;
      if (t == u) {
        covariance(t, u) += 10000.0 * std::pow(bin.ey[t] / bin.y[t], 2);
        covariance(t, u) += diagnosticFloorPercent * diagnosticFloorPercent;
      }
    }
  }
  return true;
}

FitResult fitBins(const std::vector<ResponseBin> &bins,
                  const std::array<double, kNFitFlavors> &anchor,
                  double diagnosticFloorPercent) {
  FitResult result;
  TMatrixD normal(kNPar, kNPar);
  TVectorD rhs(kNPar);
  normal.Zero();
  rhs.Zero();

  for (const ResponseBin &bin : bins) {
    if (!(bin.pt >= 90.0 && bin.pt <= 300.0))
      continue;
    TVectorD target;
    TMatrixD covariance;
    if (!observationBlock(bin, diagnosticFloorPercent, target, covariance))
      continue;
    Bool_t inverted = false;
    TDecompSVD covarianceDecomposition(covariance);
    const TMatrixD inverse = covarianceDecomposition.Invert(inverted);
    if (!inverted)
      continue;

    double design[kNTags][kNPar] = {{0}};
    const double shape = std::pow(bin.pt / 100.0, -0.3) - 1.0;
    for (int t = 0; t < kNTags; ++t) {
      double contrast[kNFitFlavors] = {0};
      for (int f = 0; f < kNFitFlavors; ++f)
        contrast[f] = bin.alpha[t][f] - bin.inclusiveAlpha[f];
      for (int f = 0; f < kNFitFlavors; ++f) {
        design[t][2 * f] = contrast[f];
        design[t][2 * f + 1] = contrast[f] * shape;
      }
    }
    for (int p = 0; p < kNPar; ++p) {
      for (int t = 0; t < kNTags; ++t)
        for (int u = 0; u < kNTags; ++u)
          rhs(p) += design[t][p] * inverse(t, u) * target(u);
      for (int q = 0; q < kNPar; ++q)
        for (int t = 0; t < kNTags; ++t)
          for (int u = 0; u < kNTags; ++u)
            normal(p, q) += design[t][p] * inverse(t, u) * design[u][q];
    }
    result.rowsSample[bin.sample] += kNTags;
  }

  const int observations = result.rowsSample[0] + result.rowsSample[1];
  if (observations <= kNPar - 2)
    return result;

  // The double ratios determine only flavor contrasts.  Center both common
  // modes on an explicit equal-Z/gamma inclusive MC mixture, so its weighted
  // residual response is zero at every pT in this relative-response model.
  double interceptConstraint[kNPar] = {anchor[0], 0, anchor[1], 0,
                                       anchor[2], 0, anchor[3], 0};
  double slopeConstraint[kNPar] = {0, anchor[0], 0, anchor[1],
                                   0, anchor[2], 0, anchor[3]};
  constexpr int kNConstraints = 2;
  TMatrixD kkt(kNPar + kNConstraints, kNPar + kNConstraints);
  TVectorD kktRhs(kNPar + kNConstraints);
  kkt.Zero();
  kktRhs.Zero();
  for (int p = 0; p < kNPar; ++p) {
    kktRhs(p) = rhs(p);
    for (int q = 0; q < kNPar; ++q)
      kkt(p, q) = normal(p, q);
    kkt(p, kNPar) = interceptConstraint[p];
    kkt(kNPar, p) = interceptConstraint[p];
    kkt(p, kNPar + 1) = slopeConstraint[p];
    kkt(kNPar + 1, p) = slopeConstraint[p];
  }

  Bool_t solved = false;
  TDecompSVD decomposition(kkt);
  const TVectorD solution = decomposition.Solve(kktRhs, solved);
  Bool_t inverted = false;
  const TMatrixD inverseKkt = decomposition.Invert(inverted);
  if (!solved || !inverted)
    return result;
  for (int p = 0; p < kNPar; ++p) {
    result.p(p) = solution(p);
    for (int q = 0; q < kNPar; ++q)
      result.cov(p, q) = inverseKkt(p, q);
  }
  const TVectorD singular = decomposition.GetSig();
  if (singular.GetNrows() > 0 && singular(singular.GetNrows() - 1) > 0)
    result.condition = singular(0) / singular(singular.GetNrows() - 1);

  for (const ResponseBin &bin : bins) {
    if (!(bin.pt >= 90.0 && bin.pt <= 300.0))
      continue;
    TVectorD target;
    TMatrixD covariance;
    if (!observationBlock(bin, diagnosticFloorPercent, target, covariance))
      continue;
    Bool_t covarianceOk = false;
    TDecompSVD covarianceDecomposition(covariance);
    const TMatrixD inverse =
        covarianceDecomposition.Invert(covarianceOk);
    if (!covarianceOk)
      continue;
    TVectorD residual(kNTags);
    const double shape = std::pow(bin.pt / 100.0, -0.3) - 1.0;
    for (int t = 0; t < kNTags; ++t) {
      double prediction = 0;
      const double contrast[kNFitFlavors] = {
          bin.alpha[t][0] - bin.inclusiveAlpha[0],
          bin.alpha[t][1] - bin.inclusiveAlpha[1],
          bin.alpha[t][2] - bin.inclusiveAlpha[2],
          bin.alpha[t][3] - bin.inclusiveAlpha[3]};
      for (int f = 0; f < kNFitFlavors; ++f)
        prediction +=
            contrast[f] * (result.p(2 * f) + result.p(2 * f + 1) * shape);
      residual(t) = target(t) - prediction;
    }
    double chi2 = 0;
    for (int t = 0; t < kNTags; ++t)
      for (int u = 0; u < kNTags; ++u)
        chi2 += residual(t) * inverse(t, u) * residual(u);
    result.chi2 += chi2;
    result.chi2Sample[bin.sample] += chi2;
  }
  result.ndf = observations - (kNPar - kNConstraints);
  result.covarianceScale =
      result.ndf > 0 ? std::max(1.0, result.chi2 / result.ndf) : 1.0;
  result.cov *= result.covarianceScale;
  result.ok = true;
  return result;
}

double flavorValue(int flavor, double pt, const TVectorD &p) {
  return p(2 * flavor) +
         p(2 * flavor + 1) * (std::pow(pt / 100.0, -0.3) - 1.0);
}

void flavorVector(int flavor, double pt, double vector[kNPar]) {
  for (int p = 0; p < kNPar; ++p)
    vector[p] = 0;
  vector[2 * flavor] = 1;
  vector[2 * flavor + 1] = std::pow(pt / 100.0, -0.3) - 1.0;
}

double flavorError(int flavor, double pt, const TMatrixD &covariance) {
  double vector[kNPar] = {0};
  flavorVector(flavor, pt, vector);
  double variance = 0;
  for (int i = 0; i < kNPar; ++i)
    for (int j = 0; j < kNPar; ++j)
      variance += vector[i] * covariance(i, j) * vector[j];
  return variance > 0 ? std::sqrt(variance) : 0;
}

std::vector<FlavorPoint>
makeFlavorPoints(const std::vector<ResponseBin> &bins,
                 const std::array<double, kNFitFlavors> &anchor,
                 double diagnosticFloorPercent, double covarianceScale) {
  std::vector<FlavorPoint> points;
  for (const ResponseBin &bin : bins) {
    TVectorD target;
    TMatrixD covariance;
    if (!observationBlock(bin, diagnosticFloorPercent, target, covariance))
      continue;
    Bool_t covarianceOk = false;
    TDecompSVD covarianceDecomposition(covariance);
    const TMatrixD covarianceInverse =
        covarianceDecomposition.Invert(covarianceOk);
    if (!covarianceOk)
      continue;
    TMatrixD matrix(kNTags, kNFitFlavors);
    for (int t = 0; t < kNTags; ++t)
      for (int f = 0; f < kNFitFlavors; ++f)
        matrix(t, f) = bin.alpha[t][f] - bin.inclusiveAlpha[f];
    TMatrixD normal(kNFitFlavors, kNFitFlavors);
    TVectorD rhs(kNFitFlavors);
    normal.Zero();
    rhs.Zero();
    for (int f = 0; f < kNFitFlavors; ++f) {
      for (int t = 0; t < kNTags; ++t)
        for (int u = 0; u < kNTags; ++u)
          rhs(f) += matrix(t, f) * covarianceInverse(t, u) * target(u);
      for (int g = 0; g < kNFitFlavors; ++g)
        for (int t = 0; t < kNTags; ++t)
          for (int u = 0; u < kNTags; ++u)
            normal(f, g) +=
                matrix(t, f) * covarianceInverse(t, u) * matrix(u, g);
    }
    TMatrixD kkt(kNFitFlavors + 1, kNFitFlavors + 1);
    TVectorD kktRhs(kNFitFlavors + 1);
    kkt.Zero();
    kktRhs.Zero();
    for (int f = 0; f < kNFitFlavors; ++f) {
      kktRhs(f) = rhs(f);
      for (int g = 0; g < kNFitFlavors; ++g)
        kkt(f, g) = normal(f, g);
      kkt(f, kNFitFlavors) = anchor[f];
      kkt(kNFitFlavors, f) = anchor[f];
    }
    TDecompSVD decomposition(kkt);
    Bool_t solved = false;
    const TVectorD solution = decomposition.Solve(kktRhs, solved);
    Bool_t inverted = false;
    const TMatrixD inverse = decomposition.Invert(inverted);
    if (!solved || !inverted)
      continue;
    const TVectorD singular = decomposition.GetSig();
    double condition = 0;
    if (singular.GetNrows() > 0 && singular(singular.GetNrows() - 1) > 0)
      condition = singular(0) / singular(singular.GetNrows() - 1);

    for (int f = 0; f < kNFitFlavors; ++f) {
      const double error = inverse(f, f) > 0
                               ? std::sqrt(inverse(f, f) * covarianceScale)
                               : 0;
      if (!std::isfinite(solution(f)) || !std::isfinite(error))
        continue;
      FlavorPoint point;
      point.sample = bin.sample;
      point.flavor = f;
      point.pt = bin.pt;
      point.ept = 0.5 * (bin.ptHigh - bin.ptLow);
      point.value = solution(f);
      point.error = error;
      point.condition = condition;
      points.push_back(point);
    }
  }
  return points;
}

TGraphErrors *makeCurve(int flavor, const FitResult &fit, double xmin,
                        double xmax, int points, bool errors) {
  TGraphErrors *graph = new TGraphErrors(points);
  const double logMin = std::log(xmin);
  const double logMax = std::log(xmax);
  for (int i = 0; i < points; ++i) {
    const double fraction = points > 1 ? double(i) / (points - 1) : 0;
    const double pt = std::exp(logMin + fraction * (logMax - logMin));
    graph->SetPoint(i, pt, flavorValue(flavor, pt, fit.p));
    graph->SetPointError(
        i, 0, errors ? flavorError(flavor, pt, fit.cov) : 0);
  }
  return graph;
}

void writeOutputs(int year, bool useDataFractions,
                  const std::vector<ResponseBin> &bins,
                  const std::vector<FitRow> &rows,
                  const std::vector<FlavorPoint> &points,
                  const FitResult &fit,
                  const std::array<double, kNFitFlavors> &anchor,
                  double diagnosticFloorPercent,
                  const std::vector<TGraphErrors *> &bands,
                  const std::vector<TGraphErrors *> &curves,
                  const std::vector<std::vector<TGraphErrors *>> &dataGraphs) {
  gSystem->mkdir("results", true);
  const TString variant = useDataFractions ? "_dataFractions" : "";
  const TString stem = Form("results/JetFlavorJointFit_%d%s", year,
                            variant.Data());
  TFile output(stem + ".root", "RECREATE");

  TH1D parameters("fit_parameters", ";parameter;value (%)", kNPar, 0, kNPar);
  const char *parameterLabels[kNPar] = {
      "uds@100", "uds slope", "g@100", "g slope",
      "c@100",   "c slope",   "b@100", "b slope"};
  for (int p = 0; p < kNPar; ++p) {
    parameters.GetXaxis()->SetBinLabel(p + 1, parameterLabels[p]);
    parameters.SetBinContent(p + 1, fit.p(p));
    parameters.SetBinError(p + 1, std::sqrt(std::max(0.0, fit.cov(p, p))));
  }
  parameters.Write();

  TH2D covariance("fit_covariance", ";parameter;parameter", kNPar, 0, kNPar,
                  kNPar, 0, kNPar);
  for (int i = 0; i < kNPar; ++i) {
    covariance.GetXaxis()->SetBinLabel(i + 1, parameterLabels[i]);
    covariance.GetYaxis()->SetBinLabel(i + 1, parameterLabels[i]);
    for (int j = 0; j < kNPar; ++j)
      covariance.SetBinContent(i + 1, j + 1, fit.cov(i, j));
  }
  covariance.Write();

  TH1D anchorHistogram("reference_anchor", ";flavor;response share", 4, 0, 4);
  for (int f = 0; f < kNFitFlavors; ++f) {
    anchorHistogram.GetXaxis()->SetBinLabel(f + 1, kFlavorLabels[f]);
    anchorHistogram.SetBinContent(f + 1, anchor[f]);
  }
  anchorHistogram.Write();

  for (int f = 0; f < kNFitFlavors; ++f) {
    bands[f]->SetName(Form("fit_band_%s", kFlavorLabels[f]));
    curves[f]->SetName(Form("fit_curve_%s", kFlavorLabels[f]));
    bands[f]->Write();
    curves[f]->Write();
    for (int sample = 0; sample < 2; ++sample) {
      dataGraphs[sample][f]->SetName(
          Form("points_%s_%s", sample == 0 ? "zjet" : "gamjet",
               kFlavorLabels[f]));
      dataGraphs[sample][f]->Write();
    }
  }

  TTree rowTree("fit_rows", "Tagged balance inputs to the joint fit");
  int sample = 0, tag = 0;
  double pt = 0, ratio = 0, ratioError = 0, base = 0;
  double coefficients[kNFitFlavors] = {0};
  rowTree.Branch("sample", &sample, "sample/I");
  rowTree.Branch("tag", &tag, "tag/I");
  rowTree.Branch("pt", &pt, "pt/D");
  rowTree.Branch("ratio", &ratio, "ratio/D");
  rowTree.Branch("ratioError", &ratioError, "ratioError/D");
  rowTree.Branch("base", &base, "base/D");
  rowTree.Branch("coefficients", coefficients,
                 Form("coefficients[%d]/D", kNFitFlavors));
  for (const FitRow &row : rows) {
    sample = row.sample;
    tag = row.tag;
    pt = row.pt;
    ratio = row.y;
    ratioError = row.ey;
    base = row.base;
    for (int f = 0; f < kNFitFlavors; ++f)
      coefficients[f] = row.a[f];
    rowTree.Fill();
  }
  rowTree.Write();

  TNamed metadata(
      "fit_scope",
      Form("year=%d; eta=0-1.3; jet=AK4PUPPI; samples=Zjet,gamjet; "
           "fitRangeGeV=90-300; dataFractions=%s; "
           "unmatched=explicit_fixed_unit_response; "
           "diagnosticFloorPercent=%.3f; anchor=equal_channel_inclusive_MC_100GeV; "
           "uncertainty=stored_response_errors_shared_inclusive_scaled_by_chi2ndf",
           year, useDataFractions ? "nonnegative_count_fit" : "MC",
           diagnosticFloorPercent));
  metadata.Write();
  output.Close();

  std::ofstream summary((stem + ".txt").Data());
  summary << "Joint Z+jet + photon+jet four-flavor response fit\n"
          << "year: " << year << "\n"
          << "jets: AK4 PUPPI, |eta| < 1.3\n"
          << "composition: "
          << (useDataFractions ? "non-negative data count fit"
                               : "nominal MC fractions")
          << "\n"
          << "observable: tagged balance divided by same-channel inclusive balance\n"
          << "fit range: 90-300 GeV\n"
          << "unmatched: explicit in response shares, residual response fixed to 1\n"
          << "diagnostic diagonal floor: " << diagnosticFloorPercent << " %\n"
          << "anchor weights (uds,g,c,b): " << anchor[0] << ", " << anchor[1]
          << ", " << anchor[2] << ", " << anchor[3] << "\n"
          << "anchor convention: weighted intercept and slope are zero (reference "
             "mixture preserved at every pT)\n"
          << "uncertainties: stored tagged/inclusive response errors with a shared "
             "inclusive-denominator covariance; fit covariance multiplied by "
             "max(1, chi2/ndf)\n\n"
          << std::setprecision(7);
  for (int p = 0; p < kNPar; ++p)
    summary << parameterLabels[p] << " = " << fit.p(p) << " +/- "
            << std::sqrt(std::max(0.0, fit.cov(p, p))) << " %\n";
  summary << "\nchi2/ndf = " << fit.chi2 << " / " << fit.ndf << " = "
          << (fit.ndf > 0 ? fit.chi2 / fit.ndf : 0) << "\n"
          << "Z+jet chi2/rows = " << fit.chi2Sample[0] << " / "
          << fit.rowsSample[0] << "\n"
          << "photon+jet chi2/rows = " << fit.chi2Sample[1] << " / "
          << fit.rowsSample[1] << "\n"
          << "weighted design condition = " << fit.condition << "\n"
          << "covariance scale = " << fit.covarianceScale << "\n"
          << "response bins = " << bins.size() << ", available rows = " << rows.size()
          << ", unfolded points = " << points.size() << "\n\n"
          << "Model: every flavor uses delta_f(pt)=a_f+b_f*((pt/100)^-0.3-1).\n"
          << "This is a proof of principle, not a JEC payload: tag-SF/mistag, "
             "background, MC-template, tag/inclusive overlap, and other correlated "
             "experimental uncertainties are not included.\n";

  std::ofstream csv((stem + "_points.csv").Data());
  csv << "sample,flavor,pt,pt_halfwidth,response_shift_percent,error_percent,"
         "matrix_condition\n";
  for (const FlavorPoint &point : points)
    csv << (point.sample == 0 ? "zjet" : "gamjet") << ','
        << kFlavorLabels[point.flavor] << ',' << point.pt << ',' << point.ept
        << ',' << point.value << ',' << point.error << ',' << point.condition
        << '\n';
}

void drawResult(int year, bool useDataFractions,
                const std::vector<ResponseBin> &bins,
                const std::vector<FitRow> &rows,
                const std::vector<FlavorPoint> &points,
                const FitResult &fit,
                const std::array<double, kNFitFlavors> &anchor,
                double diagnosticFloorPercent) {
  std::vector<TGraphErrors *> bands(kNFitFlavors);
  std::vector<TGraphErrors *> curves(kNFitFlavors);
  std::vector<std::vector<TGraphErrors *>> dataGraphs(
      2, std::vector<TGraphErrors *>(kNFitFlavors, nullptr));

  double yExtent = 4.0;
  for (int f = 0; f < kNFitFlavors; ++f) {
    bands[f] = makeCurve(f, fit, 90, 300, 100, true);
    curves[f] = makeCurve(f, fit, 90, 300, 100, false);
    bands[f]->SetFillColorAlpha(kColors[f], 0.14);
    bands[f]->SetLineColor(kColors[f]);
    bands[f]->SetLineWidth(0);
    curves[f]->SetLineColor(kColors[f]);
    curves[f]->SetLineWidth(3);
    for (int i = 0; i < bands[f]->GetN(); ++i) {
      double x = 0, y = 0;
      bands[f]->GetPoint(i, x, y);
      yExtent = std::max(yExtent, std::abs(y) + bands[f]->GetErrorY(i));
    }
  }
  yExtent = std::min(12.0, std::ceil(1.20 * yExtent));

  for (int sample = 0; sample < 2; ++sample) {
    for (int f = 0; f < kNFitFlavors; ++f) {
      TGraphErrors *graph = new TGraphErrors();
      for (const FlavorPoint &point : points) {
        if (point.sample != sample || point.flavor != f)
          continue;
        if (!(point.pt >= 90.0 && point.pt <= 300.0 && point.error < 20.0))
          continue;
        const int n = graph->GetN();
        graph->SetPoint(n, point.pt, point.value);
        graph->SetPointError(n, 0, point.error);
      }
      graph->SetMarkerStyle(sample == 0 ? kOpenCircle : kFullSquare);
      graph->SetMarkerSize(sample == 0 ? 0.85 : 0.70);
      graph->SetMarkerColor(kColors[f]);
      graph->SetLineColor(kColors[f]);
      dataGraphs[sample][f] = graph;
    }
  }

  lumi_136TeV = Form("%d", year);
  writeExtraText = true;
  extraText = "Work in progress";
  TH1D *frame = new TH1D(Form("frame_%d", year),
                         ";p_{T} (GeV);Flavor response relative to reference (%)",
                         100, 80, 350);
  frame->SetMinimum(-yExtent);
  frame->SetMaximum(yExtent);
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  TCanvas *canvas =
      tdrCanvas(Form("cJetFlavor%d", year), frame, 8, 0, kRectangular);
  canvas->SetLogx();

  TLine *zero = new TLine(80, 0, 350, 0);
  zero->SetLineColor(kGray + 1);
  zero->SetLineStyle(kDotted);
  zero->Draw("SAME");

  for (int f = 0; f < kNFitFlavors; ++f)
    bands[f]->Draw("3 SAME");
  for (int f = 0; f < kNFitFlavors; ++f)
    curves[f]->Draw("L SAME");
  for (int sample = 0; sample < 2; ++sample)
    for (int f = 0; f < kNFitFlavors; ++f)
      dataGraphs[sample][f]->Draw("PZ SAME");

  TLegend *flavorLegend = tdrLeg(0.14, 0.63, 0.31, 0.88);
  flavorLegend->SetTextSize(0.037);
  for (int f = 0; f < kNFitFlavors; ++f)
    flavorLegend->AddEntry(curves[f], kFlavorLabels[f], "LF");
  flavorLegend->Draw();

  TGraphErrors *zStyle = new TGraphErrors(1);
  TGraphErrors *gammaStyle = new TGraphErrors(1);
  zStyle->SetMarkerStyle(kOpenCircle);
  zStyle->SetMarkerColor(kBlack);
  gammaStyle->SetMarkerStyle(kFullSquare);
  gammaStyle->SetMarkerColor(kBlack);
  TLegend *sampleLegend = tdrLeg(0.32, 0.72, 0.53, 0.88);
  sampleLegend->SetTextSize(0.037);
  sampleLegend->AddEntry(zStyle, "Z+jet", "P");
  sampleLegend->AddEntry(gammaStyle, "#gamma+jet", "P");
  sampleLegend->AddEntry(curves[0], "joint fit", "L");
  sampleLegend->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.036);
  text.DrawLatex(0.55, 0.86, "AK4 PUPPI, |#eta| < 1.3");
  text.DrawLatex(0.55, 0.81, "4-flavor proof, 90 < p_{T} < 300 GeV");
  text.SetTextSize(0.031);
  text.DrawLatex(0.55, 0.755,
                 Form("#chi^{2}/ndf = %.1f/%d", fit.chi2, fit.ndf));
  text.DrawLatex(0.55, 0.715,
                 useDataFractions ? "fractions constrained by data counts"
                                  : "nominal MC flavor fractions");
  text.DrawLatex(0.55, 0.675,
                 Form("stored errors + %.1f%% floor", diagnosticFloorPercent));
  text.DrawLatex(0.55, 0.635, "reference-mixture mode anchored");

  canvas->RedrawAxis();
  gSystem->mkdir("plots", true);
  const TString variant = useDataFractions ? "_dataFractions" : "";
  canvas->SaveAs(
      Form("plots/JetFlavorJointFit_%d%s.pdf", year, variant.Data()));
  canvas->SaveAs(
      Form("plots/JetFlavorJointFit_%d%s.png", year, variant.Data()));

  writeOutputs(year, useDataFractions, bins, rows, points, fit, anchor,
               diagnosticFloorPercent, bands, curves, dataGraphs);
}

}  // namespace jointflavor

void JetFlavor(int year = 2024, bool useDataFractions = false,
               double diagnosticFloorPercent = 0.2) {
  using namespace jointflavor;
  const TString inputName = Form("data/jecdata%dFLAVOR.root", year);
  TFile input(inputName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Could not open " << inputName << std::endl;
    return;
  }

  const std::vector<ResponseBin> nominalBins = loadInputs(input, false);
  const std::array<double, kNFitFlavors> anchor =
      makeAnchorWeights(nominalBins);
  const std::vector<ResponseBin> bins =
      useDataFractions ? loadInputs(input, true) : nominalBins;
  const std::vector<FitRow> rows = makeRows(bins);
  std::cout << "Loaded " << nominalBins.size() << " nominal response bins and "
            << rows.size() << " valid tagged double-ratio rows" << std::endl;
  const FitResult fit =
      fitBins(bins, anchor, diagnosticFloorPercent);
  if (!fit.ok) {
    std::cerr << "Joint fit failed for " << year << " (" << rows.size()
              << " usable tagged response rows)" << std::endl;
    return;
  }
  const std::vector<FlavorPoint> points = makeFlavorPoints(
      bins, anchor, diagnosticFloorPercent, fit.covarianceScale);

  std::cout << "\nJoint four-flavor fit, " << year << "\n"
            << "  rows: Z=" << fit.rowsSample[0]
            << ", gamma=" << fit.rowsSample[1] << "\n"
            << "  chi2/ndf = " << fit.chi2 << "/" << fit.ndf << " = "
            << (fit.ndf > 0 ? fit.chi2 / fit.ndf : 0) << "\n"
            << "  uds: " << fit.p(0) << " + " << fit.p(1)
            << "*((pt/100)^-0.3-1) %\n"
            << "  g:   " << fit.p(2) << " + " << fit.p(3)
            << "*((pt/100)^-0.3-1) %\n"
            << "  c:   " << fit.p(4) << " + " << fit.p(5)
            << "*((pt/100)^-0.3-1) %\n"
            << "  b:   " << fit.p(6) << " + " << fit.p(7)
            << "*((pt/100)^-0.3-1) %\n"
            << "  anchor (uds,g,c,b): " << anchor[0] << ", " << anchor[1]
            << ", " << anchor[2] << ", " << anchor[3] << "\n"
            << "  covariance scale: " << fit.covarianceScale << std::endl;

  drawResult(year, useDataFractions, bins, rows, points, fit, anchor,
             diagnosticFloorPercent);
}
