// Four-flavor inference from the v173 dijet FlavorMatrix profiles.
//
// Run with:
//   root -l -b -q 'DijetFlavorMatrix.C+()'
//
// The producer stores pair IDs as 10*tag+probe.  Reconstructed and true jet
// flavors are grouped as uds, g(+undefined), c, b.  Data tagging transitions
// are the minimum-KL/IPF projection of the MC joint distribution subject to
// the observed data reco-pair fractions and the MC truth-pair fractions.
// Relative responses are fitted with uds fixed to one; both tag/probe
// orientations are used, so only response ratios (not an absolute scale) are
// inferred from dijets.

#include "tdrstyle_mod22.C"

#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dijetfm {

constexpr int kNFlavor = 4;
constexpr int kNPair = 16;
constexpr double kResponseN = 1.00;
constexpr double kResponseU = 0.92;

const char *kFlavorName[kNFlavor] = {"uds", "g", "c", "b"};
const char *kFlavorLabel[kNFlavor] = {"uds", "g + unmatched", "c", "b"};
const int kColor[kNFlavor] = {kBlack, kRed + 1, kGreen + 2, kBlue + 1};
const int kMarker[kNFlavor] = {kFullCircle, kFullSquare, kFullTriangleUp,
                               kFullDiamond};
const char *kAxes[] = {"ab", "ad", "tc", "pf"};
const char *kComponents[] = {"m0", "m2", "mn", "mu"};

struct Summary {
  double entries = 0.;
  double effectiveEntries = 0.;
  double mean = 0.;
  double error = 0.;
  bool valid = false;
};

struct HDMSummary {
  double entries = 0.;
  double effectiveEntries = 0.;
  double value = 0.;
  double error = 0.;
  bool valid = false;
};

struct Profiles2D {
  TProfile2D *m0 = nullptr;
  TProfile2D *mn = nullptr;
  TProfile2D *mu = nullptr;
};

struct Profiles3D {
  TProfile3D *m0 = nullptr;
  TProfile3D *mn = nullptr;
  TProfile3D *mu = nullptr;
};

struct Result {
  bool valid = false;
  bool ipfConverged = false;
  int ipfIterations = 0;
  double ipfError = 0.;
  double low = 0.;
  double high = 0.;
  double totalData = 0.;
  double totalMC = 0.;
  double dataReco[kNPair] = {0.};
  double truthPair[kNPair] = {0.};
  double jointMC[kNPair][kNPair] = {{0.}};
  double jointData[kNPair][kNPair] = {{0.}};
  double efficiencyMC[kNFlavor][kNFlavor] = {{0.}};
  double efficiencyData[kNFlavor][kNFlavor] = {{0.}};
  double purityMC[kNFlavor][kNFlavor] = {{0.}};
  double purityData[kNFlavor][kNFlavor] = {{0.}};
  double response[kNFlavor] = {1., 1., 1., 1.};
  double responseError[kNFlavor] = {0., 0., 0., 0.};
  double truthBalanceMC[kNPair] = {0.};
  double truthBalanceData[kNPair] = {0.};
  int responseRows = 0;
  int responseRank = 0;
  int responseFallbacks = 0;
  double chi2 = 0.;
  int ndof = 0;
};

int groupFlavor(int id) {
  if (id >= 1 && id <= 3) return 0;
  if (id == 0 || id == 6) return 1;
  if (id == 4) return 2;
  if (id == 5) return 3;
  return -1;
}

int pairIndex(int tag, int probe) { return tag * kNFlavor + probe; }
int pairCode(int tagId, int probeId) { return 10 * tagId + probeId; }

std::string pairLabel(int pair) {
  return std::string(kFlavorName[pair / kNFlavor]) + "#rightarrow" +
         kFlavorName[pair % kNFlavor];
}

std::vector<int> rawIdsForGroup(int group) {
  if (group == 0) return {1, 2, 3};
  if (group == 1) return {0, 6};
  if (group == 2) return {4};
  if (group == 3) return {5};
  return {};
}

std::vector<int> rawCodesForPair(int groupedPair) {
  std::vector<int> result;
  const int tag = groupedPair / kNFlavor;
  const int probe = groupedPair % kNFlavor;
  for (int tagId : rawIdsForGroup(tag))
    for (int probeId : rawIdsForGroup(probe))
      result.push_back(pairCode(tagId, probeId));
  return result;
}

template <typename Profile>
double binEffectiveEntries(Profile *profile, int bin) {
  if (!profile) return 0.;
  const double sumw = profile->GetBinEntries(bin);
  const TArrayD *sumw2 = profile->GetBinSumw2();
  if (!sumw2 || bin >= sumw2->GetSize() || !(sumw2->At(bin) > 0.))
    return std::fabs(sumw);
  return sumw * sumw / sumw2->At(bin);
}

std::pair<int, int> xBinRange(const TAxis *axis, double low, double high) {
  const double epsilon = 1.e-8;
  const int first = std::max(1, axis->FindFixBin(low + epsilon));
  const int last = std::min(axis->GetNbins(), axis->FindFixBin(high - epsilon));
  return {first, last};
}

Summary summarize2D(TProfile2D *profile, int groupedPair, double low,
                    double high) {
  Summary result;
  if (!profile) return result;
  const auto range = xBinRange(profile->GetXaxis(), low, high);
  double weighted = 0.;
  double variance = 0.;
  double effective = 0.;
  for (int code : rawCodesForPair(groupedPair)) {
    const int ybin = profile->GetYaxis()->FindFixBin(code);
    for (int xbin = range.first; xbin <= range.second; ++xbin) {
      const int bin = profile->GetBin(xbin, ybin);
      const double entries = profile->GetBinEntries(bin);
      if (!std::isfinite(entries) || std::fabs(entries) < 1.e-15) continue;
      const double value = profile->GetBinContent(bin);
      const double error = profile->GetBinError(bin);
      result.entries += entries;
      weighted += entries * value;
      if (std::isfinite(error) && error >= 0.)
        variance += entries * entries * error * error;
      effective += std::max(0., binEffectiveEntries(profile, bin));
    }
  }
  if (std::fabs(result.entries) < 1.e-12) return result;
  result.mean = weighted / result.entries;
  result.error = std::sqrt(std::max(0., variance)) / std::fabs(result.entries);
  result.effectiveEntries = effective;
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

Summary summarize3D(TProfile3D *profile, int recoPair, int truthPair,
                    double low, double high) {
  Summary result;
  if (!profile) return result;
  const auto range = xBinRange(profile->GetXaxis(), low, high);
  double weighted = 0.;
  double variance = 0.;
  double effective = 0.;
  for (int recoCode : rawCodesForPair(recoPair)) {
    const int ybin = profile->GetYaxis()->FindFixBin(recoCode);
    for (int truthCode : rawCodesForPair(truthPair)) {
      const int zbin = profile->GetZaxis()->FindFixBin(truthCode);
      for (int xbin = range.first; xbin <= range.second; ++xbin) {
        const int bin = profile->GetBin(xbin, ybin, zbin);
        const double entries = profile->GetBinEntries(bin);
        if (!std::isfinite(entries) || std::fabs(entries) < 1.e-15) continue;
        const double value = profile->GetBinContent(bin);
        const double error = profile->GetBinError(bin);
        result.entries += entries;
        weighted += entries * value;
        if (std::isfinite(error) && error >= 0.)
          variance += entries * entries * error * error;
        effective += std::max(0., binEffectiveEntries(profile, bin));
      }
    }
  }
  if (std::fabs(result.entries) < 1.e-12) return result;
  result.mean = weighted / result.entries;
  result.error = std::sqrt(std::max(0., variance)) / std::fabs(result.entries);
  result.effectiveEntries = effective;
  result.valid = std::isfinite(result.mean) && std::isfinite(result.error);
  return result;
}

HDMSummary makeHDM(const Summary &m0, const Summary &mn, const Summary &mu) {
  HDMSummary result;
  if (!m0.valid || !mn.valid || !mu.valid) return result;
  const double numerator = m0.mean - mn.mean - mu.mean;
  const double denominator = 1. - mn.mean / kResponseN - mu.mean / kResponseU;
  if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
      std::fabs(denominator) < 1.e-9)
    return result;
  result.value = numerator / denominator;
  const double d0 = 1. / denominator;
  const double dn = (-denominator + numerator / kResponseN) /
                    (denominator * denominator);
  const double du = (-denominator + numerator / kResponseU) /
                    (denominator * denominator);
  result.error = std::sqrt(std::pow(d0 * m0.error, 2) +
                           std::pow(dn * mn.error, 2) +
                           std::pow(du * mu.error, 2));
  result.entries = m0.entries;
  result.effectiveEntries = m0.effectiveEntries;
  result.valid = std::isfinite(result.value) && std::isfinite(result.error);
  return result;
}

HDMSummary summarizeHDM(const Profiles2D &profiles, int pair, double low,
                        double high) {
  return makeHDM(summarize2D(profiles.m0, pair, low, high),
                 summarize2D(profiles.mn, pair, low, high),
                 summarize2D(profiles.mu, pair, low, high));
}

HDMSummary summarizeHDM(const Profiles3D &profiles, int recoPair,
                        int truthPair, double low, double high) {
  return makeHDM(summarize3D(profiles.m0, recoPair, truthPair, low, high),
                 summarize3D(profiles.mn, recoPair, truthPair, low, high),
                 summarize3D(profiles.mu, recoPair, truthPair, low, high));
}

Profiles2D load2D(TFile &file, const std::string &prefix,
                  const std::string &axis) {
  Profiles2D result;
  result.m0 = dynamic_cast<TProfile2D *>(file.Get(
      (prefix + "/p2m0" + axis + "_flavormatrix").c_str()));
  result.mn = dynamic_cast<TProfile2D *>(file.Get(
      (prefix + "/p2mn" + axis + "_flavormatrix").c_str()));
  result.mu = dynamic_cast<TProfile2D *>(file.Get(
      (prefix + "/p2mu" + axis + "_flavormatrix").c_str()));
  return result;
}

Profiles3D load3D(TFile &file, const std::string &prefix,
                  const std::string &axis) {
  Profiles3D result;
  result.m0 = dynamic_cast<TProfile3D *>(
      file.Get((prefix + "/p3m0" + axis + "_true").c_str()));
  result.mn = dynamic_cast<TProfile3D *>(
      file.Get((prefix + "/p3mn" + axis + "_true").c_str()));
  result.mu = dynamic_cast<TProfile3D *>(
      file.Get((prefix + "/p3mu" + axis + "_true").c_str()));
  return result;
}

bool complete(const Profiles2D &profiles) {
  return profiles.m0 && profiles.mn && profiles.mu;
}
bool complete(const Profiles3D &profiles) {
  return profiles.m0 && profiles.mn && profiles.mu;
}

double hdmToRatio(double hdm) {
  const double denominator = 3. - hdm;
  if (!(std::fabs(denominator) > 1.e-10))
    return std::numeric_limits<double>::quiet_NaN();
  const double ratio = (hdm + 1.) / denominator;
  return ratio > 0. ? ratio : std::numeric_limits<double>::quiet_NaN();
}

double ratioToHDM(double ratio) {
  if (!(ratio > 0.) || !std::isfinite(ratio))
    return std::numeric_limits<double>::quiet_NaN();
  return 1. + 2. * (ratio - 1.) / (ratio + 1.);
}

bool runIPF(Result &result) {
  const double tiny = 1.e-12;
  double inferred[kNPair][kNPair] = {{0.}};
  double normalization = 0.;
  for (int reco = 0; reco < kNPair; ++reco)
    for (int truth = 0; truth < kNPair; ++truth) {
      inferred[reco][truth] = std::max(result.jointMC[reco][truth], tiny);
      normalization += inferred[reco][truth];
    }
  for (int reco = 0; reco < kNPair; ++reco)
    for (int truth = 0; truth < kNPair; ++truth)
      inferred[reco][truth] /= normalization;

  for (int iteration = 1; iteration <= 10000; ++iteration) {
    for (int reco = 0; reco < kNPair; ++reco) {
      double sum = 0.;
      for (int truth = 0; truth < kNPair; ++truth) sum += inferred[reco][truth];
      if (sum > 0.)
        for (int truth = 0; truth < kNPair; ++truth)
          inferred[reco][truth] *= result.dataReco[reco] / sum;
    }
    for (int truth = 0; truth < kNPair; ++truth) {
      double sum = 0.;
      for (int reco = 0; reco < kNPair; ++reco) sum += inferred[reco][truth];
      if (sum > 0.)
        for (int reco = 0; reco < kNPair; ++reco)
          inferred[reco][truth] *= result.truthPair[truth] / sum;
    }
    double maximum = 0.;
    for (int reco = 0; reco < kNPair; ++reco) {
      double sum = 0.;
      for (int truth = 0; truth < kNPair; ++truth) sum += inferred[reco][truth];
      maximum = std::max(maximum, std::fabs(sum - result.dataReco[reco]));
    }
    for (int truth = 0; truth < kNPair; ++truth) {
      double sum = 0.;
      for (int reco = 0; reco < kNPair; ++reco) sum += inferred[reco][truth];
      maximum = std::max(maximum, std::fabs(sum - result.truthPair[truth]));
    }
    if (maximum < 1.e-10) {
      result.ipfConverged = true;
      result.ipfIterations = iteration;
      result.ipfError = maximum;
      break;
    }
  }
  if (!result.ipfConverged) return false;
  for (int reco = 0; reco < kNPair; ++reco)
    for (int truth = 0; truth < kNPair; ++truth)
      result.jointData[reco][truth] = inferred[reco][truth];
  return true;
}

void deriveSingleJetMatrices(Result &result) {
  double singleMC[kNFlavor][kNFlavor] = {{0.}};
  double singleData[kNFlavor][kNFlavor] = {{0.}};
  for (int recoPair = 0; recoPair < kNPair; ++recoPair) {
    const int recoTag = recoPair / kNFlavor;
    const int recoProbe = recoPair % kNFlavor;
    for (int truthPair = 0; truthPair < kNPair; ++truthPair) {
      const int truthTag = truthPair / kNFlavor;
      const int truthProbe = truthPair % kNFlavor;
      singleMC[recoTag][truthTag] += 0.5 * result.jointMC[recoPair][truthPair];
      singleMC[recoProbe][truthProbe] += 0.5 * result.jointMC[recoPair][truthPair];
      singleData[recoTag][truthTag] +=
          0.5 * result.jointData[recoPair][truthPair];
      singleData[recoProbe][truthProbe] +=
          0.5 * result.jointData[recoPair][truthPair];
    }
  }
  for (int truth = 0; truth < kNFlavor; ++truth) {
    double mcColumn = 0.;
    double dataColumn = 0.;
    for (int reco = 0; reco < kNFlavor; ++reco) {
      mcColumn += singleMC[reco][truth];
      dataColumn += singleData[reco][truth];
    }
    for (int reco = 0; reco < kNFlavor; ++reco) {
      result.efficiencyMC[reco][truth] =
          mcColumn > 0. ? singleMC[reco][truth] / mcColumn : 0.;
      result.efficiencyData[reco][truth] =
          dataColumn > 0. ? singleData[reco][truth] / dataColumn : 0.;
    }
  }
  for (int reco = 0; reco < kNFlavor; ++reco) {
    double mcRow = 0.;
    double dataRow = 0.;
    for (int truth = 0; truth < kNFlavor; ++truth) {
      mcRow += singleMC[reco][truth];
      dataRow += singleData[reco][truth];
    }
    for (int truth = 0; truth < kNFlavor; ++truth) {
      result.purityMC[reco][truth] =
          mcRow > 0. ? singleMC[reco][truth] / mcRow : 0.;
      result.purityData[reco][truth] =
          dataRow > 0. ? singleData[reco][truth] / dataRow : 0.;
    }
  }
}

bool solveResponse(Result &result, const HDMSummary dataHDM[kNPair],
                   HDMSummary mcHDM[kNPair][kNPair]) {
  double fallback[kNPair] = {0.};
  for (int truth = 0; truth < kNPair; ++truth) {
    double numerator = 0.;
    double denominator = 0.;
    for (int reco = 0; reco < kNPair; ++reco)
      if (mcHDM[reco][truth].valid) {
        const double weight = result.jointMC[reco][truth];
        numerator += weight * mcHDM[reco][truth].value;
        denominator += weight;
      }
    fallback[truth] = denominator > 0. ? numerator / denominator : 1.;
  }
  for (int reco = 0; reco < kNPair; ++reco)
    for (int truth = 0; truth < kNPair; ++truth)
      if (!mcHDM[reco][truth].valid && result.jointData[reco][truth] > 1.e-8) {
        mcHDM[reco][truth].value = fallback[truth];
        mcHDM[reco][truth].error = 0.;
        mcHDM[reco][truth].valid = true;
        ++result.responseFallbacks;
      }

  std::vector<int> rows;
  for (int reco = 0; reco < kNPair; ++reco)
    if (result.dataReco[reco] > 1.e-7 && dataHDM[reco].valid &&
        dataHDM[reco].effectiveEntries >= 25.)
      rows.push_back(reco);
  result.responseRows = rows.size();
  if (rows.size() < 3) return false;

  double delta[3] = {0., 0., 0.};
  TMatrixD lastHessian(3, 3);
  int rank = 0;
  for (int iteration = 0; iteration < 20; ++iteration) {
    TMatrixD hessian(3, 3);
    TVectorD gradient(3);
    hessian.Zero();
    gradient.Zero();
    for (int reco : rows) {
      const double rowFraction = result.dataReco[reco];
      if (!(rowFraction > 0.)) continue;
      double prediction = 0.;
      double jacobian[3] = {0., 0., 0.};
      for (int truth = 0; truth < kNPair; ++truth) {
        const double purity = result.jointData[reco][truth] / rowFraction;
        if (!(purity > 0.) || !mcHDM[reco][truth].valid) continue;
        const int truthTag = truth / kNFlavor;
        const int truthProbe = truth % kNFlavor;
        const double qmc = hdmToRatio(mcHDM[reco][truth].value);
        if (!std::isfinite(qmc)) continue;
        const double tagDelta = truthTag == 0 ? 0. : delta[truthTag - 1];
        const double probeDelta = truthProbe == 0 ? 0. : delta[truthProbe - 1];
        const double q = qmc * std::exp(probeDelta - tagDelta);
        const double value = ratioToHDM(q);
        prediction += purity * value;
        const double derivative = 4. * q / ((q + 1.) * (q + 1.));
        for (int parameter = 0; parameter < 3; ++parameter) {
          const int flavor = parameter + 1;
          const double sign = (truthProbe == flavor ? 1. : 0.) -
                              (truthTag == flavor ? 1. : 0.);
          jacobian[parameter] += purity * derivative * sign;
        }
      }
      const double sigma = std::max(2.e-4, dataHDM[reco].error);
      const double weight = 1. / (sigma * sigma);
      const double residual = dataHDM[reco].value - prediction;
      for (int first = 0; first < 3; ++first) {
        gradient[first] += jacobian[first] * weight * residual;
        for (int second = 0; second < 3; ++second)
          hessian(first, second) +=
              jacobian[first] * weight * jacobian[second];
      }
    }
    TDecompSVD decomposition(hessian);
    const TVectorD singular = decomposition.GetSig();
    rank = 0;
    double largest = 0.;
    for (int index = 0; index < singular.GetNrows(); ++index)
      largest = std::max(largest, std::fabs(singular[index]));
    for (int index = 0; index < singular.GetNrows(); ++index)
      if (std::fabs(singular[index]) > std::max(1.e-12, largest * 1.e-10))
        ++rank;
    Bool_t solved = false;
    const TVectorD step = decomposition.Solve(gradient, solved);
    if (!solved || rank < 3) return false;
    double maximum = 0.;
    for (int parameter = 0; parameter < 3; ++parameter) {
      const double bounded = std::max(-0.15, std::min(0.15, step[parameter]));
      delta[parameter] += bounded;
      maximum = std::max(maximum, std::fabs(bounded));
    }
    lastHessian.ResizeTo(hessian);
    lastHessian = hessian;
    if (maximum < 1.e-9) break;
  }
  result.responseRank = rank;
  result.response[0] = 1.;
  for (int parameter = 0; parameter < 3; ++parameter)
    result.response[parameter + 1] = std::exp(delta[parameter]);

  TDecompSVD covarianceSVD(lastHessian);
  Bool_t inverted = false;
  const TMatrixD covariance = covarianceSVD.Invert(inverted);
  if (inverted)
    for (int parameter = 0; parameter < 3; ++parameter)
      result.responseError[parameter + 1] =
          result.response[parameter + 1] *
          std::sqrt(std::max(0., covariance(parameter, parameter)));

  result.chi2 = 0.;
  for (int reco : rows) {
    const double rowFraction = result.dataReco[reco];
    double prediction = 0.;
    for (int truth = 0; truth < kNPair; ++truth) {
      const double purity = result.jointData[reco][truth] / rowFraction;
      if (!(purity > 0.) || !mcHDM[reco][truth].valid) continue;
      const int truthTag = truth / kNFlavor;
      const int truthProbe = truth % kNFlavor;
      const double qmc = hdmToRatio(mcHDM[reco][truth].value);
      if (!std::isfinite(qmc)) continue;
      const double q = qmc * result.response[truthProbe] /
                       result.response[truthTag];
      prediction += purity * ratioToHDM(q);
    }
    const double sigma = std::max(2.e-4, dataHDM[reco].error);
    result.chi2 += std::pow((dataHDM[reco].value - prediction) / sigma, 2);
  }
  result.ndof = static_cast<int>(rows.size()) - 3;
  for (int truth = 0; truth < kNPair; ++truth) {
    if (!(result.truthPair[truth] > 0.)) continue;
    double weightSum = 0.;
    double mcMean = 0.;
    double dataMean = 0.;
    const int truthTag = truth / kNFlavor;
    const int truthProbe = truth % kNFlavor;
    for (int reco = 0; reco < kNPair; ++reco) {
      const double weight = result.jointData[reco][truth];
      if (!(weight > 0.) || !mcHDM[reco][truth].valid) continue;
      const double qmc = hdmToRatio(mcHDM[reco][truth].value);
      if (!std::isfinite(qmc)) continue;
      const double qdata = qmc * result.response[truthProbe] /
                           result.response[truthTag];
      mcMean += weight * mcHDM[reco][truth].value;
      dataMean += weight * ratioToHDM(qdata);
      weightSum += weight;
    }
    if (weightSum > 0.) {
      result.truthBalanceMC[truth] = mcMean / weightSum;
      result.truthBalanceData[truth] = dataMean / weightSum;
    }
  }
  return true;
}

Result analyzeRange(const Profiles2D &data, const Profiles3D &mc,
                    double low, double high) {
  Result result;
  result.low = low;
  result.high = high;
  HDMSummary dataHDM[kNPair];
  HDMSummary mcHDM[kNPair][kNPair];
  double mcRaw[kNPair][kNPair] = {{0.}};
  double dataRaw[kNPair] = {0.};
  for (int reco = 0; reco < kNPair; ++reco) {
    const Summary count = summarize2D(data.m0, reco, low, high);
    dataRaw[reco] = count.entries;
    if (dataRaw[reco] < -1.e-8)
      throw std::runtime_error("Negative data population in FlavorMatrix");
    dataRaw[reco] = std::max(0., dataRaw[reco]);
    result.totalData += dataRaw[reco];
    dataHDM[reco] = summarizeHDM(data, reco, low, high);
    for (int truth = 0; truth < kNPair; ++truth) {
      const Summary mcCount = summarize3D(mc.m0, reco, truth, low, high);
      mcRaw[reco][truth] = mcCount.entries;
      if (mcRaw[reco][truth] < -1.e-7)
        throw std::runtime_error("Negative grouped MC population; signed-weight "
                                 "IPF is not implemented");
      mcRaw[reco][truth] = std::max(0., mcRaw[reco][truth]);
      result.totalMC += mcRaw[reco][truth];
      mcHDM[reco][truth] = summarizeHDM(mc, reco, truth, low, high);
    }
  }
  if (!(result.totalData > 0.) || !(result.totalMC > 0.)) return result;
  for (int reco = 0; reco < kNPair; ++reco) {
    result.dataReco[reco] = dataRaw[reco] / result.totalData;
    for (int truth = 0; truth < kNPair; ++truth) {
      result.jointMC[reco][truth] = mcRaw[reco][truth] / result.totalMC;
      result.truthPair[truth] += result.jointMC[reco][truth];
    }
  }
  if (!runIPF(result)) return result;
  deriveSingleJetMatrices(result);
  result.valid = solveResponse(result, dataHDM, mcHDM);
  return result;
}

void labelMatrix(TH2D *histogram) {
  for (int index = 0; index < kNFlavor; ++index) {
    histogram->GetXaxis()->SetBinLabel(index + 1, kFlavorName[index]);
    histogram->GetYaxis()->SetBinLabel(index + 1, kFlavorName[index]);
  }
  histogram->GetXaxis()->LabelsOption("h");
}

void drawTripleMatrix(const Result &result, bool efficiency,
                      const std::string &outputDirectory) {
  const char *kind = efficiency ? "efficiency" : "purity";
  std::unique_ptr<TCanvas> canvas(
      new TCanvas(Form("c_%s", kind), "", 1500, 480));
  canvas->Divide(3, 1, 0.001, 0.001);
  const char *titles[] = {"Simulation", "Data inferred", "Data / simulation"};
  std::array<std::unique_ptr<TH2D>, 3> histograms;
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3f");
  for (int panel = 0; panel < 3; ++panel) {
    canvas->cd(panel + 1);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.14);
    histograms[panel].reset(new TH2D(
        Form("h_%s_%d", kind, panel),
        Form(";%s category;True flavor;%s", efficiency ? "Reco" : "Reco",
             panel == 2 ? "ratio" : (efficiency ? "efficiency" : "purity")),
        4, 0, 4, 4, 0, 4));
    TH2D *histogram = histograms[panel].get();
    labelMatrix(histogram);
    for (int reco = 0; reco < kNFlavor; ++reco)
      for (int truth = 0; truth < kNFlavor; ++truth) {
        const double mc = efficiency ? result.efficiencyMC[reco][truth]
                                     : result.purityMC[reco][truth];
        const double data = efficiency ? result.efficiencyData[reco][truth]
                                       : result.purityData[reco][truth];
        histogram->SetBinContent(
            reco + 1, truth + 1,
            panel == 0 ? mc : panel == 1 ? data : (mc > 0. ? data / mc : 0.));
      }
    histogram->SetMinimum(panel == 2 ? 0.3 : 0.0);
    histogram->SetMaximum(panel == 2 ? 3.3 : 1.0);
    histogram->SetMarkerSize(1.45);
    histogram->SetTitle("");
    histogram->Draw("COLZ TEXT");
    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.032);
    note.DrawLatex(0.14, 0.94, "CMS Work in progress");
    note.SetTextSize(0.040);
    note.DrawLatex(0.42, 0.94, titles[panel]);
  }
  canvas->SaveAs((outputDirectory + "/" + kind + "_matrices.png").c_str());
  canvas->SaveAs((outputDirectory + "/" + kind + "_matrices.pdf").c_str());
  const char *suffix[] = {"mc", "data_inferred", "ratio_data_over_mc"};
  for (int panel = 0; panel < 3; ++panel) {
    std::unique_ptr<TCanvas> single(
        new TCanvas(Form("c_%s_single_%d", kind, panel), "", 560, 520));
    single->SetRightMargin(0.17);
    single->SetLeftMargin(0.15);
    single->SetTopMargin(0.11);
    std::unique_ptr<TH2D> histogram(
        dynamic_cast<TH2D *>(histograms[panel]->Clone(
            Form("h_%s_single_%d", kind, panel))));
    histogram->SetTitle("");
    histogram->SetMarkerSize(1.75);
    histogram->GetXaxis()->SetLabelSize(0.052);
    histogram->GetYaxis()->SetLabelSize(0.052);
    histogram->GetZaxis()->SetLabelSize(0.045);
    histogram->GetXaxis()->SetTitleSize(0.055);
    histogram->GetYaxis()->SetTitleSize(0.055);
    histogram->GetXaxis()->SetTitleOffset(0.90);
    histogram->GetYaxis()->SetTitleOffset(1.05);
    histogram->Draw("COLZ TEXT");
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.040);
    title.DrawLatex(0.16, 0.95, titles[panel]);
    single->SaveAs(
        (outputDirectory + "/" + kind + "_" + suffix[panel] + ".png")
            .c_str());
  }
}

void drawResponse(const std::vector<Result> &results,
                  const std::string &axis,
                  const std::string &outputDirectory) {
  std::unique_ptr<TH1D> frame(new TH1D(
      Form("response_frame_%s", axis.c_str()),
      ";p_{T}^{dijet} (GeV);Relative response data / MC (uds = 1)", 100, 55,
      3600));
  frame->SetMinimum(0.94);
  frame->SetMaximum(1.06);
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(
      tdrCanvas(Form("c_response_%s", axis.c_str()), frame.get(), 8, 0,
                kRectangular));
  canvas->SetLogx();
  TLine unity(60, 1., 3500, 1.);
  unity.SetLineStyle(kDotted);
  unity.Draw("SAME");
  std::array<std::unique_ptr<TGraphErrors>, kNFlavor> graphs;
  for (int flavor = 0; flavor < kNFlavor; ++flavor) {
    graphs[flavor].reset(new TGraphErrors());
    graphs[flavor]->SetName(Form("g_response_%s_%s", axis.c_str(),
                                 kFlavorName[flavor]));
    for (const Result &result : results) {
      if (!result.valid) continue;
      const double x = std::sqrt(result.low * result.high);
      const int point = graphs[flavor]->GetN();
      graphs[flavor]->SetPoint(point, x, result.response[flavor]);
      graphs[flavor]->SetPointError(point, 0., result.responseError[flavor]);
    }
    graphs[flavor]->SetLineColor(kColor[flavor]);
    graphs[flavor]->SetMarkerColor(kColor[flavor]);
    graphs[flavor]->SetMarkerStyle(kMarker[flavor]);
    graphs[flavor]->SetLineWidth(2);
    graphs[flavor]->Draw("PZL SAME");
  }
  TLegend *legend = tdrLeg(0.16, 0.68, 0.34, 0.88);
  for (int flavor = 0; flavor < kNFlavor; ++flavor)
    legend->AddEntry(graphs[flavor].get(), kFlavorLabel[flavor], "PL");
  legend->Draw();
  TLatex note;
  note.SetNDC();
  note.SetTextFont(42);
  note.SetTextSize(0.033);
  note.DrawLatex(0.52, 0.86, "AK4 PUPPI, both jets |#eta| < 1.3");
  note.DrawLatex(0.52, 0.81, Form("HDM, %s coordinate", axis.c_str()));
  note.DrawLatex(0.52, 0.76, "pair-level IPF; statistical errors only");
  canvas->SaveAs((outputDirectory + "/response_vs_pt_" + axis + ".png").c_str());
  canvas->SaveAs((outputDirectory + "/response_vs_pt_" + axis + ".pdf").c_str());
}

void drawDiagonalRatio(const std::vector<Result> &results, bool efficiency,
                       const std::string &outputDirectory) {
  const char *kind = efficiency ? "efficiency" : "purity";
  std::unique_ptr<TH1D> frame(new TH1D(
      Form("%s_ratio_frame", kind),
      Form(";p_{T}^{dijet} (GeV);Diagonal %s data / MC", kind), 100, 55,
      3600));
  frame->SetMinimum(0.75);
  frame->SetMaximum(1.25);
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(
      tdrCanvas(Form("c_%s_ratio", kind), frame.get(), 8, 0, kRectangular));
  canvas->SetLogx();
  TLine unity(60, 1., 3500, 1.);
  unity.SetLineStyle(kDotted);
  unity.Draw("SAME");
  std::array<std::unique_ptr<TGraphErrors>, kNFlavor> graphs;
  for (int flavor = 0; flavor < kNFlavor; ++flavor) {
    graphs[flavor].reset(new TGraphErrors());
    for (const Result &result : results) {
      if (!result.valid) continue;
      const double mc = efficiency ? result.efficiencyMC[flavor][flavor]
                                   : result.purityMC[flavor][flavor];
      const double data = efficiency ? result.efficiencyData[flavor][flavor]
                                     : result.purityData[flavor][flavor];
      if (!(mc > 0.)) continue;
      const int point = graphs[flavor]->GetN();
      graphs[flavor]->SetPoint(point, std::sqrt(result.low * result.high),
                               data / mc);
    }
    graphs[flavor]->SetLineColor(kColor[flavor]);
    graphs[flavor]->SetMarkerColor(kColor[flavor]);
    graphs[flavor]->SetMarkerStyle(kMarker[flavor]);
    graphs[flavor]->SetLineWidth(2);
    graphs[flavor]->Draw("PL SAME");
  }
  TLegend *legend = tdrLeg(0.16, 0.68, 0.36, 0.88);
  for (int flavor = 0; flavor < kNFlavor; ++flavor)
    legend->AddEntry(graphs[flavor].get(), kFlavorLabel[flavor], "PL");
  legend->Draw();
  TLatex note;
  note.SetNDC();
  note.SetTextFont(42);
  note.SetTextSize(0.033);
  note.DrawLatex(0.52, 0.86, "AK4 PUPPI, |#eta| < 1.3");
  note.DrawLatex(0.52, 0.81, "bisector-average p_{T}; pair-level IPF");
  canvas->SaveAs((outputDirectory + "/" + kind + "_diagonal_ratio_vs_pt.png").c_str());
  canvas->SaveAs((outputDirectory + "/" + kind + "_diagonal_ratio_vs_pt.pdf").c_str());
}

void drawAxisTaggingComparison(const std::array<Result, 4> &results,
                               const std::string &outputDirectory) {
  std::unique_ptr<TCanvas> canvas(
      new TCanvas("c_axis_tagging_comparison", "", 1100, 520));
  canvas->Divide(2, 1, 0.002, 0.002);
  std::array<std::unique_ptr<TH2D>, 2> histograms;
  const char *quantity[] = {"Diagonal efficiency data / MC",
                            "Diagonal purity data / MC"};
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3f");
  for (int panel = 0; panel < 2; ++panel) {
    canvas->cd(panel + 1);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.14);
    histograms[panel].reset(new TH2D(
        Form("h_axis_tagging_%d", panel),
        ";p_{T} coordinate;Flavor", 4, 0, 4, 4, 0, 4));
    TH2D *histogram = histograms[panel].get();
    for (int axis = 0; axis < 4; ++axis)
      histogram->GetXaxis()->SetBinLabel(axis + 1, kAxes[axis]);
    for (int flavor = 0; flavor < kNFlavor; ++flavor) {
      histogram->GetYaxis()->SetBinLabel(flavor + 1, kFlavorName[flavor]);
      for (int axis = 0; axis < 4; ++axis) {
        const double mc = panel == 0
                              ? results[axis].efficiencyMC[flavor][flavor]
                              : results[axis].purityMC[flavor][flavor];
        const double data = panel == 0
                                ? results[axis].efficiencyData[flavor][flavor]
                                : results[axis].purityData[flavor][flavor];
        histogram->SetBinContent(axis + 1, flavor + 1,
                                 mc > 0. ? data / mc : 0.);
      }
    }
    histogram->SetMinimum(0.70);
    histogram->SetMaximum(1.25);
    histogram->SetMarkerSize(1.75);
    histogram->GetXaxis()->SetLabelSize(0.052);
    histogram->GetYaxis()->SetLabelSize(0.052);
    histogram->GetZaxis()->SetLabelSize(0.045);
    histogram->GetXaxis()->SetTitleSize(0.055);
    histogram->GetYaxis()->SetTitleSize(0.055);
    histogram->SetTitle(quantity[panel]);
    histogram->Draw("COLZ TEXT");
  }
  canvas->SaveAs((outputDirectory + "/axis_binning_comparison.png").c_str());
  canvas->SaveAs((outputDirectory + "/axis_binning_comparison.pdf").c_str());
}

void drawTruthBalanceMatrices(const Result &result,
                              const std::string &outputDirectory) {
  std::unique_ptr<TCanvas> canvas(
      new TCanvas("c_truth_balance_matrices", "", 1500, 480));
  canvas->Divide(3, 1, 0.001, 0.001);
  const char *titles[] = {"MC truth-pair HDM", "Inferred data truth-pair HDM",
                          "HDM data / MC"};
  std::array<std::unique_ptr<TH2D>, 3> histograms;
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3f");
  for (int panel = 0; panel < 3; ++panel) {
    canvas->cd(panel + 1);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.14);
    histograms[panel].reset(new TH2D(
        Form("h_truth_balance_%d", panel),
        ";True tag flavor;True probe flavor", 4, 0, 4, 4, 0, 4));
    TH2D *histogram = histograms[panel].get();
    for (int flavor = 0; flavor < kNFlavor; ++flavor) {
      histogram->GetXaxis()->SetBinLabel(flavor + 1, kFlavorName[flavor]);
      histogram->GetYaxis()->SetBinLabel(flavor + 1, kFlavorName[flavor]);
    }
    for (int tag = 0; tag < kNFlavor; ++tag)
      for (int probe = 0; probe < kNFlavor; ++probe) {
        const int pair = pairIndex(tag, probe);
        const double mc = result.truthBalanceMC[pair];
        const double data = result.truthBalanceData[pair];
        histogram->SetBinContent(tag + 1, probe + 1,
                                 panel == 0 ? mc : panel == 1 ? data
                                                              : (mc != 0. ? data / mc : 0.));
      }
    histogram->SetMinimum(panel == 2 ? 0.94 : 0.90);
    histogram->SetMaximum(panel == 2 ? 1.06 : 1.10);
    histogram->SetMarkerSize(1.45);
    histogram->Draw("COLZ TEXT");
    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.032);
    note.DrawLatex(0.14, 0.94, "CMS Work in progress");
    note.SetTextSize(0.037);
    note.DrawLatex(0.40, 0.94, titles[panel]);
  }
  canvas->SaveAs((outputDirectory + "/truth_pair_response_matrices.png").c_str());
  canvas->SaveAs((outputDirectory + "/truth_pair_response_matrices.pdf").c_str());
  const char *suffix[] = {"mc", "data_inferred", "ratio_data_over_mc"};
  for (int panel = 0; panel < 3; ++panel) {
    std::unique_ptr<TCanvas> single(new TCanvas(
        Form("c_truth_balance_single_%d", panel), "", 560, 520));
    single->SetRightMargin(0.17);
    single->SetLeftMargin(0.15);
    single->SetTopMargin(0.11);
    std::unique_ptr<TH2D> histogram(
        dynamic_cast<TH2D *>(histograms[panel]->Clone(
            Form("h_truth_balance_single_%d", panel))));
    histogram->SetTitle("");
    histogram->SetMarkerSize(1.75);
    histogram->GetXaxis()->SetLabelSize(0.052);
    histogram->GetYaxis()->SetLabelSize(0.052);
    histogram->GetZaxis()->SetLabelSize(0.045);
    histogram->GetXaxis()->SetTitleSize(0.055);
    histogram->GetYaxis()->SetTitleSize(0.055);
    histogram->GetXaxis()->SetTitleOffset(0.90);
    histogram->GetYaxis()->SetTitleOffset(1.05);
    histogram->Draw("COLZ TEXT");
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.036);
    title.DrawLatex(0.16, 0.95, titles[panel]);
    single->SaveAs((outputDirectory + "/truth_pair_response_" +
                    suffix[panel] + ".png")
                       .c_str());
  }
}

std::vector<std::string> triggerDirectories(TFile &file) {
  std::vector<std::string> result;
  TIter next(file.GetListOfKeys());
  while (TKey *key = dynamic_cast<TKey *>(next())) {
    const std::string name = key->GetName();
    if (name.rfind("HLT_", 0) == 0) result.push_back(name);
  }
  std::sort(result.begin(), result.end());
  return result;
}

bool sameXBin(TProfile2D *combined, TProfile2D *source, int xbin) {
  if (!combined || !source || combined->GetNbinsY() != source->GetNbinsY())
    return false;
  for (int ybin = 1; ybin <= combined->GetNbinsY(); ++ybin) {
    const int combinedBin = combined->GetBin(xbin, ybin);
    const int sourceBin = source->GetBin(xbin, ybin);
    const double first = combined->GetBinEntries(combinedBin);
    const double second = source->GetBinEntries(sourceBin);
    const double tolerance = 1.e-10 * std::max(1., std::fabs(first));
    if (std::fabs(first - second) > tolerance) return false;
    if (std::fabs(first) > 1.e-12) {
      const double firstWeighted = first * combined->GetBinContent(combinedBin);
      const double secondWeighted = second * source->GetBinContent(sourceBin);
      const double weightedTolerance =
          1.e-10 * std::max(1., std::fabs(firstWeighted));
      if (std::fabs(firstWeighted - secondWeighted) > weightedTolerance)
        return false;
    }
  }
  return true;
}

void auditTriggerMerge(TFile &combinedFile, TFile &sourceFile,
                       const std::string &outputDirectory) {
  std::ofstream output(outputDirectory + "/trigger_merge_audit.tsv");
  output << "axis\tcomponent\tpt_low\tpt_high\tnonzero\texact_sources\n";
  const std::vector<std::string> triggers = triggerDirectories(sourceFile);
  for (const char *axis : kAxes)
    for (const char *component : kComponents) {
      const std::string name =
          std::string("p2") + component + axis + "_flavormatrix";
      TProfile2D *combined = dynamic_cast<TProfile2D *>(
          combinedFile.Get(("GluonJets/FlavorMatrix/" + name).c_str()));
      if (!combined) continue;
      for (int xbin = 1; xbin <= combined->GetNbinsX(); ++xbin) {
        double total = 0.;
        for (int ybin = 1; ybin <= combined->GetNbinsY(); ++ybin)
          total += combined->GetBinEntries(combined->GetBin(xbin, ybin));
        std::vector<std::string> matches;
        if (std::fabs(total) > 1.e-12)
          for (const std::string &trigger : triggers) {
            TProfile2D *source = dynamic_cast<TProfile2D *>(sourceFile.Get(
                (trigger + "/GluonJets/FlavorMatrix/" + name).c_str()));
            if (sameXBin(combined, source, xbin)) matches.push_back(trigger);
          }
        output << axis << '\t' << component << '\t'
               << combined->GetXaxis()->GetBinLowEdge(xbin) << '\t'
               << combined->GetXaxis()->GetBinUpEdge(xbin) << '\t'
               << (std::fabs(total) > 1.e-12 ? 1 : 0) << '\t';
        for (std::size_t index = 0; index < matches.size(); ++index) {
          if (index) output << ',';
          output << matches[index];
        }
        output << '\n';
      }
    }
}

void auditProfileClosure(TFile &mcFile, const std::string &outputDirectory) {
  std::ofstream output(outputDirectory + "/profile3d_to_2d_closure.tsv");
  output << "axis\tcomponent\tchecked_bins\tmax_count_relative"
            "\tmax_weighted_mean_absolute\tfailures\n";
  for (const char *axis : kAxes)
    for (const char *component : kComponents) {
      const std::string p2name = std::string("GluonJets/FlavorMatrix/p2") +
                                 component + axis + "_flavormatrix";
      const std::string p3name = std::string("GluonJets/FlavorMatrix/p3") +
                                 component + axis + "_true";
      TProfile2D *p2 = dynamic_cast<TProfile2D *>(mcFile.Get(p2name.c_str()));
      TProfile3D *p3 = dynamic_cast<TProfile3D *>(mcFile.Get(p3name.c_str()));
      if (!p2 || !p3) continue;
      double maxCountRelative = 0.;
      double maxMeanAbsolute = 0.;
      int checked = 0;
      int failures = 0;
      for (int xbin = 1; xbin <= p2->GetNbinsX(); ++xbin)
        for (int ybin = 1; ybin <= p2->GetNbinsY(); ++ybin) {
          const int p2bin = p2->GetBin(xbin, ybin);
          const double count2 = p2->GetBinEntries(p2bin);
          if (std::fabs(count2) < 1.e-12) continue;
          double count3 = 0.;
          double weighted3 = 0.;
          for (int zbin = 1; zbin <= p3->GetNbinsZ(); ++zbin) {
            const int p3bin = p3->GetBin(xbin, ybin, zbin);
            const double entries = p3->GetBinEntries(p3bin);
            count3 += entries;
            weighted3 += entries * p3->GetBinContent(p3bin);
          }
          const double countRelative =
              std::fabs(count3 - count2) / std::max(1., std::fabs(count2));
          const double mean3 = count3 != 0. ? weighted3 / count3 : 0.;
          const double meanAbsolute = std::fabs(mean3 - p2->GetBinContent(p2bin));
          maxCountRelative = std::max(maxCountRelative, countRelative);
          maxMeanAbsolute = std::max(maxMeanAbsolute, meanAbsolute);
          if (countRelative > 1.e-8 || meanAbsolute > 1.e-10) ++failures;
          ++checked;
        }
      output << axis << '\t' << component << '\t' << checked << '\t'
             << std::setprecision(12) << maxCountRelative << '\t'
             << maxMeanAbsolute << '\t' << failures << '\n';
    }
}

bool sameXBin(TProfile3D *combined, TProfile3D *source, int xbin) {
  if (!combined || !source || combined->GetNbinsY() != source->GetNbinsY() ||
      combined->GetNbinsZ() != source->GetNbinsZ())
    return false;
  for (int ybin = 1; ybin <= combined->GetNbinsY(); ++ybin)
    for (int zbin = 1; zbin <= combined->GetNbinsZ(); ++zbin) {
      const int combinedBin = combined->GetBin(xbin, ybin, zbin);
      const int sourceBin = source->GetBin(xbin, ybin, zbin);
      const double first = combined->GetBinEntries(combinedBin);
      const double second = source->GetBinEntries(sourceBin);
      const double tolerance = 1.e-10 * std::max(1., std::fabs(first));
      if (std::fabs(first - second) > tolerance) return false;
      if (std::fabs(first) > 1.e-12) {
        const double firstWeighted = first * combined->GetBinContent(combinedBin);
        const double secondWeighted = second * source->GetBinContent(sourceBin);
        const double weightedTolerance =
            1.e-10 * std::max(1., std::fabs(firstWeighted));
        if (std::fabs(firstWeighted - secondWeighted) > weightedTolerance)
          return false;
      }
    }
  return true;
}

void auditMCCombine(TFile &combinedFile, TFile &sourceFile,
                    const std::string &outputDirectory) {
  std::ofstream output(outputDirectory + "/mc_combine_audit.tsv");
  output << "axis\tcomponent\tdimension\tpopulated_pt_bins"
            "\texact_pt_bins\n";
  for (const char *axis : kAxes)
    for (const char *component : kComponents) {
      const std::string p2name =
          std::string("p2") + component + axis + "_flavormatrix";
      TProfile2D *combined2 = dynamic_cast<TProfile2D *>(combinedFile.Get(
          ("GluonJets/FlavorMatrix/" + p2name).c_str()));
      TProfile2D *source2 = dynamic_cast<TProfile2D *>(sourceFile.Get(
          ("HLT_MC/GluonJets/FlavorMatrix/" + p2name).c_str()));
      int populated2 = 0;
      int exact2 = 0;
      if (combined2 && source2)
        for (int xbin = 1; xbin <= combined2->GetNbinsX(); ++xbin) {
          double total = 0.;
          for (int ybin = 1; ybin <= combined2->GetNbinsY(); ++ybin)
            total += combined2->GetBinEntries(combined2->GetBin(xbin, ybin));
          if (std::fabs(total) < 1.e-12) continue;
          ++populated2;
          if (sameXBin(combined2, source2, xbin)) ++exact2;
        }
      output << axis << '\t' << component << "\t2\t" << populated2 << '\t'
             << exact2 << '\n';

      const std::string p3name =
          std::string("p3") + component + axis + "_true";
      TProfile3D *combined3 = dynamic_cast<TProfile3D *>(combinedFile.Get(
          ("GluonJets/FlavorMatrix/" + p3name).c_str()));
      TProfile3D *source3 = dynamic_cast<TProfile3D *>(sourceFile.Get(
          ("HLT_MC/GluonJets/FlavorMatrix/" + p3name).c_str()));
      int populated3 = 0;
      int exact3 = 0;
      if (combined3 && source3)
        for (int xbin = 1; xbin <= combined3->GetNbinsX(); ++xbin) {
          double total = 0.;
          for (int ybin = 1; ybin <= combined3->GetNbinsY(); ++ybin)
            for (int zbin = 1; zbin <= combined3->GetNbinsZ(); ++zbin)
              total += combined3->GetBinEntries(
                  combined3->GetBin(xbin, ybin, zbin));
          if (std::fabs(total) < 1.e-12) continue;
          ++populated3;
          if (sameXBin(combined3, source3, xbin)) ++exact3;
        }
      output << axis << '\t' << component << "\t3\t" << populated3 << '\t'
             << exact3 << '\n';
    }
}

void writeTables(const std::vector<Result> &results,
                 const Result &integrated,
                 const std::string &outputDirectory) {
  std::ofstream response(outputDirectory + "/response_vs_pt.tsv");
  response << "pt_low\tpt_high\tflavor\tdata_over_mc_relative_to_uds"
              "\tstat_error\tchi2\tndof\trows\trank\tipf_iterations\tvalid\n";
  for (const Result &result : results)
    for (int flavor = 0; flavor < kNFlavor; ++flavor)
      response << result.low << '\t' << result.high << '\t'
               << kFlavorName[flavor] << '\t' << result.response[flavor]
               << '\t' << result.responseError[flavor] << '\t' << result.chi2
               << '\t' << result.ndof << '\t' << result.responseRows << '\t'
               << result.responseRank << '\t' << result.ipfIterations << '\t'
               << result.valid << '\n';

  std::ofstream matrices(outputDirectory + "/integrated_matrices.tsv");
  matrices << "quantity\treco\ttruth\tmc\tdata_inferred\tratio\n";
  for (int reco = 0; reco < kNFlavor; ++reco)
    for (int truth = 0; truth < kNFlavor; ++truth) {
      const double emc = integrated.efficiencyMC[reco][truth];
      const double edata = integrated.efficiencyData[reco][truth];
      matrices << "efficiency\t" << kFlavorName[reco] << '\t'
               << kFlavorName[truth] << '\t' << emc << '\t' << edata << '\t'
               << (emc > 0. ? edata / emc : 0.) << '\n';
      const double pmc = integrated.purityMC[reco][truth];
      const double pdata = integrated.purityData[reco][truth];
      matrices << "purity\t" << kFlavorName[reco] << '\t'
               << kFlavorName[truth] << '\t' << pmc << '\t' << pdata << '\t'
               << (pmc > 0. ? pdata / pmc : 0.) << '\n';
    }

  std::ofstream tagging(outputDirectory + "/tagging_vs_pt.tsv");
  tagging << "pt_low\tpt_high\treco\ttruth\tefficiency_mc"
             "\tefficiency_data\tefficiency_sf\tpurity_mc\tpurity_data"
             "\tpurity_ratio\tipf_iterations\n";
  for (const Result &result : results)
    for (int reco = 0; reco < kNFlavor; ++reco)
      for (int truth = 0; truth < kNFlavor; ++truth) {
        const double emc = result.efficiencyMC[reco][truth];
        const double edata = result.efficiencyData[reco][truth];
        const double pmc = result.purityMC[reco][truth];
        const double pdata = result.purityData[reco][truth];
        tagging << result.low << '\t' << result.high << '\t'
                << kFlavorName[reco] << '\t' << kFlavorName[truth] << '\t'
                << emc << '\t' << edata << '\t'
                << (emc > 0. ? edata / emc : 0.) << '\t' << pmc << '\t'
                << pdata << '\t' << (pmc > 0. ? pdata / pmc : 0.) << '\t'
                << result.ipfIterations << '\n';
      }
}

void writeRoot(const std::vector<Result> &results, const Result &integrated,
               const Result &responseExample,
               const std::string &outputDirectory,
               const std::string &dataName, const std::string &mcName) {
  TFile output((outputDirectory + "/dijet_flavor_matrix.root").c_str(),
               "RECREATE");
  TH2D efficiencyMC("h2_efficiency_mc", ";Reco category;True flavor", 4, 0,
                    4, 4, 0, 4);
  TH2D efficiencyData("h2_efficiency_data_inferred",
                      ";Reco category;True flavor", 4, 0, 4, 4, 0, 4);
  TH2D purityMC("h2_purity_mc", ";Reco category;True flavor", 4, 0, 4, 4,
                0, 4);
  TH2D purityData("h2_purity_data_inferred", ";Reco category;True flavor", 4,
                  0, 4, 4, 0, 4);
  labelMatrix(&efficiencyMC);
  labelMatrix(&efficiencyData);
  labelMatrix(&purityMC);
  labelMatrix(&purityData);
  for (int reco = 0; reco < kNFlavor; ++reco)
    for (int truth = 0; truth < kNFlavor; ++truth) {
      efficiencyMC.SetBinContent(reco + 1, truth + 1,
                                 integrated.efficiencyMC[reco][truth]);
      efficiencyData.SetBinContent(reco + 1, truth + 1,
                                   integrated.efficiencyData[reco][truth]);
      purityMC.SetBinContent(reco + 1, truth + 1,
                             integrated.purityMC[reco][truth]);
      purityData.SetBinContent(reco + 1, truth + 1,
                               integrated.purityData[reco][truth]);
    }
  efficiencyMC.Write();
  efficiencyData.Write();
  purityMC.Write();
  purityData.Write();
  TH2D truthBalanceMC("h2_truth_pair_hdm_mc",
                      ";True tag flavor;True probe flavor", 4, 0, 4, 4, 0, 4);
  TH2D truthBalanceData("h2_truth_pair_hdm_data_inferred",
                        ";True tag flavor;True probe flavor", 4, 0, 4, 4, 0,
                        4);
  labelMatrix(&truthBalanceMC);
  labelMatrix(&truthBalanceData);
  for (int tag = 0; tag < kNFlavor; ++tag)
    for (int probe = 0; probe < kNFlavor; ++probe) {
      const int pair = pairIndex(tag, probe);
      truthBalanceMC.SetBinContent(tag + 1, probe + 1,
                                   responseExample.truthBalanceMC[pair]);
      truthBalanceData.SetBinContent(tag + 1, probe + 1,
                                     responseExample.truthBalanceData[pair]);
    }
  truthBalanceMC.Write();
  truthBalanceData.Write();
  for (int flavor = 0; flavor < kNFlavor; ++flavor) {
    TGraphErrors graph;
    graph.SetName(Form("g_response_relative_vs_pt_%s", kFlavorName[flavor]));
    for (const Result &result : results) {
      if (!result.valid) continue;
      const int point = graph.GetN();
      graph.SetPoint(point, std::sqrt(result.low * result.high),
                     result.response[flavor]);
      graph.SetPointError(point, 0., result.responseError[flavor]);
    }
    graph.Write();
  }
  TNamed provenance(
      "provenance",
      Form("data=%s; mc=%s; jet=AK4PUPPI; eta=barrel; axis=ab; "
           "pair_encoding=10*tag+probe; groups=uds,g+undefined,c,b; "
           "tagging=pair-level_minimum-KL_IPF_fixed_MC_truth; "
           "response=HDM_m0_mn_mu_Rn1_Ru0.92_relative_uds; "
           "uncertainty=conditional_profile_stat_only",
           dataName.c_str(), mcName.c_str()));
  provenance.Write();
  TParameter<double>("integrated_pt_low", integrated.low).Write();
  TParameter<double>("integrated_pt_high", integrated.high).Write();
  TParameter<double>("response_example_pt_low", responseExample.low).Write();
  TParameter<double>("response_example_pt_high", responseExample.high).Write();
  output.Close();
}

}  // namespace dijetfm

void DijetFlavorMatrix(
    const char *dataCombined = "data/jmenano_data_cmb_2026BD_JME_v173.root",
    const char *mcCombined =
        "data/jmenano_mc_cmb_Summer24MG_JMENANO_v173_v2.root",
    const char *dataOutput = "data/jmenano_data_out_2026BD_JME_v173.root",
    const char *mcOutput =
        "data/jmenano_mc_out_Summer24MG_JMENANO_v173_v2.root",
    const char *resultDirectory = "results/dijetFlavorMatrix",
    const char *plotDirectory = "plots/dijetFlavorMatrix") {
  using namespace dijetfm;
  std::unique_ptr<TFile> data(TFile::Open(dataCombined, "READ"));
  std::unique_ptr<TFile> mc(TFile::Open(mcCombined, "READ"));
  std::unique_ptr<TFile> dataOut(TFile::Open(dataOutput, "READ"));
  std::unique_ptr<TFile> mcOut(TFile::Open(mcOutput, "READ"));
  if (!data || data->IsZombie() || !mc || mc->IsZombie() || !dataOut ||
      dataOut->IsZombie() || !mcOut || mcOut->IsZombie())
    throw std::runtime_error("Could not open one or more v173 dijet inputs");
  gSystem->mkdir(resultDirectory, true);
  gSystem->mkdir(plotDirectory, true);

  auditTriggerMerge(*data, *dataOut, resultDirectory);
  auditProfileClosure(*mc, resultDirectory);
  auditMCCombine(*mc, *mcOut, resultDirectory);

  const Profiles2D dataAB = load2D(*data, "GluonJets/FlavorMatrix", "ab");
  const Profiles3D mcAB = load3D(*mc, "GluonJets/FlavorMatrix", "ab");
  if (!complete(dataAB) || !complete(mcAB))
    throw std::runtime_error("Missing central ab FlavorMatrix profiles");

  const double edges[] = {60,  85,  125, 180, 250, 350, 500, 600,
                          800, 1000, 1200, 1500, 1800, 2100, 2700, 3500};
  const int nRanges = sizeof(edges) / sizeof(edges[0]) - 1;
  std::vector<Result> results;
  for (int range = 0; range < nRanges; ++range) {
    try {
      Result result = analyzeRange(dataAB, mcAB, edges[range], edges[range + 1]);
      if (result.ipfConverged) results.push_back(result);
    } catch (const std::exception &error) {
      std::cerr << "Skipping pT " << edges[range] << "-" << edges[range + 1]
                << ": " << error.what() << std::endl;
    }
  }
  Result integrated = analyzeRange(dataAB, mcAB, 60., 2000.);
  if (!integrated.ipfConverged)
    throw std::runtime_error("Integrated pair-level IPF failed");

  writeExtraText = true;
  extraText = "Work in progress";
  lumi_136TeV = "2026 B-D";
  gStyle->SetPaintTextFormat(".3f");
  drawTripleMatrix(integrated, true, plotDirectory);
  drawTripleMatrix(integrated, false, plotDirectory);
  drawResponse(results, "ab", plotDirectory);
  drawDiagonalRatio(results, true, plotDirectory);
  drawDiagonalRatio(results, false, plotDirectory);
  const Result *responseExample = nullptr;
  for (const Result &result : results)
    if (result.valid && result.low == 350. && result.high == 500.) {
      drawTruthBalanceMatrices(result, plotDirectory);
      responseExample = &result;
      break;
    }
  if (!responseExample)
    throw std::runtime_error("Missing the 350-500 GeV response example");

  std::array<Result, 4> axisResults;
  std::ofstream axes(std::string(resultDirectory) + "/axis_comparison.tsv");
  axes << "axis\tflavor\tefficiency_mc\tefficiency_data\tefficiency_sf"
          "\tpurity_mc\tpurity_data\tpurity_ratio\n";
  for (int axisIndex = 0; axisIndex < 4; ++axisIndex) {
    const char *axis = kAxes[axisIndex];
    const Profiles2D dataProfiles =
        load2D(*data, "GluonJets/FlavorMatrix", axis);
    const Profiles3D mcProfiles = load3D(*mc, "GluonJets/FlavorMatrix", axis);
    if (!complete(dataProfiles) || !complete(mcProfiles)) continue;
    axisResults[axisIndex] = analyzeRange(dataProfiles, mcProfiles, 60., 2000.);
    for (int flavor = 0; flavor < kNFlavor; ++flavor) {
      const double emc = axisResults[axisIndex].efficiencyMC[flavor][flavor];
      const double edata =
          axisResults[axisIndex].efficiencyData[flavor][flavor];
      const double pmc = axisResults[axisIndex].purityMC[flavor][flavor];
      const double pdata = axisResults[axisIndex].purityData[flavor][flavor];
      axes << axis << '\t' << kFlavorName[flavor] << '\t' << emc << '\t'
           << edata << '\t' << (emc > 0. ? edata / emc : 0.) << '\t' << pmc
           << '\t' << pdata << '\t' << (pmc > 0. ? pdata / pmc : 0.)
           << '\n';
    }
  }
  axes.close();
  drawAxisTaggingComparison(axisResults, plotDirectory);

  writeTables(results, integrated, resultDirectory);
  writeRoot(results, integrated, *responseExample, resultDirectory,
            dataCombined, mcCombined);
  std::cout << "Dijet FlavorMatrix inference complete: " << results.size()
            << " pT ranges; integrated tagging IPF iterations="
            << integrated.ipfIterations << std::endl;
}
