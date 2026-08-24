// Reconstructed dijet q/g-tag diagnostics from the primitive GluonJets output.
//
// Run with, for example:
//   root -l -b -q 'DijetReco.C(2024)'
//
// These plots are deliberately not unfolded into quark/gluon response.  They
// expose the four reconstructed tag-pair balances and fractions, plus the
// gq-qg data-minus-MC contrast needed before a purity/SF likelihood is built.

#include "tdrstyle_mod22.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

namespace dijetreco {

constexpr int kNPairs = 4;
constexpr int kNComponents = 4;
const char *kPairs[kNPairs] = {"qq", "qg", "gq", "gg"};
const char *kPairLabels[kNPairs] = {"qq", "qg", "gq", "gg"};
const int kPairColors[kNPairs] = {kBlack, kGreen + 2, kRed + 1, kBlue + 1};
const char *kComponents[kNComponents] = {"m0", "m2", "mn", "mu"};
const char *kComponentLabels[kNComponents] = {
    "MPF0", "MPF2", "neutral", "unclustered"};
const int kComponentColors[kNComponents] = {kBlack, kBlue + 1, kGreen + 2,
                                             kOrange + 7};
const int kComponentMarkers[kNComponents] = {kFullCircle, kOpenSquare,
                                              kOpenDiamond, kOpenTriangleUp};

struct Inputs {
  TH1D *dataResponse[kNComponents][kNPairs] = {{nullptr}};
  TH1D *mcResponse[kNComponents][kNPairs] = {{nullptr}};
  TH1D *dataCounts[kNPairs] = {nullptr};
  TH1D *mcCounts[kNPairs] = {nullptr};
};

TString dataEra(int year) {
  if (year == 2024)
    return "2024CDEFGHI";
  if (year == 2025)
    return "2025CDEFG";
  if (year == 2026)
    return "2026BD";
  return "";
}

TH1D *projectProfile(TFile &file, const TString &path, const TString &name) {
  TProfile2D *profile = dynamic_cast<TProfile2D *>(file.Get(path));
  if (!profile || !(profile->GetEntries() > 0)) {
    std::cerr << "Missing or empty profile " << path << " in "
              << file.GetName() << std::endl;
    return nullptr;
  }
  TProfile *profileY = profile->ProfileY(name + "_py");
  TH1D *histogram = profileY ? profileY->ProjectionX(name) : nullptr;
  if (histogram)
    histogram->SetDirectory(nullptr);
  return histogram;
}

TH1D *cloneHistogram(TFile &file, const TString &path, const TString &name) {
  TH1D *source = dynamic_cast<TH1D *>(file.Get(path));
  if (!source) {
    std::cerr << "Missing histogram " << path << " in " << file.GetName()
              << std::endl;
    return nullptr;
  }
  TH1D *clone = dynamic_cast<TH1D *>(source->Clone(name));
  clone->SetDirectory(nullptr);
  return clone;
}

bool loadInputs(TFile &data, TFile &mc, Inputs &inputs) {
  bool ok = true;
  for (int component = 0; component < kNComponents; ++component) {
    for (int pair = 0; pair < kNPairs; ++pair) {
      inputs.dataResponse[component][pair] = projectProfile(
          data,
          Form("GluonJets/tight/p2%sab_%s", kComponents[component],
               kPairs[pair]),
          Form("data_%s_%s", kComponents[component], kPairs[pair]));
      inputs.mcResponse[component][pair] = projectProfile(
          mc,
          Form("HLT_MC/GluonJets/tight/p2%sab_%s", kComponents[component],
               kPairs[pair]),
          Form("mc_%s_%s", kComponents[component], kPairs[pair]));
      ok = ok && inputs.dataResponse[component][pair] &&
           inputs.mcResponse[component][pair];
    }
  }
  for (int pair = 0; pair < kNPairs; ++pair) {
    inputs.dataCounts[pair] = cloneHistogram(
        data, Form("GluonJets/tight/h_tagprobeab_%s", kPairs[pair]),
        Form("data_counts_%s", kPairs[pair]));
    inputs.mcCounts[pair] = cloneHistogram(
        mc, Form("HLT_MC/GluonJets/tight/h_tagprobeab_%s", kPairs[pair]),
        Form("mc_counts_%s", kPairs[pair]));
    ok = ok && inputs.dataCounts[pair] && inputs.mcCounts[pair];
  }
  return ok;
}

void styleFrame(TH1D *frame) {
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
}

void addRecoLabel(double x = 0.50, double y = 0.86) {
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.033);
  text.DrawLatex(x, y, "AK4 PUPPI, probe |#eta| < 1.3");
  text.DrawLatex(x, y - 0.045, "tight reco tag pairs, bisector average");
  text.SetTextSize(0.030);
  text.DrawLatex(x, y - 0.090, "reco diagnostic - not flavor unfolded");
}

void drawBalance(int year, const Inputs &inputs) {
  TH1D *frame = new TH1D(Form("balance_frame_%d", year),
                         ";p_{T} (GeV);MPF0 dijet balance", 100, 60, 1500);
  frame->SetMinimum(0.85);
  frame->SetMaximum(1.15);
  styleFrame(frame);
  TCanvas *canvas =
      tdrCanvas(Form("cDijetBalance%d", year), frame, 8, 0, kRectangular);
  canvas->SetLogx();

  for (int pair = 0; pair < kNPairs; ++pair) {
    TH1D *mc = inputs.mcResponse[0][pair];
    TH1D *data = inputs.dataResponse[0][pair];
    mc->GetXaxis()->SetRangeUser(60, 1500);
    data->GetXaxis()->SetRangeUser(60, 1500);
    mc->SetLineColor(kPairColors[pair]);
    mc->SetLineWidth(2);
    mc->SetLineStyle(kSolid);
    data->SetMarkerColor(kPairColors[pair]);
    data->SetLineColor(kPairColors[pair]);
    data->SetMarkerStyle(kFullCircle);
    data->SetMarkerSize(0.65);
    mc->Draw("HIST SAME");
    data->Draw("PZ SAME");
  }

  TLegend *pairLegend = tdrLeg(0.15, 0.65, 0.27, 0.88);
  pairLegend->SetTextSize(0.037);
  for (int pair = 0; pair < kNPairs; ++pair)
    pairLegend->AddEntry(inputs.dataResponse[0][pair], kPairLabels[pair], "PL");
  pairLegend->Draw();

  TH1D *dataStyle = dynamic_cast<TH1D *>(inputs.dataResponse[0][0]->Clone());
  TH1D *mcStyle = dynamic_cast<TH1D *>(inputs.mcResponse[0][0]->Clone());
  dataStyle->SetMarkerColor(kBlack);
  dataStyle->SetLineColor(kBlack);
  mcStyle->SetLineColor(kBlack);
  TLegend *styleLegend = tdrLeg(0.29, 0.75, 0.42, 0.88);
  styleLegend->SetTextSize(0.037);
  styleLegend->AddEntry(dataStyle, "data", "P");
  styleLegend->AddEntry(mcStyle, "MC", "L");
  styleLegend->Draw();
  addRecoLabel();

  canvas->RedrawAxis();
  gSystem->mkdir("plots", true);
  canvas->SaveAs(Form("plots/DijetRecoBalance_%d.pdf", year));
  canvas->SaveAs(Form("plots/DijetRecoBalance_%d.png", year));
}

void drawFractions(int year, const Inputs &inputs,
                   std::array<TH1D *, kNPairs> &dataFractions,
                   std::array<TH1D *, kNPairs> &mcFractions) {
  TH1D *dataTotal = dynamic_cast<TH1D *>(inputs.dataCounts[0]->Clone(
      Form("data_pair_total_%d", year)));
  TH1D *mcTotal = dynamic_cast<TH1D *>(inputs.mcCounts[0]->Clone(
      Form("mc_pair_total_%d", year)));
  dataTotal->SetDirectory(nullptr);
  mcTotal->SetDirectory(nullptr);
  for (int pair = 1; pair < kNPairs; ++pair) {
    dataTotal->Add(inputs.dataCounts[pair]);
    mcTotal->Add(inputs.mcCounts[pair]);
  }
  for (int pair = 0; pair < kNPairs; ++pair) {
    dataFractions[pair] = dynamic_cast<TH1D *>(inputs.dataCounts[pair]->Clone(
        Form("data_fraction_%s", kPairs[pair])));
    mcFractions[pair] = dynamic_cast<TH1D *>(inputs.mcCounts[pair]->Clone(
        Form("mc_fraction_%s", kPairs[pair])));
    dataFractions[pair]->SetDirectory(nullptr);
    mcFractions[pair]->SetDirectory(nullptr);
    dataFractions[pair]->Divide(dataTotal);
    mcFractions[pair]->Divide(mcTotal);
  }

  TH1D *frame = new TH1D(Form("fraction_frame_%d", year),
                         ";p_{T} (GeV);Reco tag-pair fraction", 100, 60,
                         1500);
  frame->SetMinimum(0.0);
  frame->SetMaximum(0.55);
  styleFrame(frame);
  TCanvas *canvas =
      tdrCanvas(Form("cDijetFractions%d", year), frame, 8, 0, kRectangular);
  canvas->SetLogx();
  for (int pair = 0; pair < kNPairs; ++pair) {
    mcFractions[pair]->SetLineColor(kPairColors[pair]);
    mcFractions[pair]->SetLineWidth(2);
    dataFractions[pair]->SetMarkerColor(kPairColors[pair]);
    dataFractions[pair]->SetLineColor(kPairColors[pair]);
    dataFractions[pair]->SetMarkerStyle(kFullCircle);
    dataFractions[pair]->SetMarkerSize(0.65);
    mcFractions[pair]->GetXaxis()->SetRangeUser(60, 1500);
    dataFractions[pair]->GetXaxis()->SetRangeUser(60, 1500);
    mcFractions[pair]->Draw("HIST SAME");
    dataFractions[pair]->Draw("PZ SAME");
  }
  TLegend *pairLegend = tdrLeg(0.15, 0.65, 0.27, 0.88);
  pairLegend->SetTextSize(0.037);
  for (int pair = 0; pair < kNPairs; ++pair)
    pairLegend->AddEntry(dataFractions[pair], kPairLabels[pair], "PL");
  pairLegend->Draw();
  TH1D *dataStyle = dynamic_cast<TH1D *>(dataFractions[0]->Clone());
  TH1D *mcStyle = dynamic_cast<TH1D *>(mcFractions[0]->Clone());
  dataStyle->SetMarkerColor(kBlack);
  dataStyle->SetLineColor(kBlack);
  mcStyle->SetLineColor(kBlack);
  TLegend *styleLegend = tdrLeg(0.29, 0.75, 0.42, 0.88);
  styleLegend->SetTextSize(0.037);
  styleLegend->AddEntry(dataStyle, "data", "P");
  styleLegend->AddEntry(mcStyle, "MC", "L");
  styleLegend->Draw();
  addRecoLabel();
  canvas->RedrawAxis();
  canvas->SaveAs(Form("plots/DijetRecoFractions_%d.pdf", year));
  canvas->SaveAs(Form("plots/DijetRecoFractions_%d.png", year));
}

std::array<TH1D *, kNComponents> makeContrasts(const Inputs &inputs) {
  std::array<TH1D *, kNComponents> contrasts = {nullptr};
  for (int component = 0; component < kNComponents; ++component) {
    TH1D *data = dynamic_cast<TH1D *>(
        inputs.dataResponse[component][2]->Clone(
            Form("data_gq_minus_qg_%s", kComponents[component])));
    TH1D *mc = dynamic_cast<TH1D *>(
        inputs.mcResponse[component][2]->Clone(
            Form("mc_gq_minus_qg_%s", kComponents[component])));
    data->SetDirectory(nullptr);
    mc->SetDirectory(nullptr);
    data->Add(inputs.dataResponse[component][1], -1.0);
    mc->Add(inputs.mcResponse[component][1], -1.0);
    contrasts[component] = dynamic_cast<TH1D *>(data->Clone(
        Form("contrast_%s", kComponents[component])));
    contrasts[component]->SetDirectory(nullptr);
    contrasts[component]->Add(mc, -1.0);
    contrasts[component]->Scale(100.0);
  }
  return contrasts;
}

void drawContrast(int year,
                  const std::array<TH1D *, kNComponents> &contrasts) {
  double extent = 3.0;
  for (TH1D *histogram : contrasts) {
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
      const double pt = histogram->GetBinCenter(bin);
      if (!(pt >= 60 && pt <= 1500) || !(histogram->GetBinError(bin) > 0))
        continue;
      extent = std::max(extent, std::abs(histogram->GetBinContent(bin)) +
                                    histogram->GetBinError(bin));
    }
  }
  extent = std::min(8.0, std::ceil(1.15 * extent));
  TH1D *frame = new TH1D(
      Form("contrast_frame_%d", year),
      ";p_{T} (GeV);Reco gq-qg contrast (%)", 100, 60, 1500);
  frame->SetMinimum(-extent);
  frame->SetMaximum(extent);
  styleFrame(frame);
  TCanvas *canvas =
      tdrCanvas(Form("cDijetContrast%d", year), frame, 8, 0, kRectangular);
  canvas->SetLogx();
  TLine *zero = new TLine(60, 0, 1500, 0);
  zero->SetLineStyle(kDotted);
  zero->SetLineColor(kGray + 1);
  zero->Draw("SAME");
  for (int component = 0; component < kNComponents; ++component) {
    TH1D *histogram = contrasts[component];
    histogram->SetMarkerColor(kComponentColors[component]);
    histogram->SetLineColor(kComponentColors[component]);
    histogram->SetMarkerStyle(kComponentMarkers[component]);
    histogram->SetMarkerSize(0.75);
    histogram->GetXaxis()->SetRangeUser(60, 1500);
    histogram->Draw("PZ SAME");
  }
  TLegend *legend = tdrLeg(0.15, 0.65, 0.31, 0.88);
  legend->SetTextSize(0.036);
  for (int component = 0; component < kNComponents; ++component)
    legend->AddEntry(contrasts[component], kComponentLabels[component], "PL");
  legend->Draw();
  addRecoLabel(0.50, 0.86);
  canvas->RedrawAxis();
  canvas->SaveAs(Form("plots/DijetRecoContrast_%d.pdf", year));
  canvas->SaveAs(Form("plots/DijetRecoContrast_%d.png", year));
}

void writeOutputs(int year, const TString &dataName, const TString &mcName,
                  const Inputs &inputs,
                  const std::array<TH1D *, kNPairs> &dataFractions,
                  const std::array<TH1D *, kNPairs> &mcFractions,
                  const std::array<TH1D *, kNComponents> &contrasts) {
  gSystem->mkdir("results", true);
  TFile output(Form("results/DijetReco_%d.root", year), "RECREATE");
  for (int component = 0; component < kNComponents; ++component) {
    for (int pair = 0; pair < kNPairs; ++pair) {
      inputs.dataResponse[component][pair]->Write();
      inputs.mcResponse[component][pair]->Write();
    }
    contrasts[component]->Write();
  }
  for (int pair = 0; pair < kNPairs; ++pair) {
    inputs.dataCounts[pair]->Write();
    inputs.mcCounts[pair]->Write();
    dataFractions[pair]->Write();
    mcFractions[pair]->Write();
  }
  TNamed provenance(
      "provenance",
      Form("year=%d; jet=AK4PUPPI; data=%s; mc=%s; "
           "directory=GluonJets/tight; coordinate=ab; "
           "tagger=PNetQvG_0.45_with_UParT_heavy_veto; "
           "status=reco_category_diagnostic_not_unfolded",
           year, dataName.Data(), mcName.Data()));
  provenance.Write();
  output.Close();

  std::ofstream summary(Form("results/DijetReco_%d.txt", year));
  summary << "Dijet reconstructed tag-category diagnostic\n"
          << "year: " << year << "\n"
          << "data: " << dataName << "\n"
          << "MC: " << mcName << "\n"
          << "selection: GluonJets/tight, bisector-average pT, probe |eta|<1.3\n"
          << "tagging: PNet QvG 0.45 after UParT heavy-flavor veto\n"
          << "interpretation: reco tag categories only; no purity/SF unfolding\n\n"
          << "MPF0 contrast 100*[(gq-qg)data-(gq-qg)MC]:\n"
          << std::setprecision(6);
  const double referencePt[] = {92.5, 195.0, 450.0, 900.0};
  for (double pt : referencePt) {
    const int bin = contrasts[0]->FindBin(pt);
    summary << "  pT~" << contrasts[0]->GetBinCenter(bin) << " GeV: "
            << contrasts[0]->GetBinContent(bin) << " +/- "
            << contrasts[0]->GetBinError(bin) << " %\n";
  }

  std::ofstream csv(Form("results/DijetRecoContrast_%d.csv", year));
  csv << "component,pt,value_percent,error_percent\n";
  for (int component = 0; component < kNComponents; ++component) {
    TH1D *histogram = contrasts[component];
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
      if (!(histogram->GetBinError(bin) > 0))
        continue;
      csv << kComponents[component] << ',' << histogram->GetBinCenter(bin)
          << ',' << histogram->GetBinContent(bin) << ','
          << histogram->GetBinError(bin) << '\n';
    }
  }
}

}  // namespace dijetreco

void DijetReco(int year = 2024, const char *jecsysBase = "../jecsys3") {
  using namespace dijetreco;
  const TString era = dataEra(year);
  if (era.IsNull()) {
    std::cerr << "Supported years are 2024, 2025, and 2026" << std::endl;
    return;
  }
  const TString dataName =
      Form("%s/rootfiles/Prompt/Jet_v170/jmenano_data_cmb_%s_JME_v170.root",
           jecsysBase, era.Data());
  const TString mcName =
      Form("%s/rootfiles/Prompt/Jet_v171/"
           "jmenano_mc_out_Summer24MG_JMENANO_v171.root",
           jecsysBase);
  TFile data(dataName, "READ");
  TFile mc(mcName, "READ");
  if (data.IsZombie() || mc.IsZombie()) {
    std::cerr << "Could not open dijet inputs:\n  " << dataName << "\n  "
              << mcName << std::endl;
    return;
  }
  Inputs inputs;
  if (!loadInputs(data, mc, inputs)) {
    std::cerr << "Dijet input schema is incomplete" << std::endl;
    return;
  }

  lumi_136TeV = Form("%d + Summer24 MC", year);
  writeExtraText = true;
  extraText = "Work in progress";
  gSystem->mkdir("plots", true);
  std::array<TH1D *, kNPairs> dataFractions = {nullptr};
  std::array<TH1D *, kNPairs> mcFractions = {nullptr};
  const std::array<TH1D *, kNComponents> contrasts = makeContrasts(inputs);
  drawBalance(year, inputs);
  drawFractions(year, inputs, dataFractions, mcFractions);
  drawContrast(year, contrasts);
  writeOutputs(year, dataName, mcName, inputs, dataFractions, mcFractions,
               contrasts);
  std::cout << "Wrote dijet reco diagnostics for " << year << std::endl;
}
