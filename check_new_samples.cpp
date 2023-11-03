#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include <cmath>



void nano_analysis(){
 /*
 TChain *signal = new TChain();
 signal->Add("/vols/cms/mc3909/bparkProductionAll_V1p0/HiddenValley_vector_m_2_ctau_10_xiO_1_xiL_1_privateMC_11X_NANOAODSIM_v1p0_generationSync*.root/Events");
 ROOT::RDataFrame df(*signal);
  */

 ROOT::RDataFrame df_scenarioA("Events", "/vols/cms/jleonhol/samples/scenarioA_mpi_10_mA_3p33_ctau_100/nano.root");

 //auto h_muon_dxy = df.Histo1D({"hist_muon_dxy", "; Muon dxy (cm); Number of muons", 200, -25.0, 25.0}, "Muon_dxy");
 auto h_muon_dxy_A = df_scenarioA.Histo1D({"hist_muon_dxy", "; Muon dxy (cm); Number of muons", 200, -25.0, 25.0}, "Muon_dxy");
 
 TChain *signal_B1 = new TChain();
 signal_B1->Add("/vols/cms/jleonhol/samples/ul_pu/scenarioA_mpi_4_mA_1p33_ctau_10/nano.root/Events");
 //signal_B1->Add("/vols/cms/mc3909/bparkProductionV3/HiddenValley_vector_m_2_ctau_10_xiO_1_xiL_1_privateMC_11X_NANOAODSIM_v3_generationForBParking/*.root/Events");
 ROOT::RDataFrame df_scenarioB1(*signal_B1);

 auto h_muon_dxy_B1 = df_scenarioB1.Histo1D({"hist_muon_dxy", "; Muon dxy (cm); Number of muons", 200, -25.0, 25.0}, "Muon_dxy");

 TFile *f = new TFile("/vols/cms/khl216/darkshower_analysis_check_new_models/makeclass_code/output/ul_pu_new/outputSV_UL_scenarioA_mpi_4_mA_1p33_ctau_10.root");
 TH2D *gen_d3d_reco_d3d = (TH2D*) f->Get("gen_d3d_reco_d3d");
 TH2D *gen_dxy_reco_dxy = (TH2D*) f->Get("gen_dxy_reco_dxy");
 TH2D *gen_d3d_reco_d3d_SV = (TH2D*) f->Get("gen_d3d_reco_d3d_SV");
 TH2D *gen_dxy_reco_dxy_SV = (TH2D*) f->Get("gen_dxy_reco_dxy_SV");

 gStyle->SetOptStat(0); gStyle->SetTextFont(42);

 //TCanvas* c1 = new TCanvas("", "", 800, 700);
 //c1->SetLogy();
 //h_muon_dxy_B1->DrawClone("hist");

 TCanvas* c2 = new TCanvas("", "", 800, 700);
 gen_d3d_reco_d3d->DrawClone("COLZ");

 TCanvas* c3 = new TCanvas("", "", 800, 700);
 gen_dxy_reco_dxy->DrawClone("COLZ");

 TCanvas* c6 = new TCanvas("", "", 800, 700);
 gen_d3d_reco_d3d_SV->DrawClone("COLZ");

 TCanvas* c7 = new TCanvas("", "", 800, 700);
 gen_dxy_reco_dxy_SV->DrawClone("COLZ");

 /*
 TChain *signal1 = new TChain();
 TChain *signal2 = new TChain();
 TChain *signal3 = new TChain();
 TChain *signal4 = new TChain();

 signal1->Add("/vols/cms/jleonhol/samples/pu/scenarioA_mpi_4_mA_1p33_ctau_10/nano.root/Events");
 signal2->Add("/vols/cms/jleonhol/samples/pu/scenarioB1_mpi_4_mA_1p33_ctau_10/nano.root/Events");
 signal3->Add("/vols/cms/jleonhol/samples/pu_new/scenarioB2_mpi_4_mA_2p10_ctau_10/nano.root/Events");
 signal4->Add("/vols/cms/jleonhol/samples/pu_new/scenarioC_mpi_10_mA_8p00_ctau_10/nano.root/Events");

 ROOT::RDataFrame df_signal_original1(*signal1);
 ROOT::RDataFrame df_signal_original2(*signal2);
 ROOT::RDataFrame df_signal_original3(*signal3);
 ROOT::RDataFrame df_signal_original4(*signal4);

 Double_t no_of_entries_signal1_original = df_signal_original1.Count().GetValue();
 Double_t no_of_entries_signal2_original = df_signal_original2.Count().GetValue();
 Double_t no_of_entries_signal3_original = df_signal_original3.Count().GetValue();
 Double_t no_of_entries_signal4_original = df_signal_original4.Count().GetValue();

 TFile *f1 = new TFile("/vols/cms/khl216/darkshower_analysis_check_new_models/makeclass_code/output/pu_with_gen_muon_histo/outputSV_scenarioA_mpi_4_mA_1p33_ctau_10.root");
 TFile *f2 = new TFile("/vols/cms/khl216/darkshower_analysis_check_new_models/makeclass_code/output/pu_with_gen_muon_histo/outputSV_scenarioB1_mpi_4_mA_1p33_ctau_10.root");
 TFile *f3 = new TFile("/vols/cms/khl216/darkshower_analysis_check_new_models/makeclass_code/output/pu_with_gen_muon_histo/outputSV_scenarioB2_mpi_4_mA_2p10_ctau_10.root");
 TFile *f4 = new TFile("/vols/cms/khl216/darkshower_analysis_check_new_models/makeclass_code/output/pu_with_gen_muon_histo/outputSV_scenarioC_mpi_10_mA_8p00_ctau_10.root");

 TH1D *gen_muon_pT1 = (TH1D*) f1->Get("gen_muon_pT");
 TH1D *gen_muon_pT2 = (TH1D*) f2->Get("gen_muon_pT");
 TH1D *gen_muon_pT3 = (TH1D*) f3->Get("gen_muon_pT"); 
 TH1D *gen_muon_pT4 = (TH1D*) f4->Get("gen_muon_pT");

 TH1D *gen_muon_vertex_dxy1 = (TH1D*) f1->Get("gen_muon_vertex_dxy"); 
 TH1D *gen_muon_vertex_dxy2 = (TH1D*) f2->Get("gen_muon_vertex_dxy");
 TH1D *gen_muon_vertex_dxy3 = (TH1D*) f3->Get("gen_muon_vertex_dxy");
 TH1D *gen_muon_vertex_dxy4 = (TH1D*) f4->Get("gen_muon_vertex_dxy"); 

 Double_t branching_ratio = 0.01;
 Double_t xs = 43.9;
 Double_t lumi = 33.6*1000;

 gen_muon_pT1->Scale(lumi*xs*branching_ratio/no_of_entries_signal1_original);  
 gen_muon_pT2->Scale(lumi*xs*branching_ratio/no_of_entries_signal2_original);
 gen_muon_pT3->Scale(lumi*xs*branching_ratio/no_of_entries_signal3_original);
 gen_muon_pT4->Scale(lumi*xs*branching_ratio/no_of_entries_signal4_original); 

 gen_muon_vertex_dxy1->Scale(lumi*xs*branching_ratio/no_of_entries_signal1_original);
 gen_muon_vertex_dxy2->Scale(lumi*xs*branching_ratio/no_of_entries_signal2_original);
 gen_muon_vertex_dxy3->Scale(lumi*xs*branching_ratio/no_of_entries_signal3_original);
 gen_muon_vertex_dxy4->Scale(lumi*xs*branching_ratio/no_of_entries_signal4_original);

 TCanvas* c5 = new TCanvas("", "", 800, 700);
 c5->SetLogy();
 gen_muon_pT1->SetLineColor(kRed);
 gen_muon_pT2->SetLineColor(kBlue+1);
 gen_muon_pT3->SetLineColor(kCyan+1);
 gen_muon_pT4->SetLineColor(kGreen+1);

 gen_muon_pT2->GetXaxis()->SetTitle("pT (GeV)");
 gen_muon_pT2->GetYaxis()->SetTitle("Number of muons");

 //h_muonsvdxy_signal1->DrawClone("hist");
 gen_muon_pT2->DrawClone("hist Same");
 gen_muon_pT3->DrawClone("hist Same");
 gen_muon_pT4->DrawClone("hist Same");

 TLegend* out_legend5 = new TLegend(0.48, 0.795, 0.88, 0.875);
 //out_legend3->AddEntry(h_muonsvdxy_signal1->GetName(), "Scenario A, m_{#pi}=10 GeV, m_{A}=3.33 GeV, c#tau=10 mm)", "l");
 out_legend5->AddEntry(gen_muon_pT2, "Scenario B1, m_{#pi}=4 GeV, m_{A}=1.33 GeV, c#tau=10 mm)", "l");
 out_legend5->AddEntry(gen_muon_pT3, "Scenario B2, m_{#pi}=4 GeV, m_{A}=2.10 GeV, c#tau=10 mm)", "l");
 out_legend5->AddEntry(gen_muon_pT4, "Scenario C, m_{#pi}=10 GeV, m_{A}=8.00 GeV, c#tau=10 mm)", "l");
 out_legend5->Draw("Same");

 TCanvas* c4 = new TCanvas("", "", 800, 700);
 c4->SetLogy();
 gen_muon_vertex_dxy1->SetLineColor(kRed);
 gen_muon_vertex_dxy2->SetLineColor(kBlue+1);
 gen_muon_vertex_dxy3->SetLineColor(kCyan+1);
 gen_muon_vertex_dxy4->SetLineColor(kGreen+1);

 gen_muon_vertex_dxy2->GetXaxis()->SetTitle("dxy (cm)");
 gen_muon_vertex_dxy2->GetYaxis()->SetTitle("Number of vertices");

 //h_muonsvdxy_signal1->DrawClone("hist");
 gen_muon_vertex_dxy2->DrawClone("hist Same");
 gen_muon_vertex_dxy3->DrawClone("hist Same");
 gen_muon_vertex_dxy4->DrawClone("hist Same");

 TLegend* out_legend4 = new TLegend(0.48, 0.795, 0.88, 0.875);
 //out_legend3->AddEntry(h_muonsvdxy_signal1->GetName(), "Scenario A, m_{#pi}=10 GeV, m_{A}=3.33 GeV, c#tau=10 mm)", "l");
 out_legend4->AddEntry(gen_muon_vertex_dxy2, "Scenario B1, m_{#pi}=4 GeV, m_{A}=1.33 GeV, c#tau=10 mm)", "l");
 out_legend4->AddEntry(gen_muon_vertex_dxy3, "Scenario B2, m_{#pi}=4 GeV, m_{A}=2.10 GeV, c#tau=10 mm)", "l");
 out_legend4->AddEntry(gen_muon_vertex_dxy4, "Scenario C, m_{#pi}=10 GeV, m_{A}=8.00 GeV, c#tau=10 mm)", "l");
 out_legend4->Draw("Same");
 */

}
