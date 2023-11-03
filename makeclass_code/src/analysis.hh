#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TLatex.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "THStack.h"
#include "TString.h"
#include "TH1.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
#include "RooAbsRealLValue.h"
#include "RooAbsPdf.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooNLLVar.h"
#include "RooSimultaneous.h"
#include "RooExponential.h"
#include "RooGlobalFunc.h"
#include "RooCBShape.h"
#include "RooFormula.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include "vector"
#include "../aux/Events.hh"
#include "header.hh"
#include <cmath>

#include "TMath.h"
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

using namespace std;
using namespace RooFit;
using std::vector;
#include <sys/stat.h>

Events *p;

TH1D* ctau_distribution = new TH1D("ctau_distribution","ctau_distribution", 50, 0, 5);
TH1D* d3d_distribution = new TH1D("d3d_distribution", "d3d_distribution", 100, 0, 25);
TH1D* pdgid_distribution = new TH1D("pdgid_distribution", "pdgID distribution", 51, -25.5, 25.5);
TH1D* pdgid_original_distribution = new TH1D("pdgid_original_distribution", "pdgID distribution", 51, -25.5, 25.5);


TH1D* genmuonsv_d3d_distribution = new TH1D("genmuonsv_d3d_distribution", "Gen muon SV d3d distribution", 100, 0, 25);
TH1D* recomuonsv_d3d_distribution = new TH1D("recomuonsv_d3d_distribution", "Reco muon SV d3d distribution", 100, 0, 25);

TH1D* genmuonsv_dxy_distribution = new TH1D("genmuonsv_dxy_distribution", "Gen muon SV dxy distribution", 100, 0, 25);
TH1D* recomuonsv_dxy_distribution = new TH1D("recomuonsv_dxy_distribution", "Reco muon SV dxy distribution", 100, 0, 25);

TH1D* genminusreco_d3d_distribution = new TH1D("genminusreco_d3d_distribution", "d3d diff; Gen muonSV d3d - reco muonSV d3d;", 200, -25, 25);
TH1D* genminusreco_dxy_distribution = new TH1D("genminusreco_dxy_distribution", "dxy diff; Gen muonSV dxy - reco muonSV dxy;", 200, -25, 25);

TH2D* gen_d3d_reco_d3d = new TH2D("gen_d3d_reco_d3d", ";gen muon SV d3d; reco muon SV d3d ", 100, 0, 25, 100, 0, 25);
TH2D* gen_dxy_reco_dxy = new TH2D("gen_dxy_reco_dxy", ";gen muon SV dxy; reco muon SV dxy ", 100, 0, 25, 100, 0, 25);

TH1D* gen_muon_pT = new TH1D("gen_muon_pT", "Gen muon pT distribution", 100, 0, 100);
TH1D* gen_muon_vertex_dxy = new TH1D("gen_muon_vertex_dxy", "Gen muon vertex dxy distribution", 100, 0, 100);

TH2D* gen_d3d_reco_d3d_SV = new TH2D("gen_d3d_reco_d3d_SV", ";gen muon vertex d3d (cm); reco SV d3d (cm) ", 100, 0, 25, 100, 0, 25);
TH2D* gen_dxy_reco_dxy_SV = new TH2D("gen_dxy_reco_dxy_SV", ";gen muon vertex dxy (cm); reco SV dxy (cm)", 100, 0, 25, 100, 0, 25);

void get_event(int i) {
  if ( p->LoadTree(i) < 0) { 
    cout<<"\nProblem in LoadTree."
        <<"\nEntry: "<<i<<endl;
    exit(0);
  }
  p->fChain->GetEntry(i);
}


