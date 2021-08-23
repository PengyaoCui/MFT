#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#endif

void SetHistoStyle(TH1F *h, int index=1);
void SetStyle(Bool_t graypalette=kFALSE);
const Int_t fillColors[]  = {kRed-10,kBlue-9,kCyan+2,kOrange-9, kGreen-8,kBlue-9,kOrange-9, kRed-3,kRed-3,kOrange-9,kOrange-9,  kMagenta-9, kOrange-9,kCyan-8,kYellow-7,kGreen-8}; // for syst bands
const Int_t markerColors[]= {kRed+1,kBlack,kAzure+2,kYellow+2, kGreen+3};
const Int_t colors[]      = {kAzure+1, kViolet-8,kRed-7,kOrange-3,kYellow+2,kRed+1,kBlack,kAzure+2,kBlue,kGreen+3,kCyan+2};
// const Int_t colors[]      = {kOrange+1,kBlack,kRed+1,kAzure+2, kBlue-9,kRed-10, kCyan+2,kCyan+2,kYellow+2,kBlue+1, kCyan+2, kMagenta+1, kCyan+2};
const Int_t markers[]     = {kFullCircle,kOpenCircle,kFullCircle,kFullDiamond, kFullSquare,kOpenSquare,kOpenTriangleUp,kOpenTriangleDown,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
const double sizes[]      = {0.4, 0.5,0.7,0.5,0.7,0.5,0.5,0.7};

bool save = true;

const int nFiles = 1;

const char *path [nFiles]= {
  "/home/cuipengyao/MFT/code/run/RCut_0.0100"};
  //"/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons_LEta_moreEv/reqfz/RCut_0.0045",
  //"/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons_LEta_moreEv/reqfz/RCut_0.0050",
  //"/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons_LEta_moreEv/reqfz/RCut_0.0055",
  //"/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons_LEta_moreEv/reqfz/RCut_0.0058"};
  // "/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons/reqfzLTF_loss1/RCut_0.0058_verif"};
 // "/Users/francisco/WORKDIR/O2/LTFOptimisation//pythia/reqfzLTFandCA_loss1/RCut_0.0100"};
  //"/Users/francisco/WORKDIR/O2/LTFOptimisation/1000ev_LTFRandCAx4/RCut_0.0100_RCA_0.04",
 // "/Users/francisco/WORKDIR/O2/LTFOptimisation/10muons/reqfzLTFandCA_ELoss1/RCut_0.0100"};
// const double rVal = 0.0100;
// const double caVal = 0.04;
const char *refPath = "/home/cuipengyao/MFT/code/macro/RCut_0.0100";

const char *legText [nFiles]= {"R^{Cylinder} = 0.0100"};
  //"R^{Cylinder} = 0.0100", "R_{O}^{Cone} = 0.0045", "R_{O}^{Cone} = 0.0050", "R_{O}^{Cone} = 0.0055", "R_{O}^{Cone} = 0.0058"};//,"R^{LTF}=f(z), R^{CA}=4*max"};//,"R^{LTF} + R^{CA}=f(z)" };
const char *fileName = "MinRCone";

//_____________________________________________________________________________
void CompToAE()
{
  SetStyle();
	TFile *refFile = new TFile(Form("%s/MFTEffCheck.root",refPath));
  TFile *files[nFiles];
	TH1F *hEffpT[nFiles], *hEffp[nFiles], *hEffEta[nFiles], *hEffVz[nFiles];
	TH1F *hRatiopT[nFiles], *hRatiop[nFiles], *hRatioEta[nFiles], *hRatioVz[nFiles];
  TH1F *hLTFEffpT[nFiles], *hLTFEffp[nFiles], *hLTFEffEta[nFiles], *hLTFEffVz[nFiles];
  TH1F *hLTFPercpT[nFiles], *hLTFPercp[nFiles], *hLTFPercEta[nFiles], *hLTFPercVz[nFiles];
  TH1F *hLTFRatiopT[nFiles], *hLTFRatiop[nFiles], *hLTFRatioEta[nFiles], *hLTFRatioVz[nFiles];
  TH1F *hEffRatiopT[nFiles], *hEffRatiop[nFiles], *hEffRatioEta[nFiles], *hEffRatioVz[nFiles];
  TH1F *hPercRatiopT[nFiles], *hPercRatiop[nFiles], *hPercRatioEta[nFiles], *hPercRatioVz[nFiles];


  TH1F* hrefpT= static_cast<TH1F*>(refFile->Get("MFT_LTF_Efficiency_pT"));
  TH1F* hrefp= static_cast<TH1F*>(refFile->Get("MFT_LTF_Efficiency_p"));
  TH1F* hrefEta= static_cast<TH1F*>(refFile->Get("MFT_LTF_Efficiency_eta"));
  TH1F* hrefVz= static_cast<TH1F*>(refFile->Get("MFT_LTF_Efficiency_Vz"));


  TH1F* hLTFTrackspT[nFiles],*hLTFTracksp[nFiles],*hLTFTracksEta[nFiles],*hLTFTracksVz[nFiles];
  TH1F *hMCTrackspT[nFiles], *hMCTracksp[nFiles],*hMCTracksEta[nFiles], *hMCTracksVz[nFiles]; 
  TH1F *hTrackablespT[nFiles],*hTrackablesp[nFiles],*hTrackablesEta[nFiles],*hTrackablesVz[nFiles];
  TH1F *hAcceptancepT[nFiles],*hAcceptancep[nFiles],*hAcceptanceEta[nFiles],*hAcceptanceVz[nFiles];
  TH1F *hAcceptanceEffpT[nFiles],*hAcceptanceEffp[nFiles],*hAcceptanceEffEta[nFiles],*hAcceptanceEffVz[nFiles];
	for (int i = 0; i < nFiles; ++i){ 

		files[i]=new TFile(Form("%s/MFTEffCheck.root",path[i]));
		if(!files[i]) {cout << "File " << Form("%s/MFTEffCheck.root",path[i]) <<" not found" <<endl; return;}
		//histos : efficiency vs pt, efficiency vs p, efficiency vs eta, efficiency vs vz
		hEffpT[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_pT"));
		hEffp[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_p"));
		hEffEta[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_eta"));
		hEffVz[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_Vz"));
		//percentage of tracks found by LTF vs CA	
		hRatiopT[i]= static_cast<TH1F*>(files[i]->Get("RatioLTFtoCA_pt"));
		hRatiop[i]= static_cast<TH1F*>(files[i]->Get("RatioLTFtoCA_p"));
		hRatioEta[i]= static_cast<TH1F*>(files[i]->Get("RatioLTFtoCA_eta"));
		hRatioVz[i]= static_cast<TH1F*>(files[i]->Get("RatioLTFtoCA_Vz"));
    //histos : LTF efficiency vs pt, p, eta, vz
    hLTFEffpT[i]= static_cast<TH1F*>(files[i]->Get("MFT_LTF_Efficiency_pT"));
    hLTFEffp[i]= static_cast<TH1F*>(files[i]->Get("MFT_LTF_Efficiency_p"));
    hLTFEffEta[i]= static_cast<TH1F*>(files[i]->Get("MFT_LTF_Efficiency_eta"));
    hLTFEffVz[i]= static_cast<TH1F*>(files[i]->Get("MFT_LTF_Efficiency_Vz"));
    //percentage of tracks found by LTF vs total
    hLTFPercpT[i]= static_cast<TH1F*>(files[i]->Get("PercLTF_pt"));
    hLTFPercp[i]= static_cast<TH1F*>(files[i]->Get("PercLTF_p"));
    hLTFPercEta[i]= static_cast<TH1F*>(files[i]->Get("PercLTF_eta"));
    hLTFPercVz[i]= static_cast<TH1F*>(files[i]->Get("PercLTF_Vz"));


    hLTFRatiopT[i] = static_cast<TH1F*>(hLTFEffpT[i]->Clone());
    hLTFRatiopT[i]->Divide(hrefpT);
    hLTFRatioEta[i] = static_cast<TH1F*>(hLTFEffEta[i]->Clone());
    hLTFRatioEta[i]->Divide(hrefEta);
    hLTFRatioVz[i] = static_cast<TH1F*>(hLTFEffVz[i]->Clone());
    hLTFRatioVz[i]->Divide(hrefVz);
    hLTFRatiop[i] = static_cast<TH1F*>(hLTFEffp[i]->Clone());
    hLTFRatiop[i]->Divide(hrefp);


    hEffRatiopT[i] = static_cast<TH1F*>(hEffpT[i]->Clone());
    hEffRatiopT[i]->Divide(hEffpT[0]);
    hEffRatioEta[i] = static_cast<TH1F*>(hEffEta[i]->Clone());
    hEffRatioEta[i]->Divide(hEffEta[0]);
    hEffRatioVz[i] = static_cast<TH1F*>(hEffVz[i]->Clone());
    hEffRatioVz[i]->Divide(hEffVz[0]);
    hEffRatiop[i] = static_cast<TH1F*>(hEffp[i]->Clone());
    hEffRatiop[i]->Divide(hEffp[0]);

    hPercRatiopT[i] = static_cast<TH1F*>(hLTFPercpT[i]->Clone());
    hPercRatiopT[i]->Divide(hLTFPercpT[0]);
    hPercRatioEta[i] = static_cast<TH1F*>(hLTFPercEta[i]->Clone());
    hPercRatioEta[i]->Divide(hLTFPercEta[0]);
    hPercRatioVz[i] = static_cast<TH1F*>(hLTFPercVz[i]->Clone());
    hPercRatioVz[i]->Divide(hLTFPercVz[0]);
    hPercRatiop[i] = static_cast<TH1F*>(hLTFPercp[i]->Clone());
    hPercRatiop[i]->Divide(hLTFPercp[0]);


    hLTFTrackspT[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_pT"));
    hLTFTracksp[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_p"));
    hLTFTracksEta[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_eta"));
    hLTFTracksVz[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_Vz"));
    hMCTrackspT[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_pT"));
    hMCTracksp[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_p"));
    hMCTracksEta[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_eta"));
    hMCTracksVz[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_Vz"));

    hTrackablespT[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_pT"));
    hTrackablesp[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_p"));
    hTrackablesEta[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_eta"));
    hTrackablesVz[i]= static_cast<TH1F*>(files[i]->Get("MFT_Efficiency_Vz"));

    hAcceptancepT[i] = static_cast<TH1F*>(hTrackablespT[i]->Clone());
    hAcceptancepT[i]->Divide(hMCTrackspT[i]);
    hAcceptanceEffpT[i] = static_cast<TH1F*>(hLTFPercpT[i]->Clone());
    hAcceptanceEffpT[i]->Divide(hMCTrackspT[i]);

    hAcceptancep[i] = static_cast<TH1F*>(hTrackablesp[i]->Clone());
    hAcceptancep[i]->Divide(hMCTracksp[i]);
    hAcceptanceEffp[i] = static_cast<TH1F*>(hLTFPercp[i]->Clone());
    hAcceptanceEffp[i]->Divide(hMCTracksp[i]);

    hAcceptanceEta[i] = static_cast<TH1F*>(hTrackablesEta[i]->Clone());
    hAcceptanceEta[i]->Divide(hMCTracksEta[i]);
    hAcceptanceEffEta[i] = static_cast<TH1F*>(hLTFPercEta[i]->Clone());
    hAcceptanceEffEta[i]->Divide(hMCTracksEta[i]);

    hAcceptanceVz[i] = static_cast<TH1F*>(hTrackablesVz[i]->Clone());
    hAcceptanceVz[i]->Divide(hMCTracksVz[i]);
    hAcceptanceEffVz[i] = static_cast<TH1F*>(hLTFPercVz[i]->Clone());
    hAcceptanceEffVz[i]->Divide(hMCTracksVz[i]);


		if(!hRatioVz[i]) {cout << "Histo ratio vz " << Form("%s/MFTEffCheck.root",path[i]) <<" not found" <<endl; return;}

    	SetHistoStyle(hEffpT[i],i);
    	SetHistoStyle(hEffp[i],i);
    	SetHistoStyle(hEffEta[i],i);
    	SetHistoStyle(hEffVz[i],i);
    	SetHistoStyle(hRatiopT[i],i);
    	SetHistoStyle(hRatiop[i],i);
    	SetHistoStyle(hRatioEta[i],i);
    	SetHistoStyle(hRatioVz[i],i);
      SetHistoStyle(hLTFEffpT[i],i);
      SetHistoStyle(hLTFEffp[i],i);
      SetHistoStyle(hLTFEffEta[i],i);
      SetHistoStyle(hLTFEffVz[i],i);
      SetHistoStyle(hLTFPercpT[i],i);
      SetHistoStyle(hLTFPercp[i],i);
      SetHistoStyle(hLTFPercEta[i],i);
      SetHistoStyle(hLTFPercVz[i],i);

      SetHistoStyle(hEffRatiopT[i],i);
      SetHistoStyle(hEffRatiop[i],i);
      SetHistoStyle(hEffRatioEta[i],i);
      SetHistoStyle(hEffRatioVz[i],i);
      SetHistoStyle(hLTFRatiopT[i],i);
      SetHistoStyle(hLTFRatiop[i],i);
      SetHistoStyle(hLTFRatioEta[i],i);
      SetHistoStyle(hLTFRatioVz[i],i);

      SetHistoStyle(hPercRatiopT[i],i);
      SetHistoStyle(hPercRatiop[i],i);
      SetHistoStyle(hPercRatioEta[i],i);
      SetHistoStyle(hPercRatioVz[i],i);


      SetHistoStyle(hAcceptancepT[i],i);
      SetHistoStyle(hAcceptancep[i],i);
      SetHistoStyle(hAcceptanceEta[i],i);
      SetHistoStyle(hAcceptanceVz[i],i);
      SetHistoStyle(hAcceptanceEffpT[i],i);
      SetHistoStyle(hAcceptanceEffp[i],i);
      SetHistoStyle(hAcceptanceEffEta[i],i);
      SetHistoStyle(hAcceptanceEffVz[i],i);
	}


  TCanvas *cEff= new TCanvas("MFT Efficiency","MFT Efficiency", 1600, 1600);
  cEff->Divide(2,2,0.01,0.01);
  // TCanvas *cRatio= new TCanvas("LTF to CA Ratio","LTF to CA Ratio", 1600, 1600);
  // cRatio->Divide(2,2,0.01,0.01);
  TCanvas *cAcc= new TCanvas("LTF to CA Ratio","LTF to CA Ratio", 1600, 1600);
  cAcc->Divide(2,2,0.01,0.01);
  TCanvas *cLTFPerc= new TCanvas("LTF Percentage","LTF Percentage", 1600, 1600);
  cLTFPerc->Divide(2,2,0.01,0.01);
  TCanvas *cLTFEff= new TCanvas("LTF Efficiency","LTF Efficiency", 1600, 1600);
  cLTFEff->Divide(2,2,0.01,0.01);

  //TPad *padhEff[4];
  //TPad *padbEff[4];
  TPad *padhEffLTF[4];
  TPad *padbEffLTF[4];
  //TPad *padhPerc[4];
  //TPad *padbPerc[4];
  TLine *line[4];
  double xMin[4] = {0.,0.,-4.,-1.1};
  double xMax[4] = {2.2,11.,-2.,1.1};
  for (int i = 0; i < 4; ++i)
  {
    // cLTFEff->cd(i+1);
    //  padhEffLTF[i]= new TPad("pad1","pad1",0,0.3,1,1);
    //  padhEffLTF[i]->SetBottomMargin(0);
    //  padhEffLTF[i]->Draw();
    //  padbEffLTF[i] = new TPad("pad2","pad2",0,0,1,0.3);
    //  padbEffLTF[i]->SetTopMargin(0);
    //  padbEffLTF[i]->Draw();
     line[i] = new TLine(xMin[i],1.,xMax[i],1.);
     line[i]->SetLineColor(12);
     line[i]->SetLineStyle(7);


  }

  TLegend *legEff = new TLegend(0.6,0.24,0.8,0.51);
  legEff->SetFillColor(0);
	legEff->SetTextSize(gStyle->GetTextSize()*0.75);
	legEff->SetBorderSize(0);
  TLegend *legRatio = new TLegend(0.13,0.6,0.33,0.86);
  legRatio->SetFillColor(0);
	legRatio->SetTextSize(gStyle->GetTextSize()*0.75);
	legRatio->SetBorderSize(0);
  TLegend *legLTFEff = new TLegend(0.6,0.24,0.8,0.51);
  legLTFEff->SetFillColor(0);
  legLTFEff->SetTextSize(gStyle->GetTextSize()*0.8);
  legLTFEff->SetBorderSize(0);
  TLegend *legLTFPerc = new TLegend(0.6,0.24,0.8,0.51);
  legLTFPerc->SetFillColor(0);
  legLTFPerc->SetTextSize(gStyle->GetTextSize()*0.75);
  legLTFPerc->SetBorderSize(0);
  for (int i = 0; i < nFiles; ++i)
  {
  	cEff->cd(1);
  	legEff->AddEntry(hEffpT[i],legText[i]);
  	if(i==0) { hEffpT[i]->Draw("h"); hEffpT[i]->GetYaxis()->SetRangeUser(0.,1.1);}
  	else hEffpT[i]->Draw("same");
  	if(i==(nFiles-1)) legEff->Draw("same");
  	cEff->cd(2);
  	if(i==0) { hEffRatiopT[i]->Draw("h");  line[0]->Draw("same"); hEffRatiopT[i]->GetYaxis()->SetRangeUser(0.5,1.1);}
  	else hEffRatiopT[i]->Draw("same");
  	cEff->cd(3);
  	if(i==0) { hEffEta[i]->Draw("h"); hEffEta[i]->GetYaxis()->SetRangeUser(0.,1.1);}
  	else hEffEta[i]->Draw("same");
    if(i==(nFiles-1)) legEff->Draw("same");
  	cEff->cd(4);
  	if(i==0) { hEffRatioEta[i]->Draw("h");  line[2]->Draw("same"); hEffRatioEta[i]->GetXaxis()->SetRangeUser(-4.,-2.); hEffRatioEta[i]->GetYaxis()->SetRangeUser(0.5,1.1);}
  	else hEffRatioEta[i]->Draw("same");

    cAcc->cd(1);
    legRatio->AddEntry(hAcceptanceEffpT[i],legText[i]);
    if(i==0) { hAcceptanceEffpT[i]->Draw("h"); hAcceptanceEffpT[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hAcceptanceEffpT[i]->Draw("same");
    if(i==(nFiles-1)) legRatio->Draw("same");
    cAcc->cd(2);
    if(i==0) { hAcceptancepT[i]->Draw("h");  line[0]->Draw("same"); hAcceptancepT[i]->GetYaxis()->SetRangeUser(0.5,1.1);}
    else hAcceptancepT[i]->Draw("same");
    cAcc->cd(3);
    if(i==0) { hAcceptanceEffEta[i]->Draw("h"); hAcceptanceEffEta[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hAcceptanceEffEta[i]->Draw("same");
    if(i==(nFiles-1)) legRatio->Draw("same");
    cAcc->cd(4);
    if(i==0) { hAcceptanceEta[i]->Draw("h");  line[2]->Draw("same"); hAcceptanceEta[i]->GetXaxis()->SetRangeUser(-4.,-2.); hAcceptanceEta[i]->GetYaxis()->SetRangeUser(0.5,1.1);}
    else hAcceptanceEta[i]->Draw("same");

  	// cRatio->cd(1);
  	// // gPad->SetLogy();
  	// legRatio->AddEntry(hRatiopT[i],legText[i]);
  	// if(i==0) { hRatiopT[i]->Draw("h"); hRatiopT[i]->GetYaxis()->SetRangeUser(0.,500.);}
  	// else hRatiopT[i]->Draw("same");
  	// cRatio->cd(3);
  	// // gPad->SetLogy();
  	// if(i==0) { hRatiop[i]->Draw("h"); hRatiop[i]->GetYaxis()->SetRangeUser(0.,300.);}
  	// else hRatiop[i]->Draw("same");
  	// cRatio->cd(2);
  	// // gPad->SetLogy();
  	// if(i==0) { hRatioEta[i]->Draw("h"); hRatioEta[i]->GetYaxis()->SetRangeUser(0.,300.);}
  	// else hRatioEta[i]->Draw("same");
  	// cRatio->cd(4);
  	// // gPad->SetLogy();
  	// if(i==0) { hRatioVz[i]->Draw("h"); hRatioVz[i]->GetXaxis()->SetRangeUser(-1.,1.); hRatioVz[i]->GetYaxis()->SetRangeUser(0.,200.);}
  	// else hRatioVz[i]->Draw("same");
  	// if(i==(nFiles-1)) legRatio->Draw("same");

    cLTFPerc->cd(1);
    //gPad->SetLogy();
    legLTFPerc->AddEntry(hLTFPercpT[i],legText[i]);
    if(i==0) {hLTFPercpT[i]->Draw("h"); hLTFPercpT[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hLTFPercpT[i]->Draw("same");
    if(i==(nFiles-1)) legLTFPerc->Draw("same");
    cLTFPerc->cd(2);
    //gPad->SetLogy();
    if(i==0) {hPercRatiopT[i]->Draw("h");  line[0]->Draw("same"); hPercRatiopT[i]->GetYaxis()->SetRangeUser(0.5,3.);}
    else hPercRatiopT[i]->Draw("same");
    cLTFPerc->cd(3);
    //gPad->SetLogy();
    if(i==0) {hLTFPercEta[i]->Draw("h"); hLTFPercEta[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hLTFPercEta[i]->Draw("same");
    if(i==(nFiles-1)) legLTFPerc->Draw("same");
    cLTFPerc->cd(4);
    //gPad->SetLogy();
    if(i==0) { hPercRatioEta[i]->Draw("h"); line[2]->Draw("same"); hPercRatioEta[i]->GetXaxis()->SetRangeUser(-4.,-2.); hPercRatioEta[i]->GetYaxis()->SetRangeUser(0.5,3.);}
    else hPercRatioEta[i]->Draw("same");

    cLTFEff->cd(1);
    legLTFEff->AddEntry(hLTFEffpT[i],legText[i]);
    // padhEffLTF[0]->cd();
    if(i==0) { hLTFEffpT[i]->Draw("h"); hLTFEffpT[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hLTFEffpT[i]->Draw("same");
    if(i==(nFiles-1)) legLTFEff->Draw("same");
    // padbEffLTF[0]->cd();
    cLTFEff->cd(2);
    if(i==0) {hLTFRatiopT[i]->Draw("h"); line[0]->Draw("same");  hLTFRatiopT[i]->GetYaxis()->SetRangeUser(0.4,1.1); hLTFRatiopT[i]->GetYaxis()->SetTitle("Ratio to cylinder");}
    else hLTFRatiopT[i]->Draw("same");

    // cLTFEff->cd(2);
    // // padhEffLTF[1]->cd();
    // if(i==0) { hLTFEffp[i]->Draw("h"); hLTFEffp[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    // else hLTFEffp[i]->Draw("same");
    // // padbEffLTF[1]->cd();
    // if(i==0) {hLTFRatiop[i]->Draw("h");  line[1]->Draw("same"); hLTFRatiop[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    // else hLTFRatiop[i]->Draw("same");
    cLTFEff->cd(3);
    // padhEffLTF[2]->cd();
    if(i==0) { hLTFEffEta[i]->Draw("h"); hLTFEffEta[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    else hLTFEffEta[i]->Draw("same");
    if(i==(nFiles-1)) legLTFEff->Draw("same");
    // padbEffLTF[2]->cd();
    cLTFEff->cd(4);
    if(i==0) {hLTFRatioEta[i]->Draw("h");  line[2]->Draw("same"); hLTFRatioEta[i]->GetYaxis()->SetRangeUser(0.5,1.1); hLTFRatioEta[i]->GetXaxis()->SetRangeUser(-4.,-2.); hLTFRatioEta[i]->GetYaxis()->SetTitle("Ratio");}
    else hLTFRatioEta[i]->Draw("same");

    // cLTFEff->cd(4);
    // // padhEffLTF[3]->cd();
    // if(i==0) { hLTFEffVz[i]->Draw("h"); hLTFEffVz[i]->GetXaxis()->SetRangeUser(-1.,1.); hLTFEffVz[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    // else hLTFEffVz[i]->Draw("same");
    // // padbEffLTF[3]->cd();
    // if(i==0) {hLTFRatioVz[i]->Draw("h"); line[3]->Draw("same"); hLTFRatioVz[i]->GetYaxis()->SetRangeUser(0.,1.1);}
    // else hLTFRatioVz[i]->Draw("same");

  }
  if(save){
    cEff->SaveAs(Form("%s_Efficiency.eps",fileName));
    cLTFEff->SaveAs(Form("%s_LTFEfficiency.eps",fileName));
    cLTFPerc->SaveAs(Form("%s_LTFPerc.eps",fileName));
    // cRatio->SaveAs(Form("%s_AlgoRatios.eps",fileName));
  }
}

//__________________________________________________________________________
void SetHistoStyle(TH1F *h, int index) {

  if(!h) { cout << "histo not found " << endl; return;}

    h->SetMarkerColor(colors[index]);
    h->SetLineColor(colors[index]);
    h->SetLineWidth(1);
    h->SetMarkerSize(0.8);
    h->SetMarkerStyle(kFullCircle);
}
//__________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2.);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTextAlign(11);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

}
