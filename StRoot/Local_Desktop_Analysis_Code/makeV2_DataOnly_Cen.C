#include <stdio>
#include <iomanip>

void makeV2_DataOnly_Cen(const int flag=0)
{
  gROOT->Reset();
  const Double_t MassD = 1.865;
  const Double_t MassKs = 0.498;
  const Double_t MassPhi = 1.019;
  const Double_t MassLa = 1.1156;
  const Double_t MassXi = 1.3217;
  const Double_t MassOmega = 1.672;

  ifstream inData;
  // new data
  TFile *fin = new TFile("Systematics_D0vn_SL16d_2016-10-14.ME.root");
  TGraphErrors *gr_data_0_80 = (TGraphErrors *)(fin->Get("vnStat_0_80"));
  gr_data_0_80->RemovePoint(0);
  TGraphErrors *gr_data_0_80_sys = (TGraphErrors *)(fin->Get("vnSyst_0_80"));
  gr_data_0_80_sys->RemovePoint(0);

  cout << gr_data_0_80->GetN() << endl;
  cout << gr_data_0_80_sys->GetN() << endl;
  
  const Double_t scale_w_Ks = 1.0;
  const Int_t n_data_new = 8;
  //////////////////////////////////////////////
  // Read-in data points for 0-80% centrality
  //////////////////////////////////////////////
  Double_t nonflow[n_data_new] = { 0.032629 , 0.0336555 , 0.0336555 , 0.033947 , 0.0346236 , 0.0353009 , 0.0361988 , 0.0382869};  // full range 0-80%
  Double_t x_data_new[n_data_new], y_data_new[n_data_new], ye_data_new[n_data_new], yes_data_new[n_data_new], yesL_data_new[n_data_new];
  Double_t yeL_data_new[n_data_new], yeU_data_new[n_data_new];
  Double_t x_mTScaled_data_new[n_data_new], yScaled_data_new[n_data_new], yeScaled_data_new[n_data_new], yesScaled_data_new[n_data_new], yesLScaled_data_new[n_data_new];
  Double_t yeLScaled_data_new[n_data_new], yeUScaled_data_new[n_data_new];
  for(int i=0;i<gr_data_0_80->GetN();i++) {
    cout << i << endl;
    x_data_new[i] = gr_data_0_80->GetX()[i];
    y_data_new[i] = gr_data_0_80->GetY()[i];
    ye_data_new[i] = gr_data_0_80->GetEY()[i];
    yes_data_new[i] = gr_data_0_80_sys->GetEY()[i];
    yesL_data_new[i] = nonflow[i];

    x_mTScaled_data_new[i] = (sqrt(x_data_new[i]*x_data_new[i] + MassD*MassD) - MassD)/2;
    yScaled_data_new[i] = y_data_new[i]/2 * scale_w_Ks;
    yeScaled_data_new[i] = ye_data_new[i]/2 * scale_w_Ks;
    yesScaled_data_new[i] = yes_data_new[i]/2 * scale_w_Ks;
    yesLScaled_data_new[i] = yesL_data_new[i]/2 * scale_w_Ks;
  
    yeL_data_new[i] = sqrt(ye_data_new[i]**2+yes_data_new[i]**2+yesL_data_new[i]**2) * scale_w_Ks;
    yeU_data_new[i] = sqrt(ye_data_new[i]**2+yes_data_new[i]**2) * scale_w_Ks;

    yeLScaled_data_new[i] = yeL_data_new[i]/2 * scale_w_Ks;
    yeUScaled_data_new[i] = yeU_data_new[i]/2 * scale_w_Ks;
  }

  //////////////////////////////
  // different centrality bins
  //////////////////////////////
  const Int_t NCen = 3;
  TGraphErrors *gr_data_cen[NCen];
  TGraphErrors *gr_data_cen_sys[NCen];
  TGraphErrors *gr_data_mT_cen[NCen];
  const Char_t *CenName[NCen] = {"40_80","10_40","0_10"};
  for(int i=0;i<NCen;i++) {
    gr_data_cen[i] = (TGraphErrors *)(fin->Get(Form("vnStat_%s",CenName[i])));
    gr_data_cen[i]->RemovePoint(0);
    gr_data_cen_sys[i] = (TGraphErrors *)(fin->Get(Form("vnSyst_%s",CenName[i])));
    gr_data_cen_sys[i]->RemovePoint(0);
  }

  const Int_t n_data_cen = 8;
  Double_t nonflow_cen[NCen][n_data_cen] = {{ 0.115456 , 0.118967 , 0.118967 , 0.120065 , 0.122454 , 0.124808 , 0.128025 , 0.135539},
					    { 0.0210766 , 0.0217398 , 0.0217398 , 0.0219309 , 0.0223667 , 0.0228055 , 0.0233849 , 0.0247439},
					    { 0.0202398 , 0.0208774 , 0.0208774 , 0.0210541 , 0.0214741 , 0.021893 , 0.0224473 , 0.0237127}};

  Double_t x_data_cen[NCen][n_data_cen], y_data_cen[NCen][n_data_cen], ye_data_cen[NCen][n_data_cen], yes_data_cen[NCen][n_data_cen], yesL_data_cen[NCen][n_data_cen];
  Double_t yeL_data_cen[NCen][n_data_cen], yeU_data_cen[NCen][n_data_cen];
  Double_t x_mTScaled_data_cen[NCen][n_data_cen], yScaled_data_cen[NCen][n_data_cen], yeScaled_data_cen[NCen][n_data_cen], yesScaled_data_cen[NCen][n_data_cen], yesLScaled_data_cen[NCen][n_data_cen];
  Double_t yeLScaled_data_cen[NCen][n_data_cen], yeUScaled_data_cen[NCen][n_data_cen];
  for(int ic=0;ic<NCen;ic++) {
    for(int i=0;i<gr_data_cen[ic]->GetN();i++) {
      x_data_cen[ic][i] = gr_data_cen[ic]->GetX()[i];
      y_data_cen[ic][i] = gr_data_cen[ic]->GetY()[i];
      ye_data_cen[ic][i] = gr_data_cen[ic]->GetEY()[i];
      yes_data_cen[ic][i] = gr_data_cen_sys[ic]->GetEY()[i];
      yesL_data_cen[ic][i] = nonflow_cen[ic][i];
      
      x_mTScaled_data_cen[ic][i] = (sqrt(x_data_cen[ic][i]*x_data_cen[ic][i] + MassD*MassD) - MassD)/2;
      yScaled_data_cen[ic][i] = y_data_cen[ic][i]/2 * scale_w_Ks;
      yeScaled_data_cen[ic][i] = ye_data_cen[ic][i]/2 * scale_w_Ks;
      yesScaled_data_cen[ic][i] = yes_data_cen[ic][i]/2 * scale_w_Ks;
      yesLScaled_data_cen[ic][i] = yesL_data_cen[ic][i]/2 * scale_w_Ks;
      
      yeL_data_cen[ic][i] = sqrt(ye_data_cen[ic][i]**2+yes_data_cen[ic][i]**2+yesL_data_cen[ic][i]**2) * scale_w_Ks;
      yeU_data_cen[ic][i] = sqrt(ye_data_cen[ic][i]**2+yes_data_cen[ic][i]**2) * scale_w_Ks;
      
      yeLScaled_data_cen[ic][i] = yeL_data_cen[ic][i]/2 * scale_w_Ks;
      yeUScaled_data_cen[ic][i] = yeU_data_cen[ic][i]/2 * scale_w_Ks;      
    }
    gr_data_mT_cen[ic] = new TGraphErrors(n_data_cen, x_mTScaled_data_cen[ic], yScaled_data_cen[ic], 0, yeScaled_data_cen[ic]);
  }

  cout << " Read-in D0 data done ..." << endl;
  
 /* const Int_t n_ks_cen = 19;
  const Int_t n_la_cen = 18;

  Double_t x_ks_cen[n_ks_cen], y_ks_cen[NCen][n_ks_cen], ye_ks_cen[NCen][n_ks_cen], yes_ks_cen[NCen][n_ks_cen];
  Double_t x_la_cen[n_la_cen], y_la_cen[NCen][n_la_cen], ye_la_cen[NCen][n_la_cen], yes_la_cen[NCen][n_la_cen];

  Double_t x_mTScaled_ks_cen[n_ks_cen], yScaled_ks_cen[NCen][n_ks_cen], yeScaled_ks_cen[NCen][n_ks_cen], yesScaled_ks_cen[NCen][n_ks_cen];
  Double_t x_mTScaled_la_cen[n_la_cen], yScaled_la_cen[NCen][n_la_cen], yeScaled_la_cen[NCen][n_la_cen], yesScaled_la_cen[NCen][n_la_cen];


  inData.open("Run14/ks_v2_cen_PRC77.txt");
  for(int i=0;i<n_ks_cen;i++) {
    double a, b, c;
    inData >> x_ks_cen[i] >> a >> b >> c >> y_ks_cen[0][i] >> ye_ks_cen[0][i] >> yes_ks_cen[0][i] >> y_ks_cen[1][i] >> ye_ks_cen[1][i] >> yes_ks_cen[1][i] >> y_ks_cen[2][i] >> ye_ks_cen[2][i] >> yes_ks_cen[2][i];

    x_mTScaled_ks_cen[i] = (sqrt(x_ks_cen[i]*x_ks_cen[i]+MassKs*MassKs)-MassKs)/2;
    for(int ic=0;ic<NCen;ic++) {
      yScaled_ks_cen[ic][i] = y_ks_cen[ic][i]/2.;
      yeScaled_ks_cen[ic][i] = ye_ks_cen[ic][i]/2.;
      yesScaled_ks_cen[ic][i] = yes_ks_cen[ic][i]/2.;
    }

  }
  inData.close();

  TGraphErrors *gr_ks_cen[NCen], *gr_ks_mT_cen[NCen];
  for(int ic=0;ic<NCen;ic++) {
    gr_ks_cen[ic] = new TGraphErrors(n_ks_cen, x_ks_cen, y_ks_cen[ic], 0, ye_ks_cen[ic]);
    gr_ks_mT_cen[ic] = new TGraphErrors(n_ks_cen, x_mTScaled_ks_cen, yScaled_ks_cen[ic], 0, yeScaled_ks_cen[ic]);
    //    gr_ks_mT_cen[ic]->Print();
  }
  cout << " Read-in Ks data points done ... " << endl;
  
  inData.open("Run14/lambda_v2_cen_PRC77.txt");
  for(int i=0;i<n_la_cen;i++) {
    double a, b, c;
    inData >> x_la_cen[i] >> a >> b >> c >> y_la_cen[0][i] >> ye_la_cen[0][i] >> yes_la_cen[0][i] >> y_la_cen[1][i] >> ye_la_cen[1][i] >> yes_la_cen[1][i] >> y_la_cen[2][i] >> ye_la_cen[2][i] >> yes_la_cen[2][i];

    x_mTScaled_la_cen[i] = (sqrt(x_la_cen[i]*x_la_cen[i]+MassLa*MassLa)-MassLa)/3;
    for(int ic=0;ic<NCen;ic++) {
      yScaled_la_cen[ic][i] = y_la_cen[ic][i]/3.;
      yeScaled_la_cen[ic][i] = ye_la_cen[ic][i]/3.;
      yesScaled_la_cen[ic][i] = yes_la_cen[ic][i]/3.;
    }

  }
  inData.close();
  TGraphErrors *gr_la_cen[NCen], *gr_la_mT_cen[NCen];
  for(int ic=0;ic<NCen;ic++) {
    gr_la_cen[ic] = new TGraphErrors(n_la_cen, x_la_cen, y_la_cen[ic], 0, ye_la_cen[ic]);
    gr_la_mT_cen[ic] = new TGraphErrors(n_la_cen, x_mTScaled_la_cen, yScaled_la_cen[ic], 0, yeScaled_la_cen[ic]);
  }
  cout << " Read-in Lambda data points done ... " << endl;
  
  // Xi data points, format is a bit different
  const Int_t n_xi_cenMax = 9;
  Int_t n_xi_cen[NCen];
  Double_t x_xi_cen[NCen][n_xi_cenMax], y_xi_cen[NCen][n_xi_cenMax], ye_xi_cen[NCen][n_xi_cenMax];
  Double_t x_mTScaled_xi_cen[NCen][n_xi_cenMax], yScaled_xi_cen[NCen][n_xi_cenMax], yeScaled_xi_cen[NCen][n_xi_cenMax];
  TGraphErrors *gr_xi_cen[NCen];
  TGraphErrors *gr_xi_mT_cen[NCen];
  for(int ic=0;ic<NCen;ic++) {
    gr_xi_cen[ic] = new TGraphErrors(Form("Run14/xi_v2_%s.txt",CenName[ic]),"%lg %lg %lg");

    n_xi_cen[ic] = gr_xi_cen[ic]->GetN();
    for(int i=0;i<gr_xi_cen[ic]->GetN();i++) {
      x_xi_cen[ic][i] = gr_xi_cen[ic]->GetX()[i];
      y_xi_cen[ic][i] = gr_xi_cen[ic]->GetY()[i];
      ye_xi_cen[ic][i] = gr_xi_cen[ic]->GetEY()[i];

      x_mTScaled_xi_cen[ic][i] = (sqrt(x_xi_cen[ic][i]*x_xi_cen[ic][i]+MassXi*MassXi)-MassXi)/3;
      yScaled_xi_cen[ic][i] = y_xi_cen[ic][i]/3.;
      yeScaled_xi_cen[ic][i] = ye_xi_cen[ic][i]/3.;      
    }
    
    gr_xi_mT_cen[ic] = new TGraphErrors(n_xi_cen[ic], x_mTScaled_xi_cen[ic], yScaled_xi_cen[ic], 0, yeScaled_xi_cen[ic]);
  }
  cout << " Read-in Xi data points done ... " << endl;

  const Int_t n_phi_cenMax = 7;
  Int_t n_phi_cen[NCen];
  Double_t x_phi_cen[NCen][n_phi_cenMax], y_phi_cen[NCen][n_phi_cenMax], ye_phi_cen[NCen][n_phi_cenMax], yes_phi_cen[NCen][n_phi_cenMax];
  Double_t x_mTScaled_phi_cen[NCen][n_phi_cenMax], yScaled_phi_cen[NCen][n_phi_cenMax], yeScaled_phi_cen[NCen][n_phi_cenMax], yesScaled_phi_cen[NCen][n_phi_cenMax];
  TGraphErrors *gr_phi_cen[NCen];
  TGraphErrors *gr_phi_mT_cen[NCen];
  const Char_t *CenName_Phi[NCen] = {"40_80","10_40","0_5"};
  for(int ic=0;ic<NCen;ic++) {
    gr_phi_cen[ic] = new TGraphErrors(Form("Run14/phi_v2_%s_PRC.txt",CenName_Phi[ic]),"%lg %lg %lg");

    n_phi_cen[ic] = gr_phi_cen[ic]->GetN();
    inData.open(Form("Run14/phi_v2_%s_PRC.txt",CenName_Phi[ic]));
    for(int i=0;i<gr_phi_cen[ic]->GetN();i++) {
      x_phi_cen[ic][i] = gr_phi_cen[ic]->GetX()[i];
      y_phi_cen[ic][i] = gr_phi_cen[ic]->GetY()[i];
      ye_phi_cen[ic][i] = gr_phi_cen[ic]->GetEY()[i];

      double a, b, c, d;
      inData >> a >> b >> c >> d;
      yes_phi_cen[ic][i] = d;

      x_mTScaled_phi_cen[ic][i] = (sqrt(x_phi_cen[ic][i]*x_phi_cen[ic][i]+MassPhi*MassPhi)-MassPhi)/2.;
      yScaled_phi_cen[ic][i] = y_phi_cen[ic][i]/2.;
      yeScaled_phi_cen[ic][i] = ye_phi_cen[ic][i]/2.;
      yesScaled_phi_cen[ic][i] = yes_phi_cen[ic][i]/2.;
    }
    
    gr_phi_mT_cen[ic] = new TGraphErrors(n_phi_cen[ic], x_mTScaled_phi_cen[ic], yScaled_phi_cen[ic], 0, yeScaled_phi_cen[ic]);
  }
  cout << " Read-in Phi data points done ... " << endl;*/
  
  

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,1600,900);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTitleBorderSize(0);
  c1->SetFillColor(10);
  c1->SetFillStyle(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameFillColor(10);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  
  /*TPad* p1 = new TPad("p1","p1",0.,0.5,1.0,1.0);
  p1->SetFillColor(10);
  p1->SetFillStyle(0);
  p1->SetBorderMode(0);
  p1->SetBorderSize(0);
  p1->SetFrameFillColor(10);
  p1->SetFrameFillStyle(0);
  p1->SetFrameBorderMode(0);
  //p1->SetLogy();
  p1->SetGridx(0);
  p1->SetGridy(0);
  p1->SetLeftMargin(0.16);
  p1->SetBottomMargin(0.15);
  p1->SetTopMargin(0.02);
  p1->SetRightMargin(0.02);
  p1->Draw();
  //p1->cd();*/
  
  double x1 = 0.0;
  double x2 = 6.8;
  double y1 = -0.05;
  double y2 = 0.38;
  TH1 *h0 = new TH1D("h0","",1,x1, x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(.9);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetXaxis()->SetLabelOffset(0.005);
  h0->GetXaxis()->SetLabelSize(0.04);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(505);
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("Anisotropy Parameter, v_{2}");
  h0->GetYaxis()->SetTitleOffset(.9);
  h0->GetYaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetLabelOffset(0.005);
  h0->GetYaxis()->SetLabelSize(0.04);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetTitleFont(42);
  h0->Draw("c");

  
  double pt[2] = {1.766, 4.459};
  double pt1Bin[1] = {1.866};
  double v2[2] = {0.08329931, 0.126231279};
  double eptLow[2] = {.766, .459};
  double eptHigh[2] = {2.234, 2.341};
  double eptSyst[2] = { 0.0, 0.0};
  double eptSyst1Bin[1] = { 0.0};
  double ev2Low[2] = { 0.00736501, 0.012241649};
  double ev2High[2] = { 0.00736501, 0.012241649};
  double ev2Syst[2] = {.0124948965, .02145931743 };
  
  
  TGraphAsymmErrors *AlexData = new TGraphAsymmErrors(2, pt, v2, eptLow, eptHigh, ev2Low, ev2High);
  
  AlexData->SetMarkerColor(2);
  AlexData->SetMarkerStyle(21);
  AlexData->SetMarkerSize(2);
  AlexData->SetLineColor(2);
  AlexData->SetLineWidth(3);
  
  TGraphErrors *AlexDataSyst = new TGraphErrors(2, pt, v2, eptSyst, ev2Syst);
  
  AlexDataSyst->SetMarkerColor(2);
  AlexDataSyst->SetMarkerStyle(21);
  AlexDataSyst->SetMarkerSize(1.1);
  AlexDataSyst->SetLineColor(2);
  AlexDataSyst->SetLineWidth(33);
  
  
  
  // 8 data points
   
 
  
  double weights[6] = {0.4034, .3058, .1696, .0778, .0316, .0119};

  double avgV2[1] = {0};
  double avgV2StatErrorSquared[1] = {0};
  double avgV2SystErrorSquared[1] = {0};
  double avgV2StatError[1] = {0};
  double avgV2SystError[1] = {0};
  
  for(int i = 0; i < 6; i++){

    cout << "V2: " << y_data_cen[1][i] << endl;
    
    avgV2[0] = avgV2[0] + (weights[i]*y_data_cen[1][i]);
    
    avgV2StatErrorSquared[0] = avgV2StatErrorSquared[0] + (weights[i]*weights[i])*(ye_data_cen[1][i]*ye_data_cen[1][i]);
    avgV2SystErrorSquared[0] = avgV2SystErrorSquared[0] + (weights[i]*weights[i])*(yes_data_cen[1][i]*yes_data_cen[1][i]);
    
  }
    
  avgV2StatError[0] = TMath::Sqrt(avgV2StatErrorSquared[0]);
  avgV2SystError[0] = TMath::Sqrt(avgV2SystErrorSquared[0]);
    
  TGraphErrors *AvgV2Data = new TGraphErrors(1, pt1Bin, avgV2, eptSyst1Bin, avgV2StatError);
  
  AvgV2Data->SetMarkerColor(4);
  AvgV2Data->SetMarkerStyle(22);
  AvgV2Data->SetMarkerSize(2.5);
  AvgV2Data->SetLineColor(4);
  //AvgV2Data->SetLineWidth(15);
  
  TGraphErrors *AvgV2DataSyst = new TGraphErrors(1, pt1Bin, avgV2, eptSyst1Bin, avgV2SystError);
  
  AvgV2DataSyst->SetMarkerColor(4);
  AvgV2DataSyst->SetMarkerStyle(22);
  AvgV2DataSyst->SetMarkerSize(2.5);
  AvgV2DataSyst->SetLineColor(4);
  AvgV2DataSyst->SetLineWidth(29);
  
  
  
  /*TH1 *hAlex1 = new TH1D("","", 7, 0, 6.8);
  hAlex1->SetMarkerStyle(20);
  hAlex1->SetMarkerColor(2);
  hAlex1->SetMarkerSize(1.5);
  //hAlex->SetBinContent(3, 0.0801251);
  //hAlex->SetBinError(3, 0.00244519);
  
  hAlex1->SetBinContent(1, 0.122872);
  hAlex1->SetBinError(1, 0.0149516);
  
 // hAlex1->SetBinContent(5, 0.0789727);
  //hAlex1->SetBinError(5, 0.00444987);
  
  TH1 *hAlex2 = new TH1D("","", 1, 1, 4);
  hAlex2->SetMarkerStyle(20);
  hAlex2->SetMarkerColor(2);
  hAlex2->SetMarkerSize(1.5);
  hAlex2->SetBinContent(1, 0.0801251);
  hAlex2->SetBinError(1, 0.00244519);
  
  TH1 *hAlex3 = new TH1D("","", 1, 4, 6.8);
  hAlex3->SetMarkerStyle(20);
  hAlex3->SetMarkerColor(2);
  hAlex3->SetMarkerSize(1.5);
  hAlex3->SetBinContent(1, 0.0789727);
  hAlex3->SetBinError(1, 0.00444987);
  //////////////////////////////////////////////////////////////////
  TH1 *hAlex1Syst = new TH1D("","", 7, 0, 6.8);
  hAlex1Syst->SetMarkerStyle(20);
  hAlex1Syst->SetMarkerColor(2);
  hAlex1Syst->SetMarkerSize(1.5);
  //hAlex1Syst->SetEndErrorSize(2);
  hAlex1Syst->SetBinContent(1, 0.122872);
  hAlex1Syst->SetBinError(1, 0.02);
  
 
  TH1 *hAlex2Syst = new TH1D("","", 1, 1, 4);
  hAlex2Syst->SetMarkerStyle(20);
  hAlex2Syst->SetMarkerColor(2);
  hAlex2Syst->SetMarkerSize(1.5);
  hAlex2Syst->SetBinContent(1, 0.0801251);
  hAlex2Syst->SetBinError(1, 0.006172166);
  
  TH1 *hAlex3Syst = new TH1D("","", 1, 4, 6.8);
  hAlex3Syst->SetMarkerStyle(20);
  hAlex3Syst->SetMarkerColor(2);
  hAlex3Syst->SetMarkerSize(1.5);
  hAlex3Syst->SetBinContent(1, 0.0789727);
  hAlex3Syst->SetBinError(1, 0.004869332);*/
  
 
  
  /*TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(2);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth(2);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth(2);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(2);
  l4->Draw("same");

  TLine *l0 = new TLine(x1, 0, x2, 0);
  l0->SetLineWidth(2);
  l0->SetLineStyle(2);
  l0->Draw("same");*/

  /*for(int i=0;i<n_ks_cen;i++) {
    double x1 = x_ks_cen[i]-0.08;
    double x2 = x_ks_cen[i]+0.08;
    double y1 = y_ks_cen[1][i]-yes_ks_cen[1][i];
    double y2 = y_ks_cen[1][i]+yes_ks_cen[1][i];
    
    TLine *la = new TLine(x1, y1, x1, y1+0.003);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.003);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.003);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.003);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }*/
  
  /*gr_ks_cen[1]->Print();
  gr_ks_cen[1]->SetMarkerStyle(25);
  gr_ks_cen[1]->SetMarkerColor(1);
  gr_ks_cen[1]->SetMarkerSize(1.2);
  gr_ks_cen[1]->SetLineColor(1);
  gr_ks_cen[1]->SetLineWidth(2);
  gr_ks_cen[1]->Draw("p");*/

  
 /* for(int i=0;i<n_la_cen;i++) {
    double x1 = x_la_cen[i]-0.08;
    double x2 = x_la_cen[i]+0.08;
    double y1 = y_la_cen[1][i]-yes_la_cen[1][i];
    double y2 = y_la_cen[1][i]+yes_la_cen[1][i];
    
    TLine *la = new TLine(x1, y1, x1, y1+0.003);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.003);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.003);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.003);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }
  
  gr_la_cen[1]->Print();
  gr_la_cen[1]->SetMarkerStyle(24);
  gr_la_cen[1]->SetMarkerColor(1);
  gr_la_cen[1]->SetMarkerSize(1.2);
  gr_la_cen[1]->SetLineColor(1);
  gr_la_cen[1]->SetLineWidth(2);
  gr_la_cen[1]->Draw("p");


  gr_xi_cen[1]->Print();
  gr_xi_cen[1]->SetMarkerStyle(26);
  gr_xi_cen[1]->SetMarkerColor(1);
  gr_xi_cen[1]->SetMarkerSize(1.2);
  gr_xi_cen[1]->SetLineColor(1);
  gr_xi_cen[1]->SetLineWidth(2);
  gr_xi_cen[1]->Draw("p");
  
  gr_phi_cen[1]->Print();
  gr_phi_cen[1]->SetMarkerStyle(25);
  gr_phi_cen[1]->SetMarkerColor(1);
  gr_phi_cen[1]->SetMarkerSize(1.2);
  gr_phi_cen[1]->SetLineColor(1);
  gr_phi_cen[1]->SetLineWidth(2);
  //  gr_phi_cen[1]->Draw("p");*/

    for(int i=0;i<n_data_cen[1];i++) {
    double x1 = x_data_cen[1][i]-0.08;
    double x2 = x_data_cen[1][i]+0.08;
    double y1 = y_data_cen[1][i]-yes_data_cen[1][i];
    double y2 = y_data_cen[1][i]+yes_data_cen[1][i];

    double y3 = y_data_cen[1][i] - yesL_data_cen[1][i];
    double y4 = y_data_cen[1][i];
    TBox *box = new TBox(x1, y3, x2, y4);
    box->SetLineColor(16);
    box->SetFillColor(16);
    box->Draw("same");
    
    TLine *la = new TLine(x1, y1, x1, y1+0.003);
    la->SetLineColor(4);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.003);
    lb->SetLineColor(4);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.003);
    lc->SetLineColor(4);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.003);
    ld->SetLineColor(4);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->SetLineColor(4);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->SetLineColor(4);
    lf->Draw("same");
  }
   
   gr_data_cen[1]->SetMarkerStyle(20);
   //gr_data_cen[1]->SetMarkerColorAlpha(4, 0.35);
   gr_data_cen[1]->SetMarkerColor(4);
   
   gr_data_cen[1]->SetMarkerSize(1.5);
   gr_data_cen[1]->SetLineWidth(2);
   gr_data_cen[1]->SetLineColor(4);
   gr_data_cen[1]->Draw("p");
  
  /////MY DATA HERE////////////////////////
    for(int i = 0; i < 2; i++) {
        double x1 = pt[i]-0.08;
        double x2 = pt[i]+0.08;
        double y1 = v2[i]-ev2Syst[i];
        double y2 = v2[i]+ev2Syst[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.003);
        la->SetLineColor(2);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.003);
        lb->SetLineColor(2);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.003);
        lc->SetLineColor(2);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.003);
        ld->SetLineColor(2);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(2);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(2);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
  

  //TLatex *tex = new TLatex(5.5, 0.28, "10-40%");
  //tex->SetTextFont(42);
  //tex->SetTextSize(0.055);
  //tex->Draw();
  
  TLegend *leg = new TLegend(0.1, 0.72, 0.55, 0.9);
  leg->SetFillStyle(0);
  leg->SetLineStyle(4000);
  leg->SetLineColor(10);
  //leg->SetLineWidth(1.0);
  leg->SetBorderSize(0.0);
  leg->SetTextSize(0.045);
  leg->AddEntry(gr_data_cen[1], "D^{0} Event Plane, 10-40% ", "p");
  leg->AddEntry(AlexData, "D^{0}-Hadron Angular Corr, 20-50%", "p");
  //leg->AddEntry(AvgV2Data, "D^{0} Event Plane Avg. of 1 GeV/c < p_{t} < 4 GeV/c", "p");
  //leg->AddEntry(gr_xi_cen[1], "#Xi^{-}", "p");
  //leg->AddEntry(gr_la_cen[1], "#Lambda", "p");
  //leg->AddEntry(gr_ks_cen[1], "K_{S}", "p");
  //  leg->AddEntry(gr_phi_cen[1], "#phi", "p");
  leg->Draw();
  
  AlexData->Draw("SAME P");
  //AlexDataSyst->Draw("SAME []");
  //AvgV2Data->Draw("SAME P");
  //AvgV2DataSyst->Draw("SAME []");
  
  TString starPrelim = "STAR Preliminary";
  
  TPaveText *starPrelimTextBox = new TPaveText(0.7, 0.7, .8, .77, "NB NDC");
  starPrelimTextBox->SetFillColor(0);
  starPrelimTextBox->AddText(starPrelim);
  starPrelimTextBox->GetLine(0)->SetTextSize(.055);
  starPrelimTextBox->GetLine(0)->SetTextColor(2);
  starPrelimTextBox->Draw("SAME");

  TLatex *tex = new TLatex(4.2, 0.345, "STAR  Au+Au @ 200 GeV");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->Draw();
  //hAlex2->Draw("SAME");
  //hAlex3->Draw("SAME");
  //hAlex1Syst->Draw("SAME E1");
  //hAlex2Syst->Draw("SAME E1");
  //hAlex3Syst->Draw("SAME E1");
  
 // tex = new TLatex(0.2, 0.34, "a)");
  //tex->SetTextFont(42);
  //tex->SetTextSize(0.065);
  //tex->Draw();

  //p1->Modified();
  //c1->Update();
  //c1->cd();

  /*TPad* p2 = new TPad("p2","",0.,0.,1.0,0.5);
  p2->SetFillColor(10);
  p2->SetFillStyle(0);
  p2->SetBorderMode(0);
  p2->SetBorderSize(0);
  p2->SetFrameFillColor(10);
  p2->SetFrameFillStyle(0);
  p2->SetFrameBorderMode(0);
  //p2->SetLogy();
  p2->SetGridx(0);
  p2->SetGridy(0);
  p2->SetLeftMargin(0.16);
  p2->SetBottomMargin(0.18);
  p2->SetTopMargin(0.01);
  p2->SetRightMargin(0.02);
  p2->Draw();
  p2->cd();*/

 /* double x1 = 0.;
  double x2 = 2.8;
  double y1 = -0.025;
  double y2 = 0.18;
  TH1 *h0 = new TH1D("h0","",1,x1, x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("(m_{T} - m_{0}) / n_{q} (GeV/c^{ 2})");
  h0->GetXaxis()->SetTitleOffset(1.1);
  h0->GetXaxis()->SetTitleSize(0.065);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(505);
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("Anisotropy Parameter, v_{2} / n_{q}");
  h0->GetYaxis()->SetTitleOffset(1.1);
  h0->GetYaxis()->SetTitleSize(0.065);
  h0->GetYaxis()->SetLabelOffset(0.015);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetTitleFont(42);
  h0->Draw("c");

  TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(2);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth(2);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth(2);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(2);
  l4->Draw("same");

  TLine *l0 = new TLine(x1, 0, x2, 0);
  l0->SetLineWidth(2);
  l0->SetLineStyle(2);
  l0->Draw("same");*/

  /*for(int i=0;i<n_ks_cen;i++) {
    double x1 = x_mTScaled_ks_cen[i]-0.04;
    double x2 = x_mTScaled_ks_cen[i]+0.04;
    double y1 = yScaled_ks_cen[1][i]-yesScaled_ks_cen[1][i];
    double y2 = yScaled_ks_cen[1][i]+yesScaled_ks_cen[1][i];
    
    TLine *la = new TLine(x1, y1, x1, y1+0.0015);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.0015);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.0015);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.0015);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }
  
  gr_ks_mT_cen[1]->Print();
  gr_ks_mT_cen[1]->SetMarkerStyle(25);
  gr_ks_mT_cen[1]->SetMarkerColor(1);
  gr_ks_mT_cen[1]->SetMarkerSize(1.2);
  gr_ks_mT_cen[1]->SetLineColor(1);
  gr_ks_mT_cen[1]->SetLineWidth(2);
  gr_ks_mT_cen[1]->Draw("p");

  
  for(int i=0;i<n_la_cen;i++) {
    double x1 = x_mTScaled_la_cen[i]-0.04;
    double x2 = x_mTScaled_la_cen[i]+0.04;
    double y1 = yScaled_la_cen[1][i]-yesScaled_la_cen[1][i];
    double y2 = yScaled_la_cen[1][i]+yesScaled_la_cen[1][i];
    
    TLine *la = new TLine(x1, y1, x1, y1+0.0015);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.0015);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.0015);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.0015);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }
  
  gr_la_mT_cen[1]->Print();
  gr_la_mT_cen[1]->SetMarkerStyle(24);
  gr_la_mT_cen[1]->SetMarkerColor(1);
  gr_la_mT_cen[1]->SetMarkerSize(1.2);
  gr_la_mT_cen[1]->SetLineColor(1);
  gr_la_mT_cen[1]->SetLineWidth(2);
  gr_la_mT_cen[1]->Draw("p");


  gr_xi_mT_cen[1]->Print();
  gr_xi_mT_cen[1]->SetMarkerStyle(26);
  gr_xi_mT_cen[1]->SetMarkerColor(1);
  gr_xi_mT_cen[1]->SetMarkerSize(1.2);
  gr_xi_mT_cen[1]->SetLineColor(1);
  gr_xi_mT_cen[1]->SetLineWidth(2);
  gr_xi_mT_cen[1]->Draw("p");
  
  gr_phi_mT_cen[1]->Print();
  gr_phi_mT_cen[1]->SetMarkerStyle(25);
  gr_phi_mT_cen[1]->SetMarkerColor(1);
  gr_phi_mT_cen[1]->SetMarkerSize(1.2);
  gr_phi_mT_cen[1]->SetLineColor(1);
  gr_phi_mT_cen[1]->SetLineWidth(2);
  //  gr_phi_mT_cen[1]->Draw("p");*/

/*  for(int i=0;i<n_data_cen[1];i++) {
    double x1 = x_mTScaled_data_cen[1][i]-0.04;
    double x2 = x_mTScaled_data_cen[1][i]+0.04;
    double y1 = yScaled_data_cen[1][i]-yesScaled_data_cen[1][i];
    double y2 = yScaled_data_cen[1][i]+yesScaled_data_cen[1][i];

    double y3 = yScaled_data_cen[1][i] - yesLScaled_data_cen[1][i];
    double y4 = yScaled_data_cen[1][i];
    TBox *box = new TBox(x1, y3, x2, y4);
    box->SetLineColor(16);
    box->SetFillColor(16);
    box->Draw("same");
    
    TLine *la = new TLine(x1, y1, x1, y1+0.0015);
    la->SetLineColor(4);
    la->Draw("same");
    TLine *lb = new TLine(x2, y1, x2, y1+0.0015);
    lb->SetLineColor(4);
    lb->Draw("same");
    TLine *lc = new TLine(x1, y2, x1, y2-0.0015);
    lc->SetLineColor(4);
    lc->Draw("same");
    TLine *ld = new TLine(x2, y2, x2, y2-0.0015);
    ld->SetLineColor(4);
    ld->Draw("same");
    TLine *le = new TLine(x1, y1, x2, y1);
    le->SetLineColor(4);
    le->SetLineWidth(2);
    le->Draw("same");
    TLine *lf = new TLine(x1, y2, x2, y2);
    lf->SetLineColor(4);
    lf->SetLineWidth(2);
    lf->Draw("same");
  }*/
   
   //gr_data_mT_cen[1]->SetMarkerStyle(20);
   //gr_data_mT_cen[1]->SetMarkerColor(4);
   //gr_data_mT_cen[1]->SetMarkerSize(1.5);
   //gr_data_mT_cen[1]->SetLineWidth(2);
   //gr_data_mT_cen[1]->SetLineColor(4);
   //gr_data_mT_cen[1]->Draw("p");
  
  
  //TLegend *leg = new TLegend(0.22, 0.66, 0.5, 0.96);
  //leg->SetFillStyle(0);
  //leg->SetLineStyle(4000);
  //leg->SetLineColor(-1);
  //leg->SetLineWidth(0);
  //leg->SetTextSize(0.06);
  //leg->AddEntry(gr_data_mT_cen[1], "D^{0}", "p");
  //leg->AddEntry(gr_xi_mT_cen[1], "#Xi^{-}", "p");
  //leg->AddEntry(gr_la_mT_cen[1], "#Lambda", "p");
  //leg->AddEntry(gr_ks_mT_cen[1], "K_{S}", "p");  
  //  leg->AddEntry(gr_phi_mT_cen[1], "#phi", "p");  
  //leg->Draw();

  //TLatex *tex = new TLatex(1.2, 0.155, "STAR  Au+Au @ 200 GeV");
 // tex->SetTextFont(42);
  //tex->SetTextSize(0.065);
  //tex->Draw();

  //TLatex *tex = new TLatex(2.27, 0.13, "10-40%");
  //tex->SetTextFont(42);
  //tex->SetTextSize(0.065);
  //tex->Draw();


  //tex = new TLatex(0.08, 0.155, "b)");
  //tex->SetTextFont(42);
  //tex->SetTextSize(0.065);
  //tex->Draw();

  //p2->Modified();
  //c1->Update();
  //c1->cd();

  //c1->SaveAs("v2CompareWithData.eps");
  c1->SaveAs("v2_plot.png");

}
