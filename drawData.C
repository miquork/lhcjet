//TH1D *graphToHist(TGraphAsymmErrors
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TMultiGraph.h"
#include "TMath.h"

#include "tools.C"
#include "tdrstyle_gen15.C"

#include <vector>

using namespace std;

void drawData() {

  setTDRStyle();
  
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("lhcdata.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  //TH1D *h = new TH1D("h",";p_{T} (GeV);Data / fit",100,50,2000.);
  TH1D *h = new TH1D("h",";p_{T} (GeV);Data / F(p_{T},|#eta|)",
		     //#;N=1.2e14,#alpha=5,#beta=10,|#eta|=0)",
		       100,50,2000.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(1.5);//1.3);
  h->SetMinimum(0.5);//0.6);

  //cmsText = "CMS+ATLAS";
  cmsText = "HEPData (1406.0324, 1410.8857)";
  lumi_7TeV = "CMS 5.0 fb^{-1}, ATLAS 4.5 fb^{-1}";
  extraText = "N=1.2#times10^{14}, #alpha=5, #beta=10, |#eta|=0";
  TCanvas *c1 = tdrCanvas("c1",h,1,11,kSquare);
  // 1.1e14/fb? 1e11/pb, 1e8/nb, 1e5/mub, 1e2/mb => ~total xsec at 1 GeV?
  // Offset per unit area in QCD is 2 GeV.

  // Kinematic suppression vs eta
  // TF1 *f0 = new TF1("f0","pow(1-min(2*[0]*cosh(x),6500.)/7000.,10)",0,5)
  // Drop is steeper at higher eta, so should do 2D integration in pT and eta?

  TLegend *leg = tdrLeg(0.20,0.58,0.50,0.78);
  leg->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.40,0.17,"|y| < 0.5, anti-k_{T} R = 0.4--0.7");

  vector<string> vd;
  vd.push_back("cms_r07_y00x05");
  vd.push_back("atlas_r06_y00x05");
  vd.push_back("cms_r05_y00x05");
  vd.push_back("atlas_r04_y00x05");

  map<string, int> color;
  color["cms_r07_y00x05"] = kBlue;
  color["atlas_r06_y00x05"] = kGreen+2;
  color["cms_r05_y00x05"] = kOrange+2;
  color["atlas_r04_y00x05"] = kRed;
  
  map<string, int> marker;
  marker["cms_r07_y00x05"] = kFullCircle;
  marker["atlas_r06_y00x05"] = kOpenCircle;
  marker["cms_r05_y00x05"] = kFullCircle;
  marker["atlas_r04_y00x05"] = kOpenCircle;

  map<string, const char*> label;
  label["cms_r07_y00x05"] = "CMS R=0.7";
  label["atlas_r06_y00x05"] = "ATLAS R=0.6";
  label["cms_r05_y00x05"] = "CMS R=0.5";
  label["atlas_r04_y00x05"] = "ATLAS R=0.4";
  
  // CMS lumi 5.0/fb, ATLAS 4.5/fb, 11% difference
  // We get better agreement assuming they are the same
  // => Correct CMS up 5%, ATLAS down 5%
  double kfactor0 = 1;//1./1.05;
  map<string, double> kfactor;
  kfactor["cms_r07_y00x05"] = 1;//1.05;
  kfactor["atlas_r06_y00x05"] = 1;//0.90;
  kfactor["cms_r05_y00x05"] = 1;//1.05;
  kfactor["atlas_r04_y00x05"] = 1;//0.90;

  
  TMultiGraph *mg = new TMultiGraph();
  vector<TGraphAsymmErrors*> vg(vd.size());
  vector<TGraphAsymmErrors*> vt(vd.size());
  vector<TGraphErrors*> vr(vd.size());
  for (int i = 0; i != vd.size(); ++i) {
    
    TGraphAsymmErrors *g = (TGraphAsymmErrors*)f->Get((vd[i]+"_stat").c_str());
    assert(g);
    TGraphAsymmErrors *gt = (TGraphAsymmErrors*)f->Get(vd[i].c_str());
    assert(gt);
    if (i<=1) mg->Add(gt);
    vg[i]  = g;
    vt[i]  = gt;
    vr[i]  = (TGraphErrors*)g;
  }

  TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(x)+[3]*pow(log(x),2)"
		    "+[4]*pow(log(x),3))"
		    "*pow(1-2*x*cosh([6])/7000.,[5])",
		    60,2000.);
  f1->SetParameters(1e11,-5,0.01,0.001,0.0001,10,0.25);
  f1->FixParameter(6,0.);
  f1->FixParameter(5,10.);
  f1->FixParameter(4,0.);
  f1->FixParameter(3,0.);
  f1->FixParameter(2,0.);
  f1->FixParameter(1,-5);
  f1->FixParameter(0,1.2e14);
  mg->Fit(f1,"R");

  // 2D variant
  TF2 *f2 = new TF2("f2","[0]*pow(x,[1])"
		    "*pow(1-min(2*[3]*x*cosh(y),6990.)/7000.,[2])",
		    60,2000.,0,5.2);
  f2->SetParameters(1.2e14,-5,10,1.);//0.90);

  for (int i = 0; i != vg.size(); ++i) {
    
    TGraphAsymmErrors *g = vg[i];
    TGraphAsymmErrors *gt = vt[i];
    for (int j = 0; j != g->GetN(); ++j) {

      double x1 = g->GetX()[j] - g->GetEXlow()[j];
      double x2 = g->GetX()[j] + g->GetEXhigh()[j];
      //double fy = f1->Integral(x1, x2) / (x2 - x1)
      // / (kfactor0 * kfactor[vd[i]]);
      double fy = 2.*f2->Integral(x1, x2, 0, 0.5) / (x2 - x1)
	/ (kfactor0 * kfactor[vd[i]]);
      double y = g->GetY()[j];

      g->SetPoint(j, g->GetX()[j], y/fy);
      g->SetPointError(j, 0, 0, //g->GetEXlow()[j], g->GetEXhigh()[j],
		       g->GetEYlow()[j]/fy, g->GetEYhigh()[j]/fy);
      gt->SetPoint(j, gt->GetX()[j], y/fy);
      gt->SetPointError(j, gt->GetEXlow()[j], gt->GetEXhigh()[j],
		       gt->GetEYlow()[j]/fy, gt->GetEYhigh()[j]/fy);
    }

    g->SetLineColor(color[vd[i]]);
    g->SetMarkerColor(color[vd[i]]);
    g->SetMarkerStyle(marker[vd[i]]);
    gt->SetLineColor(color[vd[i]]);
    gt->SetMarkerColor(color[vd[i]]);
    gt->SetMarkerStyle(marker[vd[i]]);

    //if (i==0) gt->Draw("SAMEPz");
    if (i!=1 && i!=3) {
      g->Draw("SAMEP");
      leg->AddEntry(g,label[vd[i]],"PL");
    }
    //if (i==0) { gt->Draw("Pz"); g->Draw("SAMEP"); }
    //else      { gt->Draw("SAMEPz"); g->Draw("SAMEP"); }
  }

  // Estimate R=0.6 for CMS data
  if (true) {
    TGraphErrors *gr6 = new TGraphErrors();
    TGraphAsymmErrors *g7 = vg[0];
    TGraphAsymmErrors *g5 = vg[2];

    // Assuming cross section is proportional to log(R):
    double w7 = (log(0.6)-log(0.5)) / (log(0.7)-log(0.5));

    for (int i = 0; i != g7->GetN(); ++i) {
      gr6->SetPoint(i, g7->GetX()[i], w7*g7->GetY()[i]+(1-w7)*g5->GetY()[i]);
      gr6->SetPointError(i, g7->GetEXlow()[i],
			 w7*g7->GetEYlow()[i]+(1-w7)*g5->GetEYlow()[i]);
    }
    gr6->SetMarkerStyle(kFullDiamond);
    gr6->SetMarkerColor(kBlack);
    gr6->Draw("SAMEPz");
    leg->AddEntry(gr6,"CMS R=0.6 (int.)","PL");
  }

  // Estimate 1% JEC shift for ATLAS
  if (true) {
    TGraphErrors *gs6 = new TGraphErrors();
    TGraphErrors *gs6t = new TGraphErrors();
    TGraphErrors *gs4 = new TGraphErrors();
    TGraphErrors *gr5 = new TGraphErrors();
    TGraphAsymmErrors *g6 = vg[1]; assert(g6);
    TGraphAsymmErrors *gt6 = vt[1]; assert(gt6);
    TGraphAsymmErrors *g4 = vg[3]; assert(g4);
    TGraphAsymmErrors *gt4 = vt[3]; assert(gt4);

    for (int i = 0; i != gt6->GetN(); ++i) {

      double x1_0 = gt6->GetX()[i] - gt6->GetEXlow()[i];
      double x2_0 = gt6->GetX()[i] + gt6->GetEXhigh()[i];
      //double fy0 = f1->Integral(x1_0, x2_0) / (x2_0 - x1_0);
      double fy0 = 2.*f2->Integral(x1_0, x2_0, 0, 0.5) / (x2_0 - x1_0);
      double djec = +0.02;
      double x1 = (1.+djec)*x1_0;
      double x2 = (1.+djec)*x2_0;
      //double fy1 = f1->Integral(x1, x2) / (x2_0 - x1_0);
      double fy1 = 2.*f2->Integral(x1, x2, 0, 0.5) / (x2_0 - x1_0);

      double kf = 1;//1.10; // compensate JEC shift at pT=200 GeV

      gs6->SetPoint(i, gt6->GetX()[i], gt6->GetY()[i]*kf*fy1/fy0);
      gs6->SetPointError(i, gt6->GetEXlow()[i], g6->GetEYlow()[i]*kf*fy1/fy1);
      gs6t->SetPoint(i, gt6->GetX()[i], gt6->GetY()[i]*kf*fy1/fy0);
      gs6t->SetPointError(i, gt6->GetEXlow()[i], g6->GetEYlow()[i]*kf*fy1/fy1);

      gs4->SetPoint(i, gt4->GetX()[i], gt4->GetY()[i]*kf*fy1/fy0);
      gs4->SetPointError(i, gt4->GetEXlow()[i], g4->GetEYlow()[i]*kf*fy1/fy1);

      double y6 = gs6->GetY()[i];
      double y4 = gs4->GetY()[i];
      double ey6 = gs6->GetEY()[i];
      double ey4 = gs4->GetEY()[i];
      // Assuming cross section is proportional to log(R):
      double w6 = (log(0.5)-log(0.4)) / (log(0.6)-log(0.4));
      double y5 = w6*y6 + (1-w6)*y4;
      double ey5 = w6*ey6 + (1-w6)*ey4;

      gr5->SetPoint(i, gs4->GetX()[i], y5);
      gr5->SetPointError(i, gs4->GetEX()[i], ey5);

      // Show systematis around ATLAS R=0.6
      gs6t->SetPointError(i, gt6->GetEXlow()[i], gt6->GetEYlow()[i]*kf*fy1/fy1);
    }
    gs6t->SetMarkerStyle(kOpenCircle);//Diamond);
    gs6t->SetMarkerColor(kGreen+2);
    gs6t->SetLineColor(kGreen+2);
    gs6t->Draw("SAMEPz");
    leg->AddEntry(gs6t,"ATLAS R=0.6 (#Deltap_{T}=-2%)","PL");
    //leg->AddEntry(gs6,"ATLAS R=0.6","PL");
    vr[1] = gs6;

    gs4->SetMarkerStyle(kOpenCircle);//Diamond);
    gs4->SetMarkerColor(kRed);
    gs4->SetLineColor(kRed);
    gs4->Draw("SAMEPz");
    leg->AddEntry(gs4,"ATLAS R=0.4 (#Deltap_{T}=-2%)","PL");
    //leg->AddEntry(gs4,"ATLAS R=0.4","PL");
    vr[3] = gs4;

    gr5->SetMarkerStyle(kOpenCircle);//Diamond);
    gr5->SetMarkerColor(kOrange+2);
    gr5->SetLineColor(kOrange+2);
    gr5->Draw("SAMEPz");
    leg->AddEntry(gr5,"ATLAS R=0.5 (int.)","PL");
  }
		
  cout << "Fit params:" << endl;
  for (int i = 0; i != f1->GetNpar(); ++i) {
    cout << Form("p%d = %1.4g +/- %1.4g",i,
		 f1->GetParameter(i),f1->GetParError(i)) << endl;
  }
				
  gPad->SetLogx();

  TF1 *f1r = new TF1("f1r","([7]*[0]*pow(x-[6],[1]) + [4]*pow(x-[6],[5]))"
		     "*pow(1-2*(x-[6])*cosh([3])/7000.,[2]) /"
		     "([0]*pow(x,[1])"
		     "*pow(1-2*x*cosh([3])/7000.,[2]))",
		     60,2000.);
  //f1r->SetParameters(1.2e14,-5,10,0, 2.4e7,-3, 1.5*1.5, 0.9);

  f1r->SetParameters(1.2e14,-5,10,0, 3e4,-2, 1.5*1.54, 0.9);
  f1r->SetLineColor(kBlue);
  f1r->DrawClone("SAME");

  f1r->SetParameters(1.2e14,-5,10,0, 3e4,-2, 1.5*1.13, 0.85);
  f1r->SetLineColor(kGreen+2);
  f1r->DrawClone("SAME");

  f1r->SetParameters(1.2e14,-5,10,0, 3e4,-2, 1.5*0.79, 0.80);
  f1r->SetLineColor(kOrange+2);
  f1r->DrawClone("SAME");

  f1r->SetParameters(1.2e14,-5,10,0, 3e4,-2, 1.5*0.50, 0.75);
  f1r->SetLineColor(kRed);
  f1r->DrawClone("SAME");

  c1->SaveAs("pdf/drawData.pdf");


  TH1D *h2 = (TH1D*)h->Clone("h2");
  h2->SetYTitle("R(X) / R(0.7)");
  TCanvas *c2 = tdrCanvas("c2",h2,1,11,kSquare);
  
  //TGraphErrors *gr = tools::ratioGraphs((TGraphErrors*)vg[2],
  //					(TGraphErrors*)vg[0]);
  //gr->Draw("SAMEP");
  TGraphErrors *gr7 = tools::ratioGraphs(vr[0],vr[0]);
  gr7->Draw("SAMEP");
  TGraphErrors *gr6 = tools::ratioGraphs(vr[1],vr[0]);
  gr6->Draw("SAMEP");
  TGraphErrors *gr5 = tools::ratioGraphs(vr[2],vr[0]);
  gr5->Draw("SAMEP");
  TGraphErrors *gr4 = tools::ratioGraphs(vr[3],vr[0]);
  gr4->Draw("SAMEP");

  gPad->SetLogx();

  c2->SaveAs("pdf/drawData_ratio.pdf");

}
