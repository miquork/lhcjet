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

// Notes:
// - compare uncertainties vs eta vs deltaJEC
// - implement full nuisances and chi2
// - theory for better interpolation? (Marek Schonherr, LHC EW WG:
//   https://indico.cern.ch/event/730246/timetable/#20180621.detailed )

double _beta = 9;//10;
double _N0 = 1.15e14;//1.2e14;
void drawDatas(double ymin=0, double ymax=0.5,
	       double shiftjecATLAS=0, double shiftjecCMS=0);
void drawJECunc();
void drawDeltaJEC();

void drawData() {
  
  /*
  _N0=1.15e14; _beta=9.;//5;
  drawDatas(0,0.5);
  //drawDatas(0,0.5,-0.02);
  drawDatas(0,0.5,-0.01,+0.01);
  _N0=1.10e14; //_beta=8;
  drawDatas(0.5,1.0);
  drawDatas(0.5,1.0,-0.01,+0.01);
  _N0=1.05e14; //_beta=7;
  drawDatas(1.0,1.5);
  drawDatas(1.0,1.5,-0.005,+0.005);
  _N0=1.0e14; //_beta=6;
  drawDatas(1.5,2.0);
  drawDatas(1.5,2.0,-0.002,+0.002);
  //_N0=0.80e14; //_beta=5;//5.4;
  _N0=0.90e14; _beta=9.5;//5.4;
  drawDatas(2.0,2.5);
  drawDatas(2.0,2.5,+0.005,-0.005);
  //_N0=0.50e14; //_beta=4;//4.2;
  _N0=0.75e14; _beta=10;//4.2;
  drawDatas(2.5,3.0);
  drawDatas(2.5,3.0,+0.01,-0.01);

  drawDeltaJEC();
  */
  drawJECunc();
}

void drawDatas(double ymin, double ymax,
	       double shiftjecATLAS, double shiftjecCMS) {
	       

  setTDRStyle();
  
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("lhcdata.root","READ");
  TFile *fout = new TFile("lhcdataratio.root","UPDATE");
  assert(f && !f->IsZombie());

  TFile *fc = new TFile("cms/HEPData-ins1298810-v1-Table_1.root","READ");
  assert(fc && !fc->IsZombie());

  curdir->cd();

  //TH1D *h = new TH1D("h",";p_{T} (GeV);Data / fit",100,50,2000.);
  TH1D *h = new TH1D("h",";p_{T} (GeV);Data / F(p_{T},|#eta|)",
		     //#;N=1.2e14,#alpha=5,#beta=10,|#eta|=0)",
		       100,50,2000.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(2.0);//1.5);//1.3);
  h->SetMinimum(0.0);//0.5);//0.6);

  // Continuous 2D (pT,y) model of jet cross section
  TF2 *f2 = new TF2("f2","[0]*pow(x,[1])"
		    //"*pow(1-min(2*[3]*x*cosh(y),6990.)/7000.,[2])",
		    //"*pow(1-min(2*x*cosh(y),6998.)/7000.,[2])"
		    "*pow(1-min(2*x*cosh(y),6998.)/7000.,[2]-2*y)"
		    "*(2*x*cosh(y)<6998.)",
		    60,2000.,0,5.2);
  //f2->SetParameters(1.2e14,-5,10,1.);//0.90);
  f2->SetParameters(_N0,-5,_beta);//,1.);//0.90);

  //cmsText = "CMS+ATLAS";
  cmsText = "HEPData (1406.0324, 1410.8857)";
  lumi_7TeV = "CMS 5.0 fb^{-1}, ATLAS 4.5 fb^{-1}";
  extraText = Form("N=%1.2f#times10^{14}, #alpha=5, #beta=%1.2g-2|y|",
		   //, |#eta|=%1.1f",
		   _N0/1e14, _beta);//, ymin);
  TCanvas *c1 = tdrCanvas("c1",h,1,11,kSquare);
  // 1.1e14/fb? 1e11/pb, 1e8/nb, 1e5/mub, 1e2/mb => ~total xsec at 1 GeV?
  // Offset per unit area in QCD is 2 GeV.

  // Kinematic suppression vs eta
  // TF1 *f0 = new TF1("f0","pow(1-min(2*[0]*cosh(x),6500.)/7000.,10)",0,5)
  // Drop is steeper at higher eta, so should do 2D integration in pT and eta?

  //TLegend *leg = tdrLeg(0.20,0.58,0.50,0.78);
  TLegend *leg = tdrLeg(0.20,0.62,0.50,0.78);
  leg->SetTextSize(0.035);

  //TLegend *leg2 = tdrLeg(0.50, shiftjecATLAS ? 0.58 : 0.63,
  //			 0.80, shiftjecCMS ? 0.78 : 0.73);
  TLegend *leg2 = tdrLeg(0.50, shiftjecATLAS ? 0.62 : 0.66,
  			 0.80, shiftjecCMS ? 0.78 : 0.74);
  leg2->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  //tex->SetTextSize(0.035);
  tex->SetNDC();
  if (ymin==0 && ymax==0.5)
    tex->DrawLatex(0.35,0.17,"|y| < 0.5, anti-k_{T} R = 0.4--0.7");
  else
    tex->DrawLatex(0.35,0.17,Form("%1.1f<|y|<%1.1f, anti-k_{T} R = 0.4--0.7",
				  ymin, ymax));

  vector<string> vd;
  vd.push_back("cms07_r07");
  vd.push_back("atlas07_r06");
  vd.push_back("cms07_r05");
  //vd.push_back("Table 1/Graph1D_y1");
  vd.push_back("atlas07_r04");

  map<string, int> color;
  color["cms07_r07"] = kBlue;
  color["atlas07_r06"] = kGreen+2;
  color["cms07_r05"] = kOrange+2;
  //color["Table 1/Graph1D_y1"] = kOrange+2;
  color["atlas07_r04"] = kRed;
  
  map<string, int> marker;
  marker["cms07_r07"] = kFullCircle;
  marker["atlas07_r06"] = kOpenCircle;
  marker["cms07_r05"] = kFullCircle;
  //marker["Table 1/Graph1D_y1"] = kFullCircle;
  marker["atlas07_r04"] = kOpenCircle;

  map<string, const char*> label;
  label["cms07_r07"] = "CMS R=0.7";
  label["atlas07_r06"] = "ATLAS R=0.6";
  label["cms07_r05"] = "CMS R=0.5";
  //label["Table 1/Graph1D_y1"] = "CMS R=0.5";
  label["atlas07_r04"] = "ATLAS R=0.4";
  
  // CMS lumi 5.0/fb, ATLAS 4.5/fb, 11% difference
  // We get better agreement assuming they are the same
  // => Correct CMS up 5%, ATLAS down 5%
  double kfactor0 = 1;//1./1.05;
  map<string, double> kfactor;
  kfactor["cms07_r07"] = 1;//1.05;
  kfactor["atlas07_r06"] = 1;//0.90;
  kfactor["cms07_r05"] = 1;//1.05;
  //kfactor["Table 1/Graph1D_y1"] = 1;
  kfactor["atlas07_r04"] = 1;//0.90;

  
  TMultiGraph *mg = new TMultiGraph();
  vector<TGraphAsymmErrors*> vg(vd.size());
  vector<TGraphAsymmErrors*> vt(vd.size());
  vector<TGraphErrors*> vr(vd.size());
  for (int i = 0; i != vd.size(); ++i) {

    string ss = Form("%s_y%02d-%02d_%s",vd[i].c_str(),
		     int(10*ymin+0.5), int(10*ymax+0.5),
		     "stat");
    cout << ss << endl << flush;
    TGraphAsymmErrors *g = (TGraphAsymmErrors*)f->Get(ss.c_str());
    //if (!g) {
    //g = (TGraphAsymmErrors*)fc->Get((vd[i]).c_str());
    //g = (TGraphAsymmErrors*)g->Clone();
    //}
    assert(g);

    string s = Form("%s_y%02d-%02d",vd[i].c_str(),
		    int(10*ymin+0.5), int(10*ymax+0.5));
    cout << s << endl << flush;
    TGraphAsymmErrors *gt = (TGraphAsymmErrors*)f->Get(s.c_str());
    //if (!gt) {
    //gt = (TGraphAsymmErrors*)fc->Get((vd[i]).c_str());
    //}
    assert(gt);
    if (i<=1) mg->Add(gt);
    vg[i]  = g;
    vt[i]  = gt;
    vr[i]  = (TGraphErrors*)g;
  }

  /*
  TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(x)+[3]*pow(log(x),2)"
		    "+[4]*pow(log(x),3))"
		    "*pow(1-2*x*cosh([6])/7000.,[5])",
		    60,2000.);
  f1->SetParameters(_N0/1e3,-5,0.01,0.001,0.0001,_beta,ymin);
  f1->FixParameter(6,ymin);
  f1->FixParameter(5,_beta);//10.);
  f1->FixParameter(4,0.);
  f1->FixParameter(3,0.);
  f1->FixParameter(2,0.);
  f1->FixParameter(1,-5);
  f1->FixParameter(0,_N0/1e3);
  mg->Fit(f1,"R");
  */

  for (int i = 0; i != vg.size(); ++i) {
    
    TGraphAsymmErrors *g = vg[i];
    TGraphAsymmErrors *gt = vt[i];
    for (int j = 0; j != g->GetN(); ++j) {

      double x1 = g->GetX()[j] - g->GetEXlow()[j];
      double x2 = g->GetX()[j] + g->GetEXhigh()[j];
      //double fy = f1->Integral(x1, x2) / (x2 - x1)
      // / (kfactor0 * kfactor[vd[i]]);
      double fy = 2.*f2->Integral(x1, x2, ymin, ymax) / (x2 - x1)
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

    // derive JEC shift
    //TGraphAsymmErrors *gs = new TGraphAsymmErrors();
    TGraphAsymmErrors *gs = (TGraphAsymmErrors*)g->Clone();
    TGraphAsymmErrors *gst = (TGraphAsymmErrors*)gt->Clone();
    assert(g->GetN()==gt->GetN());
    for (int k = 0; k != gt->GetN(); ++k) {

      double x1_0 = gt->GetX()[k] - gt->GetEXlow()[k];
      double x2_0 = gt->GetX()[k] + gt->GetEXhigh()[k];
      double fy0 = 2.*f2->Integral(x1_0, x2_0, ymin, ymax) / (x2_0 - x1_0);
      double djec = 0;
      if (shiftjecATLAS && (i==1 || i==3)) djec = shiftjecATLAS;
      if (shiftjecCMS   && (i==0 || i==2)) djec = shiftjecCMS;
      double x1 = (1.-djec)*x1_0;
      double x2 = (1.-djec)*x2_0;
      double fy1 = 2.*f2->Integral(x1, x2, ymin, ymax) / (x2_0 - x1_0);

      gs->SetPoint(k, g->GetX()[k], g->GetY()[k]*fy1/fy0);
      gs->SetPointError(k, g->GetEXlow()[k], g->GetEXhigh()[k], 
			g->GetEYlow()[k]*fy1/fy0, g->GetEYhigh()[k]*fy1/fy0);

      gst->SetPoint(k, gt->GetX()[k], gt->GetY()[k]*fy1/fy0);
      gst->SetPointError(k, gt->GetEXlow()[k], gt->GetEXhigh()[k], 
			gt->GetEYlow()[k]*fy1/fy0, gt->GetEYhigh()[k]*fy1/fy0);
    } // for k
    g = gs;
    vg[i] = gs;
    gt = gst;
    vt[i] = gst;

    //if (i==0) gt->Draw("SAMEPz");
    //if (i!=1 && i!=3) {
    if (i==1 && !shiftjecATLAS) gt->Draw("SAMEPz");
    else g->Draw("SAMEPz");
    leg->AddEntry(g,label[vd[i]],"PL");
      //}
    //if (i==0) { gt->Draw("Pz"); g->Draw("SAMEP"); }
    //else      { gt->Draw("SAMEPz"); g->Draw("SAMEP"); }

    if (!shiftjecATLAS && !shiftjecCMS) {
      fout->cd();
      g->Write(g->GetName(),TObject::kOverwrite);
      gt->Write(gt->GetName(),TObject::kOverwrite);
      curdir->cd();
    }
  }

  // Estimate R=0.6 for CMS data
  if (true) {
    TGraphAsymmErrors *gr6 = new TGraphAsymmErrors();
    TGraphAsymmErrors *gr6t = new TGraphAsymmErrors();
    TGraphAsymmErrors *g7 = vg[0];
    TGraphAsymmErrors *g5 = vg[2];
    TGraphAsymmErrors *g7t = vt[0];
    TGraphAsymmErrors *g5t = vt[2];

    // Assuming cross section is proportional to log(R):
    double w7 = (log(0.6)-log(0.5)) / (log(0.7)-log(0.5));

    for (int i = 0; i != g7->GetN(); ++i) {
      gr6->SetPoint(i, g7->GetX()[i], w7*g7->GetY()[i]+(1-w7)*g5->GetY()[i]);
      gr6->SetPointError(i, g7->GetEXlow()[i], g7->GetEXhigh()[i],
			 w7*g7->GetEYlow()[i]+(1-w7)*g5->GetEYlow()[i],
			 w7*g7->GetEYhigh()[i]+(1-w7)*g5->GetEYhigh()[i]);

      gr6t->SetPoint(i, g7t->GetX()[i], w7*g7t->GetY()[i]+(1-w7)*g5t->GetY()[i]);
      gr6t->SetPointError(i, g7t->GetEXlow()[i], g7t->GetEXhigh()[i],
			  w7*g7t->GetEYlow()[i]+(1-w7)*g5t->GetEYlow()[i],
			  w7*g7t->GetEYhigh()[i]+(1-w7)*g5t->GetEYhigh()[i]);
    }

    TGraphAsymmErrors *gr = (shiftjecCMS ? gr6 : gr6t);
    gr->SetMarkerStyle(kFullDiamond);
    gr->SetMarkerColor(kGreen+3);
    gr->SetLineColor(kGreen+3);
    gr->Draw("SAMEPz");
    TGraphAsymmErrors *gr2 = (shiftjecCMS ? gr6t : gr6);
    gr2->SetMarkerStyle(kFullCircle);
    if (shiftjecCMS) leg2->AddEntry(gr2,Form("CMS #Deltap_{T}=%1.1f%%",
					     100.*shiftjecCMS),"P");
    leg2->AddEntry(gr,"CMS R=0.6 (int.)","PL");

    if (!shiftjecCMS) {
      fout->cd();
      gr6->Write(Form("cms07_r06_y%02d-%02d_stat",
		      int(10*ymin+0.5), int(10*ymax+0.5)),
		 TObject::kOverwrite);
      gr6t->Write(Form("cms07_r06_y%02d-%02d",
		       int(10*ymin+0.5), int(10*ymax+0.5)),
		  TObject::kOverwrite);
      curdir->cd();
    }
  } // CMS

  // Estimate R=0.5 for ATLAS data
  if (true) {
    TGraphErrors *gr5 = new TGraphErrors();
    TGraphErrors *gr5t = new TGraphErrors();
    TGraphAsymmErrors *g6 = vg[1];
    TGraphAsymmErrors *g4 = vg[3];
    TGraphAsymmErrors *g6t = vt[1];
    TGraphAsymmErrors *g4t = vt[3];

    // Assuming cross section is proportional to log(R):
    double w6 = (log(0.5)-log(0.4)) / (log(0.6)-log(0.4));

    for (int i = 0; i != g6->GetN(); ++i) {
      gr5->SetPoint(i, g6->GetX()[i], w6*g6->GetY()[i]+(1-w6)*g4->GetY()[i]);
      gr5->SetPointError(i, g6->GetEXlow()[i],
			 w6*g6->GetEYlow()[i]+(1-w6)*g4->GetEYlow()[i]);

      gr5t->SetPoint(i, g6t->GetX()[i], w6*g6t->GetY()[i]+(1-w6)*g4t->GetY()[i]);
      gr5t->SetPointError(i, g6t->GetEXlow()[i],
			 w6*g6t->GetEYlow()[i]+(1-w6)*g4t->GetEYlow()[i]);

    }
    gr5->SetMarkerStyle(kOpenDiamond);
    gr5->SetMarkerColor(kOrange+3);
    gr5->SetLineColor(kOrange+3);
    gr5->Draw("SAMEPz");
    leg2->AddEntry(gr5,"ATLAS R=0.5 (int.)","PL");
    gr5t->SetMarkerStyle(kOpenCircle);
    if (shiftjecATLAS)
      leg2->AddEntry(gr5t,Form("ATLAS #Deltap_{T}=%1.1f%%",100.*shiftjecATLAS),
		     "P");

    if (!shiftjecATLAS) {
      fout->cd();
      g6->Write(Form("atlas07_r05_y%02d-%02d_stat",
		     int(10*ymin+0.5), int(10*ymax+0.5)),
		TObject::kOverwrite);
      g6t->Write(Form("atlas07_r05_y%02d-%02d",
		      int(10*ymin+0.5), int(10*ymax+0.5)),
		 TObject::kOverwrite);
      curdir->cd();
    }
  } // ATLAS

  // Estimate 1% JEC shift for ATLAS
  if (false) {
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
      double x1 = (1.-shiftjecATLAS)*x1_0;
      double x2 = (1.-shiftjecATLAS)*x2_0;
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

      // Show systematics around ATLAS R=0.6
      gs6t->SetPointError(i, gt6->GetEXlow()[i], gt6->GetEYlow()[i]*kf*fy1/fy1);
    }
    gs6t->SetMarkerStyle(kOpenCircle);//Diamond);
    gs6t->SetMarkerColor(kGreen+2);
    gs6t->SetLineColor(kGreen+2);
    gs6t->Draw("SAMEPz");
    if (shiftjecATLAS)
      leg->AddEntry(gs6t,Form("ATLAS R=0.6 (#Deltap_{T}=%1.1f%%)",
			     100.*shiftjecATLAS),"PL");
    else 
      leg->AddEntry(gs6t,"ATLAS R=0.6","PL");
    vr[1] = gs6;

    gs4->SetMarkerStyle(kOpenCircle);//Diamond);
    gs4->SetMarkerColor(kRed);
    gs4->SetLineColor(kRed);
    gs4->Draw("SAMEPz");
    if (shiftjecATLAS)
      leg->AddEntry(gs4,Form("ATLAS R=0.4 (#Deltap_{T}=%1.1f%%)",
			     100.*shiftjecATLAS),"PL");
    else
      leg->AddEntry(gs4,"ATLAS R=0.4","PL");
    vr[3] = gs4;

    gr5->SetMarkerStyle(kOpenCircle);//Diamond);
    gr5->SetMarkerColor(kOrange+2);
    gr5->SetLineColor(kOrange+2);
    gr5->Draw("SAMEPz");
    leg->AddEntry(gr5,"ATLAS R=0.5 (int.)","PL");
  }

  /*
  cout << "Fit params:" << endl;
  for (int i = 0; i != f1->GetNpar(); ++i) {
    cout << Form("p%d = %1.4g +/- %1.4g",i,
		 f1->GetParameter(i),f1->GetParError(i)) << endl;
  }
  */
				
  gPad->SetLogx();

  /*
  TF1 *f1r = new TF1("f1r","([7]*[0]*pow(x-[6],[1]) + [4]*pow(x-[6],[5]))"
		     "*pow(1-2*(x-[6])*cosh([3])/7000.,[2]) /"
		     "([0]*pow(x,[1])"
		     "*pow(1-2*x*cosh([3])/7000.,[2]))",
		     60,2000.);
  //f1r->SetParameters(1.2e14,-5,10,0, 2.4e7,-3, 1.5*1.5, 0.9);
  f1r->SetParameters(_N0,-5,_beta,ymin, 3e4,-2*0, 1.5*1.54, 0.9);
  f1r->SetLineColor(kBlue);
  f1r->DrawClone("SAME");

  f1r->SetParameters(_N0,-5,_beta,ymin, 3e4,-2*0, 1.5*1.13, 0.85);
  f1r->SetLineColor(kGreen+2);
  f1r->DrawClone("SAME");

  f1r->SetParameters(_N0,-5,_beta,ymin, 3e4,-2*0, 1.5*0.79, 0.80);
  f1r->SetLineColor(kOrange+2);
  f1r->DrawClone("SAME");

  f1r->SetParameters(_N0,-5,_beta,ymin, 3e4,-2*0, 1.5*0.50, 0.75);
  f1r->SetLineColor(kRed);
  f1r->DrawClone("SAME");
  */
  //c1->SaveAs("pdf/drawData.pdf");
  c1->SaveAs(Form("pdf/drawData_%1.1f-%1.1f%s%s.pdf",ymin,ymax,
		  shiftjecATLAS ? "_djecATLAS" : "",
		  shiftjecCMS ? "_djecCMS" : ""));


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

  fout->Close();
}

void drawJECunc() {

  TDirectory *curdir = gDirectory;

  setTDRStyle();

  // ATLAS already differs between 37/pb and 4.5/fb, independent of radius
  TFile *f = new TFile("lhcdataratio.root","READ");
  assert(f && !f->IsZombie());

  // Continuous 2D (pT,y) model of jet cross section
  //TF2 *f2 = new TF2("f2","[0]*pow(x,[1])"
  //		    "*pow(1-min(2*x*cosh(y),6998.)/7000.,[2]-2*y)"
  //		    "*(2*x*cosh(y)<6998.)",
  //		    60,2000.,0,5.2);
  //f2->SetParameters(_N0,-5,_beta);

  TF1 *f1 = new TF2("f1","[0]*pow(x,[1])"
		    "*pow(1-min(2*x*cosh([3]),6998.)/7000.,[2]-2*[3])"
		    "*(2*x*cosh([3])<6998.)",
		    60,2000.);
  f1->SetParameters(_N0,-5,_beta,0.);

  // vsEta
  const int neta = 6;
  TGraph *gea = new TGraph(neta);
  TGraph *gec = new TGraph(neta);
  for (int ieta = 0; ieta != neta; ++ieta) {

    double y1(0.5*ieta), y2(0.5*(ieta+1));
    //double ptref(200.); // eye-balled deltaJEC
    double ptref(100.); // CMS ref plot
    f1->SetParameters(_N0,-5,_beta,0.);

    // ATLAS
    TGraphAsymmErrors *ga = (TGraphAsymmErrors*)f->Get(Form("atlas07_r06_y%02.0f-%02.0f",10.*y1,10.*y2));
    assert(ga);
    TGraphAsymmErrors *gc = (TGraphAsymmErrors*)f->Get(Form("cms07_r06_y%02.0f-%02.0f_syst",10.*y1,10.*y2));
    assert(gc);

    int ipta(-1), dpta(1e5);
    for (int i = 0; i != ga->GetN(); ++i) {
      double dx = fabs(ga->GetX()[i]-ptref);
      if (dx<dpta) { ipta = i; dpta = dx; }
    }
    int iptc(-1), dptc(1e5);
    for (int i = 0; i != gc->GetN(); ++i) {
      double dx = fabs(gc->GetX()[i]-ptref);
      if (dx<dptc) { iptc = i; dptc = dx; }
    }

    double dya = ga->GetEYhigh()[ipta]/ga->GetY()[ipta];
    double dyc = gc->GetEYhigh()[iptc]/gc->GetY()[iptc];

    double y = f1->Eval(ptref);
    double dxa = fabs(f1->GetX(y*(1+dya))/f1->GetX(y)-1);
    double dxc = fabs(f1->GetX(y*(1+dyc))/f1->GetX(y)-1);

    gea->SetPoint(ieta, 0.5*(y1+y2), 100.*dxa);
    gec->SetPoint(ieta, 0.5*(y1+y2), 100.*dxc);
  }

  TH1D *h = new TH1D("h",";|y_{jet}|;Uncertainty in #DeltaJEC (%)",100,0,3.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(6.0);
  h->SetMinimum(0.0);

  cmsText = "HEPData (1406.0324, 1410.8857)";
  lumi_7TeV = "CMS 5.0 fb^{-1}, ATLAS 4.5 fb^{-1}";
  extraText = "";
  TCanvas *c1 = tdrCanvas("c1",h,1,11,kSquare);

  TLegend *leg = tdrLeg(0.20,0.75,0.50,0.85);
  leg->SetTextSize(0.045);
  leg->AddEntry(gea,"ATLAS","PL");
  leg->AddEntry(gec,"ATLAS","PL");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.20, 0.20, Form("p_{T} = %1.0f GeV",ptref));

  gea->SetMarkerStyle(kOpenCircle);
  gea->Draw("SAMEPL");
  gec->SetMarkerStyle(kOpenCircle);
  gec->Draw("SAMEPL");

  // Notes:
  // 1) CMS uncertainty seems to include very large Time stability
  //    This was later understood to be absorbed in L2Res for whole data set
  //    (check that Time was on, and which one it was; cols 22-23, N-4,3, 12-13
  //    )
  // 2) CMS uncertainty also has large PileUp uncertainty in EC
  //    This was later understood to be partly absorbed in pT-dependent L2res
  //    (check if the pT dependence was used)
}

void drawDeltaJEC() {

  setTDRStyle();

  double x[] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75};
  double y1[] = {-1, -1, -0.5, -0.2, +0.5, +1.0};
  double y2[] = {+1, +1, +0.5, +0.2, -0.5, -1.0};
  const int n = sizeof(x)/sizeof(x[0]);
  double dy[n] = {0., 0., 0.3, 0.5, 0.7, 1.}; 

  double y1b[n], y2b[n];
  for (int i = 0; i != n; ++i) {
    y1b[i] = y1[i]+dy[i];
    y2b[i] = y2[i]+dy[i];
  }

  TH1D *h = new TH1D("h",";|y_{jet}|;Estimated #DeltaJEC (%)",100,0,3.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(3.0);
  h->SetMinimum(-2.0);

  cmsText = "HEPData (1406.0324, 1410.8857)";
  lumi_7TeV = "CMS 5.0 fb^{-1}, ATLAS 4.5 fb^{-1}";
  extraText = "";
  TCanvas *c1 = tdrCanvas("c1",h,1,11,kSquare);

  TLegend *leg = tdrLeg(0.20,0.73,0.50,0.83);
  leg->SetTextSize(0.045);
  TLegend *leg2 = tdrLeg(0.45,0.73,0.75,0.83);
  leg2->SetTextSize(0.045);

  TGraph *g1 = new TGraph(n, x, y1);
  TGraph *g2 = new TGraph(n, x, y2);
  TGraph *g1b = new TGraph(n, x, y1b);
  TGraph *g2b = new TGraph(n, x, y2b);

  g1b->SetMarkerColor(kRed);
  g1b->SetMarkerStyle(kOpenSquare);
  //g1b->Draw("SAMEP");
  g1->SetMarkerStyle(kOpenCircle);
  g1->Draw("SAMEP");

  g2b->SetMarkerColor(kRed);
  g2b->SetMarkerStyle(kFullSquare);
  //g2b->Draw("SAMEP");
  g2->SetMarkerStyle(kFullCircle);
  g2->Draw("SAMEP");

  leg->AddEntry(g1,"ATLAS","P");
  leg->AddEntry(g2,"CMS","P");
  //leg2->AddEntry(g1b,"ATLAS adj.","P");
  //leg2->AddEntry(g2b,"CMS adj.","P");


  c1->SaveAs("pdf/drawData_deltaJEC.pdf");

}
