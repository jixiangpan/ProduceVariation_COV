#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "./draw.icc"

////////////////////////////////////////////////////////////////////////////////////////////////

class TCN
{
public:
  TCN() {
    cout<<endl<<" ---> Hello TCN"<<endl<<endl;

    rand = new TRandom3(100001);

    BINS = 0;
    TOYS = 0;

    FLAG_NORM = 0;
    norm_relerr = 0;
  }

  TRandom3 *rand;

  bool FLAG_NORM;
  int BINS;
  int TOYS;
  
  TMatrixD matrix_pred;
  TMatrixD matrix_meas;
  TMatrixD matrix_syst_abscov;

  double norm_relerr;
  
  TMatrixD matrix_fake_meas;
  map<int, map<int, double> >map_toy_meas;

  void Initialization(bool flag_norm);
  
  void Set_norm_relerr(double val_relerr) { norm_relerr = val_relerr; };
  
  void ProduceVariation(int ntoys, bool flag_norm);
  
  void SetToy(int itoy) {
    matrix_fake_meas.Clear(); matrix_fake_meas.ResizeTo(1, BINS);    
    for(int ibin=1; ibin<=BINS; ibin++) matrix_fake_meas(0, ibin-1) = map_toy_meas[itoy][ibin];
  }

  void SetMeas2Use() {
    matrix_fake_meas.Clear(); matrix_fake_meas.ResizeTo(1, BINS);
    for(int ibin=1; ibin<=BINS; ibin++) matrix_fake_meas(0, ibin-1) = matrix_meas(0, ibin-1);
  }
  
  double GetChi2(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp);

  void Plotting_singlecase(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp, bool saveFIG, int index);
};

/////////////////////// ccc

void TCN::Plotting_singlecase(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp, bool saveFIG, int index)
{
  TString roostr = "";
  
  int rows = matrix_pred_temp.GetNcols();
  TMatrixD matrix_stat_cov(rows, rows);
  for(int idx=0; idx<rows; idx++) matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);
  TMatrixD matrix_total_cov(rows, rows); matrix_total_cov = matrix_syst_abscov_temp + matrix_stat_cov;

  double determinant = matrix_total_cov.Determinant();
  if( determinant==0 ) {
    cerr<<endl<<" Error: determinant of total COV is 0"<<endl<<endl; exit(1);
  }
  
  double chi2 = GetChi2(matrix_pred_temp, matrix_meas_temp, matrix_syst_abscov_temp);
  cout<<endl<<TString::Format(" ---> check chi2_normal %7.2f", chi2)<<endl;
  
  TMatrixDSym DSmatrix_cov(rows);
  for(int ibin=0; ibin<rows; ibin++)
    for(int jbin=0; jbin<rows; jbin++)
      DSmatrix_cov(ibin, jbin) = matrix_total_cov(ibin, jbin);
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();
  TMatrixD matrix_eigenvector_T = matrix_eigenvector.T(); matrix_eigenvector.T();

  TMatrixD matrix_lambda_pred = matrix_pred_temp * matrix_eigenvector;
  TMatrixD matrix_lambda_meas = matrix_meas_temp * matrix_eigenvector;
  TMatrixD matrix_delta_lambda = matrix_lambda_meas - matrix_lambda_pred;
  TMatrixD matrix_delta_lambda_T = matrix_delta_lambda.T(); matrix_delta_lambda.T();

  TMatrixD matrix_cov_lambda(rows, rows);
  for(int idx=0; idx<rows; idx++) matrix_cov_lambda(idx, idx) = matrix_eigenvalue(idx);
  TMatrixD matrix_cov_lambda_inv = matrix_cov_lambda; matrix_cov_lambda_inv.Invert();
  double chi2_lambda = (matrix_delta_lambda * matrix_cov_lambda_inv * matrix_delta_lambda_T)(0,0);
  cout<<TString::Format(" ---> check chi2_lambda %7.2f", chi2_lambda)<<endl<<endl;

  //////////////////////////////////////////////// Plotting

  int color_pred = kRed;
  int color_meas = kBlack;

  TF1 *ff_1 = new TF1("ff_1", "1", 0, 1e3);
  ff_1->SetLineColor(kGray+1); ff_1->SetLineStyle(7);

  TF1 *ff_0 = new TF1("ff_0", "0", 0, 1e3);
  ff_0->SetLineColor(kGray+1); ff_0->SetLineStyle(7);

  TH1D *h1_fake_meas = new TH1D(TString::Format("h1_fake_meas_%d", index), "", rows, 0, rows);
  TGraphErrors *gh_fake_meas = new TGraphErrors();
  TH1D *h1_pred = new TH1D("h1_pred", "", rows, 0, rows);
  TH1D *h1_meas2pred = new TH1D(TString::Format("h1_meas2pred_%d", index), "", rows, 0, rows);
  TH1D *h1_meas2pred_syst = new TH1D(TString::Format("h1_meas2pred_syst_%d", index), "", rows, 0, rows);
  
  TH1D *h1_lambda_pred = new TH1D(TString::Format("h1_lambda_pred_%d", index), "", rows, 0, rows);
  TH1D *h1_lambda_meas = new TH1D(TString::Format("h1_lambda_meas_%d", index), "", rows, 0, rows);
  TH1D *h1_lambda_absigma_dis = new TH1D(TString::Format("h1_lambda_absigma_dis_%d", index), "", 10, 0, 10);
  TH1D *h1_lambda_sigma_iii = new TH1D(TString::Format("h1_lambda_sigma_iii_%d", index), "", rows, 0, rows);

  for(int ibin=1; ibin<=rows; ibin++) {
    double pred_val = matrix_pred_temp(0, ibin-1);
    double pred_err = sqrt( matrix_total_cov(ibin-1, ibin-1) );
    double meas_val = matrix_meas_temp(0, ibin-1);
    double meas2pred_val = 0;
    if( pred_val!=0 ) {
      h1_meas2pred->SetBinContent(ibin, meas_val/pred_val);
      h1_meas2pred->SetBinError(ibin, sqrt(meas_val)/pred_val);
      h1_meas2pred_syst->SetBinContent(ibin, 1);
      h1_meas2pred_syst->SetBinError(ibin, pred_err/pred_val);
    }
    
    double pred_lambda_val = matrix_lambda_pred(0, ibin-1);
    double pred_lambda_err = sqrt( matrix_cov_lambda(ibin-1, ibin-1) );
    double meas_lambda_val = matrix_lambda_meas(0, ibin-1);
    double delta_lambda_val = matrix_delta_lambda(0, ibin-1);
    double relerr_lambda = delta_lambda_val/pred_lambda_err;

    gh_fake_meas->SetPoint( ibin-1, h1_fake_meas->GetBinCenter(ibin), meas_val );
    gh_fake_meas->SetPointError( ibin-1, h1_fake_meas->GetBinWidth(ibin)*0.5, sqrt(meas_val) );
    h1_fake_meas->SetBinContent( ibin, meas_val );
    h1_pred->SetBinContent( ibin, pred_val );
    h1_pred->SetBinError( ibin, pred_err );

    h1_lambda_pred->SetBinContent( ibin, pred_lambda_val );
    h1_lambda_pred->SetBinError( ibin, pred_lambda_err );
    h1_lambda_meas->SetBinContent( ibin, meas_lambda_val );
    
    if( fabs(relerr_lambda)>=10 ) h1_lambda_absigma_dis->Fill( 9.5 );
    else h1_lambda_absigma_dis->Fill( fabs(relerr_lambda) );

    double edge_val = 6;
    double mod_relerr_lambda = relerr_lambda;
    if( fabs(relerr_lambda)>edge_val ) {
      (relerr_lambda>0) ? (mod_relerr_lambda = edge_val) : (mod_relerr_lambda = edge_val*(-1));
    }
    h1_lambda_sigma_iii->SetBinContent( ibin, mod_relerr_lambda );
    
  }

  /////////////////////////////////////////////////////////////

  double lambda_sigma_11 = h1_lambda_absigma_dis->GetBinContent( 1 );
  double lambda_sigma_12 = h1_lambda_absigma_dis->GetBinContent( 2 );
  double lambda_sigma_23 = h1_lambda_absigma_dis->GetBinContent( 3 );
  double lambda_sigma_3p = h1_lambda_absigma_dis->Integral(4, 10);

  const int user_num = 4;
  double array_user_meas[user_num] = { lambda_sigma_11, lambda_sigma_12, lambda_sigma_23, lambda_sigma_3p };
  double array_user_pred[user_num] = {0};
  double array_user_prob[user_num] = {1-0.3173, 0.3173-0.0455, 0.0455-0.0027, 0.0027};
  for(int idx=0; idx<user_num; idx++) {
    array_user_pred[idx] = rows * array_user_prob[idx];
  }

  double chi2_user[user_num] = {0};

  for(int idx=0; idx<user_num; idx++) {
    
    double meas = array_user_meas[idx];
    double pred = array_user_pred[idx];

    double chi2_sub = 0;
    if( meas==0 ) chi2_sub = 2*pred;
    else chi2_sub = 2*( pred - meas + meas * log(meas/pred) );
    
    if( idx<1 ) chi2_user[0] += chi2_sub;
    if( idx<2 ) chi2_user[1] += chi2_sub;
    if( idx<3 ) chi2_user[2] += chi2_sub;
    if( idx<4 ) chi2_user[3] += chi2_sub;    
  }

  cout<<endl;
  double user_pValue[user_num] = {0};
  for(int idx=0; idx<user_num; idx++) {
    user_pValue[idx] = TMath::Prob( chi2_user[idx], idx+1 );
    cout<<TString::Format(" ---> check user_pValue: chi2/ndf = %3.2f/%d, pValue %8.6f", chi2_user[idx], idx+1, user_pValue[idx])<<endl;
  }
  cout<<endl;
  
  
  /////////////////////////////////////////////////////////////
  
  roostr = TString::Format("canv_h1_fake_meas_%d", index);  
  TCanvas *canv_h1_fake_meas = new TCanvas(roostr, roostr, 1000, 950);
  //func_canv_margin(canv_h1_fake_meas, 0.15, 0.1, 0.1, 0.15);

  /////////////
  canv_h1_fake_meas->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_no->Draw(); pad_top_no->cd();
  
  TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone(TString::Format("h1_pred_clone_%d", index));
  h1_pred_clone->Draw("e2"); h1_pred_clone->SetMinimum(0);
  if( h1_pred_clone->GetMaximum() < h1_fake_meas->GetMaximum() ) h1_pred_clone->SetMaximum( h1_fake_meas->GetMaximum() * 1.1 );
  h1_pred_clone->SetFillColor(kRed-10); h1_pred_clone->SetFillStyle(1001);
  h1_pred_clone->SetMarkerSize(0);
  h1_pred_clone->SetLineColor(kRed);
  h1_pred_clone->GetXaxis()->SetLabelColor(10);
  func_xy_title(h1_pred_clone, "Bin index", "Entries"); func_title_size(h1_pred_clone, 0.065, 0.065, 0.065, 0.065);
  h1_pred_clone->GetXaxis()->CenterTitle(); h1_pred_clone->GetYaxis()->CenterTitle();
  h1_pred_clone->GetYaxis()->SetTitleOffset(1.2); 

  h1_pred->Draw("hist same"); h1_pred->SetLineColor(color_pred);

  gh_fake_meas->Draw("same p");
  gh_fake_meas->SetMarkerStyle(20); gh_fake_meas->SetMarkerSize(1.2);
  gh_fake_meas->SetMarkerColor(color_meas); gh_fake_meas->SetLineColor(color_meas);

  h1_pred_clone->Draw("same axis");
    
  TLegend *lg_chi2_toy = new TLegend(0.17+0.03, 0.6, 0.4+0.04, 0.85);    
  lg_chi2_toy->AddEntry(gh_fake_meas, TString::Format("#color[%d]{Fake data}", color_meas), "lep");
  lg_chi2_toy->AddEntry(h1_pred_clone, TString::Format("#color[%d]{Prediction}", color_pred), "fl");
  lg_chi2_toy->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %4.1f/%d}", color_pred, chi2, rows), "");
  lg_chi2_toy->Draw();
  lg_chi2_toy->SetBorderSize(0); lg_chi2_toy->SetFillStyle(0); lg_chi2_toy->SetTextSize(0.08);

  /////////////
  canv_h1_fake_meas->cd();
  TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_no->Draw(); pad_bot_no->cd();
  
  h1_meas2pred_syst->Draw("e2");
  h1_meas2pred_syst->SetMinimum(0); h1_meas2pred_syst->SetMaximum(2);
  h1_meas2pred_syst->SetFillColor(kRed-10); h1_meas2pred_syst->SetFillStyle(1001); h1_meas2pred_syst->SetMarkerSize(0);
  func_title_size(h1_meas2pred_syst, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_meas2pred_syst, "Measured bin index", "Data / Pred");
  h1_meas2pred_syst->GetXaxis()->SetTickLength(0.05);  h1_meas2pred_syst->GetXaxis()->SetLabelOffset(0.005);
  h1_meas2pred_syst->GetXaxis()->CenterTitle(); h1_meas2pred_syst->GetYaxis()->CenterTitle(); 
  h1_meas2pred_syst->GetYaxis()->SetTitleOffset(0.99);
  h1_meas2pred_syst->GetYaxis()->SetNdivisions(509);

  ff_1->Draw("same");
  
  h1_meas2pred->Draw("same e1");
  h1_meas2pred->SetLineColor(color_meas);
  //h1_meas2pred->SetMarkerSytle(20); h1_meas2pred->SetMarkerSize(1.2); h1_meas2pred->SetMarkerColor(color_meas); 

  if( saveFIG ) {
    roostr = TString::Format("canv_h1_fake_meas_%d.png", index);  
    canv_h1_fake_meas->SaveAs(roostr);
  }
    
  /////////////////////////////////////////////////////////////

  roostr = TString::Format("canv_h1_lambda_pred_%d", index);
  TCanvas *canv_h1_lambda_pred = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_lambda_pred, 0.15, 0.1, 0.1, 0.15);
    
  TH1D *h1_lambda_pred_clone = (TH1D*)h1_lambda_pred->Clone("h1_lambda_pred_clone");
  h1_lambda_pred_clone->Draw("e2");
  //h1_lambda_pred_clone->SetMinimum(0);
  if( h1_lambda_pred_clone->GetMaximum()<h1_lambda_meas->GetMaximum() ) h1_lambda_pred_clone->SetMaximum (h1_lambda_meas->GetMaximum() * 1.1 );
  h1_lambda_pred_clone->SetFillColor(kRed-10); h1_lambda_pred_clone->SetFillStyle(1001);
  h1_lambda_pred_clone->SetMarkerSize(0);
  h1_lambda_pred_clone->SetLineColor(color_pred);
  func_xy_title(h1_lambda_pred_clone, "Decomposed bin index", "Entries\'");
  func_title_size(h1_lambda_pred_clone, 0.05, 0.05, 0.05, 0.05);
  h1_lambda_pred_clone->GetXaxis()->CenterTitle(); h1_lambda_pred_clone->GetYaxis()->CenterTitle();
  h1_lambda_pred_clone->GetYaxis()->SetTitleOffset(1.5); 

  h1_lambda_pred->Draw("hist same"); h1_lambda_pred->SetLineColor(color_pred);

  h1_lambda_meas->Draw("same p");
  h1_lambda_meas->SetMarkerStyle(20); h1_lambda_meas->SetMarkerSize(1.2); h1_lambda_meas->SetMarkerColor(color_meas);
  h1_lambda_meas->SetLineColor(color_meas);    

  h1_lambda_pred_clone->Draw("same axis");   

  if( saveFIG ) {
    roostr = TString::Format("canv_h1_lambda_pred_%d.png", index);
    canv_h1_lambda_pred->SaveAs(roostr);
  }
    
  /////////////////////////////////////////////////////////////

  roostr = TString::Format("canv_lambda_sigma_%d", index);
  TCanvas *canv_lambda_sigma = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_lambda_sigma, 0.15, 0.1, 0.1, 0.15);
    
  TH1D *h1_lambda_sigmax3 = new TH1D(TString::Format("h1_lambda_sigmax3_%d", index), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax3->SetBinError(ibin, 3);
  
  TH1D *h1_lambda_sigmax2 = new TH1D(TString::Format("h1_lambda_sigmax2_%d", index), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax2->SetBinError(ibin, 2);
  
  TH1D *h1_lambda_sigmax1 = new TH1D(TString::Format("h1_lambda_sigmax1_%d", index), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax1->SetBinError(ibin, 1);
  
  h1_lambda_sigmax3->Draw("e2");
  h1_lambda_sigmax3->SetMinimum(-6);
  h1_lambda_sigmax3->SetMaximum(6);
  h1_lambda_sigmax3->SetFillColor(kRed-9); h1_lambda_sigmax3->SetFillStyle(1001);
  h1_lambda_sigmax3->SetMarkerSize(0);
  h1_lambda_sigmax3->SetLineColor(color_pred);
  func_xy_title(h1_lambda_sigmax3, "Decomposed bin index", "#sigma_{i}\' value");
  func_title_size(h1_lambda_sigmax3, 0.05, 0.05, 0.05, 0.05);
  h1_lambda_sigmax3->GetXaxis()->CenterTitle(); h1_lambda_sigmax3->GetYaxis()->CenterTitle();
  h1_lambda_sigmax3->GetYaxis()->SetTitleOffset(1.2); 
  h1_lambda_sigmax3->GetYaxis()->SetNdivisions(112);

  h1_lambda_sigmax2->Draw("same e2"); h1_lambda_sigmax2->SetMarkerSize(0);
  h1_lambda_sigmax2->SetFillColor(kYellow); h1_lambda_sigmax2->SetFillStyle(1001);

  h1_lambda_sigmax1->Draw("same e2"); h1_lambda_sigmax1->SetMarkerSize(0);
  h1_lambda_sigmax1->SetFillColor(kGreen); h1_lambda_sigmax1->SetFillStyle(1001);

  h1_lambda_sigma_iii->Draw("same p");
  h1_lambda_sigma_iii->SetMarkerSize(1.4); h1_lambda_sigma_iii->SetMarkerStyle(21); h1_lambda_sigma_iii->SetMarkerColor(color_meas);
    
  ff_0->Draw("same");
  
  h1_lambda_sigmax2->Draw("same axis");
  
  TLegend *lg_lambda_sigma = new TLegend(0.35, 0.65+0.065, 0.7, 0.89+0.065);
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{Total num: %d}", kBlue, rows), "");
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| #in (1, 2]: %1.0f, expect. %3.2f}",
						kBlue, lambda_sigma_12, rows*0.2718), "");
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| #in (2, 3]: %1.0f, expect. %3.2f}",
						kBlue, lambda_sigma_23, rows*0.0428), "");
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| > 3: %1.0f, expect. %3.2f}",
						kBlue, lambda_sigma_3p, rows*0.0027), "");  
  lg_lambda_sigma->Draw();
  lg_lambda_sigma->SetBorderSize(0); lg_lambda_sigma->SetFillStyle(0); lg_lambda_sigma->SetTextSize(0.05);

  if( saveFIG ) {
    roostr = TString::Format("canv_lambda_sigma_%d.png", index);
    canv_lambda_sigma->SaveAs(roostr);
  }

  /////////////////////////////////    
    
  roostr = TString::Format("h2_lambda_transform_matrix_%d", index);
  TH2D *h2_lambda_transform_matrix = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double value = matrix_eigenvector(ibin-1, jbin-1);
      h2_lambda_transform_matrix->SetBinContent(ibin, jbin, value);
    }
  }

  roostr = TString::Format("canv_h2_lambda_transform_matrix_%d", index);
  TCanvas *canv_h2_lambda_transform_matrix = new TCanvas(roostr, roostr, 900, 850);
  func_canv_margin(canv_h2_lambda_transform_matrix, 0.15, 0.2,0.15,0.2);
  h2_lambda_transform_matrix->Draw("colz");
  func_xy_title(h2_lambda_transform_matrix, "Measured bin index", "Decomposed bin index");
  func_title_size(h2_lambda_transform_matrix, 0.05, 0.05, 0.05, 0.05);
  h2_lambda_transform_matrix->GetXaxis()->CenterTitle(); h2_lambda_transform_matrix->GetYaxis()->CenterTitle();
  if( saveFIG ) {
    roostr = TString::Format("canv_h2_lambda_transform_matrix_%d.png", index);
    canv_h2_lambda_transform_matrix->SaveAs(roostr);
  }
}

/////////////////////// ccc

double TCN::GetChi2(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp)
{
  double chi2 = 0;
  
  TMatrixD matrix_delta = matrix_pred_temp - matrix_meas_temp;
  TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T(); 

  int rows = matrix_pred_temp.GetNcols();
  
  TMatrixD matrix_stat_cov(rows, rows);
  for(int idx=0; idx<rows; idx++) matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);

  TMatrixD matrix_total_cov(rows, rows); matrix_total_cov = matrix_syst_abscov_temp + matrix_stat_cov;
  TMatrixD matrix_total_cov_inv = matrix_total_cov; matrix_total_cov_inv.Invert();
  chi2 = ( matrix_delta*matrix_total_cov_inv*matrix_delta_T )(0,0);
  
  return chi2;  
}

/////////////////////// ccc

void TCN::ProduceVariation(int ntoys, bool flag_norm)
{
  TString roostr = "";
  
  TOYS = ntoys;

  map_toy_meas.clear();

  cout<<endl<<" ---> Check, producing toys"<<endl<<endl;
    
  TMatrixD matrix_stat_abscov(BINS, BINS); for(int idx=0; idx<BINS; idx++) matrix_stat_abscov(idx, idx) = matrix_pred(0, idx);
  TMatrixD matrix_total_abscov(BINS, BINS); matrix_total_abscov = matrix_stat_abscov + matrix_syst_abscov;

  ////////////////////////////////
    
  TMatrixDSym DSmatrix_cov(BINS);
  for(int ibin=0; ibin<BINS; ibin++) {
    for(int jbin=0; jbin<BINS; jbin++) {
      DSmatrix_cov(ibin, jbin) = matrix_total_abscov(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

  ///////
  // TPrincipal principal_obj_total(BINS, "ND");
    
  for(int itoy=1; itoy<=TOYS; itoy++) {    
    TMatrixD matrix_element(BINS, 1);
    // double *array_obj_total = new double[BINS];
      
    for(int idx=0; idx<BINS; idx++) {
      if( matrix_eigenvalue(idx)>=0 ) {
	matrix_element(idx,0) = rand->Gaus( 0, sqrt( matrix_eigenvalue(idx) ) );
      }
      else {
	matrix_element(idx,0) = 0;
      }      
    }
      
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;
      
    for(int idx=0; idx<BINS; idx++) {
      double val_with_syst = matrix_variation(idx,0) + matrix_pred(0, idx);// key point
      if( val_with_syst<0 ) val_with_syst = 0;// ??? remove this requirement or not
      map_toy_meas[itoy][idx+1] = val_with_syst;
      // array_obj_total[idx] = val_with_syst;
    }// idx

    // principal_obj_total.AddRow( array_obj_total );
    // delete[] array_obj_total;
      
  }// itoy

  // //////////////////////////////////////////////////////////// obj_total
    
  // TMatrixD *matrix_abscov_obj_total = (TMatrixD *)principal_obj_total.GetCovarianceMatrix();
    
  // for(int idx=0; idx<BINS; idx++) {
  //   for(int jdx=0; jdx<BINS; jdx++) {
  //     if(idx<jdx) (*matrix_abscov_obj_total)(idx, jdx) = (*matrix_abscov_obj_total)(jdx, idx);
  //   }
  // }

  // ///////////
  // TMatrixD matrix_relcov_obj_total(BINS, BINS);
  // TMatrixD matrix_correlation_obj_total(BINS, BINS);
    
  // for(int idx=0; idx<BINS; idx++) {
  //   for(int jdx=0; jdx<BINS; jdx++) {
  //     double cov_ij = (*matrix_abscov_obj_total)(idx, jdx);
  //     double cov_ii = (*matrix_abscov_obj_total)(idx, idx);
  //     double cov_jj = (*matrix_abscov_obj_total)(jdx, jdx);        
  //     double cv_i = matrix_pred(0, idx);
  //     double cv_j = matrix_pred(0, jdx);        
  //     if( cv_i!=0 && cv_j!=0 ) matrix_relcov_obj_total(idx, jdx) = cov_ij/cv_i/cv_j;        
  //     if( cov_ii!=0 && cov_jj!=0 ) matrix_correlation_obj_total(idx, jdx) = cov_ij/sqrt(cov_ii)/sqrt(cov_jj);
  //     //if( idx==jdx ) matrix_correlation_obj_total(idx, jdx) = 1;       
  //   }
  // }
    
  // ///////////
  // roostr = "canv_matrix_relcov_obj_total";
  // TCanvas *canv_matrix_relcov_obj_total = new TCanvas(roostr, roostr, 900, 850);
  // func_canv_margin(canv_matrix_relcov_obj_total, 0.15, 0.2,0.15,0.2);
  // roostr = "h2_relcov_obj_total";
  // TH2D *h2_relcov_obj_total = new TH2D(roostr, "", BINS, 0, BINS, BINS, 0, BINS);
  // for(int idx=0; idx<BINS; idx++)
  //   for(int jdx=0; jdx<BINS; jdx++)
  //  h2_relcov_obj_total->SetBinContent( idx+1, jdx+1, matrix_relcov_obj_total(idx, jdx) );
  // h2_relcov_obj_total->Draw("colz");
  // func_xy_title(h2_relcov_obj_total, "Bin index", "Bin index");
   
  // ///////////
  // roostr = "canv_matrix_correlation_obj_total";
  // TCanvas *canv_matrix_correlation_obj_total = new TCanvas(roostr, roostr, 900, 850);
  // func_canv_margin(canv_matrix_correlation_obj_total, 0.15, 0.2,0.15,0.2);
  // roostr = "h2_correlation_obj_total";
  // TH2D *h2_correlation_obj_total = new TH2D(roostr, "", BINS, 0, BINS, BINS, 0, BINS);
  // for(int idx=0; idx<BINS; idx++)
  //   for(int jdx=0; jdx<BINS; jdx++)
  //  h2_correlation_obj_total->SetBinContent( idx+1, jdx+1, matrix_correlation_obj_total(idx, jdx) );
  // h2_correlation_obj_total->Draw("colz");
  // func_xy_title(h2_correlation_obj_total, "Bin index", "Bin index");
  // h2_correlation_obj_total->GetZaxis()->SetRangeUser(-1,1);    
  
}

/////////////////////// ccc

void TCN::Initialization(bool flag_norm)
{
  TString roostr = "";
  
  FLAG_NORM = flag_norm;
  
  if( FLAG_NORM ) {
    
  }// if( FLAG_NORM )
  else {
    
    ////////////////////////////////

    TFile *file_numu = new TFile("file_numu.root", "read");
  
    TMatrixD *matrix_gof_pred = (TMatrixD*)file_numu->Get("matrix_gof_pred");
    TMatrixD *matrix_gof_meas = (TMatrixD*)file_numu->Get("matrix_gof_meas");
    TMatrixD *matrix_gof_syst = (TMatrixD*)file_numu->Get("matrix_gof_syst");
    
    int rows = matrix_gof_syst->GetNrows();
    BINS = rows;

    matrix_pred.Clear(); matrix_pred.ResizeTo(1, rows); matrix_pred = (*matrix_gof_pred);
    matrix_meas.Clear(); matrix_meas.ResizeTo(1, rows); matrix_meas = (*matrix_gof_meas);
    matrix_syst_abscov.Clear(); matrix_syst_abscov.ResizeTo(rows, rows); matrix_syst_abscov = (*matrix_gof_syst);
  
    TH2D *h2_gof_abscov = new TH2D("h2_gof_abscov", "", rows, 0, rows, rows, 0, rows);
    TH2D *h2_gof_relcov = new TH2D("h2_gof_relcov", "", rows, 0, rows, rows, 0, rows);
    TH2D *h2_gof_correlation = new TH2D("h2_gof_correlation", "", rows, 0, rows, rows, 0, rows);

    for(int ibin=1; ibin<=rows; ibin++ ) {
      for(int jbin=1; jbin<=rows; jbin++) {
	double cov_ij = (*matrix_gof_syst)(ibin-1, jbin-1);
	double cov_ii = (*matrix_gof_syst)(ibin-1, ibin-1);
	double cov_jj = (*matrix_gof_syst)(jbin-1, jbin-1);

	double cv_i = matrix_pred(0, ibin-1);
	double cv_j = matrix_pred(0, jbin-1);

	h2_gof_abscov->SetBinContent(ibin, jbin, cov_ij);
	if( cv_i!=0 && cv_j!=0 ) h2_gof_relcov->SetBinContent(ibin, jbin, cov_ij/cv_i/cv_j);
	if( cov_ii!=0 && cov_jj!=0 ) h2_gof_correlation->SetBinContent(ibin, jbin, cov_ij/sqrt(cov_ii)/sqrt(cov_jj) );	
      }
    }
    
    // ///////////
    // roostr = "canv_h2_gof_relcov";
    // TCanvas *canv_h2_gof_relcov = new TCanvas(roostr, roostr, 900, 850);
    // func_canv_margin(canv_h2_gof_relcov, 0.15, 0.2,0.15,0.2);
    // h2_gof_relcov->Draw("colz");
    // func_xy_title(h2_gof_relcov, "Bin index", "Bin index");
    // h2_gof_relcov->GetXaxis()->CenterTitle(); h2_gof_relcov->GetYaxis()->CenterTitle();
    // canv_h2_gof_relcov->SaveAs("canv_h2_gof_relcov.png");
      
    // ///////////
    // roostr = "canv_h2_gof_correlation";
    // TCanvas *canv_h2_gof_correlation = new TCanvas(roostr, roostr, 900, 850);
    // func_canv_margin(canv_h2_gof_correlation, 0.15, 0.2,0.15,0.2);
    // h2_gof_correlation->Draw("colz");
    // h2_gof_correlation->GetZaxis()->SetRangeUser(-1,1);
    // func_xy_title(h2_gof_correlation, "Bin index", "Bin index");       
    // h2_gof_correlation->GetXaxis()->CenterTitle(); h2_gof_correlation->GetYaxis()->CenterTitle();
    // canv_h2_gof_correlation->SaveAs("canv_h2_gof_correlation.png");
      
  }// else of if( FLAG_NORM )
    
}


////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void read_obj()
{
  //////////////////////////////////////////////////////////////////////////////////////// Draw style

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  TString roostr = "";
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  bool flag_norm     = 0;  
  double norm_relerr = 0.2;
  
  TCN *testcn = new TCN();
  
  testcn->Initialization( flag_norm ); 
  testcn->Set_norm_relerr( norm_relerr );

  if( 1 ) {
    testcn->SetMeas2Use();
    double chi2_meas = testcn->GetChi2( testcn->matrix_pred, testcn->matrix_fake_meas, testcn->matrix_syst_abscov );
    cout<<endl<<" ---> check chi2_meas: "<<chi2_meas<<endl<<endl;
    
    TMatrixD matrix_meas_temp = testcn->matrix_fake_meas;    
    testcn->Plotting_singlecase(testcn->matrix_pred, matrix_meas_temp, testcn->matrix_syst_abscov, 1, 1);
  }

  if( 0 ) {
    int ntoys = 100;
    testcn->ProduceVariation( ntoys, flag_norm );
    for(int idx=1; idx<=100; idx++) {
      testcn->SetToy(idx);
      double chi2 = testcn->GetChi2( testcn->matrix_pred, testcn->matrix_fake_meas, testcn->matrix_syst_abscov );
      //cout<<TString::Format(" ---> itoy %4d, chi2 %7.2f", idx, chi2)<<endl;
    }

    int itoy = 90;
    testcn->SetToy(itoy);    
    TMatrixD matrix_meas_temp = testcn->matrix_fake_meas;    
    testcn->Plotting_singlecase(testcn->matrix_pred, matrix_meas_temp, testcn->matrix_syst_abscov, 1, 101);
  }

  //////////////////////////////////////////////////////////////////////////////
  
  /*
  int ntoys = 100000;
  testcn->ProduceVariation( ntoys, flag_norm );
  
  double low_chi2 = 0;
  double hgh_chi2 = 120;
  int bins_chi2 = 200;
  int ndf_chi2 = testcn->BINS;
  
  roostr = "h1_toy_chi2";
  TH1D *h1_toy_chi2 = new TH1D(roostr, "", bins_chi2, low_chi2, hgh_chi2);
  
  for(int idx=1; idx<=ntoys; idx++) {
    if( idx%(ntoys/10)==0 ) cout<<TString::Format(" ---> processing %8d, %8.4f", idx, idx*1./ntoys)<<endl;
    testcn->SetToy(idx);
    double chi2 = testcn->GetChi2( testcn->matrix_pred, testcn->matrix_fake_meas, testcn->matrix_syst_abscov );
    h1_toy_chi2->Fill( chi2 );
  }
  h1_toy_chi2->Scale( 1./h1_toy_chi2->Integral() );
  
  roostr = "canv_h1_toy_chi2";
  TCanvas *canv_h1_toy_chi2 = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_toy_chi2, 0.15, 0.1, 0.1, 0.15);
  h1_toy_chi2->Draw("hist");
  func_xy_title(h1_toy_chi2, "Chi2", "PDF");
  func_title_size(h1_toy_chi2, 0.05, 0.05, 0.05, 0.05);
  h1_toy_chi2->GetXaxis()->CenterTitle(); h1_toy_chi2->GetYaxis()->CenterTitle();
  h1_toy_chi2->GetYaxis()->SetTitleOffset(1.5);
  cout<<endl<<" ---> h1_toy_chi2 mean: "<<h1_toy_chi2->GetMean()<<endl<<endl;
  
  double width_chi2 = (hgh_chi2-low_chi2)/bins_chi2;
  TF1 *f1_chi2_temp = new TF1("f1_chi2_temp", TString::Format("ROOT::Math::chisquared_pdf(x,%d,0)*%f", ndf_chi2, width_chi2),
			      low_chi2, hgh_chi2);
  f1_chi2_temp->Draw("same");

  h1_toy_chi2->Draw("same axis");
  
  if( f1_chi2_temp->GetMaximum() > h1_toy_chi2->GetMaximum() ) {
    h1_toy_chi2->SetMaximum( f1_chi2_temp->GetMaximum() * 1.1 );
  }

  TLegend *lg_chi2_toy = new TLegend(0.6, 0.70, 0.83, 0.85);
  lg_chi2_toy->AddEntry(h1_toy_chi2, TString::Format("#color[%d]{Toy result}", kBlue), "l");
  lg_chi2_toy->AddEntry(f1_chi2_temp, TString::Format("#color[%d]{#chi^{2}(%d)}", kRed, ndf_chi2), "l");  
  lg_chi2_toy->Draw();
  lg_chi2_toy->SetBorderSize(0); lg_chi2_toy->SetFillStyle(0); lg_chi2_toy->SetTextSize(0.065);
  
  canv_h1_toy_chi2->SaveAs("canv_h1_toy_chi2.png");
  */
}
