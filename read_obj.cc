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

    FLAG_STATCOV_NO = 0;
  }

  TRandom3 *rand;

  bool FLAG_NORM;
  int BINS;
  int TOYS;

  bool FLAG_STATCOV_NO;
  
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
  //TMatrixD matrix_stat_cov(rows, rows);
  //for(int idx=0; idx<rows; idx++) matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);


  for(int idx=0; idx<rows; idx++) {
    double val_pred = matrix_pred_temp(0, idx);
    double val_meas = matrix_meas_temp(0, idx);
    //cout<<" ---> check "<<idx+1<<"\t"<<val_meas<<"\t"<<val_pred<<endl;
  }

  
  /// docDB 32520, when the prediction is sufficiently low
  double array_pred_protect[11] = {0, 0.461, 0.916, 1.382, 1.833, 2.298, 2.767, 3.225, 3.669, 4.141, 4.599};
  
  TMatrixD matrix_stat_cov(rows, rows);
  if( !FLAG_STATCOV_NO ) {
    for(int idx=0; idx<rows; idx++) {
      matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);// Pearson
    
      double val_meas = matrix_meas_temp(0, idx);
      double val_pred = matrix_pred_temp(0, idx);
      int int_meas = (int)(val_meas+0.1);

      if( int_meas>=1 && int_meas<=10) {
	if( val_pred<array_pred_protect[int_meas] ) {
	  double numerator = pow(val_pred-val_meas, 2);
	  double denominator = 2*( val_pred - val_meas + val_meas*log(val_meas/val_pred) );
	  matrix_stat_cov(idx, idx) = numerator/denominator;
	}
      }    
    }
  }
    
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

  TH1D *h1_percentage_cov_lambda = new TH1D(TString::Format("h1_percentage_cov_lambda_%d", index), "", rows, 0, rows);

  TH1D *h1_lambda_epsilon = new TH1D(TString::Format("h1_lambda_epsilon_%d", index), "", 20, -3, 3.5);
  
  map<int, double>map_above3sigma;
  vector<double>vec_above3sigma;

  double sum_cov_lambda = 0;
  for(int ibin=1; ibin<=rows; ibin++) {
    sum_cov_lambda += matrix_cov_lambda(ibin-1, ibin-1);
  }
  
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

    h1_percentage_cov_lambda->SetBinContent( ibin, matrix_cov_lambda(ibin-1, ibin-1)*100./sum_cov_lambda );
    
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

    if( fabs(relerr_lambda)>=3 ) {
      map_above3sigma[ibin] = fabs(relerr_lambda);
      vec_above3sigma.push_back( fabs(relerr_lambda) );
    }

    h1_lambda_epsilon->Fill( relerr_lambda );
    
  }

  cout<<endl;
  cout<<" bin above 3sigma"<<endl;
  for(auto it_map_above3sigma=map_above3sigma.begin(); it_map_above3sigma!=map_above3sigma.end(); it_map_above3sigma++) {
    int user_index = it_map_above3sigma->first;
    double user_value = it_map_above3sigma->second;
    cout<<TString::Format(" ---> %2d, %3.1f", user_index, user_value)<<endl;
  }
  cout<<endl;

  ////////////////////// Fisher's method
  /*
  double sum_chi2_lambda_epsilon = 0;
  double sum_chi2_fisher = 0;
  double sum_chi2_pearson = 0;


  const int p4_power = 20;
  double sum_p4_data = 0;
  
  for(int ibin=1; ibin<=rows; ibin++) {

    ////
    double lambda_epsilon = h1_lambda_sigma_iii->GetBinContent( ibin );
    double chi2_sub = lambda_epsilon*lambda_epsilon;
    sum_chi2_lambda_epsilon += chi2_sub;
    double pvalue_sub = TMath::Prob( chi2_sub, 1 );
    
    ////
    double chi2_sub_fisher = ( -2*log(pvalue_sub) );
    sum_chi2_fisher += chi2_sub_fisher;

    ////
    double chi2_sub_pearson = ( -2*log(1-pvalue_sub) );
    sum_chi2_pearson += chi2_sub_pearson;

    ////
    sum_p4_data += pow( lambda_epsilon, p4_power );
    
    ////
    //cout<<TString::Format(" ---> bin %3d, epsilon %8.4f, chi2 %7.4f, pvalue %8.6f, chi2_fisher %6.2f, chi2_pearson %6.2f",
    //			  ibin, lambda_epsilon, chi2_sub, pvalue_sub, chi2_sub_fisher, chi2_sub_pearson)<<endl;
  }
  
  double pvalue_normal = TMath::Prob( sum_chi2_lambda_epsilon, rows );
  double significance_normal = sqrt( TMath::ChisquareQuantile( 1-pvalue_normal, 1 ) );;
  
  double pvalue_fisher = TMath::Prob( sum_chi2_fisher, 2*rows );
  double significance_fisher = sqrt( TMath::ChisquareQuantile( 1-pvalue_fisher, 1 ) );;
  
  double pvalue_pearson = TMath::Prob( sum_chi2_pearson, 2*rows );
  double significance_pearson = sqrt( TMath::ChisquareQuantile( 1-pvalue_pearson, 1 ) );;
  
  cout<<" ---> sum_chi2_lambda_epsilon "<<sum_chi2_lambda_epsilon<<", significance "<<significance_normal<<endl;
  cout<<" ---> sum_chi2_fisher         "<<sum_chi2_fisher<<", significance "<<significance_fisher<<endl;
  cout<<" ---> sum_chi2_pearson        "<<sum_chi2_pearson<<", significance "<<significance_pearson<<endl;
  cout<<endl;


  ///////

  int ntoys = 1000000;

  vector<double>vc_sum_p4;

  for(int itoy=1; itoy<=ntoys; itoy++) {
    if( itoy%1000000==0 ) cout<<" ---> processing "<<itoy*1./ntoys<<endl;
    
    double sum_p4 = 0;
    
    for(int idx=1; idx<=rows; idx++) {
      double val_rand = rand->Gaus(0, 1);
      sum_p4 += pow(val_rand, p4_power);
    }

    vc_sum_p4.push_back( sum_p4 );
  }

  sort( vc_sum_p4.begin(), vc_sum_p4.end() );
  int size_sum_p4 = vc_sum_p4.size();

  
  double xmin_sum_p4 = 0;
  if( vc_sum_p4.at(0)>=1 ) xmin_sum_p4 = int( vc_sum_p4.at(0) ) - 1;
  double xmax_sum_p4 = int(vc_sum_p4.at( size_sum_p4-1 )) + 1;
  
  roostr = TString::Format("h1_dist_sum_p4_%d", index);
  TH1D *h1_dist_sum_p4 = new TH1D(roostr, "", 100, xmin_sum_p4, xmax_sum_p4);

  for(int idx=0; idx<size_sum_p4; idx++) {
    double val = vc_sum_p4.at( idx );
    h1_dist_sum_p4->Fill( val );
  }


  cout<<" ---> check sum_p4_data "<<sum_p4_data<<endl;
  int line_p4_pvalue = 0;
  for(int idx=size_sum_p4-1; idx>=0; idx--) {
    double val = vc_sum_p4.at( idx );
    if( val>=sum_p4_data ) line_p4_pvalue++;
    else break;
  }
  double pvalue_p4 = line_p4_pvalue*1./ntoys;
  double significance_p4 = sqrt( TMath::ChisquareQuantile( 1-pvalue_p4, 1 ) );;
  cout<<" ---> significance p4 "<<significance_p4<<endl;
  
  roostr = TString::Format("canv_h1_dist_sum_p4_%d", index);
  TCanvas *canv_h1_dist_sum_p4 = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_dist_sum_p4, 0.15, 0.1, 0.1, 0.15);
    
  h1_dist_sum_p4->Draw();
  
  roostr = TString::Format("canv_h1_dist_sum_p4_%d.png", index);
  canv_h1_dist_sum_p4->SaveAs(roostr);
  */
  /////////////////////////////////////////////////////////////
 
  double lambda_sigma_11 = h1_lambda_absigma_dis->GetBinContent( 1 );
  double lambda_sigma_12 = h1_lambda_absigma_dis->GetBinContent( 2 );
  double lambda_sigma_23 = h1_lambda_absigma_dis->GetBinContent( 3 );
  double lambda_sigma_3p = h1_lambda_absigma_dis->Integral(4, 10);

  /////////////////////////////////////////////////////////////
  
  roostr = TString::Format("canv_h1_fake_meas_%d", index);  
  TCanvas *canv_h1_fake_meas = new TCanvas(roostr, roostr, 1000, 950);

  /////////////
  canv_h1_fake_meas->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_no->Draw(); pad_top_no->cd();
  
  TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone(TString::Format("h1_pred_clone_%d", index));
  h1_pred_clone->Draw("e2"); h1_pred_clone->SetMinimum(0);
  double max_h1_pred_clone = 0;
  double max_h1_fake_meas = h1_fake_meas->GetMaximum();
  for(int ibin=1; ibin<=rows; ibin++) {
    double sub_pred = h1_pred_clone->GetBinContent(ibin) + h1_pred_clone->GetBinError(ibin);
    if( max_h1_pred_clone<sub_pred ) max_h1_pred_clone = sub_pred;
  }
  if( max_h1_pred_clone<max_h1_fake_meas ) h1_pred_clone->SetMaximum( max_h1_fake_meas*1.1 );
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

  double shift_chi2_toy_YY = 0;// -0.45
  double shift_chi2_toy_XX = 0;// 0.3
  TLegend *lg_chi2_toy = new TLegend(0.17+0.03+shift_chi2_toy_XX, 0.6+shift_chi2_toy_YY, 0.4+0.04+shift_chi2_toy_XX, 0.85+shift_chi2_toy_YY);    
  lg_chi2_toy->AddEntry(gh_fake_meas, TString::Format("#color[%d]{Data}", color_meas), "lep");
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
  func_xy_title(h1_meas2pred_syst, "Measurement bin index", "Data / Pred");
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
  
  roostr = TString::Format("canv_h1_lambda_epsilon_%d", index);  
  TCanvas *canv_h1_lambda_epsilon = new TCanvas(roostr, roostr, 1000, 950);

  h1_lambda_epsilon->Draw();
  
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
  func_xy_title(h1_lambda_pred_clone, "Decomposition bin index", "Entries\'");
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
  func_xy_title(h1_lambda_sigmax3, "Decomposition bin index", "#epsilon_{i}\' value");
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
  
  h1_lambda_sigmax3->Draw("same axis");

  double shift_x_lambda_sigma = -0.15;
  TLegend *lg_lambda_sigma = new TLegend(0.35+shift_x_lambda_sigma, 0.75, 0.7+shift_x_lambda_sigma, 0.85+0.02);
  // lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{Total num: %d}", kBlue, rows), "");
  // lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| #in (1, 2]: %1.0f, expect. %3.2f}",
  // 						kBlue, lambda_sigma_12, rows*0.2718), "");
  // lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| #in (2, 3]: %1.0f, expect. %3.2f}",
  // 						kBlue, lambda_sigma_23, rows*0.0428), "");
  // lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{|#sigma_{i}\'| > 3: %1.0f, expect. %3.2f}",
  // 						kBlue, lambda_sigma_3p, rows*0.0027), "");

  
  double pvalue_default = TMath::Prob( chi2, rows );
  double sigma_default = sqrt( TMath::ChisquareQuantile( 1-pvalue_default, 1 ) );
  double pvalue_global = 0;
  double sigma_global = 0;
  double sigma_global_AA = 0;
  double sigma_global_BB = 0;
  double sum_AA = 0;
  
  if( (int)(map_above3sigma.size())>=1 ) {    
    if( (int)(map_above3sigma.size())==1 ) {
      double chi2_local = pow(map_above3sigma.begin()->second, 2);
      double pvalue_local = TMath::Prob( chi2_local, 1 );
      pvalue_global = 1 - pow(1-pvalue_local, rows);
      sigma_global = sqrt( TMath::ChisquareQuantile( 1-pvalue_global, 1 ) );

      sum_AA = chi2_local;
    }
    else {      
      int user_vec_size = vec_above3sigma.size();      
      sum_AA = 0;
      for(int idx=0; idx<user_vec_size; idx++) {
	sum_AA += pow( vec_above3sigma.at(idx), 2 );	
      }
      double pvalue_local_AA = TMath::Prob( sum_AA, user_vec_size );      
      double coeff = TMath::Factorial(rows)/TMath::Factorial(rows-user_vec_size)/TMath::Factorial(user_vec_size);
      double pvalue_global_AA = coeff*pvalue_local_AA;
      cout<<" ---> check coeff "<<coeff<<"\t"<<pvalue_local_AA<<"\t"<<coeff*pvalue_local_AA<<"\t"<<1-pvalue_global_AA<<endl;
      sigma_global = sqrt( TMath::ChisquareQuantile( 1-pvalue_global_AA, 1 ) );      
    }
  }
  
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{%3.1f#sigma:        overall #chi^{2}/dof: %3.1f/%d}", kBlue, sigma_default, chi2, rows), "");
  if( lambda_sigma_3p>=1 ) {
    lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{%3.1f#sigma (LEE corr.): #chi^{2}/dof: %3.1f/%d (|#epsilon_{i}\'|>3)}",
						  kRed, sigma_global, sum_AA, (int)(map_above3sigma.size())), "");
  }
  else {
    lg_lambda_sigma->AddEntry("", "", "");
  }
  lg_lambda_sigma->Draw();
  lg_lambda_sigma->SetBorderSize(0); lg_lambda_sigma->SetFillStyle(0); lg_lambda_sigma->SetTextSize(0.05);

  if( saveFIG ) {
    roostr = TString::Format("canv_lambda_sigma_%d.png", index);
    canv_lambda_sigma->SaveAs(roostr);
  }

  cout<<endl<<" ---> LEE "<<sigma_global<<endl<<endl;
  
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
  func_xy_title(h2_lambda_transform_matrix, "Measurement bin index", "Decomposition bin index");
  func_title_size(h2_lambda_transform_matrix, 0.05, 0.05, 0.05, 0.05);
  h2_lambda_transform_matrix->GetXaxis()->CenterTitle(); h2_lambda_transform_matrix->GetYaxis()->CenterTitle();
  if( saveFIG ) {
    roostr = TString::Format("canv_h2_lambda_transform_matrix_%d.png", index);
    canv_h2_lambda_transform_matrix->SaveAs(roostr);
  }
     
  /////////////////////////////////////////////////////////////
  
  roostr = TString::Format("canv_comb_%d", index);  
  TCanvas *canv_comb = new TCanvas(roostr, roostr, 1000, 950);

  /////////////
  canv_comb->cd();
  TPad *pad_top_comb_no = new TPad("pad_top_comb_no", "pad_top_comb_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_comb_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_comb_no->Draw(); pad_top_comb_no->cd();
    
  h1_lambda_sigmax3->Draw("e2");
  h1_lambda_sigmax2->Draw("same e2");
  h1_lambda_sigmax1->Draw("same e2");
  h1_lambda_sigma_iii->Draw("same p");    
  ff_0->Draw("same");  
  h1_lambda_sigmax3->Draw("same axis");

  func_title_size(h1_lambda_sigmax3, 0.065, 0.065, 0.065, 0.065);
  h1_lambda_sigmax3->GetYaxis()->SetTitleOffset(1.05);
  h1_lambda_sigmax3->GetXaxis()->SetLabelColor(10);
  h1_lambda_sigmax3->GetXaxis()->SetTickLength(0.042);
  h1_lambda_sigmax3->GetYaxis()->SetTickLength(0.025);
  
  lg_lambda_sigma->Draw();
  lg_lambda_sigma->SetTextSize(0.065);
  lg_lambda_sigma->SetY1(0.7);
    
  /////////////
  canv_comb->cd();
  TPad *pad_bot_comb_no = new TPad("pad_bot_comb_no", "pad_bot_comb_no", 0, 0, 1, 0.45);
  pad_bot_comb_no->SetLogy();
  func_canv_margin(pad_bot_comb_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_comb_no->Draw(); pad_bot_comb_no->cd();

  h1_percentage_cov_lambda->Draw("hist");
  h1_percentage_cov_lambda->SetLineWidth(3);
    
  func_title_size(h1_percentage_cov_lambda, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_percentage_cov_lambda, "Decomposition bin index", "#Lambda_{i} / #sum#Lambda_{i} (%)");
  h1_percentage_cov_lambda->GetXaxis()->SetTickLength(0.05);  h1_percentage_cov_lambda->GetXaxis()->SetLabelOffset(0.005);
  h1_percentage_cov_lambda->GetXaxis()->CenterTitle(); h1_percentage_cov_lambda->GetYaxis()->CenterTitle(); 
  h1_percentage_cov_lambda->GetYaxis()->SetTitleOffset(0.8);
  h1_percentage_cov_lambda->GetYaxis()->SetNdivisions(509);
  
  h1_percentage_cov_lambda->Draw("same axis");
    
  if( saveFIG ) {
    roostr = TString::Format("canv_comb_%d.png", index);  
    canv_comb->SaveAs(roostr);
  }
    
     
}

/////////////////////// ccc

double TCN::GetChi2(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp)
{
  double chi2 = 0;
  
  TMatrixD matrix_delta = matrix_pred_temp - matrix_meas_temp;
  TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T(); 

  int rows = matrix_pred_temp.GetNcols();

  /// docDB 32520, when the prediction is sufficiently low
  double array_pred_protect[11] = {0, 0.461, 0.916, 1.382, 1.833, 2.298, 2.767, 3.225, 3.669, 4.141, 4.599};
  
  TMatrixD matrix_stat_cov(rows, rows);
  if( !FLAG_STATCOV_NO ) {
    for(int idx=0; idx<rows; idx++) {
      matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);// Pearson
    
      double val_meas = matrix_meas_temp(0, idx);
      double val_pred = matrix_pred_temp(0, idx);
      int int_meas = (int)(val_meas+0.1);

      if( int_meas>=1 && int_meas<=10) {
	if( val_pred<array_pred_protect[int_meas] ) {
	  double numerator = pow(val_pred-val_meas, 2);
	  double denominator = 2*( val_pred - val_meas + val_meas*log(val_meas/val_pred) );
	  matrix_stat_cov(idx, idx) = numerator/denominator;
	}
      }    
    }
  }

  TMatrixD matrix_total_cov(rows, rows); matrix_total_cov = matrix_syst_abscov_temp + matrix_stat_cov;
  TMatrixD matrix_total_cov_inv = matrix_total_cov; matrix_total_cov_inv.Invert();
  chi2 = ( matrix_delta*matrix_total_cov_inv*matrix_delta_T )(0,0);

  if( 1 ) {// rate analysis
    double sum_meas = 0;
    double sum_pred = 0;
    double sum_cov = 0;
    for(int idx=0; idx<rows; idx++) {
      sum_meas += matrix_meas_temp(0, idx);
      sum_pred += matrix_pred_temp(0, idx);
      for(int j=0; j<rows; j++) {
	sum_cov += matrix_total_cov( idx, j );
      }
    } 
    double chi2_rate = pow(sum_pred-sum_meas, 2)/sum_cov;   
    cout<<endl<<" ---> check chi2_rate "<<chi2_rate<<endl<<endl;
  }
  
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

    roostr = "file_numu.root";
    roostr = "file_numuPC_phi.root";
    //roostr = "file_numu_vtxZ.root";
    roostr = "file_dQdx.root";

    //roostr = "file_user_test.root";

    //roostr = "file_user_Ehad_no.root";
    //roostr = "file_user_Ehad_wi.root";
    //roostr = "file_user_Ehad70_no.root";
    //roostr = "file_user_Ehad70_wi.root";
    //roostr = "file_user_Ehad75_wi.root";
    //roostr = "file_user_Ehad80_wi.root";
    //roostr = "file_user_Ehad_PC_80.root";
	
    TFile *file_numu = new TFile(roostr, "read");
  
    TMatrixD *matrix_gof_pred = (TMatrixD*)file_numu->Get("matrix_gof_pred");
    TMatrixD *matrix_gof_meas = (TMatrixD*)file_numu->Get("matrix_gof_meas");
    TMatrixD *matrix_gof_syst = (TMatrixD*)file_numu->Get("matrix_gof_syst");
    
    int rows = matrix_gof_syst->GetNrows();
    BINS = rows;


    matrix_pred.Clear(); matrix_pred.ResizeTo(1, rows); matrix_pred = (*matrix_gof_pred);
    matrix_meas.Clear(); matrix_meas.ResizeTo(1, rows); matrix_meas = (*matrix_gof_meas);
    matrix_syst_abscov.Clear(); matrix_syst_abscov.ResizeTo(rows, rows); matrix_syst_abscov = (*matrix_gof_syst);


    // for(int idx=1; idx<=rows; idx++) {
    //   matrix_pred(0, idx-1) = 100;
    //   matrix_meas(0, idx-1) = 100;      
    //   for(int jdx=1; jdx<=rows; jdx++) {
    // 	matrix_syst_abscov(idx-1, jdx-1) = 0;
    //   }
    // }

    
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

  testcn->FLAG_STATCOV_NO = 0;
  
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
      cout<<TString::Format(" ---> itoy %4d, chi2 %7.2f", idx, chi2)<<endl;
    }

    int itoy = 90;// default
    //itoy = 22;
    testcn->SetToy(itoy);    
    TMatrixD matrix_meas_temp = testcn->matrix_fake_meas;    
    testcn->Plotting_singlecase(testcn->matrix_pred, matrix_meas_temp*0.1, testcn->matrix_syst_abscov, 1, 101);
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
