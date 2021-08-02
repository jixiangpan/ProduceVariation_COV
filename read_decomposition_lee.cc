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
    
    FLAG_INPUTFILE_COV_HAS_STAT = 0;
    
    threshold_sigma = 3;    
  }

  /////////////////////////////////////////// data member
  
  TRandom3 *rand;

  TMatrixD matrix_pred;       // inputfile
  TMatrixD matrix_meas;       // inputfile
  TMatrixD matrix_syst_abscov;// inputfile
  
  int BINS;
  int TOYS;

  bool FLAG_INPUTFILE_COV_HAS_STAT;

  TMatrixD matrix_totl_abscov_in_GetChi2;

  double threshold_sigma;
  
  TMatrixD matrix_fake_meas;
  map<int, map<int, double> >map_toy_meas;
  
  /////////////////////////////////////////// memberfunction
  
  void Initialization(TString roofile);

  void SetMeas2Expt() {
    matrix_fake_meas.Clear(); matrix_fake_meas.ResizeTo(1, BINS);
    matrix_fake_meas = matrix_meas;
  }

  double GetChi2(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp);

  void Plotting_decomposition(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp, int index, TString str_postfix);
  
  void ProduceVariation(int ntoys);
  
  void SetToy(int itoy) {
    matrix_fake_meas.Clear(); matrix_fake_meas.ResizeTo(1, BINS);    
    for(int ibin=1; ibin<=BINS; ibin++) matrix_fake_meas(0, ibin-1) = map_toy_meas[itoy][ibin];
  }

};

/////////////////////// ccc

void TCN::ProduceVariation(int ntoys)
{
  TString roostr = "";
  
  TOYS = ntoys;

  map_toy_meas.clear();

  cout<<endl<<" ---> Check, producing toys"<<endl<<endl;
    
  TMatrixD matrix_stat_abscov(BINS, BINS);
  if( !FLAG_INPUTFILE_COV_HAS_STAT ) {
    for(int idx=0; idx<BINS; idx++) matrix_stat_abscov(idx, idx) = matrix_pred(0, idx);// Pearson format
  }
  
  TMatrixD matrix_total_abscov(BINS, BINS);
  matrix_total_abscov = matrix_stat_abscov + matrix_syst_abscov;

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
      if( val_with_syst<0 ) val_with_syst = 0;// ??? remove this requirement or not. We should throw out this pesudo expt.
      map_toy_meas[itoy][idx+1] = val_with_syst;
      // array_obj_total[idx] = val_with_syst;
    }// idx

    // principal_obj_total.AddRow( array_obj_total );
    // delete[] array_obj_total;
      
  }// itoy
       
}

/////////////////////// ccc

void TCN::Plotting_decomposition(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp, int index, TString str_postfix)
{
  cout<<TString::Format(" ---------------------------------- (below) processing index %d %s", index, str_postfix.Data())<<endl<<endl;
  
  TString roostr = "";

  int rows = matrix_pred_temp.GetNcols();

  double chi2_normal = GetChi2(matrix_pred_temp, matrix_meas_temp, matrix_syst_abscov_temp);
  cout<<TString::Format(" ---> check chi2_normal %4.2f", chi2_normal)<<endl;

  ///////////////////////////////////////////////////////////////
  
  TMatrixD matrix_total_cov = matrix_totl_abscov_in_GetChi2;
  
  double matrix_determinant = matrix_total_cov.Determinant();
  if( matrix_determinant==0 ) { cerr<<endl<<" Error: determinant of matrix_total_cov is 0"<<endl<<endl; exit(1); }
    
  TMatrixDSym DSmatrix_cov(rows);
  for(int ibin=0; ibin<rows; ibin++)
    for(int jbin=0; jbin<rows; jbin++)
      DSmatrix_cov(ibin, jbin) = matrix_total_cov(ibin, jbin);
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();
  TMatrixD matrix_eigenvector_T = matrix_eigenvector.T(); matrix_eigenvector.T();
  // double check_eigenvector = 0;
  // for(int idx=0; idx<rows; idx++) check_eigenvector += pow( matrix_eigenvector(0,idx),2 );
  // cout<<endl<<" ---> check_eigenvector "<<check_eigenvector<<endl<<endl;
  
  TMatrixD matrix_lambda_pred = matrix_pred_temp * matrix_eigenvector;
  TMatrixD matrix_lambda_meas = matrix_meas_temp * matrix_eigenvector;
  TMatrixD matrix_lambda_delta = matrix_lambda_meas - matrix_lambda_pred;
  TMatrixD matrix_lambda_delta_T = matrix_lambda_delta.T(); matrix_lambda_delta.T();

  TMatrixD matrix_lambda_cov(rows, rows);
  for(int idx=0; idx<rows; idx++) matrix_lambda_cov(idx, idx) = matrix_eigenvalue(idx);
  TMatrixD matrix_lambda_cov_inv = matrix_lambda_cov; matrix_lambda_cov_inv.Invert();
  double chi2_lambda = (matrix_lambda_delta * matrix_lambda_cov_inv * matrix_lambda_delta_T)(0,0);
  double pvalue_lambda = TMath::Prob( chi2_lambda, rows );
  if( pvalue_lambda<1e-16 ) {
    cout<<" ---> WARNING (pvalue_lambda): pvalue reach the ROOT limit "<<endl;
    pvalue_lambda = 1e-16;
  }
  double significance_lambda = sqrt( TMath::ChisquareQuantile( 1-pvalue_lambda, 1 ) );
  cout<<TString::Format(" ---> check chi2_lambda %4.2f, ndf %d, pvalue %10.8f, %4.2f sigma",
                        chi2_lambda, rows, pvalue_lambda, significance_lambda)<<endl<<endl;

  /////////////////////////////////////////////////////////////// correction of look elsewhere effect

  map<int, double>map_epsilon_above_threshold_sigma;
  
  TMatrixD matrix_lambda_epsilon = matrix_lambda_delta;
  for(int idx=0; idx<rows; idx++) {
    matrix_lambda_epsilon(0, idx) = matrix_lambda_delta(0, idx)/sqrt( matrix_lambda_cov(idx, idx) );
    if( fabs(matrix_lambda_epsilon(0, idx)) >= threshold_sigma ) {
      map_epsilon_above_threshold_sigma[idx] = matrix_lambda_epsilon(0, idx);
    }
  }

  double pvalue_global = 0;
  double significance_global = 0;
  double sum_chi2_AA = 0;
  
  int size_map_epsilon_above_threshold_sigma = map_epsilon_above_threshold_sigma.size();
  cout<<TString::Format(" Decomposition bins above %3.1f sigma: %d", threshold_sigma, size_map_epsilon_above_threshold_sigma)<<endl;
  for( auto it_map=map_epsilon_above_threshold_sigma.begin(); it_map!=map_epsilon_above_threshold_sigma.end(); it_map++ ) {
    int bin_index = it_map->first; double val_epsilon = it_map->second;
    cout<<TString::Format(" ---> %3d, %5.2f", bin_index+1, val_epsilon)<<endl;
  }
  cout<<endl;

  if( size_map_epsilon_above_threshold_sigma>=1 ) {
    if( size_map_epsilon_above_threshold_sigma==1 ) {
      sum_chi2_AA = pow(map_epsilon_above_threshold_sigma.begin()->second, 2);
      double pvalue_local = TMath::Prob( pow(map_epsilon_above_threshold_sigma.begin()->second, 2), 1 );
      if( pvalue_local<1e-16 ) {
	cout<<" ---> WARNING (pvalue_local): pvalue reach the ROOT limit "<<endl;
	pvalue_local = 1e-16;
      }
      double significance_local = fabs( map_epsilon_above_threshold_sigma.begin()->second );
      pvalue_global = 1 - pow( 1-pvalue_local, rows);
      if( pvalue_global<1e-16 ) {
	cout<<" ---> WARNING (pvalue_global): pvalue reach the ROOT limit "<<endl;
	pvalue_global = 1e-16;
      }
      significance_global = sqrt( TMath::ChisquareQuantile( 1-pvalue_global, 1 ) );
      cout<<TString::Format(" LEE correction: local: %4.2f sigma, global: %4.2f sigma",
                            significance_local, significance_global)<<endl<<endl;      
    }// if( size_map_epsilon_above_threshold_sigma==1 )
    else {
      sum_chi2_AA = 0;
      for( auto it_map=map_epsilon_above_threshold_sigma.begin(); it_map!=map_epsilon_above_threshold_sigma.end(); it_map++ ) {
        int bin_index = it_map->first; double val_epsilon = it_map->second;
        sum_chi2_AA += pow( val_epsilon, 2 );
      }// for

      int user_vec_size = size_map_epsilon_above_threshold_sigma;
      double pvalue_local = TMath::Prob( sum_chi2_AA, user_vec_size );
      if( pvalue_local<1e-16 ) {
	cout<<" ---> WARNING (pvalue_local): pvalue reach the ROOT limit "<<endl;
	pvalue_local = 1e-16;
      }
      double significance_local = sqrt( TMath::ChisquareQuantile( 1-pvalue_local, user_vec_size ) );
      double coeff = TMath::Factorial(rows)/TMath::Factorial(rows-user_vec_size)/TMath::Factorial(user_vec_size);
      pvalue_global = coeff * pvalue_local;
      if( pvalue_global<1e-16 ) {
	cout<<" ---> WARNING (pvalue_global): pvalue reach the ROOT limit "<<endl;
	pvalue_global = 1e-16;
      }
      significance_global = sqrt( TMath::ChisquareQuantile( 1-pvalue_global, 1 ) );
      cout<<TString::Format(" LEE correction: local: %4.2f sigma, global: %4.2f sigma",
                            significance_local, significance_global)<<endl<<endl;  
    }// else
  }// if( size_map_epsilon_above_threshold_sigma>=1 )
    
  /////////////////////////////////////////////////////////////// PLOTTING

  int color_meas = kBlue;
  int color_pred = kRed;
  
  roostr = TString::Format("h1_meas_%d_%s", index, str_postfix.Data()); TH1D *h1_meas = new TH1D(roostr, "", rows, 0, rows);
  roostr = TString::Format("h1_pred_%d_%s", index, str_postfix.Data()); TH1D *h1_pred = new TH1D(roostr, "", rows, 0, rows);
  TGraphErrors *gh_meas = new TGraphErrors(); gh_meas->SetName( TString::Format("gh_meas_%d_%s", index, str_postfix.Data()) );
  roostr = TString::Format("h1_meas2pred_%d_%s", index, str_postfix.Data()); TH1D *h1_meas2pred = new TH1D(roostr, "", rows, 0, rows);
  roostr = TString::Format("h1_meas2pred_syst_%d_%s", index, str_postfix.Data()); TH1D *h1_meas2pred_syst = new TH1D(roostr, "", rows, 0, rows);

  roostr = TString::Format("h1_lambda_meas_%d_%s", index, str_postfix.Data()); TH1D *h1_lambda_meas = new TH1D(roostr, "", rows, 0, rows);
  roostr = TString::Format("h1_lambda_pred_%d_%s", index, str_postfix.Data()); TH1D *h1_lambda_pred = new TH1D(roostr, "", rows, 0, rows);
  roostr = TString::Format("h1_lambda_epsilon_%d_%s", index, str_postfix.Data()); TH1D *h1_lambda_epsilon = new TH1D(roostr, "", rows, 0, rows);
  roostr = TString::Format("h1_lambda_eigenvalue_%d_%s", index, str_postfix.Data()); TH1D *h1_lambda_eigenvalue = new TH1D(roostr, "", rows, 0, rows);  
  roostr = TString::Format("h2_lambda_transform_%d_%s", index, str_postfix.Data()); TH2D *h2_lambda_transform = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);

  for(int idx=0; idx<rows; idx++) {
    double meas_val = matrix_meas_temp(0, idx);
    double meas_err = sqrt( meas_val ); if( FLAG_INPUTFILE_COV_HAS_STAT ) meas_err = 0;
    double pred_val = matrix_pred_temp(0, idx);
    double pred_err = sqrt( matrix_syst_abscov_temp(idx, idx) );
    if( pred_val!=0 ) {
      h1_meas2pred->SetBinContent( idx+1, meas_val/pred_val ); h1_meas2pred->SetBinError( idx+1, meas_err/pred_val );
      h1_meas2pred_syst->SetBinContent( idx+1, 1 ); h1_meas2pred_syst->SetBinError( idx+1, pred_err/pred_val );
    }
    h1_meas->SetBinContent(idx+1, meas_val); h1_meas->SetBinError(idx+1, meas_err);
    h1_pred->SetBinContent(idx+1, pred_val); h1_pred->SetBinError(idx+1, pred_err);
    gh_meas->SetPoint(idx, h1_meas->GetBinCenter(idx+1), meas_val );
    gh_meas->SetPointError(idx, 0, meas_err );

    double meas_lambda_val = matrix_lambda_meas(0, idx);
    double pred_lambda_val = matrix_lambda_pred(0, idx);
    double pred_lambda_err = sqrt( matrix_lambda_cov(idx, idx) );
    double epsilon_val = (meas_lambda_val-pred_lambda_val)/pred_lambda_err;
    double lambda_eigenvalue = matrix_lambda_cov(idx, idx);
    
    h1_lambda_meas->SetBinContent(idx+1, meas_lambda_val);
    h1_lambda_pred->SetBinContent(idx+1, pred_lambda_val);
    h1_lambda_pred->SetBinError(idx+1, pred_lambda_err);
    h1_lambda_epsilon->SetBinContent(idx+1, epsilon_val);
    h1_lambda_eigenvalue->SetBinContent(idx+1, lambda_eigenvalue);

    for(int jdx=0; jdx<rows; jdx++) {
      h2_lambda_transform->SetBinContent(idx+1, jdx+1, matrix_eigenvector(idx,jdx));
    }// for(int jdx=0; jdx<rows; jdx++)
    
  }// for(int idx=0; idx<rows; idx++)

  //////////////////////////////////////////////////////////////////////////

  TF1 *ff_0 = new TF1(TString::Format("ff_0_%d_%s", index, str_postfix.Data()), "0", 0, 1e3);
  ff_0->SetLineColor(kGray+1); ff_0->SetLineStyle(7);
  
  TF1 *ff_1 = new TF1(TString::Format("ff_1_%d_%s", index, str_postfix.Data()), "1", 0, 1e3);
  ff_1->SetLineColor(kGray+1); ff_1->SetLineStyle(7);
  
  //////////////////////////////////////////////////////////////////////////

  roostr = TString::Format("canv_deco_normal_meas_pred_%d_%s", index, str_postfix.Data());  
  TCanvas *canv_deco_normal_meas_pred = new TCanvas(roostr, roostr, 1000, 950);

  /////////////
  canv_deco_normal_meas_pred->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_no->Draw(); pad_top_no->cd();
  
  TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone(TString::Format("h1_pred_clone_%d_%s", index, str_postfix.Data()));
  h1_pred_clone->Draw("e2"); h1_pred_clone->SetMinimum(0);
  double max_h1_pred_clone = 0; double max_h1_meas = h1_meas->GetMaximum();
  for(int ibin=1; ibin<=rows; ibin++) {
    double sub_pred = h1_pred_clone->GetBinContent(ibin) + h1_pred_clone->GetBinError(ibin);
    if( max_h1_pred_clone<sub_pred ) max_h1_pred_clone = sub_pred;
  }
  if( max_h1_pred_clone<max_h1_meas ) h1_pred_clone->SetMaximum( max_h1_meas*1.1 );
  h1_pred_clone->SetFillColor(kRed-10); h1_pred_clone->SetFillStyle(1001);
  h1_pred_clone->SetMarkerSize(0); h1_pred_clone->SetLineColor(color_pred);
  h1_pred_clone->GetXaxis()->SetLabelColor(10);
  func_xy_title(h1_pred_clone, "Measurement bin index", "Entries"); func_title_size(h1_pred_clone, 0.065, 0.065, 0.065, 0.065);
  h1_pred_clone->GetXaxis()->CenterTitle(); h1_pred_clone->GetYaxis()->CenterTitle();
  h1_pred_clone->GetYaxis()->SetTitleOffset(1.2); 

  h1_pred->Draw("same hist"); h1_pred->SetLineColor(color_pred);
  //gh_meas->Draw("same p"); gh_meas->SetMarkerStyle(20); gh_meas->SetMarkerSize(1.2);
  //gh_meas->SetMarkerColor(color_meas); gh_meas->SetLineColor(color_meas);
  h1_meas->SetLineColor(color_meas);
  if( FLAG_INPUTFILE_COV_HAS_STAT ) h1_meas->Draw("same p");
  else  h1_meas->Draw("same e1");
  h1_meas->SetMarkerStyle(20); h1_meas->SetMarkerColor(color_meas); h1_meas->SetMarkerSize(1.2);
  h1_pred_clone->Draw("same axis");
  
  double shift_normal_meas_pred_YY = 0; double shift_normal_meas_pred_XX = 0.3;
  TLegend *lg_normal_meas_pred = new TLegend(0.17+0.03+shift_normal_meas_pred_XX, 0.6+shift_normal_meas_pred_YY,
					     0.4+0.04+shift_normal_meas_pred_XX, 0.85+shift_normal_meas_pred_YY);    
  lg_normal_meas_pred->AddEntry(h1_meas, TString::Format("#color[%d]{Data}", color_meas), "lep");
  lg_normal_meas_pred->AddEntry(h1_pred_clone, TString::Format("#color[%d]{Prediction}", color_pred), "fl");
  lg_normal_meas_pred->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %4.2f/%d}", color_pred, chi2_normal, rows), "");
  lg_normal_meas_pred->Draw();
  lg_normal_meas_pred->SetBorderSize(0); lg_normal_meas_pred->SetFillStyle(0); lg_normal_meas_pred->SetTextSize(0.08);

  /////////////
  canv_deco_normal_meas_pred->cd();
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
  
  h1_meas2pred->SetLineColor(color_meas);
  if( FLAG_INPUTFILE_COV_HAS_STAT ) h1_meas2pred->Draw("same p");
  else h1_meas2pred->Draw("same e1");  
  h1_meas2pred->SetMarkerStyle(20); h1_meas2pred->SetMarkerColor(color_meas); h1_meas2pred->SetMarkerSize(1.3);
  
  h1_meas2pred_syst->Draw("same axis");

  roostr = TString::Format("canv_deco_normal_meas_pred_%d_%s.png", index, str_postfix.Data());  
  canv_deco_normal_meas_pred->SaveAs(roostr);

  //////////////////////////////////////////////////////////////////////////
  
  roostr = TString::Format("canv_deco_h2_lambda_transform_%d_%s", index, str_postfix.Data());
  TCanvas *canv_deco_h2_lambda_transform = new TCanvas(roostr, roostr, 900, 850);
  func_canv_margin(canv_deco_h2_lambda_transform, 0.15, 0.2,0.15,0.2);
  h2_lambda_transform->Draw("colz");
  func_xy_title(h2_lambda_transform, "Measurement bin index", "Decomposition bin index");
  func_title_size(h2_lambda_transform, 0.05, 0.05, 0.05, 0.05);
  h2_lambda_transform->GetZaxis()->SetLabelSize(0.05);
  h2_lambda_transform->GetXaxis()->CenterTitle(); h2_lambda_transform->GetYaxis()->CenterTitle();
  roostr = TString::Format("canv_deco_h2_lambda_transform_%d_%s.png", index, str_postfix.Data());
  canv_deco_h2_lambda_transform->SaveAs(roostr);

  //////////////////////////////////////////////////////////////////////////
  
  roostr = TString::Format("canv_deco_lambda_meas_pred_%d_%s", index, str_postfix.Data());
  TCanvas *canv_deco_lambda_meas_pred = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_deco_lambda_meas_pred, 0.15, 0.1, 0.1, 0.15);

  roostr = TString::Format("h1_lambda_pred_clone_%d_%s", index, str_postfix.Data());
  TH1D *h1_lambda_pred_clone = (TH1D*)h1_lambda_pred->Clone(roostr);
  h1_lambda_pred_clone->Draw("e2");
  h1_lambda_pred_clone->SetFillColor(kRed-10); h1_lambda_pred_clone->SetFillStyle(1001);
  h1_lambda_pred_clone->SetMarkerSize(0);
  func_xy_title(h1_lambda_pred_clone, "Decomposition bin index", "Entries\'");
  func_title_size(h1_lambda_pred_clone, 0.05, 0.05, 0.05, 0.05);
  h1_lambda_pred_clone->GetXaxis()->CenterTitle(); h1_lambda_pred_clone->GetYaxis()->CenterTitle();
  h1_lambda_pred_clone->GetYaxis()->SetTitleOffset(1.5);

  h1_lambda_pred->Draw("hist same"); h1_lambda_pred->SetLineColor(color_pred);
  double ymin_h1_lambda_pred = 1e8;
  double ymax_h1_lambda_pred = -1e8;
  for(int ibin=1; ibin<=rows; ibin++) {    
    double cv = h1_lambda_pred->GetBinContent( ibin );
    double error = h1_lambda_pred->GetBinError( ibin );
    if( ymin_h1_lambda_pred>cv-error ) ymin_h1_lambda_pred = cv-error;
    if( ymax_h1_lambda_pred<cv+error ) ymax_h1_lambda_pred = cv+error;
  }
  if( ymin_h1_lambda_pred>h1_lambda_meas->GetMinimum() && h1_lambda_meas->GetMinimum()<0 )
    h1_lambda_pred_clone->SetMinimum(h1_lambda_meas->GetMinimum()*1.1);
  if( ymin_h1_lambda_pred>h1_lambda_meas->GetMinimum() && h1_lambda_meas->GetMinimum()>=0 )
    h1_lambda_pred_clone->SetMinimum(h1_lambda_meas->GetMinimum()*0.9);
  if( h1_lambda_meas->GetMaximum()>=0 && ymax_h1_lambda_pred<h1_lambda_meas->GetMaximum() )
    h1_lambda_pred_clone->SetMaximum(h1_lambda_meas->GetMaximum()*1.1);
  if( h1_lambda_meas->GetMaximum()<0 && ymax_h1_lambda_pred<h1_lambda_meas->GetMaximum() )
    h1_lambda_pred_clone->SetMaximum(h1_lambda_meas->GetMaximum()*0.9);
  
  h1_lambda_meas->Draw("same p");
  h1_lambda_meas->SetMarkerStyle(20); h1_lambda_meas->SetMarkerSize(1.2); h1_lambda_meas->SetMarkerColor(color_meas);
  h1_lambda_meas->SetLineColor(color_meas);    

  h1_lambda_pred_clone->Draw("same axis");   

  double shift_lambda_meas_pred_YY = 0; double shift_lambda_meas_pred_XX = 0.4;
  TLegend *lg_lambda_meas_pred = new TLegend(0.17+0.03+shift_lambda_meas_pred_XX, 0.6+0.08+shift_lambda_meas_pred_YY,
					     0.4+0.04+shift_lambda_meas_pred_XX, 0.85+shift_lambda_meas_pred_YY);    
  lg_lambda_meas_pred->AddEntry(h1_meas, TString::Format("#color[%d]{Data\'}", color_meas), "p");
  lg_lambda_meas_pred->AddEntry(h1_pred_clone, TString::Format("#color[%d]{Prediction\'}", color_pred), "fl");
  lg_lambda_meas_pred->Draw();
  lg_lambda_meas_pred->SetBorderSize(0); lg_lambda_meas_pred->SetFillStyle(0); lg_lambda_meas_pred->SetTextSize(0.05);

  roostr = TString::Format("canv_deco_lambda_meas_pred_%d_%s.png", index, str_postfix.Data());
  canv_deco_lambda_meas_pred->SaveAs(roostr);
 
  //////////////////////////////////////////////////////////////////////////

  roostr = TString::Format("canv_deco_lambda_sigma_%d_%s", index, str_postfix.Data());
  TCanvas *canv_deco_lambda_sigma = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_deco_lambda_sigma, 0.15, 0.1, 0.1, 0.15);
    
  TH1D *h1_lambda_sigmax3 = new TH1D(TString::Format("h1_lambda_sigmax3_%d_%s", index, str_postfix.Data()), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax3->SetBinError(ibin, 3);  
  TH1D *h1_lambda_sigmax2 = new TH1D(TString::Format("h1_lambda_sigmax2_%d_%s", index, str_postfix.Data()), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax2->SetBinError(ibin, 2);  
  TH1D *h1_lambda_sigmax1 = new TH1D(TString::Format("h1_lambda_sigmax1_%d_%s", index, str_postfix.Data()), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_lambda_sigmax1->SetBinError(ibin, 1);
  
  h1_lambda_sigmax3->Draw("e2");
  h1_lambda_sigmax3->SetMinimum(-6); h1_lambda_sigmax3->SetMaximum(6);
  h1_lambda_sigmax3->SetFillColor(kRed-9); h1_lambda_sigmax3->SetFillStyle(1001);
  h1_lambda_sigmax3->SetMarkerSize(0);
  func_xy_title(h1_lambda_sigmax3, "Decomposition bin index", "#epsilon_{i}\' value");
  func_title_size(h1_lambda_sigmax3, 0.05, 0.05, 0.05, 0.05);
  h1_lambda_sigmax3->GetXaxis()->CenterTitle(); h1_lambda_sigmax3->GetYaxis()->CenterTitle();
  h1_lambda_sigmax3->GetYaxis()->SetTitleOffset(1.2); 
  h1_lambda_sigmax3->GetYaxis()->SetNdivisions(112);

  h1_lambda_sigmax2->Draw("same e2"); h1_lambda_sigmax2->SetMarkerSize(0);
  h1_lambda_sigmax2->SetFillColor(kYellow); h1_lambda_sigmax2->SetFillStyle(1001);

  h1_lambda_sigmax1->Draw("same e2"); h1_lambda_sigmax1->SetMarkerSize(0);
  h1_lambda_sigmax1->SetFillColor(kGreen); h1_lambda_sigmax1->SetFillStyle(1001);

  ff_0->Draw("same");
  
  h1_lambda_epsilon->Draw("same p");
  h1_lambda_epsilon->SetMarkerSize(1.4); h1_lambda_epsilon->SetMarkerStyle(21); h1_lambda_epsilon->SetMarkerColor(kBlack);

  h1_lambda_sigmax3->Draw("same axis");

  double shift_x_lambda_sigma = -0.15;
  TLegend *lg_lambda_sigma = new TLegend(0.35+shift_x_lambda_sigma, 0.75, 0.7+shift_x_lambda_sigma, 0.85+0.02);
  lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{%3.1f#sigma:        overall #chi^{2}/dof: %3.1f/%d}",
						kBlue, significance_lambda, chi2_lambda, rows), "");
  if( size_map_epsilon_above_threshold_sigma>=1 ) {
    lg_lambda_sigma->AddEntry("", TString::Format("#color[%d]{%3.1f#sigma (LEE corr.): #chi^{2}/dof: %3.1f/%d (|#epsilon_{i}\'|>3)}",
						  kRed, significance_global, sum_chi2_AA, size_map_epsilon_above_threshold_sigma), "");
  }
  
  lg_lambda_sigma->Draw();
  lg_lambda_sigma->SetBorderSize(0); lg_lambda_sigma->SetFillStyle(0); lg_lambda_sigma->SetTextSize(0.05);

  roostr = TString::Format("canv_deco_lambda_sigma_%d_%s.png", index, str_postfix.Data());
  canv_deco_lambda_sigma->SaveAs(roostr);
  
  //////////////////////////////////////////////////////////////////////////

  roostr = TString::Format("h1_lambda_eigenvalue_normalization_%d_%s", index, str_postfix.Data());
  TH1D *h1_lambda_eigenvalue_normalization = (TH1D*)h1_lambda_eigenvalue->Clone(roostr);
  double integral_h1_lambda_eigenvalue_normalization = h1_lambda_eigenvalue_normalization->Integral();
  h1_lambda_eigenvalue_normalization->Scale( 1./integral_h1_lambda_eigenvalue_normalization );

  /////////////
  roostr = TString::Format("canv_deco_lambda_lambdaA_%d_%s", index, str_postfix.Data());  
  TCanvas *canv_deco_lambda_lambdaA = new TCanvas(roostr, roostr, 1000, 950);

  /////////////
  canv_deco_lambda_lambdaA->cd();
  TPad *pad_top_lambda_lambdaA_no = new TPad("pad_top_lambda_lambdaA_no", "pad_top_lambda_lambdaA_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_lambda_lambdaA_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_lambda_lambdaA_no->Draw(); pad_top_lambda_lambdaA_no->cd();
    
  h1_lambda_sigmax3->Draw("e2");
  h1_lambda_sigmax2->Draw("same e2");
  h1_lambda_sigmax1->Draw("same e2");
  ff_0->Draw("same");  
  h1_lambda_epsilon->Draw("same p");     
  h1_lambda_sigmax3->Draw("same axis");

  func_title_size(h1_lambda_sigmax3, 0.065, 0.065, 0.065, 0.065);
  h1_lambda_sigmax3->GetYaxis()->SetTitleOffset(1.05); h1_lambda_sigmax3->GetXaxis()->SetLabelColor(10);
  h1_lambda_sigmax3->GetXaxis()->SetTickLength(0.042); h1_lambda_sigmax3->GetYaxis()->SetTickLength(0.025);
  
  lg_lambda_sigma->Draw();
  lg_lambda_sigma->SetTextSize(0.065);
  lg_lambda_sigma->SetY1(0.7);
    
  /////////////
  canv_deco_lambda_lambdaA->cd();
  TPad *pad_bot_lambda_lambdaA_no = new TPad("pad_bot_lambda_lambdaA_no", "pad_bot_lambda_lambdaA_no", 0, 0, 1, 0.45);
  pad_bot_lambda_lambdaA_no->SetLogy();
  func_canv_margin(pad_bot_lambda_lambdaA_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_lambda_lambdaA_no->Draw(); pad_bot_lambda_lambdaA_no->cd();

  h1_lambda_eigenvalue_normalization->Draw("hist");
  h1_lambda_eigenvalue_normalization->SetLineWidth(3);
  h1_lambda_eigenvalue_normalization->SetLineColor(kBlack);
  h1_lambda_eigenvalue_normalization->SetMaximum(2);
    
  func_title_size(h1_lambda_eigenvalue_normalization, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_lambda_eigenvalue_normalization, "Decomposition bin index", "#Lambda_{i} / #sum#Lambda_{i}");
  h1_lambda_eigenvalue_normalization->GetXaxis()->SetTickLength(0.05);  h1_lambda_eigenvalue_normalization->GetXaxis()->SetLabelOffset(0.005);
  h1_lambda_eigenvalue_normalization->GetXaxis()->CenterTitle(); h1_lambda_eigenvalue_normalization->GetYaxis()->CenterTitle(); 
  h1_lambda_eigenvalue_normalization->GetYaxis()->SetTitleOffset(0.8);
  h1_lambda_eigenvalue_normalization->GetYaxis()->SetNdivisions(509);
  
  roostr = TString::Format("canv_deco_lambda_lambdaA_%d_%s.png", index, str_postfix.Data());  
  canv_deco_lambda_lambdaA->SaveAs(roostr);

  
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
  if( !FLAG_INPUTFILE_COV_HAS_STAT ) {
    for(int idx=0; idx<rows; idx++) {
      matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);// Pearson format
    
      double val_meas = matrix_meas_temp(0, idx);
      double val_pred = matrix_pred_temp(0, idx);
      int int_meas = (int)(val_meas+0.1);
      
      if( int_meas>=1 && int_meas<=10) {
        if( val_pred<array_pred_protect[int_meas] ) {// protection
          double numerator = pow(val_pred-val_meas, 2);
          double denominator = 2*( val_pred - val_meas + val_meas*log(val_meas/val_pred) );
          matrix_stat_cov(idx, idx) = numerator/denominator;
        }
      }    
    }// for(int idx=0; idx<matrix_pred_temp.GetNcols(); idx++)
  }// FLAG_INPUTFILE_COV_HAS_STAT
  
  TMatrixD matrix_total_cov(rows, rows); matrix_total_cov = matrix_syst_abscov_temp + matrix_stat_cov;
  TMatrixD matrix_total_cov_inv = matrix_total_cov; matrix_total_cov_inv.Invert();
  matrix_totl_abscov_in_GetChi2.ResizeTo(rows, rows); matrix_totl_abscov_in_GetChi2 = matrix_total_cov;
  
  chi2 = ( matrix_delta*matrix_total_cov_inv*matrix_delta_T )(0,0);

  return chi2;
}

/////////////////////// ccc

void TCN::Initialization(TString roofile)
{
  TString roostr = "";
  
  TFile *file_obj = new TFile(roofile, "read");  
  TMatrixD *matrix_gof_pred = (TMatrixD*)file_obj->Get("matrix_gof_pred");
  TMatrixD *matrix_gof_meas = (TMatrixD*)file_obj->Get("matrix_gof_meas");
  TMatrixD *matrix_gof_syst = (TMatrixD*)file_obj->Get("matrix_gof_syst");

  int rows = matrix_gof_syst->GetNrows();
  BINS = rows;

  matrix_pred.Clear(); matrix_pred.ResizeTo(1, rows); matrix_pred = (*matrix_gof_pred);
  matrix_meas.Clear(); matrix_meas.ResizeTo(1, rows); matrix_meas = (*matrix_gof_meas);
  matrix_syst_abscov.Clear(); matrix_syst_abscov.ResizeTo(rows, rows); matrix_syst_abscov = (*matrix_gof_syst);
  delete matrix_gof_pred; delete matrix_gof_meas; delete matrix_gof_syst; delete file_obj;
}
  


////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// MAIN ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void read_decomposition_lee()
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
  //gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";
  
  TCN *testcn = new TCN();
  //testcn->FLAG_INPUTFILE_COV_HAS_STAT = 1;
  
  roostr = "file_numu.root";
  //roostr = "file_numuPC_phi.root";
  //roostr = "file_numu_vtxZ.root";
  //roostr = "file_dQdx.root";
  //roostr = "file_user_Ehad_no.root";
  //roostr = "file_user_Ehad_wi.root";
  //roostr = "file_user_Ehad70_no.root";
  // roostr = "file_user_Ehad70_wi.root";
  //roostr = "file_user_Ehad75_wi.root";
  //roostr = "file_user_Ehad80_wi.root";
  //roostr = "file_user_Ehad_PC_80.root";
  
  testcn->Initialization( roostr );
  
  if( 1 ) {
    testcn->SetMeas2Expt();
    testcn->threshold_sigma = 3;
    testcn->Plotting_decomposition( testcn->matrix_pred, testcn->matrix_fake_meas, testcn->matrix_syst_abscov, 1, "aa" );
  }

  if( 0 ) {
    int ntoys = 100;
    testcn->ProduceVariation( ntoys );
    
    int itoy = 90;// default
    testcn->SetToy(itoy);
    testcn->threshold_sigma = 3;
    testcn->Plotting_decomposition( testcn->matrix_pred, testcn->matrix_fake_meas*0.1, testcn->matrix_syst_abscov, 2, "aa" );
  }
  
}
