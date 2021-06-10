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
};

/////////////////////// ccc

double TCN::GetChi2(TMatrixD matrix_pred_temp, TMatrixD matrix_meas_temp, TMatrixD matrix_syst_abscov_temp)
{

  double chi2 = 0;
  
  TMatrixD matrix_delta = matrix_pred_temp - matrix_meas_temp;
  TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T(); 

  int rows = matrix_pred_temp.GetNcols();
  TMatrixD matrix_stat_cov(rows, rows);
  for(int idx=0; idx<rows; idx++) {
    matrix_stat_cov(idx, idx) = matrix_pred_temp(0, idx);
  }

  TMatrixD matrix_total_cov(rows, rows);
  matrix_total_cov = matrix_syst_abscov_temp + matrix_stat_cov;

  TMatrixD matrix_total_cov_inv = matrix_total_cov;
  matrix_total_cov_inv.Invert();

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
    
  TMatrixD matrix_stat_abscov(BINS, BINS);
  for(int idx=0; idx<BINS; idx++) matrix_stat_abscov(idx, idx) = matrix_pred(0, idx);

  TMatrixD matrix_total_abscov(BINS, BINS);
  matrix_total_abscov = matrix_stat_abscov + matrix_syst_abscov;

  // for(int idx=0; idx<BINS; idx++)
  //   cout<<" ---> check "<<idx+1<<"\t"<<matrix_pred(0, idx)<<"\t"<<sqrt( matrix_total_abscov(idx, idx) )<<endl;
    
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
  gStyle->SetEndErrorSize(2);

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
  }

  int ntoys = 10000;
  testcn->ProduceVariation( ntoys, flag_norm );
  for(int idx=1; idx<=100; idx++) {
    testcn->SetToy(idx);
    double chi2 = testcn->GetChi2( testcn->matrix_pred, testcn->matrix_fake_meas, testcn->matrix_syst_abscov );
    //cout<<TString::Format(" ---> itoy %4d, chi2 %7.2f", idx, chi2)<<endl;
  }

  if( 1 ) {
    int itoy = 90;
   
    TMatrixD matrix_toy_syst_abscov = testcn->matrix_syst_abscov ;
    int rows = matrix_toy_syst_abscov.GetNrows();
    
    TMatrixD matrix_toy_pred = testcn->matrix_pred;
    TMatrixD matrix_toy_meas(1, rows);
    
    TH1D *h1_fake_meas = new TH1D("h1_fake_meas", "", rows, 0, rows);
    TGraphErrors *gh_fake_meas = new TGraphErrors();
    TH1D *h1_pred = new TH1D("h1_pred", "", rows, 0, rows);
    
    for(int ibin=1; ibin<=rows; ibin++) {

      double scaleF_meas = 1;
      
      double val_fake_meas = testcn->map_toy_meas[itoy][ibin] * scaleF_meas;
      
      ///////
      matrix_toy_meas(0, ibin-1) = val_fake_meas;

      ///////
      h1_fake_meas->SetBinContent(ibin, val_fake_meas);
      gh_fake_meas->SetPoint( ibin-1, h1_fake_meas->GetBinCenter(ibin), val_fake_meas );
      gh_fake_meas->SetPointError( ibin-1, 0, sqrt(val_fake_meas) );      
      
      h1_pred->SetBinContent( ibin, testcn->matrix_pred(0, ibin-1) );
      h1_pred->SetBinError( ibin, sqrt(testcn->matrix_syst_abscov(ibin-1, ibin-1)) );
    }
    
    double chi2 = testcn->GetChi2( matrix_toy_pred, matrix_toy_meas, matrix_toy_syst_abscov );
    cout<<endl<<TString::Format(" ---> check itoy %4d, chi2 %7.2f", itoy, chi2)<<endl<<endl;
    
    roostr = "canv_h1_fake_meas";
    TCanvas *canv_h1_fake_meas = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_fake_meas, 0.15, 0.1, 0.1, 0.15);
    
    TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone("h1_pred_clone");
    h1_pred_clone->Draw("e2");
    h1_pred_clone->SetMinimum(0);
    if( h1_pred_clone->GetMaximum() < h1_fake_meas->GetMaximum() ) {
      h1_pred_clone->SetMaximum( h1_fake_meas->GetMaximum() * 1.1 );
    }
    h1_pred_clone->SetFillColor(kRed); h1_pred_clone->SetFillStyle(3005);
    h1_pred_clone->SetMarkerSize(0);
    h1_pred_clone->SetLineColor(kRed);
    func_xy_title(h1_pred_clone, "Bin index", "Entries");
    func_title_size(h1_pred_clone, 0.05, 0.05, 0.05, 0.05);
    h1_pred_clone->GetXaxis()->CenterTitle(); h1_pred_clone->GetYaxis()->CenterTitle();
    h1_pred_clone->GetYaxis()->SetTitleOffset(1.5); 

    h1_pred->Draw("hist same");
    h1_pred->SetLineColor(kRed);

    gh_fake_meas->Draw("same p");
    gh_fake_meas->SetMarkerStyle(20);
    gh_fake_meas->SetMarkerSize(1.1);
    gh_fake_meas->SetMarkerColor(kBlack);
    gh_fake_meas->SetLineColor(kBlack);

    h1_pred_clone->Draw("same axis");
    
    TLegend *lg_chi2_toy = new TLegend(0.17, 0.65, 0.4, 0.85);    
    lg_chi2_toy->AddEntry(gh_fake_meas, TString::Format("#color[%d]{Fake data}", kBlack), "p");
    lg_chi2_toy->AddEntry(h1_pred_clone, TString::Format("#color[%d]{Prediction}", kRed), "fl");
    lg_chi2_toy->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %4.1f/%d}", kRed, chi2, rows), "");
    lg_chi2_toy->Draw();
    lg_chi2_toy->SetBorderSize(0); lg_chi2_toy->SetFillStyle(0); lg_chi2_toy->SetTextSize(0.065);

    canv_h1_fake_meas->SaveAs("canv_h1_fake_meas.png");

    ///////////////////////////////////////////////////////////

    TMatrixD matrix_stat_abscov(rows, rows);
    for(int idx=0; idx<rows; idx++) matrix_stat_abscov(idx, idx) = matrix_toy_pred(0, idx);

    TMatrixD matrix_total_abscov(rows, rows);
    matrix_total_abscov = matrix_stat_abscov + matrix_toy_syst_abscov;

    TMatrixDSym DSmatrix_cov(rows);
    for(int ibin=0; ibin<rows; ibin++) {
      for(int jbin=0; jbin<rows; jbin++) {
	DSmatrix_cov(ibin, jbin) = matrix_total_abscov(ibin, jbin);
      }
    }
    TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
    TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
    TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

    TMatrixD matrix_eigenvector_T = matrix_eigenvector.T(); matrix_eigenvector.T();
    // TMatrixD matrix_eigenvector_IT = matrix_eigenvector * matrix_eigenvector_T;// --> unit matrix
    // TCanvas *canv_check_matrix_eigenvalue = new TCanvas("canv_check_matrix_eigenvalue", "canv_check_matrix_eigenvalue", 900, 650);
    // matrix_eigenvector_IT.Draw("colz");
    // TMatrixD matrix_lambda = matrix_eigenvector_T * matrix_total_abscov * matrix_eigenvector;
    // for(int idx=0; idx<rows; idx++) cout<<idx<<"\t"<<matrix_lambda(idx, idx)<<"\t"<<matrix_eigenvalue(idx)<<endl;

    TMatrixD matrix_delta = matrix_toy_pred - matrix_toy_meas;
    TMatrixD matrix_delta_lambda = matrix_delta * matrix_eigenvector;
    TMatrixD matrix_delta_lambda_T = matrix_delta_lambda.T(); matrix_delta_lambda.T();

    TMatrixD matrix_cov_lambda(rows, rows);
    for(int idx=0; idx<rows; idx++) matrix_cov_lambda(idx, idx) = matrix_eigenvalue(idx);
    TMatrixD matrix_cov_lambda_inv = matrix_cov_lambda;
    matrix_cov_lambda_inv.Invert();
    
    double chi2_lambda = (matrix_delta_lambda * matrix_cov_lambda_inv * matrix_delta_lambda_T)(0,0);
    cout<<endl<<TString::Format(" ---> check itoy %4d, chi2_lambda %7.2f", itoy, chi2_lambda)<<endl<<endl;

    /////////////////////////////////

    TMatrixD matrix_lambda_pred = matrix_toy_pred * matrix_eigenvector;
    TMatrixD matrix_lambda_meas = matrix_toy_meas * matrix_eigenvector;

    roostr = "h1_lambda_pred";
    TH1D *h1_lambda_pred = new TH1D(roostr, "", rows, 0, rows);
    TGraph *gh_lambda_meas = new TGraph();

    ////////
    
    int color_meas_100 = kBlack;
    int color_meas_50 = kGreen + 1;
    int color_meas_10 = kBlue;
    
    roostr = "h1_lambda_meas_100"; TH1D *h1_lambda_meas_100 = new TH1D(roostr, roostr, rows, 0, rows);
    TMatrixD matrix_lambda_meas_100(1, rows);
    for(int idx=0; idx<rows; idx++) matrix_lambda_meas_100(0, idx) = testcn->map_toy_meas[itoy][idx+1] * 1;
    matrix_lambda_meas_100 = matrix_lambda_meas_100 * matrix_eigenvector;
    for(int idx=0; idx<rows; idx++) h1_lambda_meas_100->SetBinContent( idx+1, matrix_lambda_meas_100(0, idx) );
    h1_lambda_meas_100->SetLineColor(color_meas_100);
      
    roostr = "h1_lambda_meas_50"; TH1D *h1_lambda_meas_50 = new TH1D(roostr, roostr, rows, 0, rows);
    TMatrixD matrix_lambda_meas_50(1, rows);
    for(int idx=0; idx<rows; idx++) matrix_lambda_meas_50(0, idx) = testcn->map_toy_meas[itoy][idx+1] * 0.5;
    matrix_lambda_meas_50 = matrix_lambda_meas_50 * matrix_eigenvector;
    for(int idx=0; idx<rows; idx++) h1_lambda_meas_50->SetBinContent( idx+1, matrix_lambda_meas_50(0, idx) );
    h1_lambda_meas_50->SetLineColor(color_meas_50);

    roostr = "h1_lambda_meas_10"; TH1D *h1_lambda_meas_10 = new TH1D(roostr, roostr, rows, 0, rows);
    TMatrixD matrix_lambda_meas_10(1, rows);
    for(int idx=0; idx<rows; idx++) matrix_lambda_meas_10(0, idx) = testcn->map_toy_meas[itoy][idx+1] * 0.1;
    matrix_lambda_meas_10 = matrix_lambda_meas_10 * matrix_eigenvector;
    for(int idx=0; idx<rows; idx++) h1_lambda_meas_10->SetBinContent( idx+1, matrix_lambda_meas_10(0, idx) );
    h1_lambda_meas_10->SetLineColor(color_meas_10);

    roostr = "h1_lambda_meas_diff_100"; TH1D *h1_lambda_meas_diff_100 = new TH1D(roostr, "", rows, 0, rows);
    h1_lambda_meas_diff_100->SetLineColor(color_meas_100);    
    roostr = "h1_lambda_meas_diff_50"; TH1D *h1_lambda_meas_diff_50 = new TH1D(roostr, "", rows, 0, rows);
    h1_lambda_meas_diff_50->SetLineColor(color_meas_50);    
    roostr = "h1_lambda_meas_diff_10"; TH1D *h1_lambda_meas_diff_10 = new TH1D(roostr, "", rows, 0, rows);
    h1_lambda_meas_diff_10->SetLineColor(color_meas_10);
    for(int ibin=1; ibin<=rows; ibin++) {
      double val_pred = matrix_lambda_pred(0, ibin-1);
      double val_err = sqrt( matrix_cov_lambda(ibin-1, ibin-1) );
      
      double val_meas_100 = h1_lambda_meas_100->GetBinContent(ibin);
      double val_meas_50 = h1_lambda_meas_50->GetBinContent(ibin);
      double val_meas_10 = h1_lambda_meas_10->GetBinContent(ibin);

      h1_lambda_meas_diff_100->SetBinContent(ibin, (val_meas_100-val_pred)/val_err);
      h1_lambda_meas_diff_50->SetBinContent(ibin, (val_meas_50-val_pred)/val_err);
      h1_lambda_meas_diff_10->SetBinContent(ibin, (val_meas_10-val_pred)/val_err);      
    }
    
    ////////
    
    for(int idx=0; idx<rows; idx++) {
      double val_pred = matrix_lambda_pred(0, idx);
      double val_meas = matrix_lambda_meas(0, idx);
      double val_sigma = sqrt( matrix_cov_lambda(idx, idx) );

      h1_lambda_pred->SetBinContent(idx+1, val_pred);
      h1_lambda_pred->SetBinError(idx+1, val_sigma);

      gh_lambda_meas->SetPoint( idx, h1_lambda_pred->GetBinCenter(idx+1), val_meas );
    }
    
    roostr = "canv_h1_lambda_pred";
    TCanvas *canv_h1_lambda_pred = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_lambda_pred, 0.15, 0.1, 0.1, 0.15);
    
    TH1D *h1_lambda_pred_clone = (TH1D*)h1_lambda_pred->Clone("h1_lambda_pred_clone");
    h1_lambda_pred_clone->Draw("e2");
    //h1_lambda_pred_clone->SetMinimum(0);
    h1_lambda_pred_clone->SetFillColor(kRed-10); h1_lambda_pred_clone->SetFillStyle(1001);
    h1_lambda_pred_clone->SetMarkerSize(0);
    h1_lambda_pred_clone->SetLineColor(kRed);
    func_xy_title(h1_lambda_pred_clone, "Bin index", "");
    func_title_size(h1_lambda_pred_clone, 0.05, 0.05, 0.05, 0.05);
    h1_lambda_pred_clone->GetXaxis()->CenterTitle(); h1_lambda_pred_clone->GetYaxis()->CenterTitle();
    h1_lambda_pred_clone->GetYaxis()->SetTitleOffset(1.5); 

    h1_lambda_pred->Draw("hist same");
    h1_lambda_pred->SetLineColor(kRed);

    gh_lambda_meas->Draw("same p");
    gh_lambda_meas->SetMarkerStyle(20);
    gh_lambda_meas->SetMarkerSize(1.1);
    gh_lambda_meas->SetMarkerColor(kBlack);
    gh_lambda_meas->SetLineColor(kBlack);    

    h1_lambda_pred_clone->Draw("same axis");      


    
    roostr = "canv_h1_lambda_pred2";
    TCanvas *canv_h1_lambda_pred2 = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_lambda_pred2, 0.15, 0.1, 0.1, 0.15);    
    h1_lambda_pred_clone->Draw("e2");
    h1_lambda_pred->Draw("same hist");      
    h1_lambda_meas_100->Draw("same hist");
    h1_lambda_meas_50->Draw("same hist");
    h1_lambda_meas_10->Draw("same hist");
    h1_lambda_pred_clone->Draw("same axis");
    

    
    roostr = "canv_h1_lambda_meas_diff_100";
    TCanvas *canv_h1_lambda_meas_diff_100 = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_lambda_meas_diff_100, 0.15, 0.1, 0.1, 0.15);   
    h1_lambda_meas_diff_100->Draw("hist");
    h1_lambda_meas_diff_100->SetMinimum(-6);
    h1_lambda_meas_diff_100->SetMaximum(6);    
    h1_lambda_meas_diff_100->SetMarkerSize(0);
    func_xy_title(h1_lambda_meas_diff_100, "Bin index", "#sigma_{i}\' value");
    func_title_size(h1_lambda_meas_diff_100, 0.05, 0.05, 0.05, 0.05);
    h1_lambda_meas_diff_100->GetXaxis()->CenterTitle(); h1_lambda_meas_diff_100->GetYaxis()->CenterTitle();
    h1_lambda_meas_diff_100->GetYaxis()->SetTitleOffset(1.5);

    TF1 *ff_0 = new TF1("ff_100", "0", 0, 1e3); ff_0->Draw("same"); ff_0->SetLineColor(kGray+1); ff_0->SetLineStyle(7);
    TF1 *ffp1 = new TF1("ff_100", "1", 0, 1e3); ffp1->Draw("same"); ffp1->SetLineColor(kGray+1); ffp1->SetLineStyle(7);
    TF1 *ffm1 = new TF1("ff_100", "-1", 0, 1e3); ffm1->Draw("same"); ffm1->SetLineColor(kGray+1); ffm1->SetLineStyle(7);
    TF1 *ffp3 = new TF1("ff_100", "3", 0, 1e3); ffp3->Draw("same"); ffp3->SetLineColor(kGray+1); ffp3->SetLineStyle(7);
    TF1 *ffm3 = new TF1("ff_100", "-3", 0, 1e3); ffm3->Draw("same"); ffm3->SetLineColor(kGray+1); ffm3->SetLineStyle(7);
    
    h1_lambda_meas_diff_50->Draw("same hist");
    h1_lambda_meas_diff_10->Draw("same hist");
    
    h1_lambda_meas_diff_100->Draw("same axis");
      
    /////////////////////////////////

    roostr = "h1_relerr_lambda";
    TH1D *h1_relerr_lambda = new TH1D(roostr, "", 10, 0, 10);
    
    double chi2_test = 0;
    
    for(int idx=0; idx<rows; idx++) {
      
      double val_delta = matrix_delta_lambda(0, idx);
      double val_abserr= sqrt( matrix_cov_lambda(idx, idx) );
      double val_relerr = fabs(val_delta)/val_abserr;
      double val_relerr2 = val_relerr*val_relerr;      
      chi2_test += val_relerr2;
      
      //cout<<TString::Format(" ---> ibin %3d, relerr %6.2f, relerr2 %6.2f", idx+1, val_relerr, val_relerr2)<<endl;

      if( val_relerr<=9 ) h1_relerr_lambda->Fill( val_relerr );
      else h1_relerr_lambda->Fill( 9.5 );
      
    }
    cout<<endl;
    
    roostr = "canv_h1_relerr_lambda";
    TCanvas *canv_h1_relerr_lambda = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_relerr_lambda, 0.15, 0.1, 0.1, 0.15);
    
    h1_relerr_lambda->Draw("hist text");
    h1_relerr_lambda->SetLineWidth(4);
    h1_relerr_lambda->SetMinimum(0);
    h1_relerr_lambda->SetMarkerSize(3);
    func_xy_title(h1_relerr_lambda, "#sigma_{i}\' value", "Entries");
    func_title_size(h1_relerr_lambda, 0.05, 0.05, 0.05, 0.05);
    h1_relerr_lambda->GetXaxis()->CenterTitle(); h1_relerr_lambda->GetYaxis()->CenterTitle();
    h1_relerr_lambda->GetYaxis()->SetTitleOffset(1.5); 
    h1_relerr_lambda->Draw("same axis");
    
    canv_h1_relerr_lambda->SaveAs("canv_h1_relerr_lambda.png");

    /////////////////////////////////    
    
    roostr = "h2_lambda_transform_matrix";
    TH2D *h2_lambda_transform_matrix = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      for(int jbin=1; jbin<=rows; jbin++) {
	double value = matrix_eigenvector(ibin-1, jbin-1);
	h2_lambda_transform_matrix->SetBinContent(ibin, jbin, value);
      }
    }

    roostr = "canv_h2_lambda_transform_matrix";
    TCanvas *canv_h2_lambda_transform_matrix = new TCanvas(roostr, roostr, 900, 850);
    func_canv_margin(canv_h2_lambda_transform_matrix, 0.15, 0.2,0.15,0.2);
    h2_lambda_transform_matrix->Draw("colz");
    func_xy_title(h2_lambda_transform_matrix, "Bin index", "Bin index");
    h2_lambda_transform_matrix->GetXaxis()->CenterTitle(); h2_lambda_transform_matrix->GetYaxis()->CenterTitle();
    canv_h2_lambda_transform_matrix->SaveAs("canv_h2_lambda_transform_matrix.png");
      
  }// if( 1 )
  
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
  
  // for(int ibin=1; ibin<=10; ibin++) {
  //   cout<<" ---> check "<<ibin<<"\t"<<testcn->matrix_fake_meas(0, ibin-1)<<endl;
  // }
  
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
