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

    rand = new TRandom3(0);

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

  if( flag_norm ) {

    ///////////
    TPrincipal principal_obj_syst(BINS, "ND");    
    TPrincipal principal_obj_total(BINS, "ND");    
    
    for(int itoy=1; itoy<=TOYS; itoy++) {
      double val_rand_relerr = rand->Gaus(0, norm_relerr);

      ///////////
      double *array_obj_syst = new double[BINS];
      double *array_obj_total = new double[BINS];
      
      for( int ibin=1; ibin<=BINS; ibin++ ) {
        double pred_cv = matrix_pred(0, ibin-1);
        double pred_var = pred_cv * (1+val_rand_relerr);
        if( pred_var<0 ) pred_var = 0;// ??? remove this requirement or not
        
        array_obj_syst[ibin-1] = pred_var;
        map_toy_meas[itoy][ibin] = rand->Poisson( pred_var );
        array_obj_total[ibin-1] =map_toy_meas[itoy][ibin];
        
      }// for( int ibin=1; ibin<=BINS; ibin++ )

      principal_obj_syst.AddRow( array_obj_syst );
      principal_obj_total.AddRow( array_obj_total );

      delete[] array_obj_total;
      delete[] array_obj_syst;
      
    }// for(int itoy=1; itoy<=TOYS; itoy++)

    //////////////////////////////////////////////////////////// obj_syst
    
    TMatrixD *matrix_abscov_obj_syst = (TMatrixD *)principal_obj_syst.GetCovarianceMatrix();
    
    for(int idx=0; idx<BINS; idx++) {
      for(int jdx=0; jdx<BINS; jdx++) {
        if(idx<jdx) (*matrix_abscov_obj_syst)(idx, jdx) = (*matrix_abscov_obj_syst)(jdx, idx);
      }
    }

    matrix_syst_abscov.Clear(); matrix_syst_abscov.ResizeTo(BINS, BINS);
    matrix_syst_abscov = (*matrix_abscov_obj_syst);

    ///////////
    TMatrixD matrix_relcov_obj_syst(BINS, BINS);
    TMatrixD matrix_correlation_obj_syst(BINS, BINS);
    
    for(int idx=0; idx<BINS; idx++) {
      for(int jdx=0; jdx<BINS; jdx++) {
        double cov_ij = (*matrix_abscov_obj_syst)(idx, jdx);
        double cov_ii = (*matrix_abscov_obj_syst)(idx, idx);
        double cov_jj = (*matrix_abscov_obj_syst)(jdx, jdx);
        
        double cv_i = matrix_pred(0, idx);
        double cv_j = matrix_pred(0, jdx);
        
        if( cv_i!=0 && cv_j!=0 ) matrix_relcov_obj_syst(idx, jdx) = cov_ij/cv_i/cv_j;
        
        if( cov_ii!=0 && cov_jj!=0 ) matrix_correlation_obj_syst(idx, jdx) = cov_ij/sqrt(cov_ii)/sqrt(cov_jj);
        if( idx==jdx ) matrix_correlation_obj_syst(idx, jdx) = 1;       
      }
    }
    
    // ///////////
    // roostr = "canv_matrix_relcov_obj_syst";
    // TCanvas *canv_matrix_relcov_obj_syst = new TCanvas(roostr, roostr, 900, 850);
    // func_canv_margin(canv_matrix_relcov_obj_syst, 0.15, 0.2,0.15,0.2);
    // roostr = "h2_relcov_obj_syst";
    // TH2D *h2_relcov_obj_syst = new TH2D(roostr, "", BINS, 0, BINS, BINS, 0, BINS);
    // for(int idx=0; idx<BINS; idx++)
    //   for(int jdx=0; jdx<BINS; jdx++)
    //  h2_relcov_obj_syst->SetBinContent( idx+1, jdx+1, matrix_relcov_obj_syst(idx, jdx) );
    // h2_relcov_obj_syst->Draw("colz");
    // func_xy_title(h2_relcov_obj_syst, "Bin index", "Bin index");
   
    // ///////////
    // roostr = "canv_matrix_correlation_obj_syst";
    // TCanvas *canv_matrix_correlation_obj_syst = new TCanvas(roostr, roostr, 900, 850);
    // func_canv_margin(canv_matrix_correlation_obj_syst, 0.15, 0.2,0.15,0.2);
    // roostr = "h2_correlation_obj_syst";
    // TH2D *h2_correlation_obj_syst = new TH2D(roostr, "", BINS, 0, BINS, BINS, 0, BINS);
    // for(int idx=0; idx<BINS; idx++)
    //   for(int jdx=0; jdx<BINS; jdx++)
    //  h2_correlation_obj_syst->SetBinContent( idx+1, jdx+1, matrix_correlation_obj_syst(idx, jdx) );
    // h2_correlation_obj_syst->Draw("colz");
    // func_xy_title(h2_correlation_obj_syst, "Bin index", "Bin index");

  }// if( flag_norm )
  else {

    cout<<endl<<" ---> Check, producing toys"<<endl<<endl;
    
    TMatrixD matrix_stat_abscov(BINS, BINS);
    for(int idx=0; idx<BINS; idx++) matrix_stat_abscov(idx, idx) = matrix_pred(0, idx);

    TMatrixD matrix_total_abscov(BINS, BINS);
    //matrix_total_abscov = matrix_stat_abscov + matrix_syst_abscov;
    matrix_total_abscov = matrix_syst_abscov;

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
    
  }// else of if( flag_norm )
  
}

/////////////////////// ccc

void TCN::Initialization(bool flag_norm)
{
  TString roostr = "";
  
  FLAG_NORM = flag_norm;
  
  if( FLAG_NORM ) {

    ////////////////////////////////
    
    BINS = 10;
    
    matrix_pred.Clear(); matrix_pred.ResizeTo(1, BINS);
    for(int ibin=1; ibin<=BINS; ibin++) {
      matrix_pred(0, ibin-1) = 1e6;
    }
  }// if( FLAG_NORM )
  else {
    
    ////////////////////////////////

    TFile *file_numu = new TFile("file_numu.root", "read");
  
    TMatrixD *matrix_gof_pred = (TMatrixD*)file_numu->Get("matrix_gof_pred");
    TMatrixD *matrix_gof_meas = (TMatrixD*)file_numu->Get("matrix_gof_meas");
    TMatrixD *matrix_gof_syst = (TMatrixD*)file_numu->Get("matrix_gof_syst");
    
    int rows = matrix_gof_syst->GetNrows();
    BINS = rows;

    // for(int idx=0; idx<rows; idx++)
    //   cout<<" ---> check "<<idx+1<<"\t"<<(*matrix_gof_pred)(0, idx)<<"\t"<<sqrt( (*matrix_gof_syst)(idx, idx) )<<endl;
    
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
    roostr = "canv_h2_gof_relcov";
    TCanvas *canv_h2_gof_relcov = new TCanvas(roostr, roostr, 900, 850);
    func_canv_margin(canv_h2_gof_relcov, 0.15, 0.2,0.15,0.2);
    h2_gof_relcov->Draw("colz");
    func_xy_title(h2_gof_relcov, "Bin index", "Bin index");
    h2_gof_relcov->GetXaxis()->CenterTitle(); h2_gof_relcov->GetYaxis()->CenterTitle();
    canv_h2_gof_relcov->SaveAs("canv_h2_gof_relcov.png");
      
    ///////////
    roostr = "canv_h2_gof_correlation";
    TCanvas *canv_h2_gof_correlation = new TCanvas(roostr, roostr, 900, 850);
    func_canv_margin(canv_h2_gof_correlation, 0.15, 0.2,0.15,0.2);
    h2_gof_correlation->Draw("colz");
    h2_gof_correlation->GetZaxis()->SetRangeUser(-1,1);
    func_xy_title(h2_gof_correlation, "Bin index", "Bin index");       
    h2_gof_correlation->GetXaxis()->CenterTitle(); h2_gof_correlation->GetYaxis()->CenterTitle();
    canv_h2_gof_correlation->SaveAs("canv_h2_gof_correlation.png");
      
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
  
  bool flag_norm = 0;  
  double norm_relerr = 0.2;
  
  TCN *testcn = new TCN();
  
  testcn->Initialization( flag_norm );
 
  testcn->Set_norm_relerr( norm_relerr );

  int ntoys = 100000;
  testcn->ProduceVariation( ntoys, flag_norm );

  /* 
  double low_chi2 = 0;
  double hgh_chi2 = 100;
  int bins_chi2 = 500;
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

  if( f1_chi2_temp->GetMaximum() > h1_toy_chi2->GetMaximum() ) {
    h1_toy_chi2->SetMaximum( f1_chi2_temp->GetMaximum() * 1.1 );
  }
  
  canv_h1_toy_chi2->SaveAs("canv_h1_toy_chi2.png");
  */  
}
