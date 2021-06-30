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

    FLAG_INPUTFILE_COV_HAS_STAT = 0;

    BINS = 0;
    TOYS = 0;
  }

  /////////////////////////////////////////// data member
  
  TRandom3 *rand;

  TMatrixD matrix_pred;       // inputfile
  TMatrixD matrix_meas;       // inputfile
  TMatrixD matrix_syst_abscov;// inputfile
  TMatrixD matrix_stat_abscov;
  TMatrixD matrix_totl_abscov;
  
  bool FLAG_INPUTFILE_COV_HAS_STAT;

  int BINS;
  int TOYS;
  
  TMatrixD matrix_fake_meas;
  map<int, map<int, double> >map_toy_meas;
  
  /////////////////////////////////////////// memberfunction
  
  void Initialization(TString roofile, bool flag_INPUTFILE_COV_HAS_STAT);

  void SetMeas2Expt() {
    matrix_fake_meas.Clear(); matrix_fake_meas.ResizeTo(1, BINS);
    matrix_fake_meas = matrix_meas;
  }
  
};

/////////////////////// ccc

void TCN::Initialization(TString roofile, bool flag_INPUTFILE_COV_HAS_STAT)
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

  ///////////////
  
  matrix_stat_abscov.Clear(); matrix_stat_abscov.ResizeTo(rows, rows);
  matrix_totl_abscov.Clear(); matrix_totl_abscov.ResizeTo(rows, rows);
  
  /// docDB 32520, when the prediction is sufficiently low
  double array_pred_protect[11] = {0, 0.461, 0.916, 1.382, 1.833, 2.298, 2.767, 3.225, 3.669, 4.141, 4.599};

  
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
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";

  
  TCN *testcn = new TCN();

  roostr = "file_numu.root";
  //roostr = "file_numuPC_phi.root";
  //roostr = "file_numu_vtxZ.root";
  //roostr = "file_dQdx.root";
  //roostr = "file_user_Ehad_no.root";
  //roostr = "file_user_Ehad_wi.root";
  //roostr = "file_user_Ehad70_no.root";
  //roostr = "file_user_Ehad70_wi.root";
  //roostr = "file_user_Ehad75_wi.root";
  //roostr = "file_user_Ehad80_wi.root";
  //roostr = "file_user_Ehad_PC_80.root";	
  testcn->Initialization( roostr, 0 );

  if( 1 ) {
    testcn->SetMeas2Expt();
  }

  
}
