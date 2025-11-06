// __________________________________________________________________________
/* This app is used to plot data histograms                                 */
/* The output root file computed with bkg debug mode ("DebugBkg true")      */
/* An example of how to run the analysis with this mode is available in     */
/* ConfFiles/example_configuration.txt                                      */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "Utils.h"

using namespace std;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
// * mc-files : comma-sparated list of MC files                //
// * mc-names : comma-separated list of MC names               //
// * mc-label : comma-separated list of labels for MC          //
// * output-file : file to store plots (in root, pdf... format)//
// * observable : observable used for the x axis definition    //
// * analysis-key: i.e. 1p1pim                                 //
// * beam-energy : energy                                      //
// * log-scale : y log-scale                                   //
// * scale : scaling factor, default 1                         //
// * y-max : maximum y for axis range.                         //
// * plot-ratio                                                //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting histograms..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "MC files not specified. Abort..." << std::endl;
    return 0;
  }

  
  TH1::AddDirectory(kFALSE);

  std::vector<string> mc_files, names_list, mc_label ;
  TFile* in_data = nullptr ;
  TH1D* h_data = nullptr ;
  std::string output_file = "comparisons";
  std::string observable = "ECal";
  std::string analysis_key = "1p1pim";
  double EBeam = 1 ;
  double scale = 1 ; 
  double y_max = -1 ; 
  bool add_ratio = false, is_log = false;
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("mc-files",argc,argv)) {
      string input = GetArg("mc-files",argc,argv);
      stringstream ss(input);
      while( ss.good() ){
	string substr;
	getline( ss, substr, ',' );
	mc_files.push_back( substr );
	names_list.push_back("MC_Name");
	std::cout << " --> " << mc_files.back() << std::endl;
      }
      if( mc_files.size() == 0 ) return 0;
    } else { return 0 ; }

    if( ExistArg("mc-names",argc,argv)) {
      names_list.clear();
      string input = GetArg("mc-names",argc,argv);
      stringstream ss(input);
      while( ss.good() ){
	string substr;
	getline( ss, substr, ',' );
	names_list.push_back(substr) ;
      }
      
      if( names_list.size() != mc_files.size() ) return 0;
    }

    if( ExistArg("mc-label",argc,argv)) {
      string input = GetArg("mc-label",argc,argv);
      stringstream ss(input);
      while( ss.good() ){
	string substr;
	getline( ss, substr, ',' );
	mc_label.push_back(substr);
      }

      if( mc_label.size() != mc_files.size() ) return 0;
    } else mc_label = names_list;

    for( unsigned int ni = 0 ; ni < names_list.size() ; ++ni ) { 
      std::cout << " Generator Name: " << names_list[ni] << std::endl;
    }
    
    if( ExistArg("data",argc,argv)) {
      string data = GetArg("data",argc,argv);
      in_data = new TFile(data.c_str(),"ROOT");
      if( !in_data ) { 
	std::cout << data << " file does not exist. Exiting... " <<std::endl;
	return 0 ; 
      } 
      h_data = (TH1D*)in_data->Get("Data");
      h_data->SetDirectory(0);
    
      in_data->Close();
      delete in_data;
    }    
    
    double energy  = 1 ; 
    if( ExistArg("beam-energy",argc,argv)) {
      string input = GetArg("beam-energy",argc,argv);
      energy = stod(input);
    }
    if( energy > 0  && energy < 1.5 ) EBeam = 1.161 ;
    else if( energy >1.5 && energy < 3 ) EBeam = 2.261 ;
    else if( energy >3 && energy < 5 ) EBeam = 4.461 ;

    if( ExistArg("observable", argc,argv)) { 
      observable = GetArg("observable",argc,argv); 
    }
    if( ExistArg("output-file", argc,argv)) { 
      output_file = GetArg("output-file",argc,argv); 
    }
    if( ExistArg("add-ratio", argc, argv)) { 
      add_ratio = true ; 
    }
    if( ExistArg("log-scale", argc, argv)) { 
      is_log = true ; 
    }
    if( ExistArg("scale",argc,argv)) {
      scale = stod( GetArg("scale",argc,argv) ) ; 
    }
    if( ExistArg("y-max",argc,argv)) {
      y_max = stod( GetArg("y-max",argc,argv) ) ; 
    }
  }

  std::cout << "Plotting..." << std::endl;
  // Color palette
  std::vector<int> color_list = {kBlack,kBlue-4,kOrange+1,kGreen+3,kMagenta+2,kRed-4,kSpring+7,kTeal-6,kCyan+3,kAzure-1,kBlue-8,kViolet-4,kMagenta-10,kPink-9} ;

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  TPad *pad1, *pad2 ;
  if (add_ratio) {
    pad1 = new TPad("pad1","",0,0.4,1,1);
    pad2 = new TPad("pad2","",0,0,1,0.4);
    pad1->Draw();
    pad1->cd();
    pad1->SetBottomMargin(0.015);
    pad1->SetLeftMargin(0.2);
    pad2->SetLeftMargin(0.2);
    pad1->SetRightMargin(0.01);
    pad2->SetRightMargin(0.01);
    pad2->SetBottomMargin(0.15);
  } else {
    pad1 = new TPad("pad1","",0,0,1,1);
    pad2 = new TPad("pad2","",0,0,1,1);
    pad1->Draw();
    pad1->cd();
    pad1->SetBottomMargin(0.25);
    pad1->SetTopMargin(0.1);
    pad1->SetLeftMargin(0.25);
    pad2->SetLeftMargin(0.25);
    pad1->SetRightMargin(0.01);
    pad2->SetRightMargin(0.01);
  }

  auto legend = new TLegend(0.3,0.7,0.55,0.9);
  legend->SetBorderSize(0);
  legend->SetTextFont(132);
  legend->SetTextSize(0.08);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);

  std::vector<TFile*> in_root_files ;
  std::vector<TH1D*> hists, hists_SPP;
  double ymax = 0 ;
  for( unsigned int id = 0 ; id < mc_files.size(); ++id ){
    in_root_files.push_back(new TFile((mc_files[id]).c_str(),"ROOT"));
    if( !in_root_files[id] ) { std::cout << "ERROR: the "<< mc_files[id]<<" does not exist." << std::endl; return 0;}
    else std::cout << "Reading " + mc_files[id] << std::endl;
    std::string hist_name = names_list[id]+"_1Dxsec_"+observable+"_total";
    hists.push_back( (TH1D*)in_root_files[id]->Get(hist_name.c_str()) );
    if( !hists[id] ) { std::cout << "ERROR: the histogram " << hist_name << " does not exist." <<std::endl; return 0;}
    hists[id]->SetDirectory(0);
    std::string hist_name_SPP = names_list[id]+"_1Dxsec_"+observable+"_SPP";
    hists_SPP.push_back( (TH1D*)in_root_files[id]->Get(hist_name_SPP.c_str()) );
    if( !hists_SPP[id] ) { std::cout << "ERROR: the histogram " << hist_name_SPP << " does not exist." <<std::endl; return 0;}
    hists_SPP[id]->SetDirectory(0);

    // Find maximum 
    double max = 0 ;
    for( int k = 0 ; k < hists[id]->GetNbinsX() ; ++k ) {
      if ( hists[id]->GetBinContent(k) > max ) max = hists[id]->GetBinContent(k) ;
    }
    if( max > ymax ) ymax = max;
  }
  ymax *= ( 1+0.2 ) ;

  if( y_max > 0 ) ymax = y_max; 

  for (auto f : in_root_files) {
    f->Close();
    delete f;
  }

  for( unsigned int i = 0 ; i < hists.size(); ++i ){
    StandardFormat( hists[i], "", color_list[i+1], 1, observable, is_log, ymax );
    hists[i] -> SetLineStyle(1);
    hists[i] -> SetMarkerStyle(0);
    hists[i] -> GetYaxis() -> TAxis::SetMaxDigits(3);
    hists[i] -> SetStats(0);
    hists[i] -> SetMarkerSize(1.6);
    hists[i] -> Scale( scale );

    hists[i]->GetYaxis()->SetRangeUser(0.001, ymax); 
    if( is_log ) {
      hists[i]->GetYaxis()->SetRangeUser(0.1, ymax);
    } 
    double I_error = 0 ; 
    double I = hists[i]->IntegralAndError(1, hists[i]->GetNbinsX(), I_error, "width");
    std::cout << " MC " << names_list[i] +" Integral: " << I << "+-"<< I_error << "nb/Obs"<< std::endl;

    if( i == 0 ) hists[i] -> Draw("hist err");
    else hists[i] -> Draw("hist err same");

    StandardFormat( hists_SPP[i], "", color_list[i+1], 2, observable, is_log, ymax );
    hists_SPP[i] -> SetMarkerStyle(0);
    hists_SPP[i] -> Draw("hist err same");
    legend->AddEntry(hists[i],(mc_label[i]).c_str());
    legend->AddEntry(hists_SPP[i],(mc_label[i]+" SPP").c_str());
  }

  // Plot
  if( h_data != nullptr ) {
    h_data ->Scale(scale);

    double I_error = 0 ; 
    double I = h_data->IntegralAndError(1, h_data->GetNbinsX(), I_error, "width");
    std::cout << " Data Integral: " << I << "+-"<< I_error << "nb/Obs"<< std::endl;

    h_data ->Draw("err same");
    StandardFormat( h_data, "", color_list[0], 1, observable, is_log, ymax );
    h_data -> SetMarkerStyle(8);
    h_data -> SetMarkerSize(1.7);
    legend->AddEntry(h_data,"Data");
  }

  legend->Draw();

  if( add_ratio ) {
    
    c->cd();
    pad2->Draw();
    pad2->cd();
    std::vector<TH1D*> ratios;
    for(unsigned int i = 0 ; i < hists.size() ; ++i ){
      ratios.push_back((TH1D*)hists[i]->Clone());
      ratios[i]->Scale(-1.);
      ratios[i]->Add(hists.back());
      ratios[i]->Divide(hists.back());
      ratios[i]->Scale(100);
      ratios[i]->SetStats(0);
      StandardFormat( ratios[i], "", color_list[i], 1, observable, is_log, ymax );
      ratios[i]->GetYaxis()->SetTitle("MC/Data");
      ratios[i]->SetMarkerStyle(8);
      ratios[i]->SetMarkerSize(1.6);
      ratios[i]->Scale(-1.);
      ratios[i]->GetYaxis()->SetRangeUser(-100,100);
      ratios[i]->GetXaxis()->SetTitleOffset(0.9);
      ratios[i]->GetXaxis()->SetTitleSize(0.2);
      ratios[i]->GetYaxis()->SetTitleSize(0.13);
      ratios[i]->GetXaxis()->SetLabelSize(0.2);
      ratios[i]->GetYaxis()->SetLabelSize(0.2);
      if( i == 0 ) ratios[i] -> Draw("err");
      else ratios[i] -> Draw("err same");
    }
  }
  std::string output_name = output_file ;

  c->SaveAs((output_file+".root").c_str());
  c->SaveAs((output_file+".pdf").c_str());

  return 0 ;

}
