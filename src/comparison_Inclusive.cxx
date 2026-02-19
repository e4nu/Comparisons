// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TMath.h>
#include "Utils.h"
#include <filesystem>
#include "NuHepMC/HepMC3Features.hxx"
#include "NuHepMC/Reader.hxx"
#include "NuHepMC/EventUtils.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/WriterUtils.hxx"
#include "NuHepMC/make_writer.hxx"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderFactory.h"
#include "NuHepMC/FATXUtils.hxx"
#include "TFile.h"

//////////////////////////////////////////////////////////////////////////////////////////
// comparison_1p1pi.cxx                                                                 //
// --input-hepmc3-file : input MC file in hepmc3 format                                 //
// --output-file : output MC file name, without format type. Def: comparison_Inclusive  //
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

  std::vector<std::string> input_hepmc3_files = {};
  std::string topology = "Inclusive", output_name = "comparison_Inclusive", model_name = "MC_Name" ;
  
  // process options
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-hepmc3-directory",argc,argv)) {
      for (const auto& entry : std::filesystem::directory_iterator(GetArg("input-hepmc3-directory",argc,argv))) {
        if (entry.is_regular_file() && entry.path().extension() == ".hepmc3") {
	  input_hepmc3_files.push_back(entry.path().string());  // full path
        }
      }
    } else { std::cout << "Missing input HEPMC3 file."; return 0 ; }

    if( ExistArg("output-file",argc,argv)) {
      output_name = GetArg("output-file",argc,argv);
    }
    if( ExistArg("model-name",argc,argv)) {
      model_name = GetArg("model-name",argc,argv);
    }
  } else { 
    std::cout << "Options: \n --input-hepmc3-directory : specify input directory where your hepmc files live. \n --output-file : name for output files. Default comparison_Inclusive. \n --model-name : Give a name to your model. I.e: GENIE. Default MC_Name " << std::endl;
    return 0;
  }

  if(!input_hepmc3_files.size()){
    std::cout << "found no files." << std::endl;
    return 1;
  }

  // Analyse file(s) : 
  auto rdr = NuHepMC::Reader(input_hepmc3_files[0]); // If file is nuhepmc v0.9
  // Check it doesnt fail -!!

  HepMC3::GenEvent evt;
  if ( ! rdr.read_event(evt) ) {
    std::cout << "Failed to read first event from " << input_hepmc3_files[0] << "."                
              << std::endl;                                                                                                                                    
    return 1;  
  }

  // Read beam energy and target from first event 
  evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto Beam_E = NuHepMC::Event::GetBeamParticle(evt)->momentum().e();
  auto target_pdg = NuHepMC::Event::GetTargetParticle(evt)->pid(); 
  std::cout << "Analysing electron-scattering events on " << target_pdg << " at " << Beam_E << " GeV." << std::endl;

  // Define histograms given Beam_E
  std::vector<string> observables = { "Efl", "pfl", "RecoW", "RecoQ2" };
  std::vector<string> process = {"total", "QEL", "RES", "NonRES", "MEC", "DIS", "0PP", "SPP", "MPP" } ; 

  // Define Histograms - binning will come from the data files when available.  
  std::map<string,TH1D*> histograms ; 
  for ( unsigned int i = 0 ; i < observables.size() ; ++i ) { 
    // Add breakdown 
    for ( unsigned j = 0 ; j < process.size() ; ++j ) { 
      histograms[observables[i]+process[j]] = new TH1D( (model_name+"_1Dxsec_"+observables[i]+"_"+process[j]).c_str(), process[j].c_str(),GetBinning(observables[i], Beam_E, topology).size()-1, &GetBinning(observables[i], Beam_E, topology)[0] );
    }
  }
  
  TFile * outfile = new TFile((output_name+".root").c_str(), "RECREATE");

  auto in_gen_run_info = evt.run_info();
  auto FATXAcc = FATX::MakeAccumulator(rdr.run_info());
  auto vtx_statuses = NuHepMC::GR9::ReadVertexStatusIdDefinitions(in_gen_run_info);
  auto part_statuses = NuHepMC::GR10::ReadParticleStatusIdDefinitions(in_gen_run_info);
  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // re-open the file so that you start at the beginning
  size_t nprocessed = 0 ;
  
  double has_chex = 0 ;
  double in_peak = 0;
  bool first_pass = false ; 
  for( unsigned int id = 0 ; id < input_hepmc3_files.size() ; ++id ) { 
    auto rdr = NuHepMC::Reader(input_hepmc3_files[id]);
    while (true) { // loop while there are events
      rdr.read_event(evt);
      if( rdr.failed() ) break;

      ++nprocessed;

      // Set units to GeV
      evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
    
      // Get Run Info:
      auto beampt = NuHepMC::Event::GetBeamParticle(evt);
      auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);
      auto process_id = ER3::ReadProcessID(evt);
      auto proc_ids = GR8::ReadProcessIdDefinitions(in_gen_run_info);

      if (!beampt || !tgtpt) {  // this event didn't have a beam particle or target, its an odd one
	continue;
      }

      // Read XSec
      double evw = FATXAcc->process(evt);

      // Read in-electron properties:
      auto BeamPdg = beampt->pid();
      auto Ev = beampt->momentum().e();
      auto Pxv = beampt->momentum().px();
      auto Pyv = beampt->momentum().py();
      auto Pzv = beampt->momentum().pz();
      if( BeamPdg != 11 ) { std::cout << "ERROR: Beam is not e-"<<std::endl; break; }

      // Read out-electron :
      auto primary_leptons = NuHepMC::Event::GetParticles_All(evt, NuHepMC::ParticleStatus::UndecayedPhysical, {beampt->pid()} );
      if( primary_leptons.size() == 0 ) continue ; // some NEUT events are pauli blocked but still stored. Skip them.
      auto fslep = primary_leptons.back();  

      // Smear e- according to detector Resolution
      HepMC3::GenParticlePtr fslep_reco = SmearParticles( fslep, Ev ) ;
      auto FSPrimLept = fslep_reco->pid();
      auto El = fslep_reco->momentum().e();
      auto Pxl = fslep_reco->momentum().px();
      auto Pyl = fslep_reco->momentum().py();
      auto Pzl = fslep_reco->momentum().pz();
      auto Pl = fslep_reco->momentum().p3mod();

      // Apply cut on Q2
      auto q = fslep_reco->momentum() - beampt->momentum();
      auto Q2 = -q.m2();
      if ( Q2 < GetQ2Cut( Ev ) ) continue ; 

      // Apply momentum and angle cuts for electron
      if( fslep_reco->momentum().theta() * 180 / TMath::Pi() < GetParticleMinTheta( FSPrimLept, Pl, Ev ) ) continue ;
      if( fslep_reco->momentum().theta() * 180 / TMath::Pi() > 45 /*deg*/ ) continue ;
      if( Pl < GetMinMomentumCut( FSPrimLept, Ev ) ) continue ; 

      std::string process_name = "total" ;
      
      histograms["Efl"+process_name]->Fill(El,evw);
      histograms["pfl"+process_name]->Fill(Pl,evw);
      histograms["pfl_theta"+process_name]->Fill(fslep_reco->momentum().theta()*180/TMath::Pi(),evw);
      histograms["RecoQ2"+process_name]->Fill(GetQ2(fslep_reco,beampt),evw);
      histograms["RecoW"+process_name]->Fill(GetW(fslep_reco,beampt),evw);
      
      if ( process_id >= 200 && process_id < 300 ) process_name = "QEL";
      if ( process_id >= 300 && process_id < 400 ) process_name = "MEC";
      if ( process_id >= 400 && process_id < 500 ) process_name = "RES" ;
      if ( process_id >= 600 && process_id < 700 ) {
	// Add check here if GENIE : 
	if( GetW(fslep_reco,beampt) < 1.7 ) process_name = "NonRES" ;
	else process_name = "DIS";
      }

      histograms["Efl"+process_name]->Fill(El,evw);
      histograms["pfl"+process_name]->Fill(Pl,evw);
      histograms["pfl_theta"+process_name]->Fill(fslep_reco->momentum().theta()*180/TMath::Pi(),evw);
      histograms["RecoQ2"+process_name]->Fill(GetQ2(fslep_reco,beampt),evw);
      histograms["RecoW"+process_name]->Fill(GetW(fslep_reco,beampt),evw);

    }
  }

  std::cout << " Analised " << nprocessed << " events. "<< std::endl;

  double fatx = FATXAcc->fatx() * 1E-3 ; // in nb/Atom
  double sumw = FATXAcc->sumweights();
  std::cout << fatx << " / " << sumw << std::endl;
  // Write the histogram to the file
  for ( auto it = histograms.begin(); it != histograms.end(); it++) { 
    // Normalize by bin witdh and xsection
    NormalizeHist(it->second, fatx/sumw );
    // Save
    (it->second)->Write();
  }
 
  // Close the file
  outfile->Close();

  // Clean up
  delete outfile;
       
}


