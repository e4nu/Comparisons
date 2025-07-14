// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TMath.h>
#include "NuHepMC/HepMC3Features.hxx"
//#include "NuHepMC/Reader.hxx"
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
#include <TH1D.h>
#include "Utils.h"

/////////////////////////////////////////////////////////////////////////////////////
// comparison_1p1pi.cxx                                                            //
// --input-hepmc3-file : input MC file in hepmc3 format                            //
// --output-file : output MC file name, without format type. Def: comparison_1p1pi //
/////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

  //std::string input_hepmc3_file = "/exp/genie/app/jtena/Comparisons/src/1161/output_electron_C12_1161_351.root";
  std::string input_hepmc3_file = "/pnfs/genie/persistent/users/jtenavid/e4nu_files/NuHEPMC/Carbon/4453MeV/master-routine_validation_01-eScattering/e_on_1000060120_4453MeV_0.hepmc3";
///exp/genie/app/jtena/GENIE/Generator/src/scripts/production/python/Generator/e_on_1000060120_4453MeV_0.hepmc3";
  //
  std::string topology = "1p1pim", output_name = "comparison_1p1pi", model_name = "MC_Name" ;
  
  // process options
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-hepmc3-file",argc,argv)) {
      input_hepmc3_file = GetArg("input-hepmc3-file",argc,argv);
    } //else { std::cout << " --input-hepmc3-file is not defined "; return 0 ;}
    if( ExistArg("output-file",argc,argv)) {
      output_name = GetArg("output-file",argc,argv);
    }
    if( ExistArg("model-name",argc,argv)) {
      model_name = GetArg("model-name",argc,argv);
    }
    if( ExistArg("topology",argc,argv)) {
      topology = GetArg("topology",argc,argv);
    }
  }

  auto rdr = HepMC3::deduce_reader(input_hepmc3_file);
  // auto rdr = HepMC3::Reader(input_hepmc3_file);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_file << std::endl;
    return 1;
  }

  HepMC3::GenEvent evt;
  rdr->read_event(evt);
  if (rdr->failed()) {
    std::cout << "Failed to read first event from " << input_hepmc3_file << "."
	      << std::endl;
    return 1;
  }

  // Read beam energy and target from first event 
  evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto Beam_E = NuHepMC::Event::GetBeamParticle(evt)->momentum().e();
  auto target_pdg = NuHepMC::Event::GetTargetParticle(evt)->pid(); 
  std::cout << "Analysing electron-scattering events on " << target_pdg << " at " << Beam_E << " GeV." << std::endl;

  // Define Histograms - binning will come from the data files when available.  
  TFile * outfile = new TFile((output_name+".root").c_str(), "RECREATE");
  std::vector<string> observables = { "Efl", "pfl", "pfl_theta", "proton_mom", "proton_theta", "pim_mom", "pim_theta", "ECal", "RecoW", "RecoQ2", "HadDeltaPT", "HadAlphaT" };
  std::map<string,TH1D*> histograms ; 
  for ( unsigned int i = 0 ; i < observables.size() ; ++i ) { 
    if( topology == "1p1pip" && observables[i] == "pim_mom" ) observables[i] == "pip_mom"; 
    if( topology == "1p1pip" && observables[i] == "pim_theta" ) observables[i] == "pip_theta"; 
    histograms[observables[i]] = new TH1D( (model_name+"_1Dxsec_"+observables[i]).c_str(), "",GetBinning(observables[i], Beam_E, topology).size()-1, &GetBinning(observables[i], Beam_E, topology)[0] );
  }

  auto in_gen_run_info = evt.run_info();
  //  auto FATXAcc = FATX::MakeAccumulator(rdr->run_info());
  auto vtx_statuses = NuHepMC::GR5::ReadVertexStatusIdDefinitions(in_gen_run_info);
  auto part_statuses = NuHepMC::GR6::ReadParticleStatusIdDefinitions(in_gen_run_info);
  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // re-open the file so that you start at the beginning
  rdr = HepMC3::deduce_reader(input_hepmc3_file);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_file << std::endl;
    return 1;
  }

  size_t nprocessed = 0;

  while (true) { // loop while there are events
    rdr->read_event(evt);
    if (rdr->failed()) {
      break;
    }
    ++nprocessed;
    // Set units to GeV
    evt.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

    // Get Run Info:
    auto beampt = NuHepMC::Event::GetBeamParticle(evt);
    auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);
    auto primary_vtx = Event::GetPrimaryVertex(evt);
    auto process_id = ER3::ReadProcessID(evt);
    auto proc_ids = GR4::ReadProcessIdDefinitions(in_gen_run_info);
    if (!beampt || !tgtpt) {  // this event didn't have a beam particle or target, its an odd one
      continue;
    }

    // Read in-electron properties:
    auto BeamPdg = beampt->pid();
    auto Ev = beampt->momentum().e();
    auto Pxv = beampt->momentum().px();
    auto Pyv = beampt->momentum().py();
    auto Pzv = beampt->momentum().pz();
    if( BeamPdg != 11 ) { std::cout << "ERROR: Beam is not e-"<<std::endl; break; }

    // Read out-electron :
    auto primary_leptons = NuHepMC::Event::GetParticles_All(evt, NuHepMC::ParticleStatus::UndecayedPhysical, {beampt->pid()} );
    auto fslep = primary_leptons.back();  

    // Apply cut on Q2
    auto q = fslep->momentum() - beampt->momentum();
    auto Q2 = -q.m2();
    if ( Q2 < GetQ2Cut( Ev ) ) continue ; 

    // Smear particles according to detector Resolution
    HepMC3::GenParticlePtr fslep_reco = SmearParticles( fslep, Ev ) ;
    auto FSPrimLept = fslep_reco->pid();
    auto El = fslep_reco->momentum().e();
    auto Pxl = fslep_reco->momentum().px();
    auto Pyl = fslep_reco->momentum().py();
    auto Pzl = fslep_reco->momentum().pz();
    auto Pl = fslep_reco->momentum().p3mod();
    
    // Apply momentum and angle cuts for electron
    if( fslep->momentum().theta() * 180 / TMath::Pi() < GetParticleMinTheta( FSPrimLept, Pl, Ev ) ) continue ;
    if( Pl < GetMinMomentumCut( FSPrimLept, Ev ) ) continue ; 

    // Read output particles: 
    double true_p = 0, true_pip = 0, true_pim = 0, true_gamma = 0, true_n = 0;
    double reco_p = 0, reco_pip = 0, reco_pim = 0, reco_gamma = 0, reco_n = 0;
    std::vector<HepMC3::ConstGenParticlePtr> hadrons = NuHepMC::Event::GetParticles_AllRealFinalState(evt,{});
    std::vector<HepMC3::ConstGenParticlePtr> hadrons_reco;
    HepMC3::ConstGenParticlePtr pion;
    HepMC3::ConstGenParticlePtr proton;
    for( unsigned int i = 0 ; i < hadrons.size() ; ++i ) { 
      double particle_pdg  = hadrons[i]->pid() ;
      double particle_e    = hadrons[i]->momentum().e();
      double particle_px   = hadrons[i]->momentum().px();
      double particle_py   = hadrons[i]->momentum().py();
      double particle_pz   = hadrons[i]->momentum().pz();
      double particle_mod  = hadrons[i]->momentum().p3mod();

      if( particle_pdg == 2212 ) ++true_p ;
      else if ( particle_pdg == 2112 ) ++true_n ;
      else if ( particle_pdg == -211 ) ++true_pim ;
      else if ( particle_pdg == 211 ) ++true_pip;
      else if ( particle_pdg == 22 ) ++true_gamma;
      
      // Smear Hadrons 
      HepMC3::GenParticlePtr hadron_smeared = SmearParticles( hadrons[i], Ev ) ;

      // Apply momentum and angle cuts for hadrons 
      if( hadron_smeared->momentum().theta() * TMath::RadToDeg() < GetParticleMinTheta( particle_pdg, hadron_smeared->momentum().p3mod(), Ev ) ) continue;
      if( hadron_smeared->momentum().theta() * TMath::RadToDeg() > 140 ) continue;
      if( hadron_smeared->momentum().p3mod() < GetMinMomentumCut( particle_pdg, Ev ) ) continue ;
      if( particle_pdg == 22 && !ApplyPhotRadCut(beampt->momentum(),hadron_smeared->momentum())) continue ;

      // Count detected hadrons 
      if( particle_pdg == 2212 ) {
	++reco_p ;
	proton = hadrons[i];
      }
      else if ( particle_pdg == 2112 ) ++reco_n ;
      else if ( particle_pdg == -211 ) {
	++reco_pim ;
	pion = hadrons[i];
      }
      else if ( particle_pdg == 211 ) { 
	++reco_pip;
	pion = hadrons[i];
      }
      else if ( particle_pdg == 22 ) ++reco_gamma;
      
      // Save reconstructed hadrons in new vector
      if( particle_pdg == 2212 || TMath::Abs(particle_pdg) == 211 ) hadrons_reco.push_back(hadron_smeared);
 
    }

    // Select 1p1pim (- or +) events
    if ( reco_p != 1 ) continue ; 
    if ( topology == "1p1pim" ) { 
      if ( reco_pim != 1 ) continue ; 
      if ( reco_pip != 0 ) continue ;
    } else if ( topology == "1p1pip" ) { 
      if ( reco_pim != 0 ) continue ; 
      if ( reco_pip != 1 ) continue ;
    }
    if ( reco_gamma != 0 ) continue ;
    
    // Store histogram with same binning as data. For now is only MC comparison 
    histograms["Efl"]->Fill(El);
    histograms["pfl"]->Fill(Pl);
    histograms["pfl_theta"]->Fill(fslep_reco->momentum().theta()*180/TMath::Pi());
    histograms["proton_mom"]->Fill(proton->momentum().p3mod());
    histograms["proton_theta"]->Fill(proton->momentum().theta()*180/TMath::Pi());
    histograms["pim_mom"]->Fill(pion->momentum().p3mod());
    histograms["pim_theta"]->Fill(pion->momentum().theta()*180/TMath::Pi());
    histograms["ECal"]->Fill(GetECal( fslep_reco, hadrons_reco, tgtpt->pid()));
    histograms["RecoW"]->Fill(GetW(fslep_reco,beampt));
    histograms["RecoQ2"]->Fill(GetQ2(fslep_reco,beampt));
    //    histograms["HadDeltaPT"]
   
    //double evw = FATXAcc->process(evt);
}
  std::cout << " Total processed: " << nprocessed << std::endl;
  // Write the histogram to the file
  for ( auto it = histograms.begin(); it != histograms.end(); it++) (it->second)->Write();

  // Close the file
  outfile->Close();

  // Clean up
  delete outfile;

  //  double fatx = FATXAcc->fatx(cm2ten38_PerNucleon);
  //  double fatx = FATXAcc->fatx(); // in pb/Atom
  //  double sumw = FATXAcc->sumweights();
  //  size_t nevents = FATXAcc->events();

}
