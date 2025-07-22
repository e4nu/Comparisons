// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include "Utils.h"
#include <filesystem>
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

/////////////////////////////////////////////////////////////////////////////////////
// comparison_1p1pi.cxx                                                            //
// --input-hepmc3-file : input MC file in hepmc3 format                            //
// --output-file : output MC file name, without format type. Def: comparison_1p1pi //
/////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

  //std::string input_hepmc3_file = "/exp/genie/app/jtena/Comparisons/src/1161/output_electron_C12_1161_351.root";
  std::vector<std::string> input_hepmc3_files ={ "/pnfs/genie/persistent/users/jtenavid/e4nu_files/NuHEPMC/Carbon/4453MeV/master-routine_validation_01-eScattering/e_on_1000060120_4453MeV_0.hepmc3"};
  ///exp/genie/app/jtena/GENIE/Generator/src/scripts/production/python/Generator/e_on_1000060120_4453MeV_0.hepmc3";
  //
  std::string topology = "1p1pim", output_name = "comparison_1p1pi", model_name = "MC_Name" ;
  
  // process options
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-hepmc3-directory",argc,argv)) {
      input_hepmc3_files = {};
      for (const auto& entry : std::filesystem::directory_iterator(GetArg("input-hepmc3-directory",argc,argv))) {
        if (entry.is_regular_file() && entry.path().extension() == ".hepmc3") {
	  input_hepmc3_files.push_back(entry.path().string());  // full path
        }
      }
    } 
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

  auto rdr = HepMC3::deduce_reader(input_hepmc3_files[0]);
  // auto rdr = HepMC3::Reader(input_hepmc3_file);
  if (!rdr) {
    std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_files[0] << std::endl;
    return 1;
  }

  HepMC3::GenEvent evt;
  rdr->read_event(evt);
  if (rdr->failed()) {
    std::cout << "Failed to read first event from " << input_hepmc3_files[0] << "."
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
  std::vector<string> process = {"total", "QEL", "RES", "NonRES", "MEC", "DIS"} ; 
  // In GENIE, NonRES = DIS W < 1.7 GeV. Add the specific definition for your generator if this process is included. 

  std::map<string,TH1D*> histograms ; 
  for ( unsigned int i = 0 ; i < observables.size() ; ++i ) { 
    if( topology == "1p1pip" && observables[i] == "pim_mom" ) observables[i] == "pip_mom"; 
    if( topology == "1p1pip" && observables[i] == "pim_theta" ) observables[i] == "pip_theta"; 

    // Add breakdown 
    for ( unsigned j = 0 ; j < process.size() ; ++j ) { 
      histograms[observables[i]+process[j]] = new TH1D( (model_name+"_1Dxsec_"+observables[i]+"_"+process[j]).c_str(), process[j].c_str(),GetBinning(observables[i], Beam_E, topology).size()-1, &GetBinning(observables[i], Beam_E, topology)[0] );
    }
  }

  // Adding Histograms for other studies here: 
  std::map<string,TH1D*> histograms_studies ; 
  histograms_studies["CExEnergyTransferNeutron"] = new TH1D( (model_name+"_ChExEnergyTransfer_Neutron").c_str(), "Relative Difference", 30, -10, 100 );
  histograms_studies["CExEnergyTransferPi0"] = new TH1D( (model_name+"_ChExEnergyTransfer_Pi0").c_str(), "Relative Difference", 30, -10, 100 );
  
  std::map<string,TH2D*> histograms2D ; 
  for ( unsigned j = 0 ; j < process.size() ; ++j ) { 
    histograms2D["RecoW,RecoQ2"+process[j]] = new TH2D( (model_name+"_2Dxsec_RecoW_vs_RecoQ2"+"_"+process[j]).c_str(), "", GetBinning("RecoW", Beam_E, topology).size()-1, &GetBinning("RecoW", Beam_E, topology)[0], GetBinning("RecoQ2", Beam_E, topology).size()-1, &GetBinning("RecoQ2", Beam_E, topology)[0]);
    histograms2D["RecoW,HadSystemMass"+process[j]] = new TH2D( (model_name+"_2Dxsec_RecoW_vs_MHad"+"_"+process[j]).c_str(), "", GetBinning("RecoW", Beam_E, topology).size()-1, &GetBinning("RecoW", Beam_E, topology)[0], GetBinning("HadSystemMass", Beam_E, topology).size()-1, &GetBinning("HadSystemMass", Beam_E, topology)[0]);
  }

  auto in_gen_run_info = evt.run_info();
  //  auto FATXAcc = FATX::MakeAccumulator(rdr->run_info());
  auto vtx_statuses = NuHepMC::GR5::ReadVertexStatusIdDefinitions(in_gen_run_info);
  auto part_statuses = NuHepMC::GR6::ReadParticleStatusIdDefinitions(in_gen_run_info);
  auto out_gen_run_info = std::make_shared<HepMC3::GenRunInfo>(*in_gen_run_info);

  // re-open the file so that you start at the beginning
  size_t nprocessed = 0;
  double xsec = 0 ; // in pb 

  
  double has_chex = 0 ;
  double in_peak = 0;
  for( unsigned int id = 0 ; id < input_hepmc3_files.size() ; ++id ) { 
    rdr = HepMC3::deduce_reader(input_hepmc3_files[id]);
    if (!rdr) {
      std::cout << "Failed to instantiate HepMC3::Reader from " << input_hepmc3_files[id] << std::endl;
      return 1;
    }

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

      // Read XSec 
      xsec = evt.cross_section()->xsec(); 

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
    
      bool qe = false, res = false, mec = false, nonres = false, dis = false;
      if ( process_id >= 200 && process_id < 300 ) qe = true ;
      if ( process_id >= 300 && process_id < 400 ) mec = true ;
      if ( process_id >= 400 && process_id < 500 ) res = true ;
      if ( process_id >= 600 && process_id < 700 ) {
	// Add check here if GENIE : 
	if( GetW(fslep_reco,beampt) < 1.7 ) nonres = true ; 
	else dis = true ;
      }
      // Store histogram with same binning as data. For now is only MC comparison 
      for( unsigned int pid = 0 ; pid < process.size() ; ++pid ){
	if ( process[pid] == "QEL" && !qe ) continue ;
	if ( process[pid] == "MEC" && !mec ) continue ;
	if ( process[pid] == "RES" && !res ) continue ;
	if ( process[pid] == "DIS" && !dis ) continue ;
	//if( dis ) std::cout << GetW(fslep_reco,beampt) << " **** W " << std::endl;
	
	histograms["Efl"+process[pid]]->Fill(El);
	histograms["pfl"+process[pid]]->Fill(Pl);
	histograms["pfl_theta"+process[pid]]->Fill(fslep_reco->momentum().theta()*180/TMath::Pi());
	histograms["proton_mom"+process[pid]]->Fill(proton->momentum().p3mod());
	histograms["proton_theta"+process[pid]]->Fill(proton->momentum().theta()*180/TMath::Pi());
	histograms["pim_mom"+process[pid]]->Fill(pion->momentum().p3mod());
	histograms["pim_theta"+process[pid]]->Fill(pion->momentum().theta()*180/TMath::Pi());
	histograms["ECal"+process[pid]]->Fill(GetECal( fslep_reco, hadrons_reco, tgtpt->pid()));
	histograms["RecoW"+process[pid]]->Fill(GetW(fslep_reco,beampt));
	histograms["RecoQ2"+process[pid]]->Fill(GetQ2(fslep_reco,beampt));
	histograms["HadAlphaT"+process[pid]]->Fill(DeltaAlphaT(fslep_reco, hadrons_reco));
	histograms["HadDeltaPT"+process[pid]]->Fill(DeltaPT(fslep_reco, hadrons_reco).Mag());
	
	histograms2D["RecoW,RecoQ2"+process[pid]]->Fill(GetW(fslep_reco,beampt),GetQ2(fslep_reco,beampt));
	histograms2D["RecoW,HadSystemMass"+process[pid]]->Fill(GetW(fslep_reco,beampt),HadSystemMass(hadrons_reco));
      }

      // Loop over event particles before and after FSI
      std::map<int,HepMC3::ConstGenParticlePtr> PrimaryParticlesOut ;
      std::map<int,HepMC3::ConstGenParticlePtr> SecondaryVertexIn ;
      std::map<int,HepMC3::ConstGenParticlePtr> SecondaryVertexOut ;

      int DeltaPP = 2224;
      int DeltaM = 1114;
      int Delta0 = 2114;
      int DeltaP = 2214;

      unsigned int np_vertex = 0, nn_vertex = 0, npim_vertex = 0, npip_vertex = 0;
      unsigned int np_postFSI = 0, nn_postFSI = 0, npim_postFSI = 0, npip_postFSI = 0;
      
      bool has_neutrons = false ;
      bool has_cex = false;
      bool has_resonance = false; 
      
      if( GetECal( fslep_reco, hadrons_reco, tgtpt->pid()) < 2.15 ) continue ; 
      ++in_peak;

      // Find primary particles 
      for( auto const &vtx : evt.vertices()){
	if( vtx->status() == 1 ) { 
	  for (auto const &part : vtx->particles_out()) {
	    // Count primary particles
	    if( part->status() == 26 ) { 
	      if( part->pid() == 2212 ) ++np_vertex ;
	      else if( part->pid() == 2112 ) ++nn_vertex ; 
	      else if( part->pid() == -211 ) ++npim_vertex ; 
	      else if( part->pid() == 211 ) ++npip_vertex ; 
	    }	      
	    // Store in map 
	    PrimaryParticlesOut[part->id()] = part ;
	  } 
	}
      }

      // Map with products
      for( auto const &vtx : evt.vertices()){
	if ( vtx->status() == 12 ) {
	  for (auto const &part : vtx->particles_in()) SecondaryVertexIn[part->id()] = part ;
	  for (auto const &part : vtx->particles_out()) {
	    // Count pre fsi hadrons
	    if( part->status() == 26 ) { 
	      if( part->pid() == 2212 ) ++np_vertex ;
	      else if( part->pid() == 2112 ) ++nn_vertex ; 
	      else if( part->pid() == -211 ) ++npim_vertex ; 
	      else if( part->pid() == 211 ) ++npip_vertex ; 
	    }	      

	    // Count final stable particles 
	    if( part->status() == 1 ) { 
	      if( part->pid() == 2212 ) ++np_postFSI ;
	      else if( part->pid() == 2112 ) ++nn_postFSI ; 
	      else if( part->pid() == -211 ) ++npim_postFSI ; 
	      else if( part->pid() == 211 ) ++npip_postFSI ; 
	    }
	    SecondaryVertexOut[part->id()] = part ;
	  }
	}
	       
	for (auto const &part : vtx->particles_in()) {
	  if( part->status() == 26 ) { // before FSI 
	    bool is_neutral = false ; 
	    int is_neutral_id = 0;
	    double initial_T = 0;
	    
	    // Find neutrals 
	    auto pid_neutral = part->pid();
	    if( pid_neutral == 2112 || pid_neutral == 111 ) {
	      is_neutral = true ; 
	      is_neutral_id = part->id() ;
	      initial_T = part->momentum().e()-part->momentum().m() ;
	      
	      if( vtx->particles_out().size() == 1 ) continue ; 
	      std::cout << part->pid() << " -> ";
	      double charge = 0 ;
	      double final_T = 0 ; 
	      for (auto const &part : vtx->particles_out()) {
		std::cout << part->pid() << " " ;
		if( part->pid() == 2212 || part->pid() == 211 ){
		  final_T = part->momentum().e()-part->momentum().m() ;
		  charge += 1;
		}
		if ( part->pid() == -211 ) charge -= 1;
	      }
	      std::cout << "\n";
	      std::cout << " Final charge " << charge << std::endl;
	      std::cout << " Carried relative energy %:" << (initial_T-final_T)/initial_T*100 << std::endl;
	      ++has_chex;
	      if( pid_neutral == 2112 ) histograms_studies["CExEnergyTransferNeutron"] ->Fill( (initial_T-final_T)/initial_T*100  ) ;
	      if ( pid_neutral == 111 ) histograms_studies["CExEnergyTransferPi0"] ->Fill( (initial_T-final_T)/initial_T*100  ) ;
	    }
	  }
	}

	for (auto const &part : vtx->particles_out()) {
	  //if( part->pid() > 90 && part->pid() < 93 ) std::cout << "****************************************************KS PYTHIA"<<std::endl;
	  //std::cout << " OUT particle " << part->id()  << " Pdg " << part->pid() << " status " << part->status() <<  " energy " << part->momentum().e()<<std::endl;
	}
      }

      std::cout << " Event summary. NProtons before (after) FSI: " << np_vertex << " (" << np_postFSI << ")"<<std::endl;
      std::cout << " Event summary. NNeutrons before (after) FSI: " << nn_vertex << " (" << nn_postFSI << ")"<<std::endl;
      std::cout << " Event summary. NPiM before (after) FSI: " << npim_vertex << " (" << npim_postFSI << ")"<<std::endl;
      std::cout << " Event summary. NPiP before (after) FSI: " << npip_vertex << " (" << npip_postFSI << ")"<<std::endl;

      //double evw = FATXAcc->process(evt);
    }
  }

  std::cout << " In Peak " << has_chex/in_peak*100 << std::endl;
  std::cout << " Total processed: " << nprocessed << std::endl;
  // Write the histogram to the file
  for ( auto it = histograms.begin(); it != histograms.end(); it++) { 
    // Normalize by bin witdh and xsection 
    NormalizeHist(it->second, xsec / nprocessed );
    // Save
    (it->second)->Write();
  }

  for ( auto it = histograms_studies.begin(); it != histograms_studies.end(); it++) { 
    // Normalize by number of entries
    (it->second)->Scale( 1. / in_peak );
    // Save
    (it->second)->Write();
  }

  for ( auto it = histograms2D.begin(); it != histograms2D.end(); it++) { 
    // Normalize by total number of events 
    NormalizeHist(it->second, xsec / nprocessed );
    // Save
    (it->second)->Write();
  }
 
  // Close the file
  outfile->Close();

  // Clean up
  delete outfile;

  //  double fatx = FATXAcc->fatx(cm2ten38_PerNucleon);
  //  double fatx = FATXAcc->fatx(); // in pb/Atom
  //  double sumw = FATXAcc->sumweights();
  //  size_t nevents = FATXAcc->events();

}
