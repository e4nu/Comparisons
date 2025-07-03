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
  std::string topology = "1p1pim", output_name = "myradevents", model_name = "MC_Name" ;
  
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
  TH1D * xsec_El = new TH1D( (model_name+"_1Dxsec_El").c_str(), "",GetBinning("Efl", Beam_E, topology).size()-1, &GetBinning("Efl", Beam_E, topology)[0] );
  TH1D * xsec_pl = new TH1D( (model_name+"_1Dxsec_pl").c_str(), "", GetBinning("pfl", Beam_E, topology).size()-1, &GetBinning("pfl", Beam_E, topology)[0] );
  TH1D * xsec_thetal = new TH1D( (model_name+"_1Dxsec_thetal").c_str(), "", GetBinning("pfl_theta", Beam_E, topology).size()-1, &GetBinning("pfl_theta", Beam_E, topology)[0] );
  TH1D * xsec_pp = new TH1D( (model_name+"_1Dxsec_pp").c_str(), "", GetBinning("proton_mom", Beam_E, topology).size()-1, &GetBinning("proton_mom", Beam_E, topology)[0] );
  TH1D * xsec_thetap = new TH1D( (model_name+"_1Dxsec_thetap").c_str(), "", GetBinning("proton_theta", Beam_E, topology).size()-1, &GetBinning("proton_theta", Beam_E, topology)[0] );
  TH1D * xsec_ppi, * xsec_thetapi ;
  if( topology == "1p1pim") {
    xsec_ppi = new TH1D( (model_name+"_1Dxsec_ppi").c_str(), "", GetBinning("pim_mom", Beam_E, topology).size()-1, &GetBinning("pim_mom", Beam_E, topology)[0] );
    xsec_thetapi = new TH1D( (model_name+"_1Dxsec_thetapi").c_str(), "", GetBinning("pim_theta", Beam_E, topology).size()-1, &GetBinning("pim_theta", Beam_E, topology)[0] );
  } else if( topology == "1p1pip") {
    xsec_ppi = new TH1D( (model_name+"_1Dxsec_ppi").c_str(), "", GetBinning("pip_mom", Beam_E, topology).size()-1, &GetBinning("pip_mom", Beam_E, topology)[0] );
    xsec_thetapi = new TH1D( (model_name+"_1Dxsec_thetapi").c_str(), "", GetBinning("pip_theta", Beam_E, topology).size()-1, &GetBinning("pip_theta", Beam_E, topology)[0] );
  }

  TH1D * xsec_ECal = new TH1D( (model_name+"_1Dxsec_ECal").c_str(), "", GetBinning("ECal", Beam_E, topology).size()-1, &GetBinning("ECal", Beam_E, topology)[0]);
  TH1D * xsec_W = new TH1D( (model_name+"_1Dxsec_W").c_str(), "", GetBinning("RecoW", Beam_E, topology).size()-1, &GetBinning("RecoW", Beam_E, topology)[0]);
  TH1D * xsec_Q2 = new TH1D( (model_name+"_1Dxsec_Q2").c_str(), "", GetBinning("RecoQ2", Beam_E, topology).size()-1, &GetBinning("RecoQ2", Beam_E, topology)[0] );
  TH1D * xsec_AlphaT = new TH1D( (model_name+"_1Dxsec_DeltaPT").c_str(), "", GetBinning("HadDeltaPT", Beam_E, topology).size()-1, &GetBinning("HadDeltaPT", Beam_E, topology)[0] );
  TH1D * xsec_DeltaPT = new TH1D( (model_name+"_1Dxsec_AlphaT").c_str(), "", GetBinning("HadAlphaT", Beam_E, topology).size()-1, &GetBinning("HadAlphaT", Beam_E, topology)[0] );

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
    rdr->read_event(evt);
    if (rdr->failed()) {
      break;
    }

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
      else if ( particle_pdg == -11 ) ++true_pip;
      else if ( particle_pdg == 22 ) ++true_gamma;
      
      // Smear Hadrons 
      HepMC3::GenParticlePtr hadron_smeared = SmearParticles( hadrons[i], Ev ) ;

      // Apply momentum and angle cuts for hadrons 
      if( hadron_smeared->momentum().theta() * 180 / TMath::Pi() < GetParticleMinTheta( particle_pdg, hadron_smeared->momentum().p3mod(), Ev ) ) continue;
      if( hadron_smeared->momentum().p3mod() < GetMinMomentumCut( particle_pdg, Ev ) ) continue ;
      
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
    xsec_El->Fill(El);
    xsec_pl->Fill(Pl);
    xsec_thetal->Fill(fslep_reco->momentum().theta()*180/TMath::Pi());
    xsec_pp->Fill(proton->momentum().p3mod());
    xsec_thetap->Fill(proton->momentum().theta()*180/TMath::Pi());
    xsec_ppi->Fill(pion->momentum().p3mod());
    xsec_thetapi->Fill(pion->momentum().theta()*180/TMath::Pi());
    xsec_W->Fill(GetW(fslep_reco,beampt));
    xsec_Q2->Fill(GetQ2(fslep_reco,beampt));
    xsec_ECal->Fill(GetECal( fslep_reco, hadrons_reco, tgtpt->pid()));
    //double evw = FATXAcc->process(evt);
    ++nprocessed;
  }

  // Write the histogram to the file
  xsec_El->Write();
  xsec_pl->Write();  
  xsec_thetal->Write();
  xsec_pp->Write();  
  xsec_thetap->Write();
  xsec_ppi->Write();  
  xsec_thetapi->Write();
  xsec_W->Write();
  xsec_Q2->Write();
  xsec_ECal->Write();

  // Close the file
  outfile->Close();

  // Clean up
  delete outfile;

  //  double fatx = FATXAcc->fatx(cm2ten38_PerNucleon);
  //  double fatx = FATXAcc->fatx(); // in pb/Atom
  //  double sumw = FATXAcc->sumweights();
  //  size_t nevents = FATXAcc->events();

}
