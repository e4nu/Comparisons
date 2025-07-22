// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TH1D.h>
#include <TMath.h>
#include "NuHepMC/HepMC3Features.hxx"
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
#include <TRandom3.h>
#include <TVector3.h>

using namespace NuHepMC;
static const unsigned int kPdgH   = 1000010010; 
static const unsigned int kPdgHe3 = 1000020030; 
static const unsigned int kPdgHe4 = 1000020040; 
static const unsigned int kPdgC12 = 1000060120 ; 
static const unsigned int kPdgFe56 = 1000260560 ;
static const unsigned int kPdgD = 1000010020 ; 
static const unsigned int kPdgO16 = 1000080160 ; 
static const unsigned int kPdgFreeP = 1000010010 ;
static const unsigned int kPdgFreeN = 1000000010 ;

// Binding energy in GeV
static const double kBEH = 0.008481 ; 
static const double kBEHe3 = 0.0077 ;
static const double kBEHe4 = 0.0283 ; 
static const double kBED2 = 0.00222 ; 
static const double kBEC12 = 0.09215 ; 
static const double kBEFe56 = 0.49226 ; 
static const double kBEB = 0.0762 ;
static const double kBEMn = 0.4820764 ; 

// ECalOffset
static const double kECalOffsetHe3 = 0.004 ; 
static const double kECalOffsetHe4 = 0.005 ; 
static const double kECalOffsetC12 = 0.005 ; 
static const double kECalOffsetFe56 = 0.011 ; 

std::string PartToStr(HepMC3::ConstGenParticlePtr pt) {
  if (!pt) {
    return "PARTICLE-NOTFOUND";
  }
  std::stringstream ss;

  auto mom = pt->momentum() ;

  ss << "{ id: " << pt->id() << ", pid: " << pt->pid() 
     << ", p: ( " << mom.x() << ", " << mom.y() << ", "
     << mom.z() << ", E: " << mom.e() << ") GeV }";

  return ss.str();
}


std::string GetArg(std::string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  if (strlen(&argv[2][0]) ) {
	    strcpy(argument,&argv[2][0]);
	  }
	}
      }
      argc--;
      argv++;

    }

  std::string value = std::string(argument);
  delete [] argument;
  return value ;
}


bool ExistArg(std::string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  return true ;
	}
      }
      argc--;
      argv++;

    }
  delete [] argument ;
  return false;
}

double GetParticleMass( const int pdg ) {
  static const double kProtonMass   = 0.9382720813 ; 
  static const double kNeutronMass  = 0.939565 ;
  static const double kPiMass = 0.139570 ;
  static const double kElectronMass = 0.000510998 ; 
  static const double kNucleonMass  = (kProtonMass+kNeutronMass)/2.;  //GeV
  
  if( pdg == 2212 ) return kProtonMass ;
  else if ( pdg == -211 || pdg == 211 ) return kPiMass;
  else if ( pdg == 11 ) return kElectronMass ; 
  else return 0;
}

double GetParticleResolution( const int pdg, const double Beam_E ) {
  double resolution = 0 ; 
  if ( pdg == 2212 ) resolution = 0.01 ; 
  else if ( pdg == 11 ) resolution = 0.005;
  else if ( pdg == 211 || pdg == -211 || pdg == 22 ) resolution = 0.007 ;
  if ( Beam_E < 2 ) resolution *= 3; // Is it only this value or beam_E>1.1 GeV ? 
  return resolution ; 
}

HepMC3::GenParticlePtr SmearParticles( auto particle, const double Beam_E ){
  HepMC3::GenParticlePtr smeared_particle = std::make_shared<HepMC3::GenParticle>(particle->data()); ;
  return smeared_particle;
  double res = GetParticleResolution( particle->pid(), Beam_E );
  double p = particle->momentum().p3mod();
  double M = GetParticleMass( particle->pid() ) ;
  double SmearedP = gRandom->Gaus(p,res*p);
  double SmearedE = sqrt( pow( SmearedP,2 ) + pow( M,2 ) ) ; 
  const HepMC3::FourVector smeared_mom ( SmearedP/p * particle->momentum().px(), SmearedP/p * particle->momentum().py(), SmearedP/p * particle->momentum().pz(), SmearedE ) ; 
  smeared_particle->set_momentum(smeared_mom);

  return smeared_particle;
}

double GetQ2Cut( double beam_E ) {
  if( beam_E > 1 && beam_E < 2 ) return 0.1 ; 
  else if( beam_E > 2 && beam_E < 4) return 0.4 ; 
  else if( beam_E > 4  ) return 0.8 ; 
  return 0;
}

double GetParticleMinTheta( const int pdg, const double out_mom, const double Beam_E ) { 
  double min_theta = 0 ;
  if ( pdg == 11 ) { 
    if( Beam_E < 2 ) min_theta = 17 + 7.0 / out_mom ; // deg 
    else if ( Beam_E > 2 && Beam_E < 4 ) min_theta = 16 + 10.5 / out_mom ; // deg
    else if ( Beam_E > 4 ) min_theta = 13.5 + 15 / out_mom ; // deg 
  } 
  else if ( pdg == 2212 ) min_theta = 10 ; 
  else if ( pdg == -211 ) {
    if( Beam_E < 2 ) min_theta = 17 + 4./TMath::Power(out_mom,1.);
    else if ( Beam_E > 2 ) { 
      if( out_mom<0.35) min_theta = 25.+7./TMath::Power(out_mom,1.);
      else min_theta = 16.+10./TMath::Power(out_mom,1.);
    }
  }
  else if ( pdg == 211  ) min_theta = 10 ; 
  else if ( pdg == 22   ) min_theta = 8 ; 

  return min_theta;
}

double GetMinMomentumCut( const int particle_pdg, const double Beam_E ) { 
  double min_p = 0 ;
  if( particle_pdg == 11 ) {
    if( Beam_E < 2 ) min_p = 0.4 ; 
    else if ( Beam_E > 2 && Beam_E < 4 ) min_p = 0.55 ;
    else if ( Beam_E > 4 ) min_p = 1.3 ; 
  } else if ( particle_pdg == 2212 || particle_pdg == 22 ) { 
    min_p = 0.3 ; 
  } else if ( particle_pdg == -211 ) {
    min_p = 0.15 ; 
  } else if ( particle_pdg == 211 ) { 
    min_p = 0.25 ; 
  }

  return min_p ; 
}

#include <cmath>

double AngleBetween(const HepMC3::FourVector& v1, const HepMC3::FourVector& v2) {
  // Compute dot product of spatial components
  double dot = v1.px()*v2.px() + v1.py()*v2.py() + v1.pz()*v2.pz();

  // Compute magnitudes
  double mag1 = std::sqrt(v1.px()*v1.px() + v1.py()*v1.py() + v1.pz()*v1.pz());
  double mag2 = std::sqrt(v2.px()*v2.px() + v2.py()*v2.py() + v2.pz()*v2.pz());

  // Avoid division by zero
  if (mag1 == 0 || mag2 == 0) return 0;

  double cos_theta = dot / (mag1 * mag2);
  // Clamp to avoid domain errors from precision issues
  if (cos_theta > 1.0) cos_theta = 1.0;
  if (cos_theta < -1.0) cos_theta = -1.0;

  return std::acos(cos_theta);  // in radians
}

bool ApplyPhotRadCut( const HepMC3::FourVector electron, const HepMC3::FourVector gamma ) {
  double neut_phi_mod = gamma.phi()*TMath::RadToDeg() + 30; //Add 30 degree
  if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

  double el_phi_mod = electron.phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
  if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree
  
  const double kPhotonRadCut = 40 ; 
  const double kPhotonEPhiDiffCut = 30 ;

  if( AngleBetween(gamma,electron)*TMath::RadToDeg() < kPhotonRadCut && fabs(neut_phi_mod-el_phi_mod) < kPhotonEPhiDiffCut ) return true ; 
  return false ;
}

double GetECalOffset( const unsigned int target_pdg ) {
  double ECalOffset = 0 ; 
  if ( target_pdg == kPdgHe3 ) ECalOffset = kECalOffsetHe3; 
  else if ( target_pdg == kPdgHe4 || target_pdg == kPdgC12 ) ECalOffset = kECalOffsetC12 ;
  else if ( target_pdg == kPdgFe56 ) ECalOffset = kECalOffsetFe56 ;
  return ECalOffset ; 
}

double GetBindingEnergy( const unsigned int target_pdg ) {
  double binding_energy = 0. ; 
  if ( target_pdg == kPdgH ) binding_energy = kBEH;
  else if ( target_pdg == kPdgHe3 ) binding_energy = kBEHe3 - kBED2 + GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgHe4 ) binding_energy = kBEHe4 - kBEHe3 + GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgC12 ) binding_energy = kBEC12 - kBEB + GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgFe56 ) binding_energy = kBEFe56 - kBEMn + GetECalOffset( target_pdg ) ;
  //else if ( target_pdg == kPdgCH2 ) binding_energy = kBEC12 - kBEB ;
  return binding_energy;
}

double GetEnergyTransfer( HepMC3::ConstGenParticlePtr fsl, HepMC3::ConstGenParticlePtr isl ) { 
  return isl->momentum().e() - fsl->momentum().e();
}

double GetW(HepMC3::ConstGenParticlePtr fsl, HepMC3::ConstGenParticlePtr isl){
  double mp = 0.9389 ; // GeV
  double nu = isl->momentum().e() - fsl->momentum().e();
  auto q3 = isl->momentum() - fsl->momentum(); 
  double W2 = std::pow(mp + nu, 2) - pow(q3.p3mod(),2);
  //if (W2 < 0) return 0;
  return TMath::Sqrt(W2);
  return std::pow(mp + nu, 2) ;
}

double GetQ2(HepMC3::ConstGenParticlePtr fsl, HepMC3::ConstGenParticlePtr isl){
  return -(fsl->momentum() - isl->momentum()).m2();
}

double GetECal( HepMC3::ConstGenParticlePtr fsl, std::vector<HepMC3::ConstGenParticlePtr> hadrons, const int tgt) {
  double ECal = fsl->momentum().e(); // Add energy of outgoing lepton

  for( unsigned int i = 0 ; i < hadrons.size() ; ++i ) {
    // Calculate ECal for visible particles
    if( hadrons[i]->pid() == 11 ) continue ; 
    ECal += hadrons[i]->momentum().e() ; // Add Kinetic energy of hadrons   
    if( hadrons[i]->pid() == 2212 ) ECal += GetBindingEnergy(tgt) - GetParticleMass(hadrons[i]->pid()); // Correct for proton binding energy 
  }

  return ECal;
}

TVector3 GetPT(const HepMC3::FourVector& p)
{
  TVector3 beam_dir(0, 0, 1);
  TVector3 p_vect(p.px(), p.py(), p.pz());
  double vect_parallel = p_vect.Dot(beam_dir);

  // Calculate transverse vector:
  TVector3 vect_T = p_vect - (vect_parallel * beam_dir);
  return vect_T;
}

TVector3 DeltaPT(HepMC3::ConstGenParticlePtr fsl, const std::vector<HepMC3::ConstGenParticlePtr>& hadrons)
{
  TVector3 P1_T = GetPT(fsl->momentum());

  HepMC3::FourVector tot_hadrons(0,0,0,0);
  for (const auto& had : hadrons) {
    tot_hadrons += had->momentum();
  }

  TVector3 P2_T = GetPT(tot_hadrons);

  return P1_T + P2_T;
}

double DeltaAlphaT(HepMC3::ConstGenParticlePtr fsl, const std::vector<HepMC3::ConstGenParticlePtr>& hadrons) {
  
  TVector3 P1T_dir = GetPT(fsl->momentum()).Unit();
  TVector3 DeltaPT_dir = DeltaPT(fsl, hadrons).Unit();

  return acos(-P1T_dir.Dot(DeltaPT_dir)) * TMath::RadToDeg();
}


double HadSystemMass( const std::vector<HepMC3::ConstGenParticlePtr>& hadrons ) { 
  HepMC3::FourVector tot_hadrons(0,0,0,0);
  for (const auto& had : hadrons) {
    tot_hadrons += had->momentum();
  }
  
  return tot_hadrons.p3mod();
}

std::vector<double> GetUniformBinning(unsigned int nbins, double min, double max) {
  std::vector<double> binning;
  double step = (max - min) / nbins;
  for (unsigned int i = 0; i < nbins + 1; ++i)
  {
    binning.push_back(min + i * step);
  }
  return binning;
}

std::vector<double> GetECalBinning(unsigned int nbins_tail, unsigned int nbins_peak, double min, double max, double EBeam)
{
  std::vector<double> binning;
  double temp_min = min;
  double temp_max = max;
  if (EBeam - min < max - EBeam ) temp_max = EBeam * (1 + 0.05);
  else temp_max = EBeam * (1 - 0.05);

  double step = (temp_max - temp_min) / nbins_tail;
  for (unsigned int i = 0; i < nbins_tail + 1; ++i)
  {
    binning.push_back(temp_min + i * step);
  }

  step = (max - temp_max) / nbins_peak;
  for (unsigned int i = 1; i < nbins_peak + 1; ++i)
  {
    binning.push_back(temp_max + i * step);
  }
  return binning;
}

std::vector<double> GetBinning(std::string observable, double EBeam, std::string analysis_key)
{
  std::vector<double> binning;

  if (observable == "ECal") {
    if (EBeam < 2) binning = GetECalBinning(13, 15, 0.6, EBeam + 0.15, EBeam);
    else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(13, 15, 0.6, EBeam + 0.15, EBeam);
    else if (EBeam > 4) binning = GetECalBinning(13, 15, 1.8, EBeam + 0.15, EBeam);
  } else if (observable == "Efl") {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.35, 0.9);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.5, 1.7);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 1.2, 3.8);
  } else if (observable == "DiffECal") {
    if (EBeam < 2)
    binning = GetUniformBinning(25, -0.6, 0.2);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, -0.6, 0.2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, -0.6, 0.2);
  }
  else if (observable == "pfl_theta")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 20, 50);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 20, 50);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 15, 50);
  }
  else if (observable == "pfl_phi")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  }
  else if (observable == "pfl")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.35, 0.9);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(28, 0.5, 1.9);
    else if (EBeam > 4)
    binning = GetUniformBinning(20, 1.1, 3.8);
  }
  else if (observable == "pfl_T")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.2, 0.6);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.3, 0.9);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0.5, 1.2);
  }
  else if (observable == "proton_mom")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.2, 1.1);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.2, 2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0.2, 3);
  }
  else if (observable == "proton_theta")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(23, 5, 140);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(23, 5, 140);
    else if (EBeam > 4)
    binning = GetUniformBinning(23, 5, 140);
  }
  else if (observable == "proton_phi")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(35, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(35, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(35, 0, 180);
  }
  else if (observable == "pim_mom")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.1, 0.6);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0., 1.6);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0., 2);
  }
  else if (observable == "pim_theta")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(23, 5, 140);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(23, 5, 140);
    else if (EBeam > 4)
    binning = GetUniformBinning(23, 5, 140);
  }
  else if (observable == "pip_mom")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.1, 3);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.1, 3);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0.1, 3);
  }
  else if (observable == "pip_theta")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 5, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 5, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 5, 180);
  }
  else if (observable == "RecoW")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 1, 1.5);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 1, 2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 1, 2.5);
  }
  else if (observable == "RecoXBJK")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 0.9);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 0.9);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 1);
  }
  else if (observable == "RecoQ2")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0.15, 0.45);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.3, 1.5);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0.9, 3);
  }
  else if (observable == "RecoQELEnu")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(20, 0., 0.8);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0.3, EBeam + 0.2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0.9, EBeam + 0.2);
  }
  else if (observable == "Recoq3")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, EBeam + 0.2);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, EBeam + 0.2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, EBeam + 0.2);
  }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 1);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 1);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 1);
  }
  else if (observable == "HadDeltaPTx")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, -0.6, 0.6);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, -1, 1);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, -1, 1);
  }
  else if (observable == "HadDeltaPTy")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, -0.6, 0.6);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, -1, 1);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, -1, 1);
  }
  else if (observable == "HadDeltaPhiT" || observable == "DeltaPhiT")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 80);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 80);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 80);
  }
  else if (observable == "AlphaT")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(20, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  }
  else if (observable == "HadAlphaT")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(15, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  }
  else if (observable == "RecoEnergyTransfer")
  {
    if (EBeam < 2) binning = GetUniformBinning(25, 0., 0.8);
    else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(25, 0., 2);
    else if (EBeam > 4) binning = GetUniformBinning(25, 0, 4);
  }
  else if (observable == "HadSystemMass")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 1, 1.6);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 1, 2);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 1, 2.7);
  }
  else if (observable == "MissingEnergy")
  {
    if (EBeam < 2) binning = GetECalBinning(20,60, -0.15, 1, 0.1);
    else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(20, 60, -0.15, 2, 0.1);
    else if (EBeam > 4) binning = GetECalBinning(20, 60, -0.15, 3, 0.1);
  }
  else if (observable == "MissingTransMomentum")
  {
    if (EBeam < 2)
    binning = GetECalBinning(25, 15, 0.3, 1.1, 0.9);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetECalBinning(25, 15, -0.7, 1.2, 0.9);
    else if (EBeam > 4)
    binning = GetECalBinning(25, 15, -2.5, 1.2, 0.9);
  }
  else if (observable == "CorrMissingEnergy" || observable == "CorrMissingEnergy1" || observable == "CorrMissingEnergy2" || observable == "CorrMissingEnergy3")
  {
    if (EBeam < 2)
    binning = GetECalBinning(25, 15, 0.3, 1.1, 0.9);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetECalBinning(25, 15, -0.7, 1.2, 0.9);
    else if (EBeam > 4)
    binning = GetECalBinning(25, 15, -2.5, 1.2, 0.9);
  }
  else if (observable == "MissingAngle")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 15, 180);
  }
  else if (observable == "MissingMomentum")
  {
    if (EBeam < 2) binning = GetECalBinning(20,60, -0.15, 1.5, 0.15);
    else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(20, 60, -0.15, 2.2, 0.15);
    else if (EBeam > 4) binning = GetECalBinning(20, 60, -0.15, 3.5, 0.15);
  }
  else if (observable == "InferedNucleonMom")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(30, 0, 0.8);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(30, 0, 1);
    else if (EBeam > 4)
    binning = GetUniformBinning(30, 0, 1);
  }
  else if (observable == "HadronsAngle")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(30, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(30, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(30, 0, 180);
  }
  else if (observable == "AdlerAngleThetaP")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  }
  else if (observable == "AdlerAnglePhiP")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  } else if (observable == "AdlerAngleThetaPi")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 180);
  } else if (observable == "AdlerAnglePhiPi")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(30, 20, 180);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 20, 180);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 20, 180);
  } else if (observable == "Angleqvshad")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 120);
    else if (EBeam > 2 && EBeam < 4)
    binning = GetUniformBinning(25, 0, 120);
    else if (EBeam > 4)
    binning = GetUniformBinning(25, 0, 60);
  }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
  {
    if (EBeam < 2)
    binning = GetUniformBinning(25, 0, 0.7);
  }
  else if (observable == "TrueNProtons" || observable == "TrueNNeutrons" || observable == "TrueNPiP" || observable == "TrueNPiM" || observable == "TrueNPi0" || observable == "TrueNCh" )
  {
    binning = GetUniformBinning(5, -0.5, 5.5);
  }

  if (analysis_key == "1pip")
  {
    if (observable == "ECal")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 20, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 20, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 20, 180);
    }
    else if (observable == "AdlerAnglePhiPi")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(30, 20, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 20, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 20, 180);
    }
    else if (observable == "Angleqvshad")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 70);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 50);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 30);
    }
    else if (observable == "RecoEvPion")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 2);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0.5, 3.5);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 2, 6);
    }
    else if (observable == "RecoWPion")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0.5, 2);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 1, 2);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0.5, 3.5);
    }
    else if (observable == "ElectronPT")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0.2, 0.7);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "PionPT")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 0.5);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "Angleqvshad")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 120);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 120);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 60);
    }
    else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam < 2) binning = GetUniformBinning(25, 0, 0.7);
    }
  }

  if (analysis_key == "1p1pip")
  {
    if (observable == "ECal")
    {
      if (EBeam < 2) binning = GetECalBinning(10, 10, 0.7, EBeam + 0.1, EBeam);
      else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(20, 10, 0.8, EBeam + 0.2, EBeam);
      else if (EBeam > 4) binning = GetECalBinning(20, 10, 1.2, EBeam + 0.2, EBeam);
    }
    else if (observable == "proton_mom")
    {
      if (EBeam < 2) binning = GetUniformBinning(20, 0.2, 1);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(30, 0.2, 2);
      else if (EBeam > 4) binning = GetUniformBinning(30, 0.2, 3);
    }
    else if (observable == "proton_theta")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(23, 0, 140);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(23, 0, 140);
      else if (EBeam > 4)
      binning = GetUniformBinning(23, 0, 140);
    }
    else if (observable == "pfl_theta")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 24, 48);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 20, 50);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 15, 50);
    }
    else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam < 2) binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 4) binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "HadDeltaPhiT" || observable == "DeltaPhiT")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(15, 0, 80);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 80);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 80);
    }
    else if (observable == "HadSystemMass")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(15, 1, 1.6);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 1, 2);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 1, 2.7);
    }
    else if (observable == "HadronsAngle")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(20, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(30, 0, 180);
    }
    else if (observable == "pip_mom")
    {
      if (EBeam < 2) binning = GetUniformBinning(25, 0.1, 0.6);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(25, 0.1, 1.2);
      else if (EBeam > 4) binning = GetUniformBinning(25, 0.1, 2.2);
    }
    else if (observable == "pip_theta")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(20, 10, 130);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(20, 10, 130);
      else if (EBeam > 4)
      binning = GetUniformBinning(20, 10, 130);
    }
    else if (observable == "MissingEnergy")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0.5, 1);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "CorrMissingEnergy"||observable == "CorrMissingEnergy1"||observable == "CorrMissingEnergy2"||observable == "CorrMissingEnergy3")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0.5, 1);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 1);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "MissingAngle")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(30, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(30, 0, 180);
    }
    else if (observable == "MissingMomentum")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(20, 0, 2);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 2.2);
      else if (EBeam > 4)
      binning = GetUniformBinning(30, 0, 4);
    }
    else if (observable == "InferedNucleonMom")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(30, 0, 0.8);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 1);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 1);
    }
    else if (observable == "Angleqvshad")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(30, 0, 120);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 120);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 120);
    }
    else if (observable == "HadronsAngle")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(20, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(30, 0, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(30, 0, 180);
    }
    else if (observable == "AdlerAngleThetaPi")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(10, 0, 180);
    }
    else if (observable == "AdlerAnglePhiPi")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4)
      binning = GetUniformBinning(10, 0, 180);
    }
    else if (observable == "Angleqvshad")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(25, 0, 120);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 0, 120);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 0, 60);
    }
    else if (observable == "RecoW")
    {
      if (EBeam < 2)
      binning = GetUniformBinning(15, 1.1, 1.5);
      else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(25, 1, 2);
      else if (EBeam > 4)
      binning = GetUniformBinning(25, 1., 2.5);
    }
  }

  if (binning.size() == 0)
  {
    std::cout << " ERROR: Binning for " << observable << " is null" << std::endl;
  }
  return binning;
}

void NormalizeHist(TH1D *h, double normalization_factor)
{
  // Data normalization
  h->Scale(normalization_factor);
  // h->Sumw2();///kFALSE);
  double NBins = h->GetNbinsX();

  for (int i = 1; i <= NBins; i++)
    {
      double content = h->GetBinContent(i);
      double error = h->GetBinError(i);
      double width = h->GetBinWidth(i);
      double newcontent = content / width;
      double newerror = error / width;
      h->SetBinContent(i, newcontent);
      h->SetBinError(i, newerror);
    }
}

void NormalizeHist(TH2D *h, double normalization_factor)
{
  // Data normalization
  h->Scale(normalization_factor);

  int NbinsX = h->GetNbinsX();
  int NbinsY = h->GetNbinsY();

  for (int i = 1; i <= NbinsX; i++) // Loop over X-axis bins
    {
      for (int j = 1; j <= NbinsY; j++) // Loop over Y-axis bins
	{
	  double content = h->GetBinContent(i, j);
	  double error = h->GetBinError(i, j);
	  double width_x = h->GetXaxis()->GetBinWidth(i);
	  double width_y = h->GetYaxis()->GetBinWidth(j);

	  double newcontent = content / width_x ;/// width_y;
	  double newerror = error / width_x ;/// width_y;

	  h->SetBinContent(i, j, newcontent);
	  h->SetBinError(i, j, newerror);
	}
    }
}
