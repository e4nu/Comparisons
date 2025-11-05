// Leave this at the top to enable features detected at build time in headers in
// HepMC3
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include "TStyle.h"
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

std::string GetAxisLabel(std::string observable, unsigned int id_axis, std::string units)
{
  std::string x_axis, y_axis, unit = "#mub";
  if ( units == "nb" ) unit = "nb";

  if (observable == "ECal")
    {
      x_axis = "E_{Cal} [GeV]";
      y_axis = "d#sigma/dE_{Cal} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "Efl")
    {
      x_axis = "E_{e'} [GeV]";
      y_axis = "d#sigma/dE_{e'} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "pfl_theta")
    {
      x_axis = "#theta_{e'} [deg]";
      y_axis = "d#sigma/d#theta_{e'} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "pfl_phi")
    {
      x_axis = "#phi_{e'} [deg]";
      y_axis = "d#sigma/d#phi_{e'} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "pfl")
    {
      x_axis = "p_{e'} [GeV/c]";
      y_axis = "d#sigma/dp_{e'} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "pfl_T")
    {
      x_axis = "p_{e'}^{T} [GeV/c]";
      y_axis = "d#sigma/dp_{e'}^{T} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "proton_mom")
    {
      x_axis = "p_{p} [GeV/c]";
      y_axis = "d#sigma/dp_{p} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "proton_theta")
    {
      x_axis = "#theta_{p} [deg]";
      y_axis = "d#sigma/d#theta_{p} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "proton_phi")
    {
      x_axis = "E_{Cal} [GeV]";
      y_axis = "d#sigma/dE_{Cal} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "pim_mom")
    {
      x_axis = "p_{#pi^{-}} [GeV/c]";
      y_axis = "d#sigma/dp_{#pi^{-}} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "pim_theta")
    {
      x_axis = "#theta_{#pi^{-}} [deg]";
      y_axis = "d#sigma/d#theta_{#pi^{-}} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "pip_mom")
    {
      x_axis = "p_{#pi^{+}} [GeV/c]";
      y_axis = "d#sigma/dp_{#pi^{+}} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "pip_theta")
    {
      x_axis = "#theta_{#pi^{+}} [deg]";
      y_axis = "d#sigma/d#theta_{#pi^{+}} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "RecoW")
    {
      x_axis = "W [GeV]";
      y_axis = "d#sigma/dW #left["+unit+" GeV^{-1}#right#right]";
    }
  else if (observable == "RecoQELEnu")
    {
      x_axis = "E^{QE} [GeV]";
      y_axis = "d#sigma/dE^{QE} #left["+unit+" GeV^{-1}#right#right]";
    }
  else if (observable == "RecoXBJK")
    {
      x_axis = "x_{BJK} [GeV]";
      y_axis = "d#sigma/dx_{BJK} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "RecoQ2")
    {
      x_axis = "Q^{2} [GeV^{2}]";
      y_axis = "d#sigma/dQ^{2} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "Recoq3")
    {
      x_axis = "q_{3} [GeV]";
      y_axis = "d#sigma/dq_{3} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "DeltaPT")
    {
      x_axis = "#deltap_{T} [GeV]";
      y_axis = "d#sigma/d#deltap_{T} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPT")
    {
      x_axis = "#deltap_{T} [GeV]";
      y_axis = "d#sigma/d#deltap_{T} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPTx")
    {
      x_axis = "#deltap_{Tx} [GeV]";
      y_axis = "d#sigma/d#deltap_{Tx} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "HadDeltaPTy")
    {
      x_axis = "#deltap_{Ty} [GeV]";
      y_axis = "d#sigma/d#deltap_{Ty} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "InferedNucleonMom")
    {
      x_axis = "p_{N,proxy} [GeV]";
      y_axis = "d#sigma/dp_{N,proxy} #left["+unit+" #left(GeV/c#right)^{-1}#right]";
    }
  else if (observable == "DeltaPhiT")
    {
      x_axis = "#delta#phi_{T} [deg]";
      y_axis = "d#sigma/d#delta#phi_{T} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "HadDeltaPhiT")
    {
      x_axis = "#delta#phi_{T}^{had} [deg]";
      y_axis = "d#sigma/d#delta#phi_{T}^{had} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "AlphaT")
    {
      x_axis = "#alpha_{T} [deg]";
      y_axis = "d#sigma/d#alpha_{T} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "HadAlphaT")
    {
      x_axis = "#delta#alpha_{T} [deg]";
      y_axis = "d#sigma/d#delta#alpha_{T} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "RecoEnergyTransfer")
    {
      x_axis = "#omega [GeV]";
      y_axis = "d#sigma/d#omega #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "HadSystemMass")
    {
      x_axis = "M_{had}[GeV]";
      y_axis = "d#sigma/dM_{had} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "MissingEnergy")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dE_{miss} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "MissingTransMomentum")
    {
      x_axis = "E_{miss}[GeV]";
      y_axis = "d#sigma/dp_{miss}^{T} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy")
    {
      x_axis = "E_{miss}^{corr}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy1")
    {
      x_axis = "E_{miss}^{corr}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy2")
    {
      x_axis = "E_{miss}^{corr}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "CorrMissingEnergy3")
    {
      x_axis = "E_{miss}^{corr}[GeV]";
      y_axis = "d#sigma/dE_{miss}^{corr} #left["+unit+" GeV^{-1}#right]";
    }
  else if (observable == "MissingAngle")
    {
      x_axis = "#theta_{miss}[deg]";
      y_axis = "d#sigma/d#theta_{miss} #left["+unit+" deg^{-1}#right]";
    }
  else if (observable == "MissingMomentum")
    {
      x_axis = "p_{miss}[GeV/c]";
      y_axis = "d#sigma/dp_{miss} #left["+unit+" (GeV/c)^{-1}#right]";
    }
  else if (observable == "HadronsAngle")
    {
      x_axis = "#theta_{had}[deg]";
      y_axis = "d#sigma/d#theta_{had} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAngleThetaP")
    {
      x_axis = "#theta_{p}^{*}[deg]";
      y_axis = "d#sigma/d#theta_{p}^{*} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAnglePhiP")
    {
      x_axis = "#phi_{p}^{*}[deg]";
      y_axis = "d#sigma/d#phi_{p}^{*} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAngleThetaPi")
    {
      x_axis = "#theta_{#pi}^{*}[deg]";
      y_axis = "d#sigma/d#theta_{#pi}^{*} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "AdlerAnglePhiPi")
    {
      x_axis = "#phi_{#pi}^{*}[deg]";
      y_axis = "d#sigma/d#phi_{#pi}^{*} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "Angleqvshad")
    {
      x_axis = "#theta_{#vec{q}#dot#vec{p}_{had}}[deg]";
      y_axis = "d#sigma/d#theta_{#vec{q}#dot#vec{p}_{had}} #left["+unit+" (deg)^{-1}#right]";
    }
  else if (observable == "RecoEvPion")
    {
      x_axis = "E_{rec} [GeV]";
      y_axis = "d#sigma/dE_{rec} #left["+unit+" GeV^{-1}#right#right]";
    }
  else if (observable == "RecoWPion")
    {
      x_axis = "W_{rec} [GeV]";
      y_axis = "d#sigma/dW_{rec} #left["+unit+" GeV^{-1}#right#right]";
    }
  else if (observable == "ElectronPT")
    {
      x_axis = "p_{e',T} [GeV]";
      y_axis = "d#sigma/dp_{e'T} #left["+unit+" GeV^{-1}#right#right]";
    }
  else if (observable == "PionPT")
    {
      x_axis = "p_{#pi,T} [GeV]";
      y_axis = "d#sigma/dp_{#pi,T} #left["+unit+" GeV^{-1}#right#right]";
    }

  if (id_axis == 0)
    return x_axis;
  return y_axis;
}

void StandardFormat(TH1D *prediction, std::string title, int color, int style, std::string observable, bool is_log, double y_max  )
{
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPaperSize(20, 26);
  gStyle->SetTitleFont(132, "pad");
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineStyleString(2, "[12 12]");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTextFont(132);

  prediction->SetLineColor(color);
  prediction->SetLineStyle(style);
  prediction->SetMarkerStyle(style);
  prediction->SetMarkerColor(color);
  prediction->SetLineWidth(2);

  prediction->SetTitle(title.c_str());
  // prediction -> SetTitleFont(13);
  prediction->GetXaxis()->SetTitle(GetAxisLabel(observable, 0, "nb").c_str());

  prediction->GetYaxis()->SetTitle(GetAxisLabel(observable, 1, "nb").c_str());

  prediction->GetXaxis()->CenterTitle();
  prediction->GetYaxis()->CenterTitle();

  prediction->BufferEmpty(-1);
  if (y_max == 0)
    {
      double max = -999;
      for (unsigned int k = 0; k < prediction->GetNbinsX(); ++k)
	{
	  if (prediction->GetBinContent(k) > max)
	    max = prediction->GetBinContent(k);
	}
      // y_max = (prediction -> GetMaximum()) * ( 1+0.2 );
      y_max = max * (1 + 0.2);
    }

  if( y_max == 0 ) y_max = 100; // for empty plots.

  int FontStyle = 132;
  prediction->GetXaxis()->SetTitleOffset(1.1);
  prediction->GetXaxis()->SetLabelSize(0.1);
  prediction->GetXaxis()->SetTitleSize(0.08);
  prediction->GetXaxis()->SetNdivisions(5,3,0);
  prediction->GetXaxis()->SetLabelFont(FontStyle);
  prediction->GetXaxis()->SetTitleFont(FontStyle);

  prediction->GetYaxis()->SetNdivisions(4,4,0);
  prediction->GetYaxis()->SetTitleOffset(1.4);
  prediction->GetYaxis()->SetLabelSize(0.1);
  prediction->GetYaxis()->SetTitleSize(0.09);
  prediction->GetYaxis()->SetLabelFont(43);
  prediction->GetYaxis()->SetLabelFont(FontStyle);
  prediction->GetYaxis()->SetTitleFont(FontStyle);
  if( is_log ) {
    prediction->GetYaxis()->SetRangeUser(1E-4, y_max);
  } else { prediction->GetYaxis()->SetRangeUser(0.1, y_max); }
  

  prediction->GetYaxis()->SetMaxDigits(3);
  prediction->SetTitleFont(FontStyle);

  return;
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


bool ExistArg(std::string op, int argc, char **argv)
{
  const int buf_size = 2048 * 128;
  char *argument = new char[buf_size];
  strcpy(argument, "");

  while (argc > 2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-')
	{

	  char op_cur[buf_size];
	  strcpy(op_cur, &argv[1][2]);

	  if (strcmp(op.c_str(), op_cur) == 0)
	    {
	      return true;
	    }
	}
      argc--;
      argv++;
    }
  delete[] argument;
  return false;
}

HepMC3::GenParticlePtr SmearParticles( auto particle, const double Beam_E ){
  HepMC3::GenParticlePtr smeared_particle = std::make_shared<HepMC3::GenParticle>(particle->data()); ;
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
    if ( min_theta < 15. ) min_theta = 15 ; 
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
  if (W2 < 0) return 0;
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
    if( had->pid() == 11 ) continue ; // Skip the electron is the final lepton. 
    tot_hadrons += had->momentum();
  }
  return tot_hadrons.m();
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
    if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(10, 10, 0.6, EBeam + 0.15, EBeam);
    else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(10, 10, 0.6, EBeam + 0.15, EBeam);
    else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(10, 10, 1.8, EBeam + 0.15, EBeam);
  } else if (observable == "RecoEvPion") {
    if (EBeam > 0 && EBeam < 2 )
      binning = GetUniformBinning(15, 0, 2);
    else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(15, 0.5, 3.5);
    else if (EBeam > 4 && EBeam < 5)
      binning = GetUniformBinning(15, 2, 6);
  } else if (observable == "Efl") {
    if (EBeam > 0 && EBeam < 2 )
      binning = GetUniformBinning(15, 0.35, 0.9);
    else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(15, 0.5, 1.7);
    else if (EBeam > 4 && EBeam < 5)
      binning = GetUniformBinning(15, 1.2, 3.8);
  } else if (observable == "DiffECal") {
    if (EBeam > 0 && EBeam < 2 )
      binning = GetUniformBinning(15, -0.6, 0.2);
    else if (EBeam > 2 && EBeam < 4)
      binning = GetUniformBinning(15, -0.6, 0.2);
    else if (EBeam > 4 && EBeam < 5)
      binning = GetUniformBinning(15, -0.6, 0.2);
  }
  else if (observable == "pfl_theta")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 20, 50);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 20, 50);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 15, 50);
    }
  else if (observable == "pfl_phi")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "pfl")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0.35, 0.9);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.5, 1.9);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 1.1, 3.8);
    }
  else if (observable == "pfl_T")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0.2, 0.6);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.3, 0.9);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0.5, 1.2);
    }
  else if (observable == "proton_mom")
    {
      if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0.2001, 1.0999);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0.2001, 1.9999);
      else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 0.2001, 2.9999);
    }
  else if (observable == "proton_theta")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 20, 80);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 20, 80);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 20, 80);
    }
  else if (observable == "proton_phi")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "pim_mom")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0.1, 0.6);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0., 1.6);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0., 2);
    }
  else if (observable == "pim_theta")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 5, 140);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 5, 140);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 5, 140);
    }
  else if (observable == "pip_mom")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0.1, 3);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.1, 3);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0.1, 3);
    }
  else if (observable == "pip_theta")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 5, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 5, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 5, 180);
    }
  else if (observable == "RecoW")
    {
      if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 1, 1.49);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 1, 1.9);
      else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 1, 2.49);
    }
  else if (observable == "RecoXBJK")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 0.9);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 0.9);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 1);
    }
  else if (observable == "RecoQ2")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0.15, 0.45);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.3, 1.5);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0.9, 3);
    }
  else if (observable == "RecoQELEnu")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0., 0.8);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.3, EBeam + 0.2);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0.9, EBeam + 0.2);
    }
  else if (observable == "Recoq3")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, EBeam + 0.2);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, EBeam + 0.2);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, EBeam + 0.2);
    }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0.0001, 0.999);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0.0001, 0.999);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0.0001, 0.999);
    }
  else if (observable == "HadDeltaPTx")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, -0.6, 0.6);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, -1, 1);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, -1, 1);
    }
  else if (observable == "HadDeltaPTy")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, -0.6, 0.6);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, -1, 1);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, -1, 1);
    }
  else if (observable == "HadDeltaPhiT" || observable == "DeltaPhiT")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 80);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 80);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 80);
    }
  else if (observable == "AlphaT")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "HadAlphaT")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "RecoEnergyTransfer")
    {
      if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0., 0.7);
      else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0.2, 1.2);
      else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 1.5, 3.);
    }
  else if (observable == "HadSystemMass")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 1, 1.6);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(20, 1, 2);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 1, 2.7);
    }
  else if (observable == "MissingEnergy")
    {
      if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(15,10, -0.15, 0.47, 0.1);
      else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(10, 12, -0.15, 1.7, 0.2);
      else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(8, 12, -0.15, 3.2, 0.2);
    }
  else if (observable == "MissingTransMomentum")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetECalBinning(15, 15, 0.3, 1.1, 0.9);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetECalBinning(15, 15, -0.7, 1.2, 0.9);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetECalBinning(15, 15, -2.5, 1.2, 0.9);
    }
  else if (observable == "CorrMissingEnergy" || observable == "CorrMissingEnergy1" || observable == "CorrMissingEnergy2" || observable == "CorrMissingEnergy3")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetECalBinning(15, 15, 0.3, 1.1, 0.9);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetECalBinning(15, 15, -0.7, 1.2, 0.9);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetECalBinning(15, 15, -2.5, 1.2, 0.9);
    }
  else if (observable == "MissingAngle")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(13, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "MissingMomentum")
    {
      if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(15, 10, -0.15, 1.5, 0.3);
      else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(10, 10, -0.15, 2.2, 0.35);
      else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(10, 10, -0.15, 3.5, 0.35);
    }
  else if (observable == "InferedNucleonMom")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 0.8);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 1);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 1);
    }
  else if (observable == "HadronsAngle")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "AdlerAngleThetaP")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    }
  else if (observable == "AdlerAnglePhiP")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    } else if (observable == "AdlerAngleThetaPi")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 180);
    } else if (observable == "AdlerAnglePhiPi")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(10, 20, 180);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 20, 180);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 20, 180);
    } else if (observable == "Angleqvshad")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 120);
      else if (EBeam > 2 && EBeam < 4)
	binning = GetUniformBinning(15, 0, 120);
      else if (EBeam > 4 && EBeam < 5)
	binning = GetUniformBinning(15, 0, 60);
    }
  else if (observable == "HadDeltaPT" || observable == "DeltaPT")
    {
      if (EBeam > 0 && EBeam < 2 )
	binning = GetUniformBinning(15, 0, 0.7);
    }
  else if (observable == "TrueNProtons" || observable == "TrueNNeutrons" || observable == "TrueNPiP" || observable == "TrueNPiM" || observable == "TrueNPi0" || observable == "TrueNCh" )
    {
      binning = GetUniformBinning(5, -0.5, 5.5);
    }

  if (analysis_key == "1pip")
    {
      if (observable == "ECal")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 20, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 20, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 20, 180);
	}
      else if (observable == "RecoW")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 1, 1.49);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 1, 1.9);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 1, 2.49);
	}
      else if (observable == "AdlerAnglePhiPi")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(30, 20, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 20, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 20, 180);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 70);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 50);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 30);
	}
      else if (observable == "RecoEvPion")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 2);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0.5, 3.5);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 2, 6);
	}
      else if (observable == "ElectronPT")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0.2, 0.7);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 1);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 1);
	}
      else if (observable == "PionPT")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 0.5);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 1);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 1);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 120);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 120);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 60);
	}
      else if (observable == "HadDeltaPT" || observable == "DeltaPT")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0, 0.7);
	}
    }

  if (analysis_key == "1p1pip")
    {
      if (observable == "ECal")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(10, 10, 0.6, EBeam + 0.15, EBeam);
	  else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(20, 10, 0.6, EBeam + 0.15, EBeam);
	  else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(20, 10, 1.8, EBeam + 0.15, EBeam);
	}
      else if (observable == "pfl_theta")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 20, 50);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 20, 50);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(20, 15, 50);
	}
      else if (observable == "pfl")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0.35, 0.9);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 0.5, 1.9);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(20, 1.1, 3.8);
	}
      else if (observable == "proton_mom")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(10, 0.2001, 1.0999);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0.2001, 1.9999);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 0.2001, 2.9999);
	}
      else if (observable == "proton_theta")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 80);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(20, 0, 80);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(20, 0, 80);
	}
      else if (observable == "HadDeltaPT" || observable == "DeltaPT")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0.0001, 0.999);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0.0001, 0.999);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 0.001, 0.999);
	}
      else if (observable == "HadDeltaPhiT" || observable == "DeltaPhiT")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0, 80);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0, 80);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 0, 80);
	}
      else if (observable == "HadAlphaT")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(10, 0, 180);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 0, 180);
	}
      else if (observable == "HadSystemMass")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 1, 1.6);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 1, 2);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 1, 2.7);
	}
      else if (observable == "HadronsAngle")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(20, 0, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(30, 0, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(30, 0, 180);
	}
      else if (observable == "pip_mom")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 0.1, 0.6);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 0.1, 1.2);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(20, 0.1, 2.2);
	}
      else if (observable == "pip_theta")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 5, 130);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(20, 5, 130);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(20, 5, 130);
	}
      else if (observable == "MissingEnergy")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(5,8, -0.15, 0.47, 0.1);
	  else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(10, 15, -0.15, 1.7, 0.2);
	  else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(10, 15, -0.15, 3.2, 0.2);
	}
      else if (observable == "CorrMissingEnergy"||observable == "CorrMissingEnergy1"||observable == "CorrMissingEnergy2"||observable == "CorrMissingEnergy3")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0.5, 1);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 1);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 1);
	}
      else if (observable == "MissingAngle")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(10, 0, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 180);
	}
      else if (observable == "MissingMomentum")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetECalBinning(5,8, -0.1, 1.5, 0.1);
	  else if (EBeam > 2 && EBeam < 4) binning = GetECalBinning(5, 10, -0.15, 2, 0.15);
	  else if (EBeam > 4 && EBeam < 5) binning = GetECalBinning(5, 10, -0.15, 4, 0.2);
	}
      else if (observable == "InferedNucleonMom")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(30, 0, 0.8);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(30, 0, 1);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 1);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(30, 0, 120);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(30, 0, 120);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 120);
	}
      else if (observable == "HadronsAngle")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(20, 0, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(30, 0, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(30, 0, 180);
	}
      else if (observable == "AdlerAngleThetaPi")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(10, 0, 180);
	}
      else if (observable == "AdlerAnglePhiPi")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 180);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(10, 0, 180);
	}
      else if (observable == "Angleqvshad")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0, 120);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 0, 120);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 60);
	}
      else if (observable == "RecoW")
	{
	  if (EBeam > 0 && EBeam < 2 ) binning = GetUniformBinning(15, 1, 1.49);
	  else if (EBeam > 2 && EBeam < 4) binning = GetUniformBinning(20, 1, 1.9);
	  else if (EBeam > 4 && EBeam < 5) binning = GetUniformBinning(15, 1, 2.49);
	}
      else if (observable == "RecoQ2")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0.15, 0.45);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(20, 0.3, 1.5);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(20, 0.9, 3);
	}
    }
  else if (analysis_key == "1pim")
    {
      if (observable == "ECal")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetECalBinning(15, 10, 0.6, EBeam + 0.2, EBeam);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetECalBinning(15, 10, 0.6, EBeam + 0.2, EBeam);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetECalBinning(15, 10, 1.2, EBeam + 0.2, EBeam);
	}
      if (observable == "RecoEvPion")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(50, 0.5, 2);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(50, 0.5, 3.5);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(50, 1, 7);
	}
      if (observable == "RecoWPion")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(15, 0.9, 2);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(15, 1, 2);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(15, 0, 4);
	}
      else if (observable == "ElectronPT")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(50, 0, 0.5);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(50, 0, 1);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(50, 0, 1);
	}
      else if (observable == "PionPT")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(50, 0, 0.5);
	}
      else if (observable == "RecoEvPion")
	{
	  if (EBeam > 0 && EBeam < 2 )
	    binning = GetUniformBinning(30, 0.5, 2);
	  else if (EBeam > 2 && EBeam < 4)
	    binning = GetUniformBinning(30, 0.5, 3.5);
	  else if (EBeam > 4 && EBeam < 5)
	    binning = GetUniformBinning(30, 2, 6);
	}
    }
  else if (analysis_key == "1pip") {

    if (observable == "ECal")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetECalBinning(15, 10, 0.6, EBeam + 0.2, EBeam);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetECalBinning(15, 10, 1, EBeam + 0.2, EBeam);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetECalBinning(15, 10, 1, EBeam + 0.2, EBeam);
      }
    else if (observable == "pip_mom")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(15, 0, 0.6);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(15, 0.3, 1.2);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(15, 0, 2.5);
      }
    else if (observable == "pip_theta")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(15, 0, 130);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(15, 0, 120);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(15, 0, 100);
      }
    if (observable == "RecoEvPion")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(50, 0.5, 2);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(50, 0.5, 3.5);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(50, 1, 7);
      }
    if (observable == "RecoWPion")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(50, 0.9, 2);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(50, 1, 2);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(50, 2, 4);
      }
    else if (observable == "ElectronPT")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(50, 0, 0.5);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(50, 0, 1);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(50, 0, 1);
      }
    else if (observable == "PionPT")
      {
	if (EBeam > 0 && EBeam < 2 )
	  binning = GetUniformBinning(50, 0, 0.5);
	else if (EBeam > 2 && EBeam < 4)
	  binning = GetUniformBinning(50, 0, 1);
	else if (EBeam > 4 && EBeam < 5)
	  binning = GetUniformBinning(50, 0, 1);
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

