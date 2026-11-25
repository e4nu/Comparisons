// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {

  /// @brief e4nu/CLAS6: C(e,e'p pi-) and C(e,e'p pi+) cross sections in W and E_cal at 1.159, 2.257 and 4.453 GeV
  /// Reference: [Publication link]
  /// The event selection accounts for the following: momenta are smeared with the CLAS resolutions,
  /// the CLAS momentum and angle thresholds are applied to the electron and to
  /// the hadrons, and exactly one proton and one charged pion (and no photon) above threshold must be detected.

  class CLAS_2026_I0000000 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CLAS_2026_I0000000);

    /// @name Analysis methods
    ///@{

    void init() {
      declare(FinalState(), "FS");

      // One histogram per HEPData table: W tables first, then E_cal; 1p1pi- then 1p1pi+; 1, 2, 4 GeV.
      // The binning is given explicitly, so the analysis also runs without the reference data
      for (size_t iobs = 0; iobs < 2; ++iobs)
        for (size_t itopo = 0; itopo < 2; ++itopo)
          for (size_t ie = 0; ie < 3; ++ie) {
            const size_t id = 6*iobs + 3*itopo + ie;
            book(_h[iobs][itopo][ie], mkAxisCode(id + 1, 1, 1), _edges[id]);
          }
    }


    void analyze(const Event& event) {
      // Incoming electron and nucleus
      ConstGenParticlePtr beam_e, target;
      for (ConstGenParticlePtr p : HepMCUtils::particles(event.genEvent())) {
        if (p->status() == 4 && p->pid() == PID::ELECTRON) beam_e = p;
        else if ((p->status() == 4 || p->status() == 11) && p->pid() > 1000000000) target = p;  // NuHepMC: 11
      }
      if (!beam_e) vetoEvent;
      const FourMomentum pbeam(beam_e->momentum());
      const double Ebeam = pbeam.E();
      const int ie = energyIndex(Ebeam);
      if (ie < 0) {
        MSG_WARNING("Beam energy " << Ebeam << " GeV is none of 1.159, 2.257, 4.453 GeV: event skipped");
        vetoEvent;
      }
      if (target && target->pid() != _targetPdg) {
        MSG_WARNING("Target " << target->pid() << " is not 12C: event skipped");
        vetoEvent;
      }

      // Scattered electron: the most energetic final-state electron, smeared
      const FinalState& fs = apply<FinalState>(event, "FS");
      const Particles electrons = fs.particles(Cuts::pid == PID::ELECTRON, cmpMomByE);
      if (electrons.empty()) vetoEvent;
      const FourMomentum el = smear(electrons[0].mom(), PID::ELECTRON, Ebeam);
      const double pel = el.p3().mod();

      // Electron cuts: Q2, polar angle and momentum
      const FourMomentum q = el - pbeam;
      if (-q.mass2() < q2Cut(Ebeam)) vetoEvent;
      if (el.theta()/degree < minTheta(PID::ELECTRON, pel, Ebeam)) vetoEvent;
      if (el.theta()/degree > 45.) vetoEvent;
      if (pel < minMomentum(PID::ELECTRON, Ebeam)) vetoEvent;

      // Detected hadrons: smeared, above the momentum and angle thresholds
      size_t np = 0, npim = 0, npip = 0, ngamma = 0;
      vector<pair<int, FourMomentum>> detected;  // protons and charged pions
      for (const Particle& p : fs.particles()) {
        const int pid = p.pid();
        const FourMomentum ps = smear(p.mom(), pid, Ebeam);
        const double pmod = ps.p3().mod();
        if (ps.theta()/degree < minTheta(pid, pmod, Ebeam)) continue;
        if (ps.theta()/degree >= 140.) continue;
        if (pmod < minMomentum(pid, Ebeam)) continue;
        if (pid == PID::PHOTON && !photonRadCut(pbeam, ps)) continue;

        if (pid == PID::PROTON) ++np;
        else if (pid == PID::PIMINUS) ++npim;
        else if (pid == PID::PIPLUS) ++npip;
        else if (pid == PID::PHOTON) ++ngamma;
        if (pid == PID::PROTON || pid == PID::PIMINUS || pid == PID::PIPLUS) detected.push_back({pid, ps});
      }

      // Topology: 1p1pi- or 1p1pi+, no photons
      if (np != 1 || ngamma != 0) vetoEvent;
      int itopo = -1;
      if (npim == 1 && npip == 0) itopo = 0;
      else if (npim == 0 && npip == 1) itopo = 1;
      if (itopo < 0) vetoEvent;

      // Observables
      const double nu = Ebeam - el.E();
      const double W2 = sqr(kNucleonMass + nu) - (pbeam - el).p3().mod2();
      const double W = W2 > 0 ? sqrt(W2) : 0.;
      double Ecal = el.E();
      for (const auto& h : detected) {
        Ecal += h.second.E();
        if (h.first == PID::PROTON) Ecal += kBindingC12 - kProtonMass;
      }

      _h[0][itopo][ie]->fill(W);
      _h[1][itopo][ie]->fill(Ecal);
    }


    /// Cross section per carbon nucleus, in nb / GeV
    void finalize() {
      const double sf = crossSection()/nanobarn/sumOfWeights();
      for (auto& hobs : _h)
        for (auto& htopo : hobs)
          scale(htopo, sf);
    }

    ///@}


  private:

    /// Index of the beam energy (1.159, 2.257, 4.453 GeV), within 5%, or -1
    int energyIndex(double E) const {
      for (size_t i = 0; i < _energies.size(); ++i)
        if (fabs(E - _energies[i]) < 0.05*_energies[i]) return i;
      return -1;
    }

    /// Momentum resolution (relative), three times worse at 1.1 GeV
    double resolution(int pid, double Ebeam) const {
      double res = 0.;
      if (pid == PID::PROTON) res = 0.01;
      else if (pid == PID::ELECTRON) res = 0.005;
      else if (pid == PID::PIPLUS || pid == PID::PIMINUS || pid == PID::PHOTON) res = 0.007;
      if (Ebeam < 2) res *= 3;
      return res;
    }

    double mass(int pid) const {
      if (pid == PID::PROTON) return kProtonMass;
      if (pid == PID::PIPLUS || pid == PID::PIMINUS) return 0.139570;
      if (pid == PID::ELECTRON) return 0.000510998;
      return 0.;
    }

    /// Gaussian smearing of the momentum magnitude; the direction is kept
    FourMomentum smear(const FourMomentum& mom, int pid, double Ebeam) const {
      const double p = mom.p3().mod();
      if (p <= 0) return mom;
      const double ps = randnorm(p, resolution(pid, Ebeam)*p);
      const Vector3 p3 = mom.p3() * (ps/p);
      return FourMomentum(sqrt(sqr(ps) + sqr(mass(pid))), p3.x(), p3.y(), p3.z());
    }

    double q2Cut(double Ebeam) const {
      if (Ebeam > 1 && Ebeam < 2) return 0.1;
      if (Ebeam > 2 && Ebeam < 4) return 0.4;
      if (Ebeam > 4) return 0.8;
      return 0.;
    }

    /// Minimum polar angle in degrees
    double minTheta(int pid, double p, double Ebeam) const {
      double th = 0.;
      if (pid == PID::ELECTRON) {
        if (Ebeam < 2) th = 17. + 7./p;
        else if (Ebeam > 2 && Ebeam < 4) th = 16. + 10.5/p;
        else if (Ebeam > 4) th = 13.5 + 15./p;
        if (th < 15.) th = 15.;
      }
      else if (pid == PID::PROTON) th = 10.;
      else if (pid == PID::PIMINUS) {
        if (Ebeam < 2) th = 17. + 4./p;
        else th = p < 0.35 ? 25. + 7./p : 16. + 10./p;
      }
      else if (pid == PID::PIPLUS) th = 10.;
      else if (pid == PID::PHOTON) th = 8.;
      return th;
    }

    /// Minimum momentum in GeV
    double minMomentum(int pid, double Ebeam) const {
      if (pid == PID::ELECTRON) {
        if (Ebeam < 2) return 0.4;
        if (Ebeam > 2 && Ebeam < 4) return 0.55;
        if (Ebeam > 4) return 1.3;
      }
      else if (pid == PID::PROTON || pid == PID::PHOTON) return 0.3;
      else if (pid == PID::PIMINUS) return 0.15;
      else if (pid == PID::PIPLUS) return 0.25;
      return 0.;
    }

    /// Photon within 40 deg of the electron and 30 deg of its azimuth (radiation-like photon).
    /// The electron given is the beam electron
    bool photonRadCut(const FourMomentum& electron, const FourMomentum& gamma) const {
      double gphi = gamma.phi(MINUSPI_PLUSPI)/degree + 30.;
      if (gphi < 0) gphi += 360.;
      double ephi = electron.phi(MINUSPI_PLUSPI)/degree + 30.;
      if (ephi < 0) ephi += 360.;
      return gamma.angle(electron)/degree < 40. && fabs(gphi - ephi) < 30.;
    }


    /// @name Constants (GeV)
    ///@{
    static constexpr double kProtonMass = 0.9382720813;
    static constexpr double kNucleonMass = 0.9389;
    /// E_cal correction for the knocked-out proton in 12C: B(12C) - B(11B) + offset
    static constexpr double kBindingC12 = 0.09215 - 0.0762 + 0.005;
    ///@}

    const int _targetPdg = 1000060120;
    const vector<double> _energies{1.159, 2.257, 4.453};

    /// Bin edges of the HEPData tables d01 ... d12
    const vector<vector<double>> _edges{
      /* d01 */ {1, 1.0327, 1.0653, 1.098, 1.1307, 1.1633, 1.196, 1.2287, 1.2613, 1.294, 1.3267, 1.3593, 1.392, 1.4247, 1.4573, 1.49},
      /* d02 */ {1, 1.06, 1.12, 1.18, 1.24, 1.3, 1.36, 1.42, 1.48, 1.54, 1.6, 1.66, 1.72, 1.78, 1.84, 1.9},
      /* d03 */ {1, 1.0993, 1.1987, 1.298, 1.3973, 1.4967, 1.596, 1.6953, 1.7947, 1.894, 1.9933, 2.0927, 2.192, 2.2913, 2.3907, 2.49},
      /* d04 */ {1.0653, 1.098, 1.1307, 1.1633, 1.196, 1.2287, 1.2613, 1.294, 1.3267, 1.3593, 1.392, 1.4247, 1.4573, 1.49},
      /* d05 */ {1, 1.06, 1.12, 1.18, 1.24, 1.3, 1.36, 1.42, 1.48, 1.54, 1.6, 1.66, 1.72, 1.78, 1.84, 1.9},
      /* d06 */ {1, 1.0993, 1.1987, 1.298, 1.3973, 1.4967, 1.596, 1.6953, 1.7947, 1.894, 1.9933, 2.0927, 2.192, 2.2913, 2.3907, 2.49},
      /* d07 */ {0.65029, 0.70059, 0.75088, 0.80118, 0.85147, 0.90177, 0.95206, 1.0024, 1.0527, 1.1029, 1.1238, 1.1446, 1.1654, 1.1862, 1.207, 1.2278},
      /* d08 */ {0.75479, 0.90959, 1.0644, 1.2192, 1.374, 1.5288, 1.6836, 1.8384, 1.9932, 2.1479, 2.1743, 2.2006, 2.2269, 2.2532, 2.2795, 2.3058},
      /* d09 */ {1.8, 2.0438, 2.2876, 2.5314, 2.7752, 3.019, 3.2628, 3.5066, 3.7504, 3.9942, 4.2379, 4.2753, 4.3126, 4.3499, 4.3872, 4.4245, 4.4618, 4.4991, 4.5364},
      /* d10 */ {0.75088, 0.80118, 0.85147, 0.90177, 0.95206, 1.0024, 1.0527, 1.1029, 1.1238, 1.1446, 1.1654, 1.1862, 1.207, 1.2278},
      /* d11 */ {0.90959, 1.0644, 1.2192, 1.374, 1.5288, 1.6836, 1.8384, 1.9932, 2.1479, 2.1743, 2.2006, 2.2269, 2.2532, 2.2795, 2.3058},
      /* d12 */ {1.8, 2.0438, 2.2876, 2.5314, 2.7752, 3.019, 3.2628, 3.5066, 3.7504, 3.9942, 4.2379, 4.2753, 4.3126, 4.3499, 4.3872, 4.4245, 4.4618, 4.4991, 4.5364},
    };

    /// @name Histograms [W, E_cal][1p1pi-, 1p1pi+][1.159, 2.257, 4.453 GeV]
    Histo1DPtr _h[2][2][3];

  };


  RIVET_DECLARE_PLUGIN(CLAS_2026_I0000000);

}
