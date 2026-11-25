# e4nu Comparisons (Rivet)

Rivet analysis to compare any event generator with the e4nu measurement of single charged pion
production on carbon, C(e,e'p π⁻) and C(e,e'p π⁺), with the CLAS detector. The generated events
are read in HepMC3 format (e.g. NuHepMC), so any generator in electron mode can be compared.

The analysis is `CLAS_2026_I0000000`. The name is a placeholder until the paper has an Inspire ID.

| File | Content |
|---|---|
| `CLAS_2026_I0000000.cc` | the analysis |
| `CLAS_2026_I0000000.info` | metadata |
| `CLAS_2026_I0000000.plot` | titles, axis labels and ranges of the plots |

The measured cross sections are not part of this repository. Without them the analysis runs and
produces the model predictions; to overlay the data, see [Reference data](#reference-data).

## Distributions

Differential cross sections dσ/dW and dσ/dE_cal, in nb/GeV per carbon nucleus, at beam energies of
1.159, 2.257 and 4.453 GeV. Each histogram corresponds to the HEPData table of the same number.

| | 1.159 GeV | 2.257 GeV | 4.453 GeV |
|---|---|---|---|
| W, 1p1π⁻ | d01 | d02 | d03 |
| W, 1p1π⁺ | d04 | d05 | d06 |
| E_cal, 1p1π⁻ | d07 | d08 | d09 |
| E_cal, 1p1π⁺ | d10 | d11 | d12 |

## Event selection

1. The momentum of every final-state particle is smeared with a Gaussian of the CLAS resolution
   (e⁻ 0.5%, p 1%, π± and γ 0.7%; three times larger at 1.159 GeV). The direction is unchanged.
2. Electron (the most energetic final-state e⁻): Q² above 0.1 / 0.4 / 0.8 GeV² at 1.159 / 2.257 /
   4.453 GeV, polar angle between the momentum-dependent CLAS minimum and 45°, momentum above
   0.4 / 0.55 / 1.3 GeV.
3. Hadrons are detected above the CLAS thresholds: momentum (p and γ 0.3 GeV, π⁻ 0.15 GeV,
   π⁺ 0.25 GeV), polar angle above the minimum of each particle and below 140°.
4. Exactly one detected proton, one detected π⁻ (and no π⁺) or one detected π⁺ (and no π⁻), and no
   detected photon. Photons count only if they pass the radiation cut (within 40° of the beam
   electron and 30° of its azimuth).
5. W = √((m_N + ν)² − |q|²), with m_N = 0.9389 GeV. E_cal = E_e' + E_π + T_p + 20.95 MeV, with the
   total energy of the electron and of the pion and the kinetic energy of the proton; 20.95 MeV is
   the ¹²C separation energy plus the E_cal offset.

## Requirements

[Rivet](https://rivet.hepforge.org) 4 (tested with 4.1.3) with YODA 2. Plotting needs a LaTeX
installation.

## Build

```bash
rivet-build RivetCLAS_2026_I0000000.so CLAS_2026_I0000000.cc
export RIVET_ANALYSIS_PATH=$PWD RIVET_DATA_PATH=$PWD
```

## Run

One run per beam energy. Each run only fills the histograms of its own energy, so the runs are
merged before plotting:

```bash
rivet -a CLAS_2026_I0000000 --ignore-beams --weight-nominal CV -o model_1GeV.yoda events_eC_1159MeV.hepmc3
rivet -a CLAS_2026_I0000000 --ignore-beams --weight-nominal CV -o model_2GeV.yoda events_eC_2257MeV.hepmc3
rivet -a CLAS_2026_I0000000 --ignore-beams --weight-nominal CV -o model_4GeV.yoda events_eC_4453MeV.hepmc3

rivet-merge -o model.yoda model_1GeV.yoda model_2GeV.yoda model_4GeV.yoda
rivet-mkhtml -o plots model.yoda:"Title=My model"
```

Several models can be shown together: `rivet-mkhtml -o plots genie.yoda:"Title=GENIE" achilles.yoda:"Title=Achilles"`.
Labels go through LaTeX: avoid `_` in titles, or write it as `\_`.

### Requirements on the event files

- **Process:** e⁻ on ¹²C. The beam energy is read from the incoming electron (status 4) of each event
  and must be within 5% of 1.159, 2.257 or 4.453 GeV. Events with another target are skipped.
- **Cross section:** read from the event files (`GenCrossSection`, in pb) and used together with the
  sum of the event weights, so weighted and unweighted events are both fine. It must be **per carbon
  nucleus** (NuHepMC: `NuHepMC.Units.CrossSection.TargetScale PerTargetAtom`). If it is not, give the
  right value with `rivet -x <cross section in pb>`.
- **Units:** GeV or MeV; Rivet converts them.

### Options

- `--ignore-beams`: needed for NuHepMC files. The target nucleus is not flagged as a beam particle
  (status 11 or 20), and some generators write the struck nucleon without a parent particle, so
  Rivet cannot always identify the beams. The analysis checks the electron energy and the target
  itself.
- `--weight-nominal CV`: GENIE and Achilles call their nominal weight `CV`, which Rivet does not
  recognise; without this option the output has no statistical uncertainties. Use the name of the
  nominal weight of your generator.

## Reference data

The measured cross sections are not distributed with this repository. To compare with the data,
download the YODA file of the HEPData record, save it as `CLAS_2026_I0000000.yoda.gz` in a directory
of `RIVET_DATA_PATH` and run `rivet-mkhtml`: the plots then show the data, the models and their
ratio.
