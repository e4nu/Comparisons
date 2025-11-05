# E4Nu Comparisons

This reposiory allows for quick comparison against e4nu data using neutrino event generators on electron mode. The only requirement is for the generated events to be stored in NuHEPMC format. The available data is stored in the data repository. There is a simple macro with the corresponding kinematic cuts used to compare the simulations against data.

Two main applications are available. 
- comparison_1p1pi: this app computes the (e,e'1p1pi-) or (e,e'1p1pi+) cross-section as a function of a given observable for a model of interest. The code is specific for the (e,e'1p1pi) topology
- plot_comparisons: This is a generalized code used to plot the cross-section results (i.e. from comparison_1p1pi) against data

## Build software
To build the code simply follow these instructions: 
```
  source e4nu_comparisons_gpvm_env.sh ; mkdir build ; cd build ; cmake .. ; make;
```

## comparison_1p1pi
An example to run the comparison_1p1pi code is given below: 
```
    ./comparison_1p1pi --input-hepmc3-directory /pnfs/genie/persistent/users/jtenavid/e4nu_files/NuHEPMC/Carbon/1GeV/ --output-file comparison_1p1pim_1GeV_GENIE --topology 1p1pim --model-name GENIE
```


## plot_comparisons
An example to run the plot_comparisons code is given below: 
```
  
observables_log=("ECal" "RecoW" "RecoQ2" "proton_mom" "proton_theta" "pim_mom" "pim_theta" "HadAlphaT" "HadDeltaPT")
for OBS in "${observables_log[@]}"; do
  ./plot_comparisons --mc-files comparison_1p1pim_1GeV_Genie.root,comparison_1p1pim_1GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pim  --data e4nu_1p1pim_1GeV_${OBS}.root --output-file comparisons_1p1pim_1GeV_${OBS}
  ./plot_comparisons --mc-files comparison_1p1pim_2GeV_Genie.root,comparison_1p1pim_2GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pim  --output-file comparisons_1p1pim_2GeV_${OBS}
  ./plot_comparisons --mc-files comparison_1p1pim_4GeV_Genie.root,comparison_1p1pim_4GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pim  --output-file comparisons_1p1pim_4GeV_${OBS}
done

observables=("ECal" "RecoW" "RecoQ2" "proton_mom" "proton_theta" "pip_mom" "pip_theta" "HadAlphaT" "HadDeltaPT")
for OBS in "${observables_log[@]}"; do
  ./plot_comparisons --mc-files comparison_1p1pip_1GeV_Genie.root,comparison_1p1pip_1GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pip  --data e4nu_1p1pip_1GeV_${OBS}.root --output-file comparisons_1p1pip_1GeV_${OBS}
  ./plot_comparisons --mc-files comparison_1p1pip_2GeV_Genie.root,comparison_1p1pip_2GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pip  --output-file comparisons_1p1pip_2GeV_${OBS}
  ./plot_comparisons --mc-files comparison_1p1pip_4GeV_Genie.root,comparison_1p1pip_4GeV_Achilles.root --mc-names GENIE,Achilles --observable ${OBS} --topology 1p1pip  --output-file comparisons_1p1pip_4GeV_${OBS}
done
```
This computes the comparisons for 1p1pi- and 1p1pi+ data using GENIE and Achilles.
