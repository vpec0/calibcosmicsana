# calibcosmicsana
_Modified: Aug 25, 2020_

## Intended use

This code is to be used by analyzers from the DUNE calibration working group who are interested in analysing the Sheffield cosmic muon production for DUNE FD.

The code contains a LArSoft analyser module to be used to create a common ntuple(s) for the calibration-related studies.

At the moment, the analyser creates the 1st version of a simple ntuple good for use in electron lifetime study. Contents of the ntuple should be discussed by all interested parties.

I chose to make this package separate from dunetpc for the sake of ease of use, but, of course, the product depends on `dunetpc` (and consequently on `larsoft`).

## Installation
```
# in a fresh shell

mkdir testdev # use any name you desire
cd testdev/
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_05_00 -q e19:prof
mrb newDev
. localProducts_larsoft_*/setup
mrb g calibcosmicsana%ssh://git@github.com/vpec0/calibcosmicsana.git 
# or one can use https protocol instead:
# mrb g https://github.com/vpec0/calibcosmicsana
mrbsetenv
mrb i
```

## Test/run

```
# in a fresh shell

cd testdev

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v08_60_00 -q e19:prof
. localProducts_larsoft_*/setup
mrbslp

mkdir work
cd work

MYFILE=/pnfs/dune/persistent/users/calibration/cosmic_muons/sheffield_prod/dunetpc_v08_50_00/20003100/20003100/MUSUN_dunefd_20003100_gen_g4_detsim_reco.root

lar -c CalibCosmicsAna.fcl --no-output $MYFILE
```

### Note
The `dunetpc` version above may change in the future, check in the `ups/product_deps` file which version the code currently depends on.
