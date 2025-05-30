# seti_end_to_end

This repository contains Jupyter notebooks for running seticore/seticore2 or turbo_seti from end to end, starting with .raw files
(time domain signals from actual captures or synthetic signals produced by setigen and seti_test_file_gen). 
The search program is used to search and detect drifting narrow-band tones or technosignatures which might be produced by 
an extra-terrestrial source. 

Only dedoppler mode is currently evaluated here, which assumes a single-antenna input like Green Bank Telescope (GBT), 
though parameters can be set up for other telescopes (e.g. MeerKAT or VLA) with the raw file corresponding
to simulated beam data.

Of particular interest are synthetic raw files containing a large number of drifting tones which can be used
to verify seticore's detection performance.

References: 

K. M. Houston, "Fine-Tuning the Narrowband SETI Signal Processing Pipeline", 
75th International Astronautical Congress (IAC), Milan, Italy, 14-18 October 2024.

K. M. Houston, "Fine-Tuning the Narrowband SETI Signal Processing Pipeline", 
Acta Astronautica, AA-D-25-00258, in preparation.

## Notebooks description

Notebook "01_seti_end_to_end.ipynb" does the basic end-to-end operation.  A raw file is input
to rawspec with fine fft size and the number of sti averages ("n_sti", specified as "-t" or "--ints" in rawspec) as parameters.  Rawspec
produces averaged spectrograms (matrices of energy vs time and frequency) in the form of filterbank files (.h5).
A search program (seticore/seticore2/turbo_seti/bliss2) is then called to do the search for narrowband drifting tones.  Detections or "hits"
are declared and output to a text file.  Plots of hits (e.g. signal-to-noise ratio (SNR) vs drift rate or frequency) are produced.

Notebook "02_seti_param_sweep.ipynb" runs "01_seti_end_to_end.ipynb" many times for a single raw file, 
sweeping fine fft size and n_sti over a wide range, and does composite plots.  This evaluates the fft size/n_sti trade space.
It is intended to work with special chirp test files created by 00_multichirp_raw_file_gen.ipynb (in seti_test_file_gen repository), 
which inserts a large number of chirp signals covering a wide range of drift rates.
One may observe considerable variation in the SNR of the detections as a function of drift rate,
plus variation of total processing time.  The goal is to maximize SNR for a given level of compute time.
The detection data are also stored in numpy npz files for later analysis.

Notebook "03_multi_param_sweep.ipynb" runs "02_seti_param_sweep.ipynb" over multiple raw input files.

Notebook "04_multi_param_sweep_plot.ipynb" inputs detection data in npz files and produces SNR or compute time
plots over multiple raw files, or multiple test cases.  Performance of multiple seticore2 branches can be compared.

Notebook "05_multi_raw_e2e.ipynb" does a bulk run of multiple raw capture files, such as Voyager or other non-synthetic captures, 
by running "01_seti_end_to_end.ipynb" multiple times.  Spectrograms of individual detections may be created as png files.

## About seticore2

The seticore2 fork evaluates possible improvements to the De-Doppler mode of seticore, including a
variant De-Doppler (DD) function called "fastDD" (including a cpu-only version), 
use of explicit host-gpu memory copying (avoiding the unified memory model and "managed" 
background copying), and improved spectrogram normalization/threshold setting/SNR estimation.

At present, changes have been made only to De-Doppler mode, and assume an h5 single-antenna input
analogous to TurboSETI on GBT.  Operation of cadence mode and beamforming mode (critical to 
MeerKAT and VLA) have not been tested and probably need additional fixes.

Analogous runs of bliss can be run to enable side-by-side comparisons with seticore2 and turbo_seti.

Related repositories include:

https://github.com/khouston22/seticore2 (a fork of seticore https://github.com/lacker/seticore)

https://github.com/khouston22/seti_test_file_gen,

https://github.com/UCBerkeleySETI/rawspec, and

https://github.com/khouston22/seti_detect_sim.

https://github.com/khouston22/bliss2 (a fork of bliss https://github.com/n-west/bliss)


## Setting up environment variables

Examples of environment variables that may need to be set up in ~/.profile include the following:

```
export TURBO_SETI_PATH=$HOME/Dropbox/kgit/turbo_seti
export SETICORE2_PATH=$HOME/kgit/seticore2/build
export SETICORE2_PY_PATH=$HOME/kgit/seticore2/python
export BLISS2_PATH=$HOME/kgit/bliss2/build/bliss
export BLISS2_PLOT_PATH=$HOME/kgit/bliss2/bliss/python/blissdedrift/plot_utils
export DATADIR=/datax/scratch/khouston
export RAWDIR=/datax/scratch/khouston/raw_test_files
export RAW=/datax/scratch/khouston/raw_test_files
export RAW_BACKUP_BASE_DIR=/datax/scratch/khouston/temp
export SGDIR=/datax/scratch/khouston/sg_det_files
export SG=/datax/scratch/khouston/sg_det_files
export SC2=$HOME/kgit/seticore2
export E2E=$HOME/kgit/seti_end_to_end
export BL2=$HOME/kgit/bliss2/build/bliss
```

## Setting up a conda environment

The notebooks have been successfully run in python 3.10.16.

Commands to install additional packages include:

```
pip install matplotlib
pip install astropy
pip install blimpy
pip install setigen
pip install turbo_seti
pip install nbformat
pip install nbclient
```
