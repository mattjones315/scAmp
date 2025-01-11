# scAmp
scAmp (single-cell Amplicon) is a python-based workflow for detecting and
analyzing focal amplifications from single-cell data.

It consists of three main modules (in progress):

* copy-number inference: a set of modules for inferring copy-numbers from single-cell data
* ecDNA detection: detection of extrachromsomal or chromosomal DNA amplifications
* clonal analysis: Analysis of clonal history with respect to amplifications

## Installation

You can install scAmp by cloning this directory and running `pip install .` This should install scAmp and all dependencies.

## Running scAmp from the command line

You can invoke `scamp` modules from the command line by running

`scamp [module] [arguments]`

Currently there are two command line modules you can run:

* `scamp atac-cnv`: Computes copy-numbers across genes from a scATAC fragments files.
* `scamp predict`: Predicts ecDNA status from copy-number data.

You can look at usage instructions for these modules by running `scamp [module] --help`.