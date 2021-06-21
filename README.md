# Variant interpretation using the world’s largest population database: gnomAD
This is the repository containing the supplementary code for the paper _Variant interpretation using the world’s largest population database: gnomAD_.

Description of files:

`get_data_for_figures.py` Script containing all the code needed to generate a random subset of samples and filter gnomAD v.2.1.1 exomes to rare variants of interest for those samples.

### Environment and packages used to generate data for paper:
* Python 3.8.8
* [Hail](https://hail.is/docs/0.2/index.html) version 0.2.65-367cf1874d85
* [gnomad](https://github.com/broadinstitute/gnomad_methods) - `get_data_for_figures.py` uses additional code currently not found in a release version of gnomad. You will need to add the gnomad repo to your python path. The latest commit as of June 21, 2021 is [here](https://github.com/broadinstitute/gnomad_methods/commit/7449ca6d9880b9074776156dfc38e6db4f3d2f92).
* [gnomad_qc](https://github.com/broadinstitute/gnomad_qc)

### Usage
Run the code yourself:
```sh
python3 get_data_for_figures.py

#Get help descriptions on the optional arguments:
python3 get_data_for_figures.py -h
```
