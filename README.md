# gnomad_review_hum_mut
This is the repository containing the supplementary code for the paper _Variant interpretation using the world’s largest population database: gnomAD_.

Description of files:

`get_data_for_figures.py` Script contaning all the code needed to generate a random subset of samples and filter gnomAD v.2.1.1 to rare variants of interest for those samples.

Environment and packages used:
* Python 3.8.8
* Hail version 0.2.65-367cf1874d85
* [gnomad](https://github.com/broadinstitute/gnomad_methods) - `get_data_for_figures.py` uses additional code currenlty not found in a release version of gnomad. You will need to add the gnomad repo to your python path.
