import argparse
import logging

import hail as hl
from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.grch37.gnomad import EXOME_POPS, public_release as v2_public_release

from file_utils import load_public_resources

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("filter_neurodev_variants")
logger.setLevel(logging.INFO)
gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)

def main(args):
    
    logger.info("Loading variant table, v2 exomes, v2 genomes, v3 variants, and v2 liftover...")
    ht = hl.read_table(args.sample_with_variants_path)
    original_count = ht.count()
    v2_exomes_ht = v2_public_release("exomes").ht()
    v2_genomes_ht, v3_genomes_ht, _ = load_public_resources()
    v2_exome_liftover_ht = hl.read_table(
            "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
        )
    v2_liftover_index = v2_exome_liftover_ht[ht.key]
    ht = ht.annotate(
        grch37_locus=v2_liftover_index.original_locus,
        grch37_alleles=v2_liftover_index.original_alleles,
    )
    logger.info(
        "Filtering variants not found in v2 exomes, v2 genomes or v3 samples but found in v2 exomes..."
    )
    freq_adj_idx = hl.eval(v2_exomes_ht.freq_meta.index({"group": "adj"}))
    ht = ht.annotate(
        v2_exomes_adj_freq=v2_exomes_ht[ht.grch37_locus, ht.grch37_alleles].freq[freq_adj_idx]
    )
    #Keep variants that are undefined.
    ht = ht.filter(hl.or_else(ht.v2_exomes_adj_freq.AC == 0, True))
    freq_adj_idx = hl.eval(v2_genomes_ht.freq_meta.index({"group": "adj"}))
    ht = ht.annotate(
        v2_genomes_adj_freq=v2_genomes_ht[ht.grch37_locus, ht.grch37_alleles].freq[freq_adj_idx]
    )
    # Keep variants that are undefined.
    ht = ht.filter(hl.or_else(ht.v2_genomes_adj_freq.AC == 0, True))
    freq_adj_idx = hl.eval(v3_genomes_ht.freq_meta.index(
        {"group": "adj", "subset": "non_v2"}
    ))
    ht = ht.annotate(
        v3_non_v2_adj_freq=v3_genomes_ht[ht.key].freq[
            freq_adj_idx
        ]
    )
    # Keep variants that are undefined.
    ht = ht.filter(hl.or_else(ht.v3_non_v2_adj_freq.AC == 0, True))

    logger.info(
        f"Writing Hail Table to {args.output_path_prefix}/neurodev_unique_samples_with_variants.ht"
    )
    ht = ht.checkpoint(
        f"{args.output_path_prefix}/neurodev_unique_samples_with_variants.ht",
        overwrite=args.overwrite,
    )
    logger.info(f"{original_count} variants were filtered to {ht.count()} variants.")
    if args.overwrite:
        logger.info(
            f"Writing TSV to {args.output_path_prefix}/neurodev_unique_samples_with_variants.tsv"
        )
        ht.export(
            f"{args.output_path_prefix}/neurodev_unique_samples_with_variants.tsv",
            header=True,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script filters a Hail Table of rare variants of interest from a random sample of individuals from each population present in gnomad exomes v.2.1.1 to variants not found in v2 genomes or the v3-non-v2 subset."
    )
    parser.add_argument(
        "--sample-with-variants-path",
        default="gs://gnomad-wphu/grant_proposal_01312022/samples_with_variants.ht",
        type=str,
        help="Path to Hail Table with variants from gnomad v2 exomes.",
    )
    parser.add_argument(
        "--output-path-prefix",
        default="gs://gnomad-wphu/grant_proposal_01312022",
        type=str,
        help="Google folder to store output.",
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
