import argparse
import logging

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("filter_v2_exome_variants")
logger.setLevel(logging.INFO)


def main(args):
    logger.info("Loading variant table, v2 genomes, v3 variants, and v2 liftover...")
    ht = hl.read_table(args.sample_with_variants_path)
    original_count = ht.count()
    v2_genomes_ht = hl.read_table(
        "gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht"
    )
    v3_genomes_ht = hl.read_table(
        "gs://gcp-public-data--gnomad/release/3.1.1/ht/genomes/gnomad.genomes.v3.1.1.sites.ht"
    )
    v2_exome_liftover_ht = hl.read_table(
        "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
    )
    v2_exome_liftover_ht = v2_exome_liftover_ht.key_by(
        "original_locus", "original_alleles"
    )
    v2_liftover_index = v2_exome_liftover_ht[ht.key]
    ht = ht.annotate(
        liftover_locus=v2_liftover_index.locus,
        liftover_alleles=v2_liftover_index.alleles,
    )

    logger.info(
        "Filtering variants not found in v2 genomes or v3 samples but found in v2 exomes..."
    )
    freq_adj_idx = hl.eval(v2_genomes_ht.freq_meta.index({"group": "adj"}))
    ht = ht.annotate(
        v2_genomes_adj_freq=v2_genomes_ht[ht.key].freq[freq_adj_idx]
    )
    v2_genomes_filtering_expr = (ht.v2_genomes_adj_freq.AC < 1) | hl.is_missing(
        ht.v2_genomes_adj_freq.AC
    )
    ht = ht.filter(v2_genomes_filtering_expr)
    freq_adj_idx = hl.eval(v3_genomes_ht.freq_meta.index(
        {"group": "adj", "subset": "non_v2"}
    ))
    ht = ht.annotate(
        v3_non_v2_adj_freq=v3_genomes_ht[ht.liftover_locus, ht.liftover_alleles].freq[
            freq_adj_idx
        ]
    )
    v3_non_v2_filtering_expr = (ht.v3_non_v2_adj_freq.AC < 1) | hl.is_missing(
        ht.v3_non_v2_adj_freq.AC
    )
    ht = ht.filter(v3_non_v2_filtering_expr)

    logger.info(
        f"Writing Hail Table to {args.output_path_prefix}/v2_exomes_unique_samples_with_variants.ht"
    )
    ht = ht.checkpoint(
        f"{args.output_path_prefix}/v2_exomes_unique_samples_with_variants.ht",
        overwrite=args.overwrite,
    )
    logger.info(f"{original_count} variants were filtered to {ht.count()} variants.")
    if args.overwrite:
        logger.info(
            f"Writing TSV to {args.output_path_prefix}/v2_exomes_unique_samples_with_variants.tsv"
        )
        ht.export(
            f"{args.output_path_prefix}/v2_exomes_unique_samples_with_variants.tsv",
            header=True,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script filters a Hail Table of rare variants of interest from a random sample of individuals from each population present in gnomad exomes v.2.1.1 to variants not found in v2 genomes or the v3-non-v2 subset."
    )
    parser.add_argument(
        "--sample-with-variants-path",
        default="gs://gnomad-wphu/review-hum-mut/samples_with_variants.ht",
        type=str,
        help="Path to Hail Table with variants from gnomad v2 exomes.",
    )
    parser.add_argument(
        "--output-path-prefix",
        default="gs://gnomad-wphu/review-hum-mut",  # consider where path should go
        type=str,
        help="Google folder to store output.",
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
