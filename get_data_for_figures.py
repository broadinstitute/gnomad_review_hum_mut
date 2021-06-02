import argparse
import logging
import random

import hail as hl
from gnomad.resources.grch37.gnomad import EXOME_POPS
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import filter_low_conf_regions, filter_to_adj
from gnomad.utils.vep import (
    filter_vep_to_canonical_transcripts,
    get_most_severe_consequence_for_summary,
)
EXOME_POPS = [pop.lower() for pop in EXOME_POPS]
LOF_VEP = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
}
MISSENSE_INDEL_VEP = {"missense_variant", "inframe_insertion", "inframe_deletion"}
SYNONYMOUS_VEP = {"synonymous_variant"}

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
logger = logging.getLogger("generate_figures")
logger.setLevel(logging.INFO)

def get_random_subset(meta_ht, n, pop):
    """
    Return a random subset of `n` samples from population `pop` defined by pop in the `meta_ht`.
    
    :param meta_ht: Sample metadata table
    :param n: Size of random sample sample
    :param pop: population to select from
    """
    subpop_ht = meta_ht.filter(meta_ht.pop == pop)
    subpop_samples = subpop_ht.s.collect()

    if len(subpop_samples) >= n:
        vals = random.sample(range(len(subpop_samples)), n)
        rand_samples = [subpop_samples[x] for x in vals]
        return rand_samples
    else:
        raise ValueError(f"There are fewer total samples than the requested random sample size in the population: {pop}.")

def get_hardcalls_of_samples(mt, samples):
    """
    Filter a MatrixTable to a specific set of samples found in `samples`.
    :param mt: MatrixTable to filter
    :param samples: Set of specific samples in `mt` to filter to
    """
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
    if len(samples) != mt.count_cols():
        raise ValueError(
            f"The number of samples {len(samples)} does not equal the number of selected columns from the supplied matrix table:  {mt.count_cols()}"
        )
    else:
        return mt

    
def get_random_samples_of_populations(mt, meta_ht, pops, n):
    """
    Get a random sample of `n` columns in `mt` from each population in `pops`.
    :param mt: Input MatrixTable
    :param meta_ht: Sample metadata Table
    :param pops: List of populations to select samples from
    :param n: Number of samples to select from each population
    """
    selected_samples = set([])
    for pop in pops:
        random_samples = get_random_subset(meta_ht, n, pop)
        selected_samples = selected_samples | set(random_samples)
    mt = get_hardcalls_of_samples(mt, selected_samples)
    meta_ht = meta_ht.filter(hl.literal(selected_samples).contains(meta_ht.s))
    return meta_ht, mt

def filter_hardcalls_variants_interest(mt):
    """
    Annotate variants with a `group` annotations and filter to only those in groups of interest (missense, synonymous, and pLOF).
    The following annotations are used to define the variants of interest:
        - pLOF (splice disrupting, nonsense variants, frameshift variants) that are predicted high-confidence (HC) by LOFTEE. VEP annotations: splice_acceptor, splice_donor_variant, stop_gained, frameshift_variant (that pass LOFTEE HC)
        - Missense variants and indels. VEP annotations: missense_variant, inframe_insertion, inframe_deletion
        - Synonymous variants. VEP annotations: synonymous_variant
    Must include the `most_severe_csq` and `lof` annotations on `mt`.
    :param mt: Input MatrixTable
    """
    mt = mt.annotate_rows(
        group=hl.case()
        .when(hl.literal(LOF_VEP).contains(mt.most_severe_csq) & (mt.lof == "HC"), "LoF")
        .when(hl.literal(MISSENSE_INDEL_VEP).contains(mt.most_severe_csq), "missense_indels")
        .when(hl.literal(SYNONYMOUS_VEP).contains(mt.most_severe_csq), "synonymous")
        .or_missing()
    )
    return mt.filter_rows(hl.is_defined(mt.group))


def main(args):
    hl.init(log="/gnomad_review_hum_mut.log")
    random.seed(1)
    random_samples_path = f"gs://gnomad-tmp/review-hum-mut/random_samples{'_test' if args.test else ''}.ht"
    random_samples_hardcalls_path = f"gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls{'_test' if args.test else ''}.mt"
    
    if args.use_checkpoint:
        if (file_exists(random_samples_path) and file_exists(random_samples_hardcalls_path)):
            meta_ht = hl.read_table(random_samples_path)
            mt = hl.read_matrix_table(random_samples_hardcalls_path)
        else:
            raise DataException("There is currently no checkpointed files for this dataset.")
    else:
        logger.info("Reading in gnomAD v2.1.1 exome hardcalls MatrixTable...")
        mt = hl.read_matrix_table(
            "gs://gnomad/hardcalls/hail-0.2/mt/exomes/gnomad.exomes.mt"
        )
        if args.test:
            mt = mt._filter_partitions(range(args.test_n_partitions))
        logger.info(
            "Reading in the gnomAD v2.1.1 metadata HailTable from 2018-10-11 and filtering to only release samples (releasable and pass sample QC)...")
        meta_ht = hl.read_table(
            "gs://gnomad/metadata/exomes/gnomad.exomes.metadata.2018-10-11.ht"
        )
        # Filter metadata to releasable samples
        meta_ht = meta_ht.filter(meta_ht.release)
        logger.info(
            "Reading in the gnomAD v2.1.1 release sites Hail Table to annotate with VEP, freq, popmax, and variant QC filters...")
        ht = hl.read_table(
            "gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"
        )
        ht_indexed = ht[mt.row_key]
        mt = mt.annotate_rows(
            filters=ht_indexed.filters,
            vep=ht_indexed.vep,
            freq=ht_indexed.freq,
            popmax=ht_indexed.popmax,
        )
        # Filter to only variants with a popmax allele frequency of 0.1%
        mt = mt.filter_rows(mt.popmax[0].AF < 0.001)
        meta_ht, mt = get_random_samples_of_populations(mt, meta_ht, EXOME_POPS, 100)
        meta_ht = meta_ht.checkpoint(random_samples_path, overwrite=args.overwrite)
        mt = mt.checkpoint(random_samples_hardcalls_path, overwrite=args.overwrite)
    
    logger.info("Filtering to PASS variants present in randomly sampled individuals, removing low confidence regions, and filtering VEP to canonical transcripts only...")
    mt = mt.filter_rows(
        (hl.is_defined(mt.filters) & (hl.len(mt.filters) == 0))
        & (hl.agg.any(mt.GT.is_non_ref()))
    )
    mt = filter_low_conf_regions(mt)
    mt = filter_vep_to_canonical_transcripts(mt)
    logger.info("Getting the most severe consequence from the VEP annotation of the canonical transcript...")
    most_severe_csq_summary = get_most_severe_consequence_for_summary(mt.rows())[mt.row_key]
    mt = mt.annotate_rows(
        most_severe_csq=most_severe_csq_summary.most_severe_csq,
        protein_coding=most_severe_csq_summary.protein_coding,
        lof=most_severe_csq_summary.lof,
        no_lof_flags=most_severe_csq_summary.no_lof_flags,
    )  # get_most_severe_consequence_for_summary works only on tables

    logger.info("Filtering genotypes to adj and the matrix table to variants with at least one non ref after adj filtering...")
    mt = filter_to_adj(mt)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    logger.info("Annotating and filtering the MT to only to variants with a VEP consequence of interest...")
    mt = filter_hardcalls_variants_interest(mt)

    mt = mt.annotate_rows(samples_with_variant=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect_as_set(mt.s)))
    ht = mt.rows()
    ht = ht.select(
        "samples_with_variant",
        VEP=ht.most_severe_csq,
        group=ht.group,
        variant=hl.str(ht.locus) + '-' + hl.delimit(ht.alleles, '-'),
        AC=ht.freq[0].AC,
        AF=ht.freq[0].AF,
        AN=ht.freq[0].AN,
        popmax_AC=ht.popmax[0].AC,
        popmax_AF=ht.popmax[0].AF,
        popmax_AN=ht.popmax[0].AN,
    )
    ht = ht.explode("samples_with_variant")
    random_samples_pop_map = hl.dict(hl.zip(meta_ht.s.collect(), meta_ht.pop.collect()))
    ht = ht.annotate(pop=random_samples_pop_map[ht.samples_with_variant])  # does not require a shuffle this way

    meta_ht.write(f"gs://gnomad-wphu/review-hum-mut/random_samples_metadata{'_test' if args.test else ''}.ht", overwrite=args.overwrite)
    mt.write(f"gs://gnomad-wphu/review-hum-mut/random_samples_hardcalls_filtered{'_test' if args.test else ''}.mt", overwrite=args.overwrite)
    ht.write(f"gs://gnomad-wphu/review-hum-mut/samples_with_variants{'_test' if args.test else ''}.ht", overwrite=args.overwrite)
    if args.overwrite:
        ht.export(f"gs://gnomad-wphu/review-hum-mut/samples_with_variants{'_test' if args.test else ''}.tsv", header=True)


    logger.info("Wrote out table with %s rows.", ht.count())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script generates all of the rare variants of interest from a random sample of individuals from each population present in gnomad exomes v.2.1.1"
    )
    parser.add_argument(
        "--use_checkpoint",
        action="store_true",
        help="Use previously checkpointed random sample if it exists.")
    parser.add_argument(
        "--test",
        action="store_true",
        help='subset hardcalls table to a small number of partitions for testing.'
    )
    parser.add_argument(
        "--test_n_partitions",
        default=5,
        type=int,
        help='Number of partitions to use for testing.'
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
