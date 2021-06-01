import argparse
import logging
import random

import hail as hl
from hail.expr.aggregators.aggregators import info_score
import gnomad
from gnomad import vep, filtering

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p"
)
logger = logging.getLogger("generate_figures")
logger.setLevel(logging.INFO)

def get_random_subset(n, pop, metadata):
    '''
    n: size of sample
    pop: population to select from
    metadata: the sample metadata table
    '''
    subgroup = metadata.filter(metadata.pop == pop)
    subgroup_samples = subgroup.s.collect()
    if len(subgroup_samples)>=n:
        vals = random.sample(range(len(subgroup_samples)),n)
        rand_samples = [subgroup_samples[x] for x in vals]
        subgroup = subgroup.filter(hl.set(rand_samples).contains(subgroup.s))
        return subgroup
    else:
        raise ValueError("There are fewer total samples than the requested random sample size.")

def get_hardcalls_of_samples(sample_table, hardcalls):
    '''
    sample_table: filtered metadata table filtered to specific samples
    hardcalls: hardcalls table
    '''
    selected_hardcalls = hardcalls.filter_cols(hl.set(sample_table.s.collect()).contains(hardcalls.s))
    if sample_table.count() != selected_hardcalls.count_cols():
        raise ValueError(f"The number of samples {sample_table.count()} does not equal the number of selected columns from the hardcalls table:  {selected_hardcalls.count_cols()}")
    else:
        return selected_hardcalls

    
def get_random_samples_of_populations(pops, n, metadata, hardcalls):
    '''
    pops: list of populations to select samples from
    n: number of samples to select from each population
    metadata: sample metadata table
    hardcalls: gnomad hardcalls table
    '''
    selected_samples = metadata.head(0)
    for pop in pops:
        random_samples = get_random_subset(n, pop, metadata)
        selected_samples = selected_samples.union(random_samples)
    selected_samples_hardcalls = get_hardcalls_of_samples(selected_samples,hardcalls)
    return (selected_samples,selected_samples_hardcalls)

def filter_hardcalls_variants_interest(hardcalls):
    hardcalls_missense_and_synonymous = hardcalls.filter_rows(hl.set(["missense_variant", "inframe_deletion", "inframe_insertion", "synonymous_variant"]).contains(hardcalls.vep.most_severe_consequence))
    hardcalls_LOF = hardcalls.filter_rows(hl.set(["splice_acceptor", "splice_donor_variant", "stop_gained", "frameshift_variant"]).contains(hardcalls.vep.most_severe_consequence))
    hardcalls_LOF = hardcalls_LOF.filter()


def main(args):
    hl.init()
    random.seed(1)
    v2er_hardcalls = hl.read_matrix_table("gs://gnomad/hardcalls/hail-0.2/mt/exomes/gnomad.exomes.mt")
    v2er_hardcalls = hl.split_multi(v2er_hardcalls)
    #v2er_hardcalls = v2er_hardcalls.filter_rows(v2er_hardcalls.info.AF[0]<0.001) #replace this filter with popmax AF below
    metadata = hl.read_table("gs://gnomad/metadata/exomes/gnomad.exomes.metadata.2018-10-11.ht")
    metadata = metadata.filter(metadata.release) #filter to releaseable samples
    v2er = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht/")
    v2er_indexed = v2er[v2er_hardcalls.row_key]
    v2er_hardcalls = v2er_hardcalls.annotate_rows(
        vep = v2er_indexed.vep,
        freq = v2er_indexed.freq,
        popmax = v2er_indexed.popmax
    )

    v2er_hardcalls = v2er_hardcalls.filter_rows(v2er_hardcalls.popmax[0].AF<0.001)
    populations = ["oth", "sas", "nfe", "fin", "afr", "amr", "asj","eas"]

    if args.checkpoint:
        if (hl.hadoop_exists("gs://gnomad-tmp/review-hum-mut/random_samples.ht") and hl.hadoop_exists("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")):
            random_samples = hl.read_table("gs://gnomad-tmp/review-hum-mut/random_samples.ht")
            random_samples_hardcalls = hl.read_matrix_table("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")
        else:
            random_samples, random_samples_hardcalls = get_random_samples_of_populations(populations, 100, metadata, v2er_hardcalls)
            random_samples = random_samples.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples.ht")
            random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")
    else:
        random_samples, random_samples_hardcalls = get_random_samples_of_populations(populations, 100, metadata, v2er_hardcalls)
        random_samples = random_samples.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples.ht", overwrite= True)
        random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt", overwrite=True)

    #Filter to variants of interest

    if (args.test):
        random_samples_hardcalls = random_samples_hardcalls._filter_partitions(range(2))
    
    random_samples_hardcalls = random_samples_hardcalls.filter_rows(
        (~hl.is_defined(random_samples_hardcalls.filters) | (hl.len(random_samples_hardcalls.filters) == 0))
        & (hl.agg.any(random_samples_hardcalls.GT.is_non_ref()))
    )
    random_samples_hardcalls  = filtering.filter_low_conf_regions(random_samples_hardcalls)
    random_samples_hardcalls = vep.filter_vep_to_canonical_transcripts(random_samples_hardcalls)
    most_severe_csq_summary = vep.get_most_severe_consequence_for_summary(random_samples_hardcalls.rows())[random_samples_hardcalls.row_key]
    random_samples_hardcalls = random_samples_hardcalls.annotate_rows(
        most_severe_csq = most_severe_csq_summary.most_severe_csq,
        protein_coding = most_severe_csq_summary.protein_coding,
        lof = most_severe_csq_summary.lof,
        no_lof_flags = most_severe_csq_summary.no_lof_flags
    ) #get_most_severe_consequence_for_summary works only on tables
    random_samples_hardcalls = filtering.filter_to_adj(random_samples_hardcalls)
    random_samples_hardcalls = random_samples_hardcalls.filter_rows(hl.agg.any(random_samples_hardcalls.GT.is_non_ref())) #possibly repeated above?
    random_samples_hardcalls_LOF = random_samples_hardcalls.filter_rows(
        (random_samples_hardcalls.lof == "HC") &
        (hl.set(["splice_acceptor", "splice_donor_variant", "stop_gained", "frameshift_variant"]).contains(random_samples_hardcalls.most_severe_csq))
    )
    random_samples_hardcalls_LOF = random_samples_hardcalls_LOF.annotate_rows(group = "LoF")
    random_samples_hardcalls_missense = random_samples_hardcalls.filter_rows(hl.set(["missense_variant", "inframe_deletion", "inframe_insertion"]).contains(random_samples_hardcalls.most_severe_csq))
    random_samples_hardcalls_missense = random_samples_hardcalls_missense.annotate_rows(group = "missense_indels")
    random_samples_hardcalls_synonymous = random_samples_hardcalls.filter_rows(hl.set(["synonymous_variant"]).contains(random_samples_hardcalls.most_severe_csq))
    random_samples_hardcalls_synonymous = random_samples_hardcalls_synonymous.annotate_rows(group = "synonymous")
    random_samples_hardcalls = random_samples_hardcalls_LOF.union_rows(random_samples_hardcalls_missense, random_samples_hardcalls_synonymous)
    random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls_filtered.mt", overwrite=True)

    #random_samples_hardcalls = random_samples_hardcalls.filter_rows((hl.len(random_samples_hardcalls.filters) == 0) & hl.agg.any(random_samples_hardcalls.GT.is_non_ref()))
    #random_samples_hardcalls  = vep.filter_low_conf_regions(random_samples_hardcalls)
    #random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")
    #random_samples_hardcalls = vep.filter_vep_to_canonical_transcripts(random_samples_hardcalls)
    #random_samples_hardcalls = vep.get_most_severe_consequence_for_summary(random_samples_hardcalls)
    #random_samples_hardcalls = vep.filter_to_adj(random_samples_hardcalls)
    #random_samples_hardcalls = random_samples_hardcalls.filter_rows(hl.agg.any(random_samples_hardcalls.GT.is_non_ref())) #possibly repeated above?
    #random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls-filtered.mt")

    #GT.is_non_ref() possibly repeated again?
    random_samples_hardcalls = random_samples_hardcalls.annotate_rows(samples_with_variant=hl.agg.filter(random_samples_hardcalls.GT.is_non_ref(), hl.agg.collect_as_set(random_samples_hardcalls.s)))
    ht = random_samples_hardcalls.rows()
    ht = ht.select(
        "samples_with_variant",
        VEP=ht.most_severe_csq, # From https://github.com/broadinstitute/gnomad_methods/blob/35066ffc01d63ac2d7b20e069ea6703013ae9da7/gnomad/utils/vep.py#L484
        group=ht.group, # You will need to add this annotation based on the groupings that Sanna wants
        variant=hl.str(ht.locus) + '-' + hl.delimit(ht.alleles, '-'),
        AC=ht.freq[0].AC,
        AF=ht.freq[0].AF,
        AN=ht.freq[0].AN,
        popmax_AC=ht.popmax[0].AC,
        popmax_AF=ht.popmax[0].AF,
        popmax_AN=ht.popmax[0].AN,
    )
    ht = ht.explode("samples_with_variant")
    random_samples_pop_map = hl.dict(hl.zip(random_samples.s.collect(),random_samples.pop.collect()))
    ht = ht.annotate(pop = random_samples_pop_map[ht.samples_with_variant]) #does not require a shuffle this way
    ht = ht.select(
        ht.samples_with_variant,
        ht.pop,
        ht.VEP,
        ht.group,
        ht.variant,
        ht.AC,
        ht.AF,
        ht.AN,
        ht.popmax_AC,
        ht.popmax_AF,
        ht.popmax_AN
    )
    '''
    ht_summary= ht.group_by(ht.samples_with_variant).aggregate(
        num_Lof = hl.agg.count_where(ht.group=="LoF"),
        num_missense_indel = hl.agg.count_where(ht.group=="missense_indels"),
        num_synonymous = hl.agg.count_where(ht.group=="synonymous"),
        num_total = hl.agg.count()
    )
    '''
    #ht.order_by(ht.pop,ht.samples_with_variant,ht.group)
    #ht.group_by()

    random_samples.write("gs://gnomad-wphu/review-hum-mut/random_samples_metadata.ht", overwrite=True)
    random_samples_hardcalls.write("gs://gnomad-wphu/review-hum-mut/random_samples_hardcalls_filtered.mt", overwrite=True)
    ht.write("gs://gnomad-wphu/review-hum-mut/samples_with_variants.ht", overwrite=True)
    ht.export("gs://gnomad-wphu/review-hum-mut/samples_with_variants.tsv.bgz", header=True)
    #ht_summary.write("gs://gnomad-wphu/review-hum-mut/samples_with_variants_summary.ht", overwrite= True)
    #ht_summary.export("gs://gnomad-wphu/review-hum-mut/samples_with_variants_summary.tsv.bgz", header=True)

    logger.info("Wrote out table with %s rows.", ht.count())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "This script generates all of the rare variants of interest from a random sample of individuals from each population present in gnomad exomes v.2.1.1"
    )
    parser.add_argument(
        "--checkpoint",
        action="store_true",
        help="Use previously checkpointed random sample if it exists.")
    parser.add_argument(
        "--test",
        action="store_true",
        help='subset hardcalls table to a small number of partitions for testing.'
    )
    args = parser.parse_args()

    main(args)
