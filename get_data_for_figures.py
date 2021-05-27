import random

import hail as hl
import gnomad

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
    selected_samples_hardcalls = hl.split_multi(selected_samples_hardcalls)
    return (selected_samples,selected_samples_hardcalls)

def filter_hardcalls_variants_interest(hardcalls):
    hardcalls_missense_and_synonymous = hardcalls.filter_rows(hl.set(["missense_variant", "inframe_deletion", "inframe_insertion", "synonymous_variant"]).contains(hardcalls.vep.most_severe_consequence))
    hardcalls_LOF = hardcalls.filter_rows(hl.set(["splice_acceptor", "splice_donor_variant", "stop_gained", "frameshift_variant"]).contains(hardcalls.vep.most_severe_consequence))
    hardcalls_LOF = hardcalls_LOF.filter()


def main():
    hl.init()
    random.seed(1)
    v2er_hardcalls = hl.read_matrix_table("gs://gnomad/hardcalls/hail-0.2/mt/exomes/gnomad.exomes.mt")
    v2er_hardcalls = v2er_hardcalls.filter_rows(v2er_hardcalls.info.AF[0]<0.001)
    metadata = hl.read_table("gs://gnomad/metadata/exomes/gnomad.exomes.metadata.2018-10-11.ht")
    metadata = metadata.filter(metadata.releasable_2_1) #filter to releaseable samples
    v2er = hl.read_table("gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht/")

    if (hl.hadoop_exists("gs://gnomad-tmp/review-hum-mut/random_samples.ht") and hl.hadoop_exists("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")):
        random_samples = hl.read_table("gs://gnomad-tmp/review-hum-mut/random_samples.ht")
        random_samples_hardcalls = hl.read_matrix_table("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")
    else:
        populations = ["oth", "sas", "nfe", "fin", "afr", "amr", "asj","eas"]
        random_samples, random_samples_hardcalls = get_random_samples_of_populations(populations, 100, metadata, v2er_hardcalls)
        random_samples = random_samples.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples.ht")
        random_samples_hardcalls = random_samples_hardcalls.checkpoint("gs://gnomad-tmp/review-hum-mut/random_samples_hardcalls.mt")

    #Annotate VEP and other infos
    random_samples_hardcalls = random_samples_hardcalls.annotate_rows(
                                vep = v2er[random_samples_hardcalls.locus, random_samples_hardcalls.alleles].vep,
                                popmax = v2er[random_samples_hardcalls.locus, random_samples_hardcalls.alleles].popmax,
                                freq = v2er[random_samples_hardcalls.locus, random_samples_hardcalls.alleles].freq
                            )
    #Filter to variants of interest
    
    #Filter out reference alleles
    random_samples_hardcalls = random_samples_hardcalls.annotate_rows(samples_with_variant=hl.agg.filter(random_samples_hardcalls.GT.is_non_ref(), hl.agg.collect_as_set(random_samples_hardcalls.s)))
    ht = random_samples_hardcalls.rows()
    ht = ht.select(
    "samples_with_variant",
    pop=ht.pop,
    VEP=ht.most_severe_csq, # From https://github.com/broadinstitute/gnomad_methods/blob/35066ffc01d63ac2d7b20e069ea6703013ae9da7/gnomad/utils/vep.py#L484
    group=ht.group, # You will need to add this annotation based on the groupings that Sanna wants
    variant=hl.str(ht.locus) + '-' + hl.delimit(ht.alleles, '-'),
    AC=ht.freq[0].AC,
    AF=ht.freq[0].AF,
    AN=ht.freq[0].AN
    popmax_AC=ht.popmax[0].AC,
    popmax_AF=ht.popmax[0].AF,
    popmax_AN=ht.popmax[0].AN,
)

if __name__ == "__main__":
    main()
