import hail as hl
from gnomad.resources.grch37.gnomad import public_release as v2_public_release
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.resources.resource_utils import DataException

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)

gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)

def load_public_resources():
    '''
    Return public resources for v2 genomes, v3 genomes, and v2 liftover variants.
    '''
    v2_genomes_ht = v2_public_release("genomes").ht()
    v3_genomes = v3_public_release("genomes").ht()
    v2_liftover = hl.read_table(
        "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"
    )
    return v2_genomes, v3_genomes, v2_liftover