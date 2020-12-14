import pytest
import os
from clinvarome.clinvarome_annotation import run_pipeline


@pytest.mark.integration
def test_pipeline():
    # Check if the pipeline run complete
    vcf_file = "tests/full_data/clinvar_GRCh38_2020-11.vcf.gz"
    clinvarome = "tests/full_data/clinvarome_20201031.tsv"
    compare_gene = "tests/full_data/compare-gene_total.tsv.gz"
    compare_variant = "tests/full_data/compare-variant_total.tsv.gz"
    outdir = "tests/out"
    assert (
        run_pipeline(vcf_file, clinvarome, compare_gene, compare_variant, False, outdir)
        is True
    )

    # check result
    outfile = outdir + "/clinvar_GRCh38_2020-11_clinvarome_annotation.tsv"
    pipeline_output = os.path.isfile(outfile)
    assert pipeline_output is True
    with open(outfile, "rb") as f:
        num_lines = len(f.read().splitlines())
        assert num_lines == 4581
