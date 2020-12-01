#!/usr/bin/env python

import argparse
import os
from clinvarome.clinvarome_annotation_functions import gather_annotation


def get_args():
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(
        description="Add additional annotations to the ClinVarome."
    )
    parser.add_argument(
        "--vcf",
        action="store",
        dest="vcf_file",
        help="vcf input from clinVCF",
        required=True,
    )
    parser.add_argument(
        "--clinvarome",
        action="store",
        dest="clinvarome",
        required=True,
        help="ClinVarome file from Variant Alert! to be annotated.",
    )
    parser.add_argument(
        "--compare-gene",
        action="store",
        dest="compare_gene",
        required=True,
        help="Variant Alert! compare gene concatenation gzip file of all previous ClinVar release.",
    )
    parser.add_argument(
        "--compare-variant",
        action="store",
        dest="compare_variant",
        required=True,
        help="Variant Alert! compare variant concatenation gzip file of all previous ClinVar release.",
    )
    parser.add_argument(
        "--gnomad",
        action="store_true",
        help="Add max pathogenic gnomAD filtered allele frequency for each gene, if clinVCF is annotated with gnomAD v3 by vcfanno.",
    )
    parser.add_argument(
        "--output-dir",
        action="store",
        dest="output_dir",
        required=True,
        help="Path to the output folder.",
    )

    args = parser.parse_args()
    return args


def run_pipeline(
    vcf_file: str,
    clinvarome: str,
    compare_gene: str,
    compare_variant: str,
    gnomad: bool,
    output_dir: str,
):

    clinvarome_file_output = os.path.splitext(os.path.basename(vcf_file))[0]

    clinvarome_annotation = gather_annotation(
        vcf_file,
        clinvarome,
        compare_gene,
        compare_variant,
        gnomad,
    )

    clinvarome_annotation.to_csv(
        os.path.join(
            str(output_dir), clinvarome_file_output + "_clinvarome_annotation.tsv"
        ),
        sep="\t",
        index=False,
    )

    return True


if __name__ == "__main__":
    # Parse arguments
    args = get_args()

    # Run pipeline
    run_pipeline(
        args.vcf_file,
        args.clinvarome,
        args.compare_gene,
        args.compare_variant,
        args.gnomad,
        args.output_dir,
    )
