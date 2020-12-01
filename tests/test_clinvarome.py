#############################
# Import
#############################

import pysam
import pandas as pd

from clinvarome.clinvarome_annotation_functions import (
    gather_clinical_features,
    get_clinical_dataframe,
    calcul_max_AF,
    gather_dict_gene_max_AF,
    get_AF_max_by_gene,
    mol_consequences_by_variant,
    count_type_mol_consequences,
    get_mol_consequences_dataframe,
    get_nhomalt,
    get_max_nhomalt_by_gene,
    gene_first_pathogenic_entry_date,
    gene_latest_pathogenic_entry_date,
)

CLIN_VCF = "tests/full_data/clinvar_GRCh38_2020-11.vcf.gz"
CLIN_VCF_GNOMAD = "tests/toys/clinvar_gnomad_anno.vcf"
COMPARE_GENE = "tests/toys/compare-gene.tsv.gz"
COMPARE_VARIANT = "tests/toys/compare-variant.tsv.gz"
CLINVAROME = "tests/toys/clinvarome.tsv.gz"

VCF = pysam.VariantFile(CLIN_VCF)
nb_record = 1
RECORDS = {}
for record in VCF:
    if nb_record == 0:
        break
    if (
        not ("CLN" in RECORDS.keys())
        and "CLNDISEASE" in record.info
        and "CLNFINDING" in record.info
    ):
        RECORDS["CLN"] = record
        nb_record -= 1
VCF.close()

VCF = pysam.VariantFile(CLIN_VCF_GNOMAD)
nb_record = 1
for record in VCF:
    if nb_record == 0:
        break
    if (
        "CLNDISEASE" in record.info
        and "CLNFINDING" in record.info
        and "nhomalt" in record.info
        and record.info["nhomalt"][0] > 1
    ):
        RECORDS["nhomalt"] = record
        nb_record -= 1

for k in RECORDS.keys():
    print(k)
    print(RECORDS[k])


def test_gather_clinical_features():
    gene_disease = {}
    gene_finding = {}
    gather_clinical_features(RECORDS["CLN"], gene_finding, gene_disease)
    assert gene_disease["AGRN"][0] == ["myasthenic_syndrome_congenital_8"]
    assert gene_finding["AGRN"][0] == [
        "anxiety",
        "depressivity",
        "short_attention_span",
        "abnormality_of_muscle_physiology",
    ]


def test_get_clinical_dataframe():
    gene_disease = {}
    gene_disease["AGRN"] = []
    gene_disease["AGRN"].append(["myasthenic_syndrome_congenital_8"])

    gene_finding = {}
    gene_finding["AGRN"] = []
    gene_finding["AGRN"].append(
        [
            "anxiety",
            "depressivity",
            "short_attention_span",
            "abnormality_of_muscle_physiology",
        ]
    )
    df = get_clinical_dataframe(gene_disease, gene_finding)
    print(df)
    assert df["gene_info"][0] == "AGRN"
    assert df["clinical_disease"][0] == "myasthenic_syndrome_congenital_8"
    assert (
        df["clinical_finding"][0]
        == "abnormality_of_muscle_physiology;anxiety;depressivity;short_attention_span"
    )


def test_calcul_max_AF():
    assert calcul_max_AF(0, 0) == 0
    assert calcul_max_AF(0.5, 0) == 0
    assert calcul_max_AF(0, 0.5) == 2
    assert calcul_max_AF(0.75, 0.5) == 4


def test_gather_dict_gene_max_AF():
    gene_AF_pois_dict = {}
    gather_dict_gene_max_AF(RECORDS["CLN"], gene_AF_pois_dict)
    gather_dict_gene_max_AF(RECORDS["CLN"], gene_AF_pois_dict)
    assert gene_AF_pois_dict["AGRN"] == [0, 0]


def test_get_AF_max_by_gene():
    gene_AF_dict = {"GENE": [0.5, 0.3, 0.75]}
    df = get_AF_max_by_gene(gene_AF_dict)
    assert df["gene_info"][0] == "GENE"
    assert df["FAF"][0] == 0.75


def test_mol_consequences_by_variant():
    gene_var_dict = {}
    mol_consequences_by_variant(RECORDS["CLN"], gene_var_dict)
    assert gene_var_dict["AGRN"] == ["missense_inframe"]


# WARNING the list have an extra layer  !!! can corrupdate  count_type_meca and get_MC_dataframe
# mecanism   missense_inframe  stop_fs_splice
# gene_info
# GENE                      3               1
def test_count_type_mol_consequences():
    gene_var_dict = {
        "GENE": [
            "missense_inframe",
            "missense_inframe",
            "missense_inframe",
            "stop_fs_splice",
        ]
    }
    gene_meca_count = count_type_mol_consequences(gene_var_dict)
    print(gene_meca_count)
    assert gene_meca_count["GENE"] == [[[1, "stop_fs_splice"], [3, "missense_inframe"]]]
    # For the dataframe
    df = get_mol_consequences_dataframe(gene_var_dict)
    assert df["stop_fs_splice"]["GENE"] == 1
    assert df["missense_inframe"]["GENE"] == 3


def test_get_nhomalt():
    gene_nhomalt = {}
    get_nhomalt(RECORDS["nhomalt"], gene_nhomalt)
    assert gene_nhomalt["CCDC28B"] == [22]


def test_get_max_nhomalt_by_gene():
    gene_nhomalt = {"GENE": [3, 5, 20, 10]}
    df = get_max_nhomalt_by_gene(gene_nhomalt)
    assert df["gene_info"][0] == "GENE"
    assert df["nhomalt"][0] == 20


def test_first_gene_date():
    clinvarome_df = pd.read_csv(CLINVAROME, sep="\t")
    df = gene_first_pathogenic_entry_date(clinvarome_df, COMPARE_GENE)
    assert df["gene_info"][0] == "TUBB2B"
    assert str(df["first_path_var_date"][0]) == "2017-07-03 00:00:00"


def test_last_gene_date():
    clinvarome_df = pd.read_csv(CLINVAROME, sep="\t")
    df = gene_latest_pathogenic_entry_date(clinvarome_df, COMPARE_VARIANT)
    # print(df)
    assert df["gene_info"][3] == "TUBB2B"
    assert str(df["last_pathogenic_variant"][3]) == "2018-12-31 00:00:00"


# def date_management(clinvarome_df):
# def scale_clinvarome_data(clinvarome_annotation, ARRAY_TRANSFORM):
# def give_name_to_cluster(cluster_df):

# Untested function
# Those function are tested as they require a large un data - the result of this unction will be test during the functional test
# clusterize_clinvarome => require full data set better assessed during functionalTest
# clinvarome_annotation => merging of multiple table, test during functional test should show probelm at this step
# number_of_variants => add additional feature to the annotaed_clinvarome data frame, test during functional test should show problem at this step
