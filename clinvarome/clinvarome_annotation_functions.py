# ClinVarome annotation functions
# Gather all genes annotations : gene, gene_id,
# (AF, FAF,) diseases, clinical features, mecanismes counts, nhomalt.
# Give score for genes according their confidence criteria
import gzip
import pandas as pd
import numpy as np
from datetime import datetime
import pysam
from scipy.stats import poisson
from pyselib.logger.logger import get_logger
from sklearn.preprocessing import QuantileTransformer
from sklearn.cluster import AgglomerativeClustering
from clinvarome.utils.dictionary import (
    EFFECTS,
    MC_CATEGORIES,
    MC_SHORT,
    ARRAY_TRANSFORM,
    CLUSTER_NAMES,
)

logger = get_logger(__name__)

# Clinical features


def gather_clinical_features(record, gene_finding, gene_disease):
    """
    update gene_finding and gene_disease dictionary using information from a VCF record
    """
    geneinfo = record.info["GENEINFO"].split("|")[0].split(":")[0]
    if "CLNDISEASE" in record.info:
        clndisease = record.info["CLNDISEASE"][0].split("|")
        gene_disease.setdefault(geneinfo, [])
        gene_disease[geneinfo].append(clndisease)
    if "CLNFINDING" in record.info:
        clnfinding = record.info["CLNFINDING"][0].split("|")
        gene_finding.setdefault(geneinfo, [])
        gene_finding[geneinfo].append(clnfinding)


def get_clinical_dataframe(gene_disease, gene_finding):
    """
    Process dictionary output from gather_clinical_features function
    into a dataframe
    """
    for key, value in gene_disease.items():
        flat_list = [j for i in value for j in i]
        gene_disease[key] = ";".join(sorted(list(set(flat_list))))
    gene_disease_df = pd.DataFrame(
        gene_disease.items(), columns=["gene_info", "clinical_disease"]
    )

    for key, value in gene_finding.items():
        flat_list = [j for i in value for j in i]
        gene_finding[key] = ";".join(sorted(list(set(flat_list))))
    gene_finding_df = pd.DataFrame(
        gene_finding.items(), columns=["gene_info", "clinical_finding"]
    )

    gene_features = gene_disease_df.merge(gene_finding_df, how="outer")
    return gene_features


# FAF


def calcul_max_AF(AC, AN):
    """
    For a given AC and AN, calcul the maximum AF: the
    upper bound of the Poisson 95 % CI.
    """
    if (AC == 0) and (AN != 0):
        max_AF_pois = 1 / AN
    elif (AC != 0) and (AN != 0):
        max_AC_pois = poisson.ppf(0.95, AC)
        max_AF_pois = float(max_AC_pois / AN)
    else:
        max_AF_pois = 0
    return max_AF_pois


def gather_dict_gene_max_AF(record, gene_AF_pois_dict):
    """
    Update the maximum FAF of a gene using information in a VCF record
    """
    ls_AC = []
    ls_AN = []
    ls_AF_pois = []

    geneinfo = record.info["GENEINFO"].split("|")[0].split(":")[0]
    gene_AF_pois_dict.setdefault(geneinfo, [])
    if "AC_afr" in record.info:
        AC_afr = record.info["AC_afr"]
        AC_amr = record.info["AC_amr"]
        AC_nfe = record.info["AC_nfe"]
        AC_eas = record.info["AC_eas"]
        AN_afr = record.info["AN_afr"]
        AN_amr = record.info["AN_amr"]
        AN_nfe = record.info["AN_nfe"]
        AN_eas = record.info["AN_eas"]
        ls_AC = [AC_afr, AC_amr, AC_nfe, AC_eas]
        ls_AN = [AN_afr, AN_amr, AN_nfe, AN_eas]
        for k in range(0, len(ls_AC)):
            ls_AF_pois.append(calcul_max_AF(ls_AC[k], ls_AN[k]))
        max_af_pois = max(ls_AF_pois)
        gene_AF_pois_dict[geneinfo].append(max_af_pois)
    else:
        gene_AF_pois_dict[geneinfo].append(0)


def get_AF_max_by_gene(gene_AF_dict):
    """For a given gene, return the maximum FAF (among its variants)
    and get a dataframe."""
    gene_AF_max = {}
    for key, values in gene_AF_dict.items():
        gene_max_AF = max(values)
        gene_AF_max.setdefault(key, [])
        gene_AF_max[key].append(gene_max_AF)
    print(gene_AF_max)
    gene_anno_pois = pd.DataFrame.from_dict(
        gene_AF_max, orient="index", columns=["FAF"]
    )
    gene_anno_pois = gene_anno_pois.reset_index()
    gene_anno_pois = gene_anno_pois.rename(columns={"index": "gene_info"})
    print(gene_anno_pois)
    return gene_anno_pois


# Molecular consequence counts
def meca_by_variant(record, gene_var_dict):
    """
    Parse molecular consequences available for a variant and
    return the highest predicted effect
    """
    geneinfo = record.info["GENEINFO"].split("|")[0].split(":")[0]
    gene_var_dict.setdefault(geneinfo, [])
    if "MC" in record.info:
        mc = record.info["MC"]
        mc_only = [i.split("|")[1] for i in mc]
        min_value = min([v for k, v in EFFECTS.items() if k in mc_only])

        for key, value in EFFECTS.items():
            if min_value == value:
                gene_var_dict[geneinfo].append(MC_CATEGORIES[key])
                break
    else:
        gene_var_dict[geneinfo].append("Not_provided")


def count_type_meca(gene_var_dict):
    """
    Count occurence of molecular consequence from pathogenic
    variant for each gene
    """
    gene_meca_count = {}
    for key, values in gene_var_dict.items():
        list_MC = []
        for k in MC_SHORT.keys():
            if k in values:
                count = values.count(k)
                list_MC.append([count, k])
        gene_meca_count.setdefault(key, [])
        gene_meca_count[key].append(list_MC)
    return gene_meca_count


def get_MC_dataframe(gene_var_dict):
    """
    Format molecular consequences occurences by gene dictionary into dataframe.
    """
    gene_meca_count = count_type_meca(gene_var_dict)
    df_tot = pd.DataFrame()
    for key, values in gene_meca_count.items():
        for k in range(len(values[0])):
            mecanism_dict = {}
            mecanism_dict[key] = values[0][k]
            df = pd.DataFrame.from_dict(
                mecanism_dict, orient="index", columns=["count", "mecanism"]
            )
            df_tot = df_tot.append(df)
    df_tot.index.name = "gene_info"
    df_tot_piv = pd.pivot_table(
        df_tot,
        values="count",
        index="gene_info",
        columns=["mecanism"],
        fill_value=0,
    )

    return df_tot_piv


# nhomalt annotation


def get_nhomalt(record, gene_nhomalt):
    """
    Return count of homozygous allele in gnomad for a pathogenic variant.
    """
    if "nhomalt" in record.info:
        nhomalt = record.info["nhomalt"][0]
        geneinfo = record.info["GENEINFO"].split("|")[0].split(":")[0]
        gene_nhomalt.setdefault(geneinfo, [])
        gene_nhomalt[geneinfo].append(nhomalt)
    return gene_nhomalt


def get_max_nhomalt_by_gene(gene_nhomalt):
    """
    Get the maximum count of homozygous pathogenic allele in gnomad by gene.
    Return a dataframe.
    """
    gene_nhomalt_max = {}
    for key, values in gene_nhomalt.items():
        nhomalt_max = max(values)
        gene_nhomalt_max.setdefault(key, [])
        gene_nhomalt_max[key].append(nhomalt_max)
    gene_nhomalt_max_df = pd.DataFrame.from_dict(
        gene_nhomalt_max, orient="index", columns=["nhomalt"]
    )
    gene_nhomalt_max_df = gene_nhomalt_max_df.reset_index()
    gene_nhomalt_max_df = gene_nhomalt_max_df.rename(columns={"index": "gene_info"})
    return gene_nhomalt_max_df


# Gene date
def first_gene_date(clinvarome_df, compare_gene):
    """
    Return the first occurence of pathogenic variant for a gene in ClinVar.
    """
    compare_gene_df = pd.read_csv(compare_gene, sep="\t", compression="gzip")
    compare_gene_df = compare_gene_df[
        compare_gene_df["pathogenic_class_status"] == "NEW_PATHOGENICITY"
    ]
    compare_gene_df = compare_gene_df.sort_values(by="name_clinvar_new", ascending=True)
    compare_gene_df.drop_duplicates("gene_info", inplace=True)
    clinvar_gene = clinvarome_df[["gene_info"]].merge(
        compare_gene_df[["gene_info", "name_clinvar_new"]],
        on="gene_info",
        how="outer",
    )
    clinvar_gene.fillna(value="20170703", inplace=True)
    clinvar_gene["name_clinvar_new"] = pd.to_datetime(
        clinvar_gene["name_clinvar_new"], format="%Y%m%d"
    )
    clinvar_gene.rename(
        columns={"name_clinvar_new": "first_path_var_date"}, inplace=True
    )
    return clinvar_gene


def last_gene_date(clinvarome_df, compare_variant):
    """
    Return the last occurence of pathogenic variant for a gene in ClinVar.
    """
    compare_variant_df = pd.read_csv(compare_variant, sep="\t", compression="gzip")
    filter_patho = (compare_variant_df["breaking_change"] == "major") & (
        (compare_variant_df["new_classification"] == "Pathogenic")
        | (compare_variant_df["new_classification"] == "Likely_pathogenic")
        | (compare_variant_df["new_classification"] == "Pathogenic/Likely_pathogenic")
    )
    compare_variant_df = compare_variant_df[filter_patho]
    clinvarome_compare = clinvarome_df[["gene_info"]].merge(
        compare_variant_df[["gene_info", "name_clinvar_new"]],
        on="gene_info",
        how="outer",
    )
    clinvarome_compare.fillna(value="20170703", inplace=True)
    clinvarome_compare["name_clinvar_new"] = pd.to_datetime(
        clinvarome_compare["name_clinvar_new"], format="%Y%m%d"
    )

    clinvarome_compare = clinvarome_compare.sort_values(
        by="name_clinvar_new", ascending=False
    )
    clinvarome_compare.drop_duplicates("gene_info", inplace=True)

    clinvarome_compare.rename(
        columns={"name_clinvar_new": "last_pathogenic_variant"}, inplace=True
    )
    return clinvarome_compare


# Merge results
def merge_dataframe(clinvarome_df, list_df):
    """
    Add all ClinVarome additional informations into
    the original ClinVarome available in Variant Alert!
    """
    for df in range(len(list_df)):
        clinvarome_df = clinvarome_df.merge(list_df[df], on="gene_info", how="outer")
        clinvarome_df.dropna(subset=["gene_info_id"], inplace=True)
    return clinvarome_df


# number of variants for one gene


def number_of_variants(annoted_clinvarome):
    """
    take the clinvarome and count the number of
    molecular consequences reflecting the number of variants
    by genes.
    """
    annoted_clinvarome["variant_number"] = (
        annoted_clinvarome["missense_inframe"]
        + annoted_clinvarome["other"]
        + annoted_clinvarome["stop_fs_splice"]
        + annoted_clinvarome["Not_provided"]
    )
    return annoted_clinvarome


# calcul date between the first and the last submission but also age of genes


def date_management(clinvarome_df):
    """
    Add the age of gene and the number of months between
    the first and the last submission.

    If date between first date is negative,
    then probably a variant changes of gene name.

    Gene name warning column:
    If there is a gene name warning, please verify
    if the gene is OK.
    """
    # Time output
    clinvarome_df["name_clinvar_new"] = pd.to_datetime(
        clinvarome_df["name_clinvar_new"], format="%Y%m%d"
    )
    # Get time between first and last pathogenic variant
    clinvarome_df["first_path_var_date"] = pd.to_datetime(
        clinvarome_df["first_path_var_date"], format="%Y-%m-%d"
    )
    clinvarome_df["last_pathogenic_variant"] = pd.to_datetime(
        clinvarome_df["last_pathogenic_variant"], format="%Y-%m-%d"
    )
    clinvarome_df["date_between_first_date"] = (
        clinvarome_df["last_pathogenic_variant"] - clinvarome_df["first_path_var_date"]
    ).dt.days
    clinvarome_df["date_between_first_date"] = (
        clinvarome_df["date_between_first_date"] // 30
    )
    # set gene_name_check warning column
    clinvarome_df["gene_name_check"] = np.where(
        clinvarome_df["date_between_first_date"] < 0, "warning", "pass"
    )
    # set default value if negative value
    clinvarome_df.loc[
        clinvarome_df["date_between_first_date"] < 0, "last_pathogenic_variant"
    ] = pd.NaT
    clinvarome_df.loc[
        clinvarome_df["date_between_first_date"] < 0, "date_between_first_date"
    ] = 0
    clinvarome_df["diff_date_today"] = (
        (
            pd.to_datetime("2019-12-01", format="%Y-%m-%d")
            - clinvarome_df["first_path_var_date"]
        ).dt.days
    ) // 30

    return clinvarome_df


# Automatic clustering confidence score


def scale_clinvarome_data(clinvarome_annotation, ARRAY_TRANSFORM):
    """
    Quantile scale the 4 features needed for clustering.
    """
    clinvarome_cluster = clinvarome_annotation[
        [
            "date_between_first_date",
            "variant_number",
            "highest_review_confidence",
            "highest_pathogenic_class",
        ]
    ]
    clinvarome_cluster = clinvarome_cluster.replace(ARRAY_TRANSFORM)
    clinvarome_cluster_array = np.array(clinvarome_cluster.fillna(0))
    scaler = QuantileTransformer(n_quantiles=10, random_state=0)
    scale_clinvarome = scaler.fit_transform(clinvarome_cluster_array)
    scale_clinvarome_df = pd.DataFrame(
        scale_clinvarome, columns=clinvarome_cluster.columns
    )
    return scale_clinvarome_df


def give_name_to_cluster(cluster_df):
    """
    Rename cluster to 3-level stars clinical validity of gene classification
    """
    renamed_df = cluster_df.copy()
    for cluster, clinvarome_stars in CLUSTER_NAMES.items():
        renamed_df[renamed_df == cluster] = clinvarome_stars
    return renamed_df


def clusterize_clinvarome(clinvarome_annotation, ARRAY_TRANSFORM):
    """
    Apply AgglomerativeClustering to data.
    Return 3-level stars clinical validity of gene classification.
    """
    scale_clinvarome_df = scale_clinvarome_data(clinvarome_annotation, ARRAY_TRANSFORM)
    cluster = AgglomerativeClustering(
        n_clusters=4, affinity="euclidean", linkage="ward"
    )
    cluster.fit_predict(scale_clinvarome_df)
    cluster_name = pd.DataFrame(cluster.labels_, columns=["cluster_name"])
    renamed_cluster = give_name_to_cluster(cluster_name)
    clinvarome_clusterized = pd.concat([clinvarome_annotation, renamed_cluster], axis=1)
    return clinvarome_clusterized


# Main function


def gather_annotation(vcf_file, clinvarome, compare_gene, compare_variant, gnomad):
    """
    Open the clinvarome from Variant-Alert!. This file will be annoted
    by features from compare_gene, compare_variant and a ClinVar VCF file.
    We also add classifications of genes according to their confidence score.
    """
    gene_disease = {}
    gene_finding = {}
    gene_MC = {}
    gene_AF = {}
    var_nhomalt = {}
    vcf = pysam.VariantFile(vcf_file)
    logger.info("Parse " + vcf_file + " INFO field, stock the data into dictionaries")
    for record in vcf:
        if (
            (record.info["CLNSIG"][0] == "Pathogenic")
            or (record.info["CLNSIG"][0] == "Pathogenic/Likely_pathogenic")
            or (record.info["CLNSIG"][0] == "Likely_pathogenic")
        ):
            gather_clinical_features(record, gene_finding, gene_disease)
            meca_by_variant(record, gene_MC)
            if gnomad:
                gather_dict_gene_max_AF(record, gene_AF)
                get_nhomalt(record, var_nhomalt)
    logger.info("Data processing into DataFrame")
    clinvarome_df = pd.read_csv(clinvarome, sep="\t")
    gene_features_df = get_clinical_dataframe(gene_disease, gene_finding)
    gene_MC_df = get_MC_dataframe(gene_MC)
    first_date_df = first_gene_date(clinvarome_df, compare_gene)
    last_date_df = last_gene_date(clinvarome_df, compare_variant)
    if gnomad:
        gene_AF_df = get_AF_max_by_gene(gene_AF)
        gene_nhomalt_df = get_max_nhomalt_by_gene(var_nhomalt)
        list_df = [
            gene_features_df,
            gene_MC_df,
            gene_AF_df,
            gene_nhomalt_df,
            first_date_df,
            last_date_df,
        ]
    else:
        list_df = [
            gene_features_df,
            gene_MC_df,
            first_date_df,
            last_date_df,
        ]
    clinvarome_annotation = merge_dataframe(clinvarome_df, list_df)
    clinvarome_annotation_var_numb = number_of_variants(clinvarome_annotation)
    clinvarome_annotation_var_numb_date = date_management(
        clinvarome_annotation_var_numb
    )
    clinvarome_annotation_total = clinvarome_annotation_var_numb_date.rename(
        columns={
            "pathogenic_class_new": "highest_pathogenic_class",
            "pathogenic_class_best_review_confidence_new": "highest_pathogenic_class_associate_review_confidence",
            "review_confidence_best_pathogenic_class_new": "highest_review_confidence_associate_pathogenic_class",
            "review_confidence_new": "highest_review_confidence",
            "name_clinvar_new": "clinvarome_date",
            "FAF": "gnomad_max_FAF",
            "nhomalt": "gnomad_max_nhomalt",
        }
    )
    logger.info("Clinical validity gene score attribution")
    clinvarome_annotation_score = clusterize_clinvarome(
        clinvarome_annotation_total, ARRAY_TRANSFORM
    )

    return clinvarome_annotation_score
