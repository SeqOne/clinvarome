#!/usr/bin/env python3

EFFECTS = {
    "frameshift_variant": 1,
    "stop_lost": 2,
    "nonsense": 3,
    "splice_donor_variant": 4,
    "splice_acceptor_variant": 5,
    "inframe_deletion": 6,
    "inframe_insertion": 7,
    "inframe_indel": 8,
    "missense_variant": 9,
    "synonymous_variant": 10,
    "5_prime_UTR_variant": 11,
    "3_prime_UTR_variant": 12,
    "genic_upstream_transcript_variant": 13,
    "genic_downstream_transcript_variant": 14,
    "upstream_transcript_variant": 15,
    "downstream_transcript_variant": 16,
    "non-coding_transcript_variant": 17,
    "initiator_codon_variant": 18,
    "initiatior_codon_variant": 19,
    "intron_variant": 20,
    "no_sequence_alteration": 21,
}

MC_CATEGORIES = {
    "frameshift_variant": "stop_fs_splice",
    "splice_donor_variant": "stop_fs_splice",
    "splice_acceptor_variant": "stop_fs_splice",
    "stop_lost": "stop_fs_splice",
    "nonsense": "stop_fs_splice",
    "synonymous_variant": "stop_fs_splice",
    "missense_variant": "missense_inframe",
    "inframe_deletion": "missense_inframe",
    "inframe_insertion": "missense_inframe",
    "inframe_indel": "missense_inframe",
    "5_prime_UTR_variant": "other",
    "3_prime_UTR_variant": "other",
    "genic_upstream_transcript_variant": "other",
    "upstream_transcript_variant": "other",
    "genic_downstream_transcript_variant": "other",
    "downstream_transcript_variant": "other",
    "non-coding_transcript_variant": "other",
    "initiator_codon_variant": "other",
    "initiatior_codon_variant": "other",
    "intron_variant": "other",
    "no_sequence_alteration": "other",
    "Not_provided": "Not_provided",
}


MC_SHORT = {}
for value in MC_CATEGORIES.values():
    MC_SHORT[value] = 1


ARRAY_TRANSFORM = {
    "no_assertion_criteria_provided": 1,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
    "Pathogenic": 1,
    "Pathogenic/Likely_pathogenic": 1,
    "Likely_pathogenic": 0,
}

CLUSTER_NAMES = {0: "3 stars", 1: "2 stars", 2: "0 star", 3: "1 star"}
