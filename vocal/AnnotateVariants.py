"""
Annotate a set of variants

@hrichard
"""
from itertools import groupby
from operator import itemgetter
import warnings

from data_loader_tmp import AASTOP
from data_loader_tmp import codons_stop
from data_loader_tmp import get_aa_dict
from data_loader_tmp import mutation_pattern
import numpy as np
import pandas as pd

debug = False

variantsdf_column = [
    "variant_start_aligned",
    "variant_end_aligned",
    "variant_type",
    "variant_size",
    "aa_ref",
    "aa_variant",
    "aa_pos_ref_start",
    "aa_pos_ref_end",
    "aa_pos_query_start",
    "aa_pos_query_end",
    "nt_pos_ref_start",
    "nt_pos_ref_end",
    "nt_pos_query_start",
    "nt_pos_query_end",
    "nt_pattern",
]


def runs_of_ones(lst):
    """
    list[int] -> list[Tuple(int,int, int)]
    Input: a list of values with 0 or 1
    Returns: a list of 3-tuples describing the runs of 1: (start, end+1, size)
    example
    runs_of_ones([ 1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0]) -> [(0, 5, 5), (12, 16, 4)]
    """
    lout = []
    for k, v in groupby(enumerate(lst), key=itemgetter(1)):
        if k:
            v = list(v)
            lout.append((v[0][0], v[-1][0] + 1, v[-1][0] - v[0][0] + 1))
    return lout


def translateseq(nt_seq):
    """
    str -> str
    ----------
    translate the sequence nt_seq, taking into account gap character as well as
    possible frameshifts
    The sequence can be of any length, however a message is printed if 1 or 2 nt characters
    have not been translated into an amino acid
    If the sequence contains "N" characters, the corresponding codons are translated
    into the "X" letter, irrespective of the possible translations
    ----------------------
    the characters after a stop codon are always gaps '-'
    translateseq("ATGGAGCTAGTCTAAGTGTAGGCT") -> MELV*---
    translateseq("ATG---GAG") -> "M-E"
    translateseq("ATG--GAG-CTA") -> "MEL"
    translateseq("CCN") -> "X" (this should give "P" for proline)
    """
    aa_seq = ""
    codon = ""
    nt_codon = 0
    nt_gap = 0
    stop_seen = False
    ##Very pedestrian version for now
    for nt in nt_seq:
        if nt == "-":
            nt_gap += 1
        else:
            codon += nt
            nt_codon += 1
        if nt_gap == 3:
            aa_seq += "-"
            nt_gap = 0
        if nt_codon == 3:
            ### There are only - after the stop letter *
            if stop_seen:
                aa_seq += "-"
            elif codon in codons_stop:
                aa_seq += AASTOP
                stop_seen = True
            elif "N" in codon:
                aa_seq += "X"
            else:
                aa_seq += get_aa_dict()[codon]
            codon = ""
            nt_codon = 0
            nt_gap = 0
    if nt_codon != 0:
        # print("In tranlateseq")
        # print(nt_seq)
        warnings.warn(f"Sequence length was not a multiple of 3, remaining: {codon}")
    return aa_seq


def pos_stop(aaseq):
    """
    str -> int
    returns the number of AA between the last stop to the end of the sequence
    pos_stop("ASFC*--") -> 2
    pos_stop("ASRF*") -> 0
    pos_stop("RFS") -> -1
    """
    pos = -1
    if AASTOP in aaseq:
        pos = aaseq[::-1].index(AASTOP)
    return pos


def induced_truncation_len(nt_ref_seq, nt_query_seq):
    """
    str * str -> Tuple(int, int, int)
    returns the length of the nearest stop to the end
    on the translated version of the sequence
    and the position after the stop in the ref and the query nt sequence
    ###TODO: the gap characters are not necessarily taken in to account
    """
    aa_ref = translateseq(nt_ref_seq)
    aa_query = translateseq(nt_query_seq)
    pos_stop_ref = pos_stop(aa_ref)
    pos_stop_query = pos_stop(aa_query)
    trunc_len = pos_stop_query - pos_stop_ref
    return (
        pos_stop_query - pos_stop_ref,
        3 * (len(aa_query) - pos_stop_ref),
        3 * (len(aa_query) - pos_stop_query),
    )


def correct_truncations(
    nt_ref_seq, nt_query_seq, max_allowed_truncation=5, max_tryout=5
):
    """
    str * str -> str
    Translates nt_ref_seq and nt_query_seq and checks if a truncation
    was introduced earlier in the translated nt_query_seq
    If this is the case, corrects the faulty nucleotides in query_seq in the following
    order:
        - search for a single mutation on the created stop
        - search for deletions or insertions upstream of the stop that are not a multiple of 3
        and correct them in order
    If there are multiple out of frame deletions, they will be corrected as long as there is
    a stop induced

    TODO: there may be a problem when the gap is just on the stop codon.
    """
    ref_orig, query_orig = nt_ref_seq, nt_query_seq
    trunc_len, r_stop_pos, q_stop_pos = induced_truncation_len(nt_ref_seq, nt_query_seq)
    if debug:
        print(
            f"Checking for truncations: diff of pos: {trunc_len}\nQuery stop: {q_stop_pos}, Reference stop: {r_stop_pos}"
        )
    tryout = 0
    while trunc_len > max_allowed_truncation and tryout < max_tryout:
        if debug:
            print(
                f"Found an early truncation at position {q_stop_pos} in Query (length: {trunc_len})"
            )
        ##we search for the first deletion and insertion in 5' which is not of length multiple of 3
        deletions = runs_of_ones([int(s == "-") for s in nt_query_seq])
        deletions.reverse()
        insertions = runs_of_ones([int(s == "-") for s in nt_ref_seq])
        insertions.reverse()
        if debug:
            print(f"found deletions:{deletions} and insertions:{insertions}")
        ##Correct the sequence
        if len(insertions) > 0:
            for s, e, l in insertions:
                if e <= r_stop_pos and l % 3 != 0:
                    if debug:
                        print(
                            f"Found an out of frame insertion at position {s}-{e}, correcting the ref and query sequences"
                        )
                    nt_query_seq = nt_query_seq[:s] + nt_query_seq[e:]
                    nt_ref_seq = nt_ref_seq[:s] + nt_ref_seq[e:]
        elif nt_query_seq[q_stop_pos : q_stop_pos + 3] in codons_stop:
            correction_at_stop = nt_ref_seq[q_stop_pos : q_stop_pos + 3]
            if debug:
                print(f"This is an induced stop, correcting with {correction_at_stop}")
            nt_query_seq = (
                nt_query_seq[:q_stop_pos]
                + correction_at_stop
                + nt_query_seq[q_stop_pos + 3 :]
            )
        else:
            for s, e, l in deletions:
                if e <= q_stop_pos and l % 3 != 0:
                    if debug:
                        print(
                            f"Found an out of frame deletion at position {s}-{e}, correcting the query sequence"
                        )
                    nt_query_seq = nt_query_seq[:s] + nt_ref_seq[s:e] + nt_query_seq[e:]
                    break
        trunc_len, r_stop_pos, q_stop_pos = induced_truncation_len(
            nt_ref_seq, nt_query_seq
        )
        tryout += 1
    if tryout == max_tryout:
        print("In correct_truncations")
        print(ref_orig)
        print(query_orig)
        print("**************")
        print(nt_ref_seq)
        print(nt_query_seq)
        print(
            f"************** Truncation of size {trunc_len}, at positions {r_stop_pos}, and {q_stop_pos}"
        )
        print(translateseq(nt_ref_seq))
        print(translateseq(nt_query_seq))
        print("**************")
        warnings.warn("did not manage to correct the induced truncation here")

    return nt_ref_seq, nt_query_seq


def check_gap_boundaries(ref, query, max_deletion_size=10):
    """
    str * str * int -> str * str
    from two sequences ref and query aligned (of same length), detects the positions with
    gaps, and attempt to shift the letters around those gaps to see if the number of mismatches
    decreases
    Only considers deletions of less than max_deletion_size AA.
    """
    if len(ref) != len(query):
        warnings.warn(
            f"sequences are not of same length when checking gap boundaries, diff: {len(ref) - len(query)}"
        )
    deletions = runs_of_ones([int(s == "-") for s in query])
    ##The only time when we can improve the alignment next to a deletion will be if one of the
    ##boundary letter is a mismatch and can be paired up at the other end of the deletion
    for s, e, l in deletions:
        if l > max_deletion_size:
            warnings.warn(f"Skipping correction of deletion: length:{l}")
            continue
        i = 1
        # print(f"deletion {s}-{e} (length {l})")
        # print(f"{ref[s-1:e+1]}")
        # print(f"{query[s-1:e+1]}")
        # print("*****")
        # print(i, ref[s-i], query[s-i], ref[e-i])
        while i <= l and ref[s - i] != query[s - i] and query[s - i] == ref[e - i]:
            i += 1
        if i > 1:
            lb = i - 1
            # print(f"found a block of size {lb} at the beginning")
            # let's exchange the letters positions on the size lb bloc
            query = query[: s - lb] + query[s:e] + query[s - lb : s] + query[e:]
            continue
        i = 0
        # print("around the deletion", i, ref[e+i], query[e+i], ref[s+i])
        while i <= l and ref[e + i] != query[e + i] and query[e + i] == ref[s + i]:
            i += 1
        if i > 0:
            # The other exchange
            # print(f"found a block of size {i} at the end")
            query = query[:s] + query[e : e + i] + query[s:e] + query[e + i :]
    insertions = runs_of_ones([int(s == "-") for s in ref])
    ##TODO insertions
    if len(insertions) > 0:
        warnings.warn(
            "Insertions detected, the check_gap_boundaries function is not taking those into account"
        )

    return (ref, query)


def diff_type(r, q):
    """
    str * str -> str
    return a str giving the type of difference between 2 letters from ref an query
    r == "-" -> "I"nsertion
    q == "-" -> "D"eletion
    q == "*" -> "T"runcation (stop codon)
    else "M"utation
    """
    if r == "-":
        return "I"
    elif q == "-":
        return "D"
    elif q == AASTOP:
        return "T"
    else:
        return "M"


def get_raw_coords(seq):
    """
    str -> np.array(int)
    -------------------
    give a nparray of the coordinates of all positions in seq when not taking gaps into account
    get_raw_coord("GV--A") --> [1,2,3,3,3]
    """
    return np.array([1] + [int(not c == "-") for c in seq[:-1]]).cumsum()


def list_of_diffs(rseq, qseq, exception_char="X"):
    """
    str * str -> list(Tuple(int, str, str, str) )
    ---------------------------------------------
    returns all the positions where a character is different between the two sequences
    numbering of positions starts at 1.

    """
    ldiffs = [
        (i + 1, diff_type(r, q), r, q)
        for i, (r, q) in enumerate(zip(rseq, qseq))
        if r != q and q != exception_char and r != exception_char
    ]
    return ldiffs


def compare_aa_seq(aa_ref, aa_query, optimize_alignment=True, verbose=False):
    """
    str * str * Bool -> DataFrame
    ----------------------
    compare two aligned Amino acid sequences and return a DataFrame summarizing
    all the variations between the sequences
    The parameter optimize_alignment allows to test if the alignment can be
    improved by shifting amino acids next to the deletions

    """
    if optimize_alignment:
        aa_ref, aa_query = check_gap_boundaries(aa_ref, aa_query)
        if verbose:
            print("AA sequences after correction")
            print(aa_ref)
            print(aa_query)
    diffs = list_of_diffs(aa_ref, aa_query)
    if len(diffs) == 0:
        df_variants = pd.DataFrame(columns=variantsdf_column)
        return df_variants
    # print(f"Found {len(diffs)} differences")
    ##  merge the runs of insertion or deletion
    ##  3 cases :
    ##    - basic mutation
    ##    - deletion 1+ Amino acids
    ##    - insertion 1+ amino acids
    ##    - truncation at the stop
    groups = [1]
    stop_seen = False
    for cur, next in zip(diffs[:-1], diffs[1:]):
        if stop_seen:
            groups.append(groups[-1])
            continue
        if (
            cur[1] in ("D", "I")
            and cur[1] == next[1]
            and cur[0] == next[0] - 1
            or cur[1] == "T"
        ):
            groups.append(groups[-1])
        else:
            if next[1] == "T":
                stop_seen = True
            groups.append(groups[-1] + 1)

    variants_df_long = pd.DataFrame(
        diffs, columns=["aa_pos_align", "variant_type", "aa_ref", "aa_variant"]
    )
    variants_df_long["IndexMutations"] = groups
    variants_df = variants_df_long.groupby("IndexMutations").agg(
        variant_start_aligned=pd.NamedAgg(column="aa_pos_align", aggfunc="min"),
        variant_end_aligned=pd.NamedAgg(column="aa_pos_align", aggfunc="max"),
        variant_type=pd.NamedAgg(column="variant_type", aggfunc=lambda x: x.iloc[0]),
        aa_ref=pd.NamedAgg(column="aa_ref", aggfunc="".join),
        aa_variant=pd.NamedAgg(column="aa_variant", aggfunc="".join),
    )
    variants_df["variant_size"] = (
        variants_df["variant_end_aligned"] - variants_df["variant_start_aligned"] + 1
    )
    return variants_df


def alignedCDSvariants(nt_ref, nt_query, verbose=False):
    """
    str * str -> DataFrame
    columns: nt_ref, nt_variant, nt_pos_ref_start, nt_pos_ref_end, aa_ref, aa_variant, aa_pos_ref_start, aa_pos_ref_end, variant_type, variant_size
    ----------------------
    Input :
        two aligned CDS sequences ('-' is the gap character)
        they have to be of the same length
    Return :
        The set of mutations that are changing the translated sequence
        The stop codon is translated into * and the following codons into -

    Examples:
        r: "ATGGATCTAGGTGGAGGCTAA"  (translated MDLGGG*)
        q: "ATGGAGCTAGGTGTGGGCTAA" (translated MELGVG*)
        alignedCDSvariants(r, q) returns infos about D2E, and G5V

        r: "CAGCTGGCCCGACCCTGCATG------TCAGCTTAA"  (translated QLARPC--MSA*)
        q: "CACCTG------GGCCGACCCTGCATGTCAGCTTAA"  (translated HL--GRPCMSA*)
        alignedCDSvariants(r, q) returns Q1H, AR3-4del, P5G, C6R, M7P, CM8-9ins

        r: "GCCCGACCCTGCATGTCAGCTTAA"  (translated ARPCMSA*)
        q: "GGCCGACCCTAAATGTCAGCTTAA"  (translated GRP*----)
        alignedCDSvariants(r, q) returns A1G and a premature termination at position 4
    """
    # translate ref and query
    if len(nt_ref) != len(nt_query):
        raise ValueError("Query and reference are not of the same length")
    ###First take care of possible in frame deletions and correct them
    nt_ref, nt_query = correct_truncations(nt_ref, nt_query)

    nt_ref_coord = get_raw_coords(nt_ref)
    nt_query_coord = get_raw_coords(nt_query)
    aa_ref = translateseq(nt_ref)
    aa_query = translateseq(nt_query)
    aa_ref_coord = get_raw_coords(aa_ref)
    aa_query_coord = get_raw_coords(aa_query)
    if verbose:
        print("=== Translated sequences")
        print(aa_ref)
        print(aa_query)
    variants_df = compare_aa_seq(aa_ref, aa_query, verbose=verbose)
    # if there are no variants in the sequence
    if variants_df.empty:
        return variants_df
    ##Now add the variants in nt sequence together with their coordinate in ref space
    ## (e.g. without gap characters)
    ##rather do that with 2 functions
    variants_df["aa_pos_ref_start"] = aa_ref_coord[
        variants_df["variant_start_aligned"] - 1
    ]
    variants_df["aa_pos_ref_end"] = aa_ref_coord[variants_df["variant_end_aligned"] - 1]
    variants_df["aa_pos_query_start"] = aa_query_coord[
        variants_df["variant_start_aligned"] - 1
    ]
    variants_df["aa_pos_query_end"] = aa_query_coord[
        variants_df["variant_end_aligned"] - 1
    ]
    ###We find the corresponding positions on the nt sequence
    ### This will not work in the case of out of frame indels. we should work this out from the
    ### the list of differences in the nt sequence.
    variants_df["nt_pos_ref_start"] = (variants_df["aa_pos_ref_start"] - 1) * 3 + 1
    variants_df["nt_pos_ref_end"] = variants_df["aa_pos_ref_end"] * 3
    variants_df.loc[lambda df: df.variant_type == "I", "nt_pos_ref_end"] = 0
    variants_df["nt_pos_query_start"] = (variants_df["aa_pos_query_start"] - 1) * 3 + 1
    variants_df["nt_pos_query_end"] = variants_df["aa_pos_query_end"] * 3
    variants_df.loc[lambda df: df.variant_type == "D", "nt_pos_query_end"] = 0
    # print(variants_df)
    ###Now we list the differences here
    nt_ref_raw = nt_ref.replace("-", "")
    nt_query_raw = nt_query.replace("-", "")
    nt_mut_patterns = []
    aa_mut_patterns = []
    for row in variants_df.itertuples(index=False):
        rseq, qseq = "", ""
        rseq = nt_ref_raw[(row.nt_pos_ref_start - 1) : row.nt_pos_ref_end]
        qseq = nt_query_raw[(row.nt_pos_query_start - 1) : row.nt_pos_query_end]
        if row.variant_type == "D":
            # print("Logging Deletion type variation")
            nt_variations = [
                row.nt_pos_ref_start,
                row.nt_pos_query_start,
                "D",
                rseq,
                "-" * row.variant_size * 3,
            ]
            nt_patt = f"del:{nt_variations[0]}:{len(nt_variations[3])}"
        elif row.variant_type == "I":
            # print("Logging Insertion type variation")
            nt_variations = [
                row.nt_pos_ref_start,
                row.nt_pos_query_start,
                "I",
                "-" * row.variant_size * 3,
                qseq,
            ]
            nt_patt = f"ins:{nt_variations[0]}:{len(nt_variations[4])}"
        else:
            # print("Logging Mutation or Truncation type variation")
            nt_variations = [
                (row.nt_pos_ref_start - 1 + p, row.nt_pos_query_start - 1 + p, t, r, q)
                for (p, t, r, q) in list_of_diffs(rseq, qseq)
            ]
            nt_patt = ",".join(
                [f"{nt_r}{rpos}{nt_q}" for rpos, qpos, t, nt_r, nt_q in nt_variations]
            )
        # print(nt_variations)
        ##We log the pattern as an example for the moment
        nt_mut_patterns.append(nt_patt)
        aa_mut_patterns.append(
            mutation_pattern(
                row.variant_type,
                row.aa_pos_ref_start,
                row.aa_pos_ref_end,
                row.aa_ref,
                row.aa_variant,
            )
        )

    variants_df["nt_pattern"] = nt_mut_patterns
    variants_df["aa_pattern"] = aa_mut_patterns
    return variants_df


if __name__ == "__main__":
    ##example 1
    nt_ref1 = "ATGGATCTAGGTGGAGGCTAA"
    nt_query1 = "ATGGAGCTAGGTGTGGGCTAA"
    print(nt_ref1)
    print(nt_query1)
    df_e1 = alignedCDSvariants(nt_ref1, nt_query1)
    print(df_e1)
    nt_ref2 = "CAGCTGGCCCGACCCTGCATG------TCAGCTTAA"
    nt_query2 = "CACCTG------GGCCGACCCTGCATGTCAGCTTAA"
    print(nt_ref2)
    print(nt_query2)
    df_e2 = alignedCDSvariants(nt_ref2, nt_query2)
    df_e2.to_csv("example_mutations.csv")
    print(df_e2)
    nt_ref3 = "GCCCGACCCTGCATGTCAGCTTAA"
    nt_query3 = "GGCCGACCCTAAATGTCAGCTTAA"
    print(nt_ref3)
    print(nt_query3)
    df_e3 = alignedCDSvariants(nt_ref3, nt_query3)
    print(df_e3)
