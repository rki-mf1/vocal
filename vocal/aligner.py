"""

Historical pairwise global aligner using the needleman_wunsch algorithm and
helper functions for aligning sequences from the BioPython pairwise2 library

@StefanFrankBio @hrichard
"""

import time

from Bio import pairwise2

"""

Historical aligner part


"""


def initialize_matrix(len1, len2, gap_penalty):
    matrix = [0 for i in range((len1 + 1) * (len2 + 1))]
    for i in range(1, len1 + 1):
        matrix[i * (len2 + 1)] = gap_penalty * i
    for j in range(1, len2 + 1):
        matrix[j] = gap_penalty * j
    return matrix


def needleman_wunsch(string1, string2, matrix, scoring_scheme):
    """
    Fills the DP matrix for NW alignment of string1 to string2

    """
    len1 = len(string1)
    len2 = len(string2)
    traceback_matrix = [1 for i in range((len1 + 1) * (len2 + 1))]
    for i in range(1, len2 + 1):
        traceback_matrix[i] = 2
    for j in range(1, len2 + 1):
        for i in range(1, len1 + 1):
            cursor = i * (len2 + 1) + j
            if string1[i - 1] == string2[j - 1]:
                match_score = matrix[cursor - len2 - 2] + scoring_scheme[0]
            else:
                match_score = matrix[cursor - len2 - 2] + scoring_scheme[1]
            insert_score = matrix[cursor - len2 - 1] + scoring_scheme[2]
            delete_score = matrix[cursor - 1] + scoring_scheme[2]
            possible_scores = [match_score, insert_score, delete_score]
            matrix[cursor] = max(possible_scores)
            traceback_matrix[cursor] = possible_scores.index(max(possible_scores))
    return traceback_matrix


def traceback(traceback_matrix, string1, string2):
    ## Here we need to do semi global: take the max over the last column
    ##
    traceback_matrix = traceback_matrix[::-1]
    cursor = 0
    len2 = len(string2)
    aligned1 = []
    aligned2 = []
    position1 = len(string1) - 1
    position2 = len(string2) - 1
    while not position1 == position2 == -1:
        if traceback_matrix[cursor] == 0:
            aligned1.append(string1[position1])
            aligned2.append(string2[position2])
            position1 -= 1
            position2 -= 1
            cursor += len2 + 2
        elif traceback_matrix[cursor] == 1:
            aligned1.append(string1[position1])
            aligned2.append("-")
            position1 -= 1
            cursor += len2 + 1
        else:
            aligned1.append("-")
            aligned2.append(string2[position2])
            position2 -= 1
            cursor += 1
    aligned1 = aligned1[::-1]
    aligned2 = aligned2[::-1]
    start_offset = 1
    while aligned1[0] == "-":
        del aligned1[0]
        del aligned2[0]
        start_offset += 1
    return (start_offset, "".join(aligned1), "".join(aligned2))


def align(string1, string2, scoring_scheme=[1, -1, -2]):
    """
    str * str * list[int, int, int] -> Tuple(int, str, str)

    Align string1 onto string2 in a semi global way (with
    string1 included in string2).
    The scoring scheme is given as a 3 elements list [match, mismatch, gap]

    returns a Tuple with:
     - position in seq2 of the first aligned base of seq1 (start at 1)
     - aligned seq 1
     - aligned seq 2
    """
    start_time = time.time()
    print("Init Matrix: ", end="")
    matrix = initialize_matrix(len(string1), len(string2), scoring_scheme[2])
    t_minit = time.time()
    print("{:.3f} seconds".format(t_minit - start_time))
    print("NW algo: ", end="")
    traceback_matrix = needleman_wunsch(string1, string2, matrix, scoring_scheme)
    t_NW = time.time()
    print("{:.3f} seconds".format(t_NW - t_minit))
    print("Tracing back", end="")
    t_traceback = time.time()
    print(" {:.3f} seconds".format(t_traceback - t_NW))
    return traceback(traceback_matrix, string1, string2)


"""

Biopython part

"""


def Biopairwise_align(
    ref_seq, query_seq, match=2, mismatch=-1, gap_open=-2, gap_extend=-0.5
):
    alignments = pairwise2.align.localms(
        ref_seq, query_seq, match, mismatch, gap_open, gap_extend
    )
    al_start, al_end = alignments[0].start, alignments[0].end
    ref_al = alignments[0].seqA[al_start:al_end]
    query_al = alignments[0].seqB[al_start:al_end]
    offset = al_start + 1
    return (offset, ref_al, query_al)


"""

Helper function to simply return the alignment of the two sequences

"""
