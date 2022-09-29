"""
Functions to process PSL alignments fieldnames

@hrichard
"""

from Bio import SeqIO
from data_loader_tmp import GAP_CHAR
import pandas as pd


PSL_LABELS = [
    "Matches",
    "Mismatches",
    "RepMatch",
    "Ns",
    "QGapCount",
    "QGapBases",
    "TGapCount",
    "TGapBases",
    "Strand",
    "qName",
    "qSize",
    "qStart",
    "qEnd",
    "tName",
    "tSize",
    "tStart",
    "tEnd",
    "blockCount",
    "blockSizes",
    "qStarts",
    "tStarts",
]

PSL_TYPES = {
    "Matches": int,
    "Mismatches": int,
    "RepMatch": int,
    "Ns": int,
    "QGapCount": int,
    "QGapBases": int,
    "TGapCount": int,
    "TGapBases": int,
    "Strand": str,
    "qName": str,
    "qSize": int,
    "qStart": int,
    "qEnd": int,
    "tName": str,
    "tSize": int,
    "tStart": int,
    "tEnd": int,
    "blockCount": int,
    "blockSizes": str,
    "qStarts": str,
    "tStarts": str,
}

PSL_LABELS_PRESIDENT = [
    "Matches",
    "Mismatches",
    "RepMatch",
    "Ns",
    "QGapCount",
    "QGapBases",
    "TGapCount",
    "TGapBases",
    "Strand",
    "QName",
    "QSize",
    "QStart",
    "QEnd",
    "TName",
    "TSize",
    "TStart",
    "TEnd",
    "BlockCount",
    "BlockSizes",
    "QStarts",
    "TStarts",
]

pl_lower = [x.lower() for x in PSL_LABELS]
plp_lower = [x.lower() for x in PSL_LABELS_PRESIDENT]

d_PSL_prescorr = dict(
    [
        ("PSL_" + PSL_LABELS_PRESIDENT[i], PSL_LABELS[pl_lower.index(x)])
        for i, x in enumerate(plp_lower)
    ]
)


def readPSL(filein, index="qName"):
    """
    str -> DataFrame
    Reads in a PSL alignment and return the corresponding pandas data frame
    Using index as the index of the dataframe
    we verify that the index is unique when reading and keep the
    alignment with the highest number of Matches
    """
    df_PSL = pd.read_table(
        filein, sep="\t", header=None, names=PSL_LABELS, skiprows=5, index_col=index
    )
    ##This one does not work... maybe because of columns type to str for the Matches?
    # df_PSL = df_PSL.sort_values('Matches',
    #                             ascending=False).groupby(df_PSL,
    #                                                      sort=False).first().reset_index().set_index('qName')
    df_PSL = df_PSL.groupby(df_PSL.index).first()

    ###
    # df_PSL.to_csv("toto.tsv", sep = "\t")
    return df_PSL


def readPSLfromPresident(filein, index="qName", PSL_prefix="PSL_"):
    """
    str -> DataFrame
    Reads the PSL part of a president output by getting the columns that are
    with the PSL prefix
    """
    df_PSL_president = pd.read_table(
        filein, sep="\t", header=0, usecols=d_PSL_prescorr.keys()
    )
    df_PSL_president.columns = [d_PSL_prescorr[x] for x in df_PSL_president.columns]
    df_PSL_president = df_PSL_president.dropna().astype(PSL_TYPES)
    df_PSL_president = df_PSL_president.set_index(index)
    return df_PSL_president


def PSLpretty(target, query, PSLrow):
    """
    str * str * namedTuple -> Tuple(str, str, list[Tuple(int, int)], list(Tuple[int, int]))
    -----------------------------------------
    Simplified version of the PSLpretty code from kent for BLAT
    reads in two sequences with their PSL alignment information and
    output the two aligned sequences with gap characters where needed
    and two lists of pairs that inform on the new coordinates in the target and in the query

    Hypothesis: there are only alignments on the positive strand
                The PSL is well formed
    TODO: add a sequence position mapping function. it can just be the
    positions with gaps and the corresponding cumulated number of gaps at those positions
    for instance the following aligned sequence
    01234   56  789
    ACTGT---GT--CGC
    012345678901234
    will results in the following list of positions
    [(0,0), (5,8), (7,12)]
    TODO: check for PSL parameter consistency, add negative strand alignment
    """
    b_Sizes = [int(x) for x in PSLrow.blockSizes.split(",") if len(x) > 0]
    q_Starts = [int(x) for x in PSLrow.qStarts.split(",") if len(x) > 0]
    t_Starts = [int(x) for x in PSLrow.tStarts.split(",") if len(x) > 0]
    t_bloc_end, q_bloc_end = t_Starts[0], q_Starts[0]
    al_target, al_query = "", ""
    cum_query_pos = [(0, -q_Starts[0])]
    cum_target_pos = [(0, -t_Starts[0])]
    ngq_tot, ngt_tot = 0, 0
    for t_bloc_start, q_bloc_start, bloc_size in zip(t_Starts, q_Starts, b_Sizes):
        # print the gaps from last bloc
        ngaps_query = t_bloc_start - t_bloc_end
        ngaps_target = q_bloc_start - q_bloc_end
        # Note: there can be gaps on both target and query, and there multiples possibilities for gap placement
        ## we pad the target gap first and then the query gaps
        al_target += ngaps_target * GAP_CHAR + target[t_bloc_end:t_bloc_start]
        al_query += query[q_bloc_end:q_bloc_start] + ngaps_query * GAP_CHAR
        # if ngaps_query*ngaps_target>0:
        #    print(f"Double gap: ({ngaps_query}={query[q_bloc_end:q_bloc_start]}, {ngaps_target}={target[t_bloc_end:t_bloc_start]})")
        #
        t_bloc_end, q_bloc_end = t_bloc_start + bloc_size, q_bloc_start + bloc_size
        al_target += target[t_bloc_start:t_bloc_end]
        al_query += query[q_bloc_start:q_bloc_end]
        if ngaps_query > 0:
            ngq_tot += ngaps_query
            q_corresp = (q_bloc_start, q_bloc_start + ngq_tot - q_Starts[0])
            cum_query_pos.append(q_corresp)
        if ngaps_target > 0:
            ngt_tot += ngaps_target
            t_corresp = (t_bloc_start, t_bloc_start + ngt_tot - t_Starts[0])
            cum_target_pos.append(t_corresp)
    ## We take out gap aligning to gap
    ## Commented out to be able to work with the cumulative positions
    # al_target, al_query = ["".join(x) for x in zip(*[(r,q)
    #                                                  for r,q in zip(al_target, al_query)
    #                                                         if r != GAP_CHAR or q!=GAP_CHAR
    #                                                  ]
    #                                                )]
    if t_bloc_end != PSLrow.tEnd:
        print(
            f"Error the end of target is different from field: {t_bloc_end} and {PSLrow.tEnd}"
        )
    if q_bloc_end != PSLrow.qEnd:
        print(
            f"Error the end of target is different from field: {q_bloc_end} and {PSLrow.qEnd}"
        )
    # Use that for check?
    # t_bloc_start = PSLrow.tStart
    # q_bloc_start = PSLrow.qStart
    # check last END is the same as qEND and tEND
    return (al_target, al_query, cum_target_pos, cum_query_pos)


def posconverter(pos, corresp):
    """
    int * list[(int, int)] -> int
    Converts a position pos from a list of pairs
    TODO: do the same function for a sorted list
    """
    convpos = [cp + pos - p for p, cp in corresp if p <= pos].pop()
    return convpos


def main():
    print("Testing PSL output")
    t_wuhan = SeqIO.read("psl_tests/NC_045512.2.fasta", "fasta")
    q_test = SeqIO.read("psl_tests/pblat_seq_test.fa", "fasta")
    df_PSL = pd.read_table("psl_tests/pblat_row_test.tsv", sep="\t", header=None)
    df_PSL.columns = PSL_LABELS
    print(df_PSL)
    print(t_wuhan)
    print(q_test)
    ###Seems to be working
    for row in df_PSL.itertuples(index=False):
        a, b = PSLpretty(str(t_wuhan.seq), str(q_test.seq), row)
        print(a)
        print(b)


if __name__ == "__main__":
    main()
