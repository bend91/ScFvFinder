import pandas as pd
import subprocess
import logging
import os
import configparser
import numpy as np


data_path = os.path.join(os.path.dirname(__file__), 'data_files')



aa_fp = os.path.join(f"{data_path}", "tidy_aa.csv")
tidy_df = pd.read_csv(aa_fp)


def parse_config_file(fp=None):
    if fp is None:
        config_fp = os.path.join(os.path.dirname(__file__), "./config.ini")
    else:
        config_fp = fp
    config = configparser.ConfigParser()
    config.read(config_fp)
    return config


def check_translation(seq):
    frames = [1, 2, 3]
    valid_frames = []
    for frame in frames:
        codon_seq = [seq[x:x+3] for x in range(frame-1, len(seq)-2, 3)]
        translated_seq = "".join([tidy_df.loc[tidy_df["codon"] == x, "amino_acid"].values[0] for x in codon_seq])
        if "*" not in translated_seq:
            valid_frames.append(frame)
    return valid_frames


def translate_sequence(seq, frame=1, return_codons=False):
    seq = seq.upper()
    codon_seq = [seq[x:x+3] for x in range(frame-1, len(seq)-2, 3)]
    translated_seq = "".join([tidy_df.loc[tidy_df["codon"] == x, "amino_acid"].values[0] if np.isin(list(x), ["A", "C", "T", "G"]).all() else "X" for x in codon_seq])
    if return_codons:
        return translated_seq, codon_seq
    return translated_seq


def find_frame(dna_seq, ref_aa):
    if len(dna_seq) < 3:
        return None
    for frame in range(1, 4):
        if translate_sequence(dna_seq, frame) not in ref_aa:
            return frame
    print("Frame couldn't be determined, check the reference sequence")
    return None


def compare_sequences(seq1, seq2):
    diffs = []
    seqlen = min(len(seq1), len(seq2))
    seq1 = seq1[:seqlen]
    seq2 = seq2[:seqlen]
    for i, char in enumerate(seq1):
        if char != seq2[i]:
            diffs.append(i)
    return diffs


def write_process_to_file(process, file_path):
    with open(file_path, "wb") as f:
        f.write(process.stdout)


def process_command_output(process, name):
    if process.returncode == 0:
        logging.info(f"{name} successful")
        return True
    else:
        logging.error(f"{name} failed")
        logging.error("Stdout: ", process.stdout.decode("utf-8"), "\n")
        logging.error("StdErr: ", process.stderr.decode("utf-8"), "\n")
        raise Exception(f"{name} failed")

def process_fasta_file(fasta_file, id_list):
    header = None
    sequence = b""
    collect_sequence = False
    with open(fasta_file, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if header and collect_sequence:
                    yield header, sequence
                header = line[1:].decode("utf-8").strip()
                sequence = ""
                if header in id_list:
                    collect_sequence = True
                else:
                    collect_sequence = False
            elif collect_sequence:
                sequence += line.decode("utf-8").strip()
        if header and collect_sequence:
            yield header, sequence


def sequence_alignment(query_seq: str, ref_seq: str, similarity_matrix, d: float):
    """
    Sequence alignment using the Needleman-Wunsch algorithm
    """
    if isinstance(similarity_matrix, pd.DataFrame):
        ref_mtx = similarity_matrix
    elif isinstance(similarity_matrix, str):
        ref_mtx = pd.read_csv(similarity_matrix, index_col=0)
    query_mtx = pd.DataFrame(np.zeros((len(ref_seq), len(query_seq))), index=list(ref_seq), columns=list(query_seq))
    match_mtx = pd.DataFrame(index=list(ref_seq), columns=list(query_seq))
    # d = -0.1 # Gap penalty score
    # for i, _ in enumerate(ref_seq):
    #     query_mtx.iloc[i, 0] = d * i
    # for i, _ in enumerate(query_seq):
    #     query_mtx.iloc[0, i] = d * i
    for i, char in enumerate(ref_seq, 1):
        for j, char2 in enumerate(query_seq, 1):
            match = query_mtx.iloc[i-1, j-1] + ref_mtx.loc[char, char2]
            delete = query_mtx.iloc[i-1, j] + d if j < len(query_seq) else -1
            insert = query_mtx.iloc[i, j-1] + d if i < len(ref_seq) else -1
            match_mtx.iloc[i-1, j-1] = pd.DataFrame({
                "m": match,
                "d":  delete,
                "i": insert
            }, index=[char2]).idxmax(axis=1).loc[char2]
            if i < len(ref_seq) and j < len(query_seq):
                query_mtx.iloc[i, j] = max([match, delete, insert])
    i = len(ref_seq) - 1
    j = len(query_seq) - 1
    alignment_a = ""
    alignment_b = ""
    while (i > -1 or j > -1):
        ij_pos = (i > -1 and j > -1)
        if ij_pos and match_mtx.iloc[i, j] == "m":
            alignment_a = ref_seq[i] + alignment_a
            alignment_b = query_seq[j] + alignment_b
            i -= 1
            j -= 1
        elif ij_pos and match_mtx.iloc[i, j] == "i":
            alignment_a = "-" + alignment_a
            alignment_b = query_seq[j] + alignment_b
            j -= 1
        elif ij_pos and match_mtx.iloc[i, j] == "d":
            alignment_a = ref_seq[i] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1
        elif i < 0 and j > -1:
            alignment_a = "-" + alignment_a
            alignment_b = query_seq[j] + alignment_b
            j -= 1
        elif j < 0 and i > -1:
            alignment_a = ref_seq[i] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1
        else:
            print("done")
    return alignment_a, alignment_b


def smith_waterman_align(query_seq: str, germline_seq: str) -> tuple[str, str]:
    """
    Performs local sequence alignment using the Smith-Waterman algorithm with
    scoring parameters empirically determined by abstar for antibody germline
    gene alignment:
        match reward    = +3
        mismatch penalty = -2
        gap open        = -22
        gap extend      = -1

    The high gap-open penalty strongly discourages spurious indels, meaning
    gaps introduced in the alignment are very likely to represent real
    single-nucleotide sequencing errors rather than alignment artefacts.
    This makes the function suitable as a pre-annotation correction step to
    resolve frameshifts before passing sequences to IgBLAST.

    Reference:
        Briney B, Burton DR. (2018). Massively scalable genetic analysis of
        antibody repertoires. bioRxiv 447813.
        https://doi.org/10.1101/447813

    Args:
        query_seq:    the raw query nucleotide sequence (e.g. from Nanopore read)
        germline_seq: the germline reference nucleotide sequence

    Returns:
        A tuple of (aligned_query, aligned_germline) strings, with '-'
        representing gaps. Gaps in aligned_query indicate a deletion in the
        read (likely a sequencing error); gaps in aligned_germline indicate
        an insertion.
    """
    query_seq = query_seq.upper()
    germline_seq = germline_seq.upper()

    MATCH    =  3
    MISMATCH = -2
    GAP_OPEN = -22
    GAP_EXT  = -1

    n_rows = len(germline_seq) + 1  # germline along rows
    n_cols = len(query_seq) + 1     # query along columns

    # Score and traceback matrices
    score_mtx    = [[0.0] * n_cols for _ in range(n_rows)]
    traceback    = [[None] * n_cols for _ in range(n_rows)]

    # Gap score matrices (affine gap penalty)
    # H_gap: best score ending with a gap in the germline (insertion in query)
    # V_gap: best score ending with a gap in the query (deletion in query)
    H_gap = [[float('-inf')] * n_cols for _ in range(n_rows)]
    V_gap = [[float('-inf')] * n_cols for _ in range(n_rows)]

    best_score = 0.0
    best_pos   = (0, 0)

    for i in range(1, n_rows):
        for j in range(1, n_cols):
            # Match/mismatch
            match_score = (MATCH if germline_seq[i - 1] == query_seq[j - 1]
                           else MISMATCH)
            diag = score_mtx[i - 1][j - 1] + match_score

            # Affine gap penalties
            # Gap in germline = insertion in query
            H_gap[i][j] = max(
                score_mtx[i][j - 1] + GAP_OPEN,
                H_gap[i][j - 1]     + GAP_EXT
            )
            # Gap in query = deletion in query (likely sequencing error)
            V_gap[i][j] = max(
                score_mtx[i - 1][j] + GAP_OPEN,
                V_gap[i - 1][j]     + GAP_EXT
            )

            best = max(0.0, diag, H_gap[i][j], V_gap[i][j])
            score_mtx[i][j] = best

            if best == 0.0:
                traceback[i][j] = "stop"
            elif best == diag:
                traceback[i][j] = "diag"
            elif best == H_gap[i][j]:
                traceback[i][j] = "left"
            else:
                traceback[i][j] = "up"

            if best >= best_score:
                best_score = best
                best_pos   = (i, j)

    # Traceback from best scoring cell
    aligned_query    = ""
    aligned_germline = ""
    i, j = best_pos

    while i > 0 and j > 0 and traceback[i][j] != "stop":
        direction = traceback[i][j]
        if direction == "diag":
            aligned_query    = query_seq[j - 1]    + aligned_query
            aligned_germline = germline_seq[i - 1] + aligned_germline
            i -= 1
            j -= 1
        elif direction == "left":
            aligned_query    = query_seq[j - 1] + aligned_query
            aligned_germline = "-"              + aligned_germline
            j -= 1
        elif direction == "up":
            aligned_query    = "-"               + aligned_query
            aligned_germline = germline_seq[i-1] + aligned_germline
            i -= 1

    return aligned_query, aligned_germline


def correct_indels(query_seq: str, germline_seq: str) -> str:
    """
    Uses smith_waterman_align to identify and correct single-nucleotide indels
    in query_seq relative to germline_seq. Gaps in the aligned query are filled
    with the corresponding germline base; insertions in the query relative to
    germline are removed.

    Only processes gaps of length 1 to avoid over-correcting genuine
    SHM-induced indels or larger structural differences.

    Args:
        query_seq:    raw query nucleotide sequence
        germline_seq: germline reference nucleotide sequence

    Returns:
        Corrected query sequence as a string, or the original query_seq if
        no single-nucleotide indels are found.
    """
    aligned_query, aligned_germline = smith_waterman_align(query_seq, germline_seq)

    corrected = []
    i = 0
    while i < len(aligned_query):
        q_base = aligned_query[i]
        g_base = aligned_germline[i]

        if q_base == "-":
            # Deletion in query: check it's isolated (not part of a longer gap)
            prev_gap = (i > 0 and aligned_query[i - 1] == "-")
            next_gap = (i < len(aligned_query) - 1 and aligned_query[i + 1] == "-")
            if not prev_gap and not next_gap:
                corrected.append(g_base)  # restore germline base
            # else: part of a multi-base gap, leave uncorrected (skip)
        elif g_base == "-":
            # Insertion in query: check it's isolated
            prev_gap = (i > 0 and aligned_germline[i - 1] == "-")
            next_gap = (i < len(aligned_germline) - 1 and aligned_germline[i + 1] == "-")
            if not prev_gap and not next_gap:
                pass  # drop the inserted base
            else:
                corrected.append(q_base)  # keep if part of larger insertion
        else:
            corrected.append(q_base)
        i += 1

    return "".join(corrected)


