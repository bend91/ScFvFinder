import pandas as pd
import subprocess
import logging
import os
import configparser


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

"""
Needleman-Wunsch algorithm for global sequence alignment
construct a similarity matrix, e.g.:

  A  C  T  G
A 1 -1 -1 -1
C -1 1 -1 -1
T -1 -1 1 -1
G -1 -1 -1 1

Can be weighted on importance of which base pairs you care more about, e.g. T-T could be 10 and would bias towards making sure Ts are matched

Maybe just aligning aa sequence and then can get the corresponding codon from the aa index
aas = 'ACDEFGHIKLMNPQRSTVWYX*'
aa_mtx = pd.DataFrame(index=list(aas), columns=list(aas))
for aa in aas:
    for aa2 in aas:
        if aa == aa2:
            aa_mtx.loc[aa, aa2] = 1
        else:
            aa_mtx.loc[aa, aa2] = 0

query_mtx = pd.DataFrame(np.zeros((len(ref_seq), len(query_seq))), index=list(ref_seq), columns=list(query_seq))
match_mtx = pd.DataFrame(index=list(ref_seq), columns=list(query_seq))
d = 0.1 # Gap penalty score

for i, _ in enumerate(ref_seq):
    query_mtx.iloc[i, 0] = d * i

for i, _ in enumerate(query_seq):
    query_mtx.iloc[0, i] = d * i


for i, char in enumerate(ref_seq, 1):
    for j, char2 in enumerate(query_seq, 1):
        match = query_mtx.iloc[i-1, j-1] + aa_mtx.loc[char, char2]
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

The above works!



i = len(ref_seq) - 1
j = len(query_seq) - 1
alignment_a = ""
alignment_b = ""
while (i > 0 or j > 0):
    print(f"i: {i}, j: {j}")
    if (i > 0 and j > 0 and (query_mtx.iloc[i, j] == query_mtx.iloc[i-1, j-1] + aa_mtx.loc[ref_seq[i], query_seq[j]])):
        alignment_a = ref_seq[i] + alignment_a
        alignment_b = query_seq[j] + alignment_b
        i -= 1
        j -= 1
    elif i > 0 and (query_mtx.iloc[i, j] == query_mtx.iloc[i-1, j] + d):
        alignment_a = ref_seq[i] + alignment_a
        alignment_b = "-" + alignment_b
        i -= 1
    else:
        alignment_a = "-" + alignment_a
        alignment_b = query_seq[j] + alignment_b
        j -= 1

positions = [i for i, _ in enumerate(ref_seq)]
query_position = 0
alignment = []
for i in positions:
    if query_mtx.iloc[i, query_position] == 1:
        alignment.append(query_seq[query_position])
        query_position += 1
    else:
        if query_mtx.iloc[i, query_position + 1] == 1:
            alignment.append("-")
            query_position += 1
        else:
            alignment.append(query_seq[query_position])
            query_position += 1



will provide dna sequence and amino acid sequence

seq1
seq2
score = 0
offset = 0
alignment = ""
for char in seq1:

    for i, char2 in enumerate(seq2):
        score2 = 0
        if similarity_mtx[char, char2] < 0:
            for j in range(i, len(seq2)):
                score2 = score2 // j

                offset += j


            # score2.append(similarity_mtx[char, seq2[i+1]])
        else:
            alignnment += char2




"""
