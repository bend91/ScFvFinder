import pandas as pd
import subprocess
import logging
import os
import configparser


data_path = './data_files'


aa_fp = os.path.join(f"{data_path}", "tidy_aa.csv")
tidy_df = pd.read_csv(aa_fp)


def parse_config_file():
    config_fp = "./config.ini"
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


def translate_sequence(seq):
    seq = seq.upper()
    codon_seq = [seq[x:x+3] for x in range(0, len(seq)-2, 3)]
    translated_seq = "".join([tidy_df.loc[tidy_df["codon"] == x, "amino_acid"].values[0] for x in codon_seq])
    return translated_seq


def compare_sequences(seq1, seq2):
    diffs = []
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
