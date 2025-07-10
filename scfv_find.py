import os
import subprocess
from utils import *
import numpy as np
import logging
import hashlib
import difflib


# TODO - Add logging to a logfile instead of printing - done this, kept some prints for users to see progress
# TODO - check if blastn is in path, if not then ask to download  - Donne this
# TODO - Add multi-platform support, at the moment just MacOS is supported


env = os.environ.copy()


def check_tools(tool):
    try:
        subprocess.run([tool])
        return True
    except FileNotFoundError:
        return False


def install_tool(tool):
    tool_execs = {
        "blastn": "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-universal-macosx.tar.gz",
        "igblast": "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-macosx.tar.gz"
        }
    if tool not in list(tool_execs.keys()):
        print(f"{tool} not supported yet, something went wrong...")
        return None
    process = subprocess.run(["curl", tool_execs[tool], "-o", f"{data_path}/{tool}.tar.gz"], capture_output=True)
    process_command_output(process, f"Download {tool}")
    process = subprocess.run(["tar", "-xvzf", f"{data_path}/{tool}.tar.gz", "-C", f"{data_path}"], capture_output=True)
    process_command_output(process, f"Extract {tool}")
    if tool == "igblast":
        print("Downloading mouse database...")
        process = subprocess.run(["curl", "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/mouse_gl_VDJ.tar", "-o", f"{data_path}/mouse_gl_VDJ.tar"], capture_output=True)
        process_command_output(process, "Download igblast database")
        os.mkdir(f"{data_path}/{tool}/database")
        process = subprocess.run(["tar", "-xvf", f"{data_path}/mouse_gl_VDJ.tar", "-C", f"{data_path}/{tool}/database"], capture_output=True)
        process_command_output(process, "Extract igblast database")
    env["PATH"] = f"{data_path}/{tool}/bin" + os.pathsep + env["PATH"]


def get_folder_path()->str:
    """
    Validates a user entered folder path to be the one that contains the file they want
    """
    fpath = input("Enter the full path to the folder containing the fastq file you want to process: ")
    fpath = fpath.replace('\'', "").replace("\"", "")
    files = os.listdir(fpath)
    print("These files found in this directory: ")
    [print(x) for x in files]
    if input("Continue (Y/N)? ") in "Yesyes":
        return fpath
    get_folder_path()


# reference for chopper paper - https://doi.org/10.1093/bioinformatics/btad311
# TODO - move from NanoFilt to chopper

def run_nanofilt(min_length, max_length, quality, fastq, folder_path, **kwargs):
    """
    Runs NanoFilt as a command line process and writes the output to a fastq file
    """
    if fastq.endswith(".gz"):
        p1 = subprocess.Popen(["gunzip", "-c", fastq], stdout=subprocess.PIPE)
        process = subprocess.run([
                "NanoFilt",
                "-l",
                f"{min_length}",
                "--maxlength",
                f"{max_length}",
                "-q",
                f"{quality}"
            ], stdin=p1.stdout, capture_output=True)
        p1.stdout.close()
        p1.wait()
    else:
        process = subprocess.run([
                 "NanoFilt",
                "-l",
                 f"{min_length}",
                 "--maxlength",
                 f"{max_length}",
                 "-q",
                 f"{quality}",
                 fastq
         ], capture_output=True)
    if process_command_output(process, "NanoFilt"):
         write_process_to_file(process, f"{folder_path}/sample_filtered.fastq")


def run_cutadapt(fprimer, rprimer, folder_path, **kwargs):
    """
    Runs cutadapt as a command line process and then writes the output to a new fastq file
    """
    process = subprocess.run([
             "cutadapt",
             "-g",
             fprimer,
             "-a",
             rprimer,
             f"{folder_path}/sample_filtered.fastq"
     ], capture_output=True)
    if process_command_output(process, "cutadapt"):
         write_process_to_file(process, f"{folder_path}/sample_trimmed.fastq")


def run_awk(folder_path, **kwargs):
    """
    Uses awk to convert a fastq file to a fasta file
    """
    process = subprocess.run([
                 "awk",
                 'NR%4==1 {printf ">%s\\n", substr($1,2); next} NR%4==2 {print}',
                 f"{folder_path}/sample_trimmed.fastq"
         ], capture_output=True)

    if process_command_output(process, "awk"):
             write_process_to_file(process, f"{folder_path}/sample_trimmed.fasta")


def run_blast(reference, folder_path, **kwargs):
    """
    runs blast of a fasta file using the given reference and writes to a tsv file
    """
    process = subprocess.run([
            "blastn",
            "-query",
            f"{folder_path}/sample_trimmed.fasta",
            "-subject",
            reference,
            "-out",
            f"{folder_path}/blast_output.tsv",
            "-outfmt",
            "6 std sseq",
            "-evalue",
            "1e-5"
    ], capture_output=True, env=env)
    process_command_output(process, "BLAST")


def process_blast_df(min_match_length: int, min_pct_match: float, query_start: int, folder_path: str, **kwargs)->pd.DataFrame:
    """
    Processes the output of blast to identify whole sequences that match the given arguments
    :param min_match_length: the minimum length of a match to keep
    :param min_pct_match: threshold of the minimum percentage match to keep
    :param query_start: location of where the matches should start
    :param folder_path: path to the folder containing the tsv to process
    """
    fp = f"{folder_path}/blast_output.tsv"
    df1 = pd.read_csv(fp, sep="\t", header=None)
    logging.info(f"BLAST found {df1.shape[0]} matches")
    df1.columns = [
            "query_id",
            "ref_id",
            "pct_match",
            "match_length",
            "mismatches",
            "gaps",
            "query_start",
            "query_end",
            "ref_start",
            "ref_end",
            "evalue",
            "bit_score",
            "ref_seq"
            ]
    match_mask = df1["match_length"] > min_match_length
    pct_mask = df1["pct_match"] > min_pct_match
    query_start_mask = df1["query_start"] > query_start
    filtered_df1 = df1.loc[match_mask & pct_mask & query_start_mask].copy()
    filtered_df1.to_csv(f"{folder_path}/blast_output_filtered.csv", index=False)
    logging.info(f"After filtering {filtered_df1.shape[0]} sequences remain")
    return filtered_df1


def split_df(df: pd.DataFrame, folder_path: str)->pd.DataFrame:
    """
    Takes a dataframe of potential ScFv sequences and splits them into each variable chain, removing the linker
    :param df: dataframe containing the output of blast, post filtering
    :param folder_path: folder_path to which the output will be saved
    """
    filtered_ids = df["query_id"].unique()
    # Use the filtered and trimmed fasta file
    fasta_file = f"{folder_path}/sample_trimmed.fasta"
    # Identify the sequences from the fasta file that are in the filtered ids
    processed_fasta = pd.DataFrame(list(process_fasta_file(fasta_file, filtered_ids)), columns=["query_id", "sequence"])

    # Check that the resulting dataframes are the same size
    if processed_fasta.shape[0] == df.shape[0]:
        logging.info("All IDS found in fasta file")
    else:
        logging.error(f"Processed_fasta: {processed_fasta.shape[0]}\tDataframe: {df.shape[0]}")
        logging.error("Something has gone wrong with ID processing")

    # Merge the dataframes and save them
    merged_df = pd.merge(df, processed_fasta, on="query_id")
    # merged_df.to_csv(f"{folder_path}/merged_df.csv", index=False)

    # This gets the pre linker sequence and the post-linker sequence anmd saves the dataframe as merged_df.csv
    # TODO - is there a more efficient way of doing this, or at least a way of combining into a single apply - yes, passing result_type="expand" expands iterable returns into columns
    merged_df.loc[:, ["pre_linker_sequence", "post_linker_sequence"]] = merged_df.apply(lambda row: [row["sequence"][:row["query_start"]], row["sequence"][row["query_end"]:]], axis=1, result_type="expand")
    # merged_df.loc[:, "pre_linker_sequence"] = merged_df.apply(lambda x: x["sequence"][:x["query_start"]], axis=1)
    # merged_df.loc[:, "post_linker_sequence"] = merged_df.apply(lambda x: x["sequence"][x["query_end"]:], axis=1)
    merged_df.to_csv(f"{folder_path}/merged_df.csv", index=False)
    pd.DataFrame.apply
    # This makes a new fasta file which splits up the pre-linker sequence and the post-linker sequence
    print("Writing new fasta file with split sequences...")
    split_fasta_file = f"{folder_path}/split_fasta.fasta"
    if "split_fasta.fasta" not in os.listdir(folder_path):
        with open(split_fasta_file, "wb") as f:
            for qid in filtered_ids:
                f.write(f">{qid}_pre\n".encode("utf-8"))
                f.write(merged_df.loc[merged_df["query_id"] == qid, "pre_linker_sequence"].values[0].encode("utf-8") + b"\n")
                f.write(f">{qid}_post\n".encode("utf-8"))
                f.write(merged_df.loc[merged_df["query_id"] == qid, "post_linker_sequence"].values[0].encode("utf-8") + b"\n")
    logging.info("Sequences split and new fasta file made")
    return merged_df


# def check_igblast(igblast: str, data_path: str)->bool:
#     """
#     Checks whether IgBLAST is found in the given directory, if not then downloads and extracts.
#     :param igblast: string of the path to where igblast is located or where you want it to be downloaded and extracted
#     :param data_path: sttring of the path to data files, should be the parent directory of igblast
#     """
#     try:
#         os.listdir(igblast)
#         print("IGBLAST already installed")
#         return True

#     except FileNotFoundError:
#         print("Downloading igblast...")
#         process = subprocess.run(["curl","https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-macosx.tar.gz", "-o", f"{data_path}/ncbi-igblast-1.22.0-x64-macosx.tar.gz"], capture_output=True)
#         process_command_output(process, "Download igblast")
#         process = subprocess.run(["tar", "-xvzf", f"{data_path}/ncbi-igblast-1.22.0-x64-macosx.tar.gz", "-C", f"{data_path}"], capture_output=True)
#         process_command_output(process, "Extract igblast")
#         print("Downloading mouse database...")
#         process = subprocess.run(["curl", "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/mouse_gl_VDJ.tar", "-o", f"{data_path}/mouse_gl_VDJ.tar"], capture_output=True)
#         process_command_output(process, "Download igblast database")
#         os.mkdir(f"{igblast}/database")
#         process = subprocess.run(["tar", "-xvf", f"{data_path}/mouse_gl_VDJ.tar", "-C", f"{igblast}/database"], capture_output=True)
#         process_command_output(process, "Extract igblast database")
#         return False


def run_igblast(igblast: str, folder_path: str):
    """
    Runs igblast on the data in folder_path
    :param igblast: path to where igblast is located
    :param folder_path: path to the directory containing the fasta file to be analysed
    """
    process = subprocess.run([
            f"igblastn",
            "-germline_db_V",
            f"{igblast}/database/mouse_gl_V",
            "-germline_db_J",
            f"{igblast}/database/mouse_gl_J",
            "-germline_db_D",
            f"{igblast}/database/mouse_gl_D",
            "-organism",
            "mouse",
            "-query",
            f"{folder_path}/split_fasta.fasta",
            "-num_alignments_V",
            "1",
            "-num_alignments_J",
            "1",
            "-num_alignments_D",
            "1",
            "-outfmt",
            "19",
            "-extend_align5end",
            "-extend_align3end",
            "-auxiliary_data",
            f"{igblast}/optional_file/mouse_gl.aux",
            "-out",
            f"{folder_path}/split_fasta_output.tsv"
        ], capture_output=True, env=env)
    process_command_output(process, "IGBlast")


# TODO - check the stop correction
def process_igblast_results(folder_path: str, aa_mtx:pd.DataFrame):
    fp = f"{folder_path}/split_fasta_output.tsv"
    split_fasta_df = pd.read_csv(fp, sep="\t")
    print(f"{split_fasta_df.shape[0]} total number of sequences to process")
    sequence_mask = split_fasta_df["sequence"].notna()
    filtered_split_fasta_df = split_fasta_df.loc[sequence_mask]

    # Can fwr4 be omitted from this and just add a 'standard' fwr4 at the end?
    conserved_parts = ["fwr1", "fwr2", "fwr3", "fwr4"]
    variable_regions = ["cdr1", "cdr2", "cdr3"]
    # Filters for functional sequences
    #    Only well delimited V-domains (no NAs for the V gene start nor J gene end) from the IgBLAST output are considered.
    # Removing sequences with stop codons
    v_mask = filtered_split_fasta_df["v_alignment_start"].notna()
    logging.info(f"V Mask total: {v_mask.sum()}")
    j_mask = filtered_split_fasta_df["j_alignment_end"].notna()
    logging.info(f"J Mask total: {j_mask.sum()}")
    # Filtering out parts that haven't aligned or been identified in the framework or cdrs
    part_mask = filtered_split_fasta_df[conserved_parts + variable_regions].notna().all(axis=1)

    df2 = filtered_split_fasta_df.loc[sequence_mask & v_mask & j_mask & part_mask & part_len_mask]
    df2.loc["codon_fixed"] = False
    # TODO - find the frame for each fwr
    for part in conserved_parts:
        columns_to_add = [
            f"{part}_germline",
            f"{part}_germline_aa",
            f"{part}_aa_frame",
            f"{part}_germline_aa_frame",
            f"{part}_aa_correct",
            f"{part}_aa_codons"
        ]

        processing_cols = [
            f"{part}_start",
            f"{part}_end",
            "v_sequence_start",
            "v_sequence_end",
            "germline_alignment",
            "germline_alignment_aa",
            part,
        ]

        stop_fixing_cols = [
            f"{part}_germline",
            f"{part}_germline_aa_frame",
            f"{part}_aa_correct",
            f"{part}_aa_codons",
            f"{part}_aa_frame"
        ]
        print(part)
        df2 = df2.join(df2.loc[:, processing_cols].apply(lambda row: part_df_processing(row, part), axis=1, result_type="expand").add_prefix(f"{part}_"))
        stop_mask = df2[f"{part}_aa_correct"].str.contains("\*") & (df2[columns_to_add].notna().all())
        df2 = df2.join(df2.loc[stop_mask, stop_fixing_cols].apply(lambda row: stop_find_df_processing(row, part), axis=1, result_type="expand").add_prefix(f"{part}_"))
        df2.loc[:,"codon_fixed"] = df2[f"codon_fixed"] | df2[f"{part}_codon_fixed"]

    msg = f"Actual stops fixed: {df2["codon_fixed"].sum()}"
    print(msg)
    logging.info(msg)
    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    # Removing the _pre and _post from the sequence ids and putting them in their own column
    df2.loc[:, "sequence_id_clean"] = df2["sequence_id"].str.replace("_pre", "").str.replace("_post", "")
    # Only keeping sequences where there are 2 alignments (heavy and light)
    gdf = df2.groupby("sequence_id_clean").filter(lambda x: len(x) == 2)
    # Getting the unique sequences
    sequence_ids = gdf["sequence_id_clean"].unique()
    logging.info(f"Unique sequence ids: {len(sequence_ids)}")
    # Grouping by sequence ids
    gdf = gdf.groupby("sequence_id_clean")
    return gdf


def fix_codon(codon_list, fix_idx, replacement_codon):
    codon_list[fix_idx] = replacement_codon
    return codon_list


def part_df_processing(row, part):
    part_start = row[f"{part}_start"] - (row["v_sequence_start"] - 1)
    part_end = row[f"{part}_end"] - (row["v_sequence_start"] - 1)
    germline = row["germline_alignment"][int(part_start): int(part_end)]
    germline_aa = row["germline_alignment_aa"][int(part_start) // 3: int(part_end) // 3]
    aa_frame = find_translation_frame(row[part], germline_aa)
    germline_aa_frame = find_translation_frame(germline, germline_aa)
    aa_correct, aa_codons = translate_sequence(row[part], frame=aa_frame, return_codons=True)
    return = {
        f"germline": germline,
        f"germline_aa": germline_aa,
        f"aa_frame": aa_frame,
        f"germline_aa_frame": germline_aa_frame,
        f"aa_correct": aa_correct,
        f"aa_codons": aa_codon
    }



def stop_find_df_processing(row, part):
    germ_translated = translate_sequence(row[f"{part}_germline"], int(row[f"{part}_germline_aa_frame"]), return_codons=True)
    alignment_a, alignment_b = sequence_alignment(row[f"{part}_aa_correct"], germ_translated[0], aa_mtx, d=-0.1)
    stop_loc = alignment_b.find("*")
    new_part = None
    codon_fixed = False
    if stop_loc != -1:
        if alignment_a[stop_loc] not in ["X", "-", "*"]:
            stop_codon_loc = len(alignment_b[:stop_loc + 1].replace("-", "")) - 1
            replacement_codon = germ_translated[1][len(alignment_a[:stop_loc + 1].replace("-", ""))-1]
            fixed_codons = fix_codon(row[f"{part}_aa_codons"], stop_codon_loc, replacement_codon)
            new_part = row[f"{part}"][:row[f"{part}_aa_frame"] - 1] + "".join(fixed_codons)
            codon_fixed = True
    return {
        "germline_aa_correct": germ_translated[0],
        "germline_aa_codons": germ_translated[1],
        "alignment_a": alignment_a,
        "alignment_b": alignment_b,
        "new": new_part,
        "codon_fixed": codon_fixed
    }




def find_translation_frame(query_seq, ref_seq):
    """
    Finds translation name based on similarity to reference aa sequence
    """
    sequence_match = 0
    best_frame = None
    for frame in range(1, 4):
        match = SequenceMatcher(None, ref_seq, translate_sequence(query_seq, frame)).ratio()
        if match > sequence_match:
            sequence_match = match
            best_frame = frame
    return best_frame

# TODO - add in cdr consensus identification
def cdr_consensus(df):
    df.loc[:, "all_cdr"] = df["cdr1"].astype(str) + df["cdr2"].astype(str) + df["cdr3"].astype(str)
    fasta_string = df.apply(lambda row: f">{row['sequence_id']}\n{row['all_cdr']}\n", axis=1)
    with open("cdr_sequences.fasta", "w") as f:
        f.writelines(fasta_string)

    subprocess.run(["vsearch", "--cluster_fast", "cdr_sequences.fasta", "--id", "0.98", "--centroids", "cdr_centroids.fasta", "--uc", "cdr_clusters.uc"], capture_output=True, env=env)
    # vsearch --cluster_fast cdr_concat.fasta --id 0.98 --centroids cdr_centroids.fasta --uc cdr_clusters.uc

"""
vsearch output:
cdr_clusters.uc headings:
- Record type - H is a hit; S is centroid; C is cluster; N is no hit
- Cluster number
- Sequence length
- Identity
- Strand
- Query Start
- Target Start
- Alignment length
- Query label
- Target label
- Alignment

Enabling vsearch and filtering out cdrs that contain stop codons being used as centroids:
Filter centroids_raw.fasta:

    Purpose: Identify which of the chosen centroids contain stop codons.

    Action: Run the Python script (or similar) on centroids_raw.fasta to create centroids_filtered.fasta (containing only centroids without stop codons). Also, identify the IDs of the centroids you discarded.

Re-cluster Reads from "Bad" Centroid Clusters:

    Purpose: Take all the reads that clustered around a "bad" centroid (one with a stop codon) and try to re-assign them to a "good" centroid, or see if they form new "good" clusters.

    Action:

        Parse clusters_raw.uc. For each 'H' record whose Target label (field 10) corresponds to a "bad" centroid ID, extract the Query label (field 9). Collect all these query labels into a new FASTA file, say reads_from_bad_clusters.fasta.

        Now, re-cluster reads_from_bad_clusters.fasta against your centroids_filtered.fasta using --usearch_global as in Method 1. Or, if you want them to potentially form new "good" clusters, run vsearch --cluster_fast on reads_from_bad_clusters.fasta by itself (this is risky, as it could generate more "bad" centroids).

        A more robust approach here would be to merge reads_from_bad_clusters.fasta with any unclustered_reads.fasta from the first pass and then re-run cluster_fast or usearch_global against all good centroids, combining all information.

Then use MAFTT to identify consensus sequence for each cluster
"""


# TODO - optimise this to use in the main code
def fix_stops(df, conserved_parts):
    conserved_mask = np.zeros(len(df), dtype=bool)
    for part in conserved_parts:
        in_range = (df["stop_location"] > df[f"{part}_start"]) & (df["stop_location"] < df[f"{part}_end"])
        conserved_mask |= in_range
        df.loc[in_range, "part_to_change"] = part
        # part_name[in_range & (pd.isnull(part_name))] = part
        df.loc[in_range, part] = df.loc[in_range & stop_fix_mask].apply(fix_codon, axis=1)
        df.loc[in_range, "stop_fixed"] = True


# igblast_columns = ['sequence_id', 'sequence', 'sequence_aa', 'locus', 'stop_codon',
#        'vj_in_frame', 'v_frameshift', 'productive', 'rev_comp', 'complete_vdj',
#        'd_frame', 'v_call', 'd_call', 'j_call', 'sequence_alignment',
#        'germline_alignment', 'sequence_alignment_aa', 'germline_alignment_aa',
#        'v_alignment_start', 'v_alignment_end', 'd_alignment_start',
#        'd_alignment_end', 'j_alignment_start', 'j_alignment_end',
#        'v_sequence_alignment', 'v_sequence_alignment_aa',
#        'v_germline_alignment', 'v_germline_alignment_aa',
#        'd_sequence_alignment', 'd_sequence_alignment_aa',
#        'd_germline_alignment', 'd_germline_alignment_aa',
#        'j_sequence_alignment', 'j_sequence_alignment_aa',
#        'j_germline_alignment', 'j_germline_alignment_aa', 'fwr1', 'fwr1_aa',
#        'cdr1', 'cdr1_aa', 'fwr2', 'fwr2_aa', 'cdr2', 'cdr2_aa', 'fwr3',
#        'fwr3_aa', 'fwr4', 'fwr4_aa', 'cdr3', 'cdr3_aa', 'junction',
#        'junction_length', 'junction_aa', 'junction_aa_length', 'v_score',
#        'd_score', 'j_score', 'v_cigar', 'd_cigar', 'j_cigar', 'v_support',
#        'd_support', 'j_support', 'v_identity', 'd_identity', 'j_identity',
#        'v_sequence_start', 'v_sequence_end', 'v_germline_start',
#        'v_germline_end', 'd_sequence_start', 'd_sequence_end',
#        'd_germline_start', 'd_germline_end', 'j_sequence_start',
#        'j_sequence_end', 'j_germline_start', 'j_germline_end', 'fwr1_start',
#        'fwr1_end', 'cdr1_start', 'cdr1_end', 'fwr2_start', 'fwr2_end',
#        'cdr2_start', 'cdr2_end', 'fwr3_start', 'fwr3_end', 'fwr4_start',
#        'fwr4_end', 'cdr3_start', 'cdr3_end', 'np1', 'np1_length', 'np2',
#        'np2_length']


def check_fw(x)->int:
    """
    Checks the difference in the identified stop codon locations between the sequence and the germline alignment to identify the specific base that's different
    """
    stop_loc = int(x["stop_location"])
    sequence_start = int(x["v_sequence_start"])
    sequence_stop_codon = x["sequence"][stop_loc + sequence_start - 1:stop_loc + sequence_start - 1 + 3]
    germline_codon = x["germline_alignment"][stop_loc:stop_loc + 3]
    if "N" in germline_codon:
        return None
    if (len(sequence_stop_codon) < 3) or (len(germline_codon) < 3):
        return None
    diff = compare_sequences(sequence_stop_codon, germline_codon)
    if len(diff) > 0:
        return int(diff[0])
    return None


def fix_codon(x)->str:
    """
    Fixes a stop codon in the constant region of a sequence by mutating the differing base to that found in the germline sequence and returns the "fixed" region
    """
    seq = x["sequence"]
    germ_seq = x["germline_alignment"]
    start = int(x[f"{x['part_to_change']}_start"]) - 1
    end = int(x[f"{x['part_to_change']}_end"])
    stop_off = int(x["codon_change_pos"])
    stop = int(x["stop_location"]) + stop_off
    return seq[start:stop] + germ_seq[stop] + seq[stop + 1:end]


def find_scfv_sequences(grouped_df, sequence_parts: list, linker_seq: str)->pd.DataFrame:
    """
    Finds the valid variable chain sequences in the grouped dataframe that contain all the parts defined in sequence_parts
    :param gdf: Grouped Dataframe grouped by sequence id, should be the output of process_igblast_results()
    :param sequence_parts: a list of all the parts that a valid sequence should contain, probably a minimum of ["cdr1", "cdr2", "cdr3"]
    :param linker_seq: dna sequence of the linker you want to separate heavy and light chains
    """
    sequence_ids = list(grouped_df.groups.keys())
    results_df = pd.DataFrame(index=sequence_ids, columns=["pre_seq", "post_seq", "ScFv_seq", "orientation"])
    interval = max(1, len(sequence_ids) // 50)
    results = []
    for i, (seqid, tdf) in enumerate(grouped_df):
        try:
            pre_row = tdf[tdf["sequence_id"] == f"{seqid}_pre"].iloc[0]
            post_row = tdf[tdf["sequence_id"] == f"{seqid}_post"].iloc[0]

            pre_seq = "".join([pre_row[part] for part in sequence_parts])
            post_seq = "".join([post_row[part] for part in sequence_parts])

            pre_trans_frames = check_translation(pre_seq)
            if len(pre_trans_frames) == 1:
                pre_seq = pre_seq[pre_trans_frames[0]-1:]

            post_trans_frames = check_translation(post_seq)
            if len(post_trans_frames) == 1:
                post_seq = post_seq[post_trans_frames[0]-1:]

            scfv_seq = linker_seq.join([pre_seq, post_seq])
            orientation = f"{pre_row['locus']}-{post_row['locus']}"
            results.append({
              "sequence_id": seqid,
              "pre_seq": pre_seq,
              "post_seq": post_seq,
              "ScFv_seq": scfv_seq,
              "orientation": orientation,
              "translated_seq": translate_sequence(scfv_seq),
              })

            if i % interval == 0:
                print(f"{i}/{len(sequence_ids)} completed")
        except Exception as e:
            logging.warning(f"Exception with entry: {i}\n\t Sequence ID: {seqid} \n\tError Message: {e}")
    results_df = pd.DataFrame(results)
    return results_df


def check_stops(df):
    """
    Post find check to make sure no stop codons have slipped through
    """
    df.loc[:, "translated_seq"] = df.apply(lambda x: translate_sequence(x["ScFv_seq"]), axis=1)
    df.loc[:, "contains_stop"] = df["translated_seq"].str.contains("\*")
    logging.info(f"Number of final sequences with stop codons: {df['contains_stop'].sum()}")


def main():
    logging.basicConfig(
        filename="scfv_finder.log",
        level=logging.INFO,
        format="%(asctime)s:%(levelname)s:\n\t%(message)s"
        )
    config_fp = None
    if input("Use a custom config file? ") in "Yesyes":
        config_fp = input("Enter the full path to your config file: ")

    config = parse_config_file(config_fp)
    sequence_processing = "SEQUENCE_PROCESSING"
    results_processing = "RESULTS_PROCESSING"

    tools = ["blastn", "igblast"]
    print("Checking tools...")
    for tool in tools:
        if not check_tools(tool):
            print(f"{tool} not found in path")
            if input(f"Install {tool} (Y/N)?") in "Yesyes":
                install_tool(tool)
                logging.info(f"{tool} was installed")
            else:
                tool_path = input("Enter the path to the tool \"bin\" folder: ")
                env["PATH"] = tool_path + os.pathsep + env["PATH"]

    igblast = subprocess.run(["whereis", "igblastn"], capture_output=True, env=env).stdout.decode()
    igblast = igblast[igblast.find(":") + 2:igblast.rfind("/bin/")]
    env["IGDATA"] = igblast

    folder_path = config[sequence_processing].get("folder_path")
    # fastq = [x for x in os.listdir(folder_path) if ".fastq" in x]
    # reference_fasta = input("Enter the full path of the fasta file to use as reference for BLAST, this should be just the linker sequence between the two variable chains: ")
    # reference_fasta = reference_fasta.strip().replace('\'', "").replace("\"", "")
    # TODO - Move this to a config file

    variables = {
        "min_length": config[sequence_processing].getint("min_length"),
        "max_length": config[sequence_processing].getint("max_length"),
        "quality": config[sequence_processing].getint("quality"),
        "fprimer": config[sequence_processing].get("fprimer"),
        "rprimer": config[sequence_processing].get("rprimer"),
        "fastq": config[sequence_processing].get("fastq"), # Path to the fastq file
        "reference": config[sequence_processing].get("reference"), # Path to the reference fasta file
        "folder_path": folder_path
    }
    print("Running Nanofilt...")
    run_nanofilt(**variables)
    print("Nanofilt successful", "Running Cutadapt...", sep="\n")
    run_cutadapt(**variables)
    print("Cutadapt run successfully", "Running awk...", sep="\n")
    run_awk(folder_path)
    print("Awk successfully run", "Running BLAST to find sequences that have a central linker...", sep="\n")
    run_blast(**variables)
    print("BLAST run successfully\nFiltering BLAST dataset...")
    # TODO - put this in config with above variables

    filtered_df = process_blast_df(
        min_match_length=config[results_processing].getint("min_match_length"),
        min_pct_match=config[results_processing].getfloat("min_pct_match"),
        query_start=config[results_processing].getint("query_start"),
        folder_path=folder_path
      )
    print("Filtering completed", "Splitting ScFv chains...", sep="\n")
    merged_df = split_df(filtered_df, folder_path)
    print("Splitting completed", "Running IgBLAST", sep="\n")
    run_igblast(igblast, folder_path)
    print("IgBLAST Successfully run")
    aa_mtx = pd.read_csv(config[results_processing].get("aa_mtx"), index_col=0)
    gdf = process_igblast_results(folder_path, aa_mtx)
    print("IgBLAST Results filtered for complete variable sequences")

    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    linker_seq = config["SEQUENCE_PARTS"].get("linker_seq")

    results_df = find_scfv_sequences(gdf, sequence_parts, linker_seq)
    logging.info(f"{results_df.shape[0]} ScFv sequences found")
    check_stops(results_df)
    logging.info("sequences checked for stop codons")

    logging.info(f"Results saved to: {folder_path}/results_df.csv")
    logging.info(f"{results_df['ScFv_seq'].unique().shape[0]} ScFv sequences found")
    print(results_df["ScFv_seq"].unique().shape[0], "ScFv sequences found")
    # Gets the counts for each amino acid sequence
    counts = pd.DataFrame(results_df["translated_seq"].value_counts())
    # Assigns a unique id to each sequence - this is a hash of the sequence so in theory translatable across runs
    counts.loc[:, "id"] = [hashlib.sha256(x.encode("utf-8")).hexdigest() for x in counts.index]
    results_df.join(counts["id"], on="translated_seq")
    results_df.to_csv(f"{folder_path}/results_df.csv", index=False)
    counts.reset_index(inplace=True)
    counts.columns = ["translated_seq", "count", "id"]
    counts.to_csv(f"{folder_path}/count_results.csv", index=False)





if __name__ == "__main__":
        main()
        print("All done!")
