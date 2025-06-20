import os
import subprocess
from utils import *
import numpy as np
import logging


# TODO - Add logging to a logfile instead of printing - done this, kept some prints for users to see progress
# TODO - check if blastn is in path, if not then ask to download  - Donne this
# TODO - Add multi-platform support, at the moment just MacOS is supported


env = os.environ.copy()


def check_tools(tool):
    try:
        subprocess.run([tool])
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
        logging.error(f"Processed_fasta: {processed_fasta.shape[0]}", f"Dataframe: {df.shape[0]}")
        logging.error("Something has gone wrong with ID processing")

    # Merge the dataframes and save them
    merged_df = pd.merge(df, processed_fasta, on="query_id")
    # merged_df.to_csv(f"{folder_path}/merged_df.csv", index=False)

    # This gets the pre linker sequence and the post-linker sequence anmd saves the dataframe as merged_df.csv
    merged_df.loc[:, "pre_linker_sequence"] = merged_df.apply(lambda x: x["sequence"][:x["query_start"]], axis=1)
    merged_df.loc[:, "post_linker_sequence"] = merged_df.apply(lambda x: x["sequence"][x["query_end"]:], axis=1)
    merged_df.to_csv(f"{folder_path}/merged_df.csv", index=False)

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


def process_igblast_results(folder_path: str):
    fp = f"{folder_path}/split_fasta_output.tsv"
    split_fasta_df = pd.read_csv(fp, sep="\t")
    print(f"{split_fasta_df.shape[0]} total number of sequences to process")
    sequence_mask = split_fasta_df["sequence"].notna()
    filtered_split_fasta_df = split_fasta_df.loc[sequence_mask]
    # Filters for functional sequences
    #    Only well delimited V-domains (no NAs for the V gene start nor J gene end) from the IgBLAST output are considered.
    # Removing sequences with stop codons
    v_mask = filtered_split_fasta_df["v_alignment_start"].notna()
    logging.info(f"V Mask total: {v_mask.sum()}")
    j_mask = filtered_split_fasta_df["j_alignment_end"].notna()
    logging.info(f"J Mask total: {j_mask.sum()}")
    # Fixing stop codons in conserved regions that are likely sequencing errors
    stop_find = filtered_split_fasta_df["sequence_aa"].str.find("*")
    stop_mask2 = (stop_find != -1) & ((stop_find * 3) < filtered_split_fasta_df["fwr4_end"])
    stop_fix_mask3 = filtered_split_fasta_df["sequence_aa"].str.count("\*") < 2
    stop_fix_mask = stop_mask2 & stop_fix_mask3
    logging.info(f"Potential stops to fix: {stop_fix_mask.sum()}")
    filtered_split_fasta_df.loc[stop_fix_mask, "stop_location"] = stop_find * 3
    filtered_split_fasta_df.loc[stop_fix_mask, "codon_change_pos"] = filtered_split_fasta_df.loc[stop_fix_mask].apply(check_fw, axis=1)
    stop_fix_mask = stop_fix_mask & filtered_split_fasta_df["codon_change_pos"].notna()
    conserved_parts = ["fwr1", "fwr2", "fwr3", "fwr4"]
    conserved_mask = np.zeros(len(filtered_split_fasta_df), dtype=bool)
    # part_name = np.full(len(filtered_split_fasta_df), None)
    for part in conserved_parts:
        in_range = (filtered_split_fasta_df["stop_location"] > filtered_split_fasta_df[f"{part}_start"]) & (filtered_split_fasta_df["stop_location"] < filtered_split_fasta_df[f"{part}_end"])
        conserved_mask |= in_range
        filtered_split_fasta_df.loc[in_range, "part_to_change"] = part
        # part_name[in_range & (pd.isnull(part_name))] = part
        filtered_split_fasta_df.loc[in_range & stop_fix_mask, part] = filtered_split_fasta_df.loc[in_range & stop_fix_mask].apply(fix_codon, axis=1)
        filtered_split_fasta_df.loc[in_range & stop_fix_mask, "stop_fixed"] = True
    msg = f"Actual stops fixed: {conserved_mask.sum()}"
    print(msg)
    logging.info(msg)
    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    parts_mask = filtered_split_fasta_df[sequence_parts].isna().any(axis=1)
    filtered_split_fasta_df = filtered_split_fasta_df.loc[v_mask & j_mask & ~parts_mask | conserved_mask]
    print(f"{filtered_split_fasta_df.shape[0]} number of sequences left after filtering")
    # Removing the _pre and _post from the sequence ids and putting them in their own column
    filtered_split_fasta_df.loc[:, "sequence_id_clean"] = filtered_split_fasta_df["sequence_id"].str.replace("_pre", "").str.replace("_post", "")
    # Only keeping sequences where there are 2 alignments (heavy and light)
    gdf = filtered_split_fasta_df.groupby("sequence_id_clean").filter(lambda x: len(x) == 2)
    # Getting the unique sequences
    sequence_ids = gdf["sequence_id_clean"].unique()
    logging.info(f"Unique sequence ids: {len(sequence_ids)}")
    # Grouping by sequence ids
    gdf = gdf.groupby("sequence_id_clean")
    return gdf


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
    df.loc[:, "contains_stop"] = df["ScFv_seq"].str.contains("\*")
    logging.info(f"Number of final sequences with stop codons: {df['contains_stop'].sum()}")


def main():
    logging.basicConfig(
        filename="scfv_finder.log",
        level=logging.INFO,
        format="%(asctime)s:%(levelname)s:\n\t%(message)s"
        )
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

    folder_path = get_folder_path()
    fastq = [x for x in os.listdir(folder_path) if ".fastq" in x]
    reference_fasta = input("Enter the full path of the fasta file to use as reference for BLAST, this should be just the linker sequence between the two variable chains: ")
    reference_fasta = reference_fasta.strip().replace('\'', "").replace("\"", "")
    # TODO - Move this to a config file
    variables = {
        "min_length": 700, # minimum length of sequences you want to keep
        "max_length": 800, # maximum length of sequences you want to keep
        "quality": 15, # minimum quality of sequences you want to keep - fullcircle give you a graph showing quality, wouldn't recommend going lower than 10
        "fprimer": "GTCCCTGGCTCCACTGGA", # Sequence of the forward primer used to amplify the PCR fragment
        "rprimer": "GCGCTGGCGTCGTGGT", # Sequence of the reverse primer used to amplify the PCR fragment
        "fastq": f"{folder_path}/{fastq[0]}", # Path to the fastq file
        "reference": reference_fasta, # Path to the reference fasta file
        "folder_path": folder_path
    }
    print("Running Nanofilt...")
    run_nanofilt(**variables)
    print("Nanofilt successful\nRunning Cutadapt...")
    run_cutadapt(**variables)
    print("Cutadapt run successfully\nRunning awk...")
    run_awk(folder_path)
    print("Awk successfully run\nRunning BLAST to find sequences that have a central linker...")
    run_blast(**variables)
    print("BLAST run successfully\nFiltering BLAST dataset...")
    # TODO - put this in config with above variables
    filtered_df = process_blast_df(
        min_match_length=40,
        min_pct_match=95,
        query_start=300,
        folder_path=folder_path
      )
    print("Filtering completed\nSplitting ScFv chains...")
    merged_df = split_df(filtered_df, folder_path)
    print("Splitting completed\nPreparing IGBLAST...")
    # igblast = f"{data_path}/ncbi-igblast-1.22.0"
    # if not check_igblast(igblast, data_path):
    #     print("IGBLAST successfully downloaded")
    print("Running IgBLAST...")

    run_igblast(igblast, folder_path)
    print("IgBLAST Successfully run")
    gdf = process_igblast_results(folder_path)
    print("IgBLAST Results filtered for complete variable sequences")

    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    linker_seq = merged_df["ref_seq"].unique()[0]

    results_df = find_scfv_sequences(gdf, sequence_parts, linker_seq)
    logging.info("ScFv sequences found")
    check_stops(results_df)
    logging.info("sequences checked for stop codons")
    results_df.to_csv(f"{folder_path}/results_df.csv", index=False)
    logging.info(f"Results saved to: {folder_path}/results_df.csv")
    print(results_df.shape[0])
    logging.info(f"{results_df['ScFv_seq'].unique().shape[0]} ScFv sequences found")
    print(results_df["ScFv_seq"].unique().shape[0], "ScFv sequences found")




if __name__ == "__main__":
        main()
        print("All done!")
