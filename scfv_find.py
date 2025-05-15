import os
import pprint
import subprocess
from utils import *
import numpy as np


# This is the path to the folder that has all the files in, change it if you need to
# folder_path = "/Users/rivanishah/Desktop/Alk_ScFv/Alk5-8"
env = os.environ.copy()




def get_folder_path():
    """
    """
    fpath = input("Enter the full path to the folder containing the fastq file you want to process: ")
    files = os.listdir(fpath)
    print("These files found in this directory: ")
    [print(x) for x in files]
    if input("Continue (Y/N)? ") in "Yesyes":
        return fpath
    get_folder_path()





def run_nanofilt(min_length, max_length, quality, fastq, folder_path, **kwargs):
    process = subprocess.run([
             "NanoFilt",
            "-l",
             f"{min_length}",
             "--maxlength",
             f"{max_length}",
             "-q",
             f"{quality}",
             fastq,
     ], capture_output=True)
    if process_command_output(process, "NanoFilt"):
         write_process_to_file(process, f"{folder_path}/sample_filtered.fastq")


def run_cutadapt(fprimer, rprimer, folder_path, **kwargs):
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
    process = subprocess.run([
                 "awk",
                 'NR%4==1 {printf ">%s\\n", substr($1,2); next} NR%4==2 {print}',
                 f"{folder_path}/sample_trimmed.fastq"
         ], capture_output=True)

    if process_command_output(process, "awk"):
             write_process_to_file(process, f"{folder_path}/sample_trimmed.fasta")


def run_blast(reference, folder_path, **kwargs):
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
    ], capture_output=True)
    process_command_output(process, "BLAST")


def process_blast_df(min_match_length, min_pct_match, query_start, folder_path, **kwargs):
    fp = f"{folder_path}/blast_output.tsv"
    df1 = pd.read_csv(fp, sep="\t", header=None)
    print(f"BLAST found {df1.shape[0]} matches")
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
    print(f"After filtering {filtered_df1.shape[0]} sequences remain")
    return filtered_df1


def split_df(df, folder_path):
    filtered_ids = df["query_id"].unique()
    # Use the filtered and trimmed fasta file
    fasta_file = f"{folder_path}/sample_trimmed.fasta"
    # Identify the sequences from the fasta file that are in the filtered ids
    processed_fasta = pd.DataFrame(list(process_fasta_file(fasta_file, filtered_ids)), columns=["query_id", "sequence"])

    # Check that the resulting dataframes are the same size
    if processed_fasta.shape[0] == df.shape[0]:
        print("All IDS found in fasta file")
    else:
        print("Something has gone wrong with ID processing")

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
    print("Sequences split and new fasta file made")
    return merged_df


def check_igblast(igblast, data_path):
    try:
        os.listdir(igblast)
        print("IGBLAST already installed")
        return True

    except FileNotFoundError:
        print("Downloading igblast...")
        process = subprocess.run(["curl","https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-macosx.tar.gz", "-o", f"{data_path}/ncbi-igblast-1.22.0-x64-macosx.tar.gz"], capture_output=True)
        process_command_output(process, "Download igblast")
        process = subprocess.run(["tar", "-xvzf", f"{data_path}/ncbi-igblast-1.22.0-x64-macosx.tar.gz", "-C", f"{data_path}"], capture_output=True)
        process_command_output(process, "Extract igblast")
        print("Downloading mouse database...")
        process = subprocess.run(["curl", "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/mouse_gl_VDJ.tar", "-o", f"{data_path}/mouse_gl_VDJ.tar"], capture_output=True)
        process_command_output(process, "Download igblast database")
        os.mkdir(f"{igblast}/database")
        process = subprocess.run(["tar", "-xvf", f"{data_path}/mouse_gl_VDJ.tar", "-C", f"{igblast}/database"], capture_output=True)
        process_command_output(process, "Extract igblast database")
        return False


def run_igblast(igblast, folder_path):
        process = subprocess.run([
                f"{igblast}/bin/igblastn",
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


def process_igblast_results(folder_path):
    fp = f"{folder_path}/split_fasta_output.tsv"
    split_fasta_df = pd.read_csv(fp, sep="\t")
    print(f"{split_fasta_df.shape[0]} total number of sequences to process")
    sequence_mask = split_fasta_df["sequence"].notna()
    filtered_split_fasta_df = split_fasta_df.loc[sequence_mask]
    # Filters for functional sequences
    #    Only well delimited V-domains (no NAs for the V gene start nor J gene end) from the IgBLAST output are considered.
    # Removing sequences with stop codons
    v_mask = filtered_split_fasta_df["v_alignment_start"].notna()
    print(f"{v_mask.sum()} v mask kept")
    j_mask = filtered_split_fasta_df["j_alignment_end"].notna()
    print(f"{j_mask.sum()} j mask kept")
    # TODO - put the fixing stops in sequence here?
    stop_find = filtered_split_fasta_df["sequence_aa"].str.find("*")
    stop_mask2 = (stop_find != -1) & ((stop_find * 3) < filtered_split_fasta_df["fwr4_end"])
    stop_fix_mask3 = filtered_split_fasta_df["sequence_aa"].str.count("\*") < 2
    # stop_fix_mask4 = stop_find.notna()
    # stop_fix_mask4 =
    stop_fix_mask = stop_mask2 & stop_fix_mask3
    print(f"{stop_fix_mask.sum()} potential Stops to fix")
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
    # df.loc[:, "part_to_change"] = part_name
    print(f"{conserved_mask.sum()} stops fixed")
    stop_fix_mask = stop_fix_mask & conserved_mask
    # print(filtered_split_fasta_df.loc[stop_fix_mask, "new_fwr3"].head())
    # print(f"{filtered_split_fasta_df.loc[stop_fix_mask].shape[0]} actual Stops to fix")
    # stop_mask = filtered_split_fasta_df["sequence_aa"].str.contains("\*")
    # print(f"{stop_mask.sum()} sequence with stops")
    # stop_find = filtered_split_fasta_df["sequence_aa"].str.find("*")
    # stop_mask2 = (stop_find != -1) & ((stop_find * 3) < filtered_split_fasta_df["fwr4_end"])
    # print(f"{stop_mask2.sum()} stops in the sequence")
    # stop_mask = stop_mask | filtered_split_fasta_df["v_sequence_alignment_aa"].isna()
    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    parts_mask = filtered_split_fasta_df[sequence_parts].isna().any(axis=1)
    # filtered_split_fasta_df = filtered_split_fasta_df.loc[v_mask & j_mask & ~stop_mask & ~stop_mask2 & ~parts_mask]
    # no_stops_mask = filtered_split_fasta_df["stop_location"].isna()
    # Filter out stop codons that are not in conserved places - "stop_location".isna()
    filtered_split_fasta_df = filtered_split_fasta_df.loc[v_mask & j_mask & ~parts_mask | conserved_mask]
    print(f"{filtered_split_fasta_df.shape[0]} number of sequences left after filtering")
    # Removing the _pre and _post from the sequence ids and putting them in their own column
    filtered_split_fasta_df.loc[:, "sequence_id_clean"] = filtered_split_fasta_df["sequence_id"].str.replace("_pre", "").str.replace("_post", "")
    # return filtered_split_fasta_df
    # Only keeping sequences where there are 2 alignments (heavy and light)
    gdf = filtered_split_fasta_df.groupby("sequence_id_clean").filter(lambda x: len(x) == 2)
    # Getting the unique sequences
    sequence_ids = gdf["sequence_id_clean"].unique()
    print(len(sequence_ids), "unique sequence ids")
    # Grouping by sequence ids
    gdf = gdf.groupby("sequence_id_clean")
    return gdf


# const_sequences = ["fwr1", "fwr2", "fwr3", "fwr4"]
# germline_sequence = "v_germline_alignment"

# columns = ['sequence_id', 'sequence', 'sequence_aa', 'locus', 'stop_codon',
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


def check_fw(x):
    stop_loc = int(x["stop_location"])
    sequence_start = int(x["v_sequence_start"])
    sequence_stop_codon = x["sequence"][stop_loc + sequence_start - 1:stop_loc + sequence_start - 1 + 3]
    germline_codon = x["germline_alignment"][stop_loc:stop_loc + 3]
    # print(sequence_stop_codon, germline_codon)
    if "N" in germline_codon:
        return None
    # print(len(sequence_stop_codon), len(germline_codon))
    if (len(sequence_stop_codon) < 3) or (len(germline_codon) < 3):
        return None
    diff = compare_sequences(sequence_stop_codon, germline_codon)
    if len(diff) > 0:
        return int(diff[0])
    return None

def fix_codon(x):
    seq = x["sequence"]
    germ_seq = x["germline_alignment"]
    start = int(x[f"{x['part_to_change']}_start"]) - 1
    end = int(x[f"{x['part_to_change']}_end"])
    stop_off = int(x["codon_change_pos"])
    stop = int(x["stop_location"]) + stop_off
    return seq[start:stop] + germ_seq[stop] + seq[stop + 1:end]
    # stop_df.loc[stop_loc, f"new_{fw}"] = [seq[int(start):int(stop)-1] + germ_seq[int(stop)] + seq[int(stop):int(end)] for seq, germ_seq, start, stop, end in stop_df.loc[stop_loc,["sequence", "germline_alignment", f"{fw}_start", f"stop_location", f"{fw}_end"]].values]


# # If stop codon is in conserved regions of the sequencee
# # Find the location of the stop codon in the sequence_alignment_aa -
# stop_mask = df["sequence_aa"].str.count("\*") < 2
# stop_df = df.loc[stop_mask]
# aa_loc = stop_df.loc["sequence_aa"].str.find("*")
# dna_loc = aa_loc * 3 # Change this to identify where in the codon the change is
# stop_df.loc[:, "stop_location"] = dna_loc
# for fw in ["fwr1", "fwr2", "fwr3", "fwr4"]:
#      stop_loc = (stop_df["stop_location"] > stop_df[f"{fw}_start"]) & (stop_df["stop_location"] < stop_df[f"{fw}_end"])
#      stop_df.loc[stop_loc, f"new_{fw}"] = [seq[int(start):int(stop)-1] + germ_seq[int(stop)] + seq[int(stop):int(end)] for seq, germ_seq, start, stop, end in stop_df.loc[stop_loc,["sequence", "germline_alignment", f"{fw}_start", f"stop_location", f"{fw}_end"]].values]
#      # stop_df.loc[stop_loc, "stop_part"] = fw
#      # stop_df.loc[stop_loc, "fw_stop_start"] = stop_df.loc[stop_loc, f"{fw}_start"]
#      # stop_df.loc[:, "fwr_stop_loc"] = stop_df["stop_location"] - stop_df["fw_stop_start"] + stop_df["v_sequence_start"]


#   fw_loc = df["stop_location"].values - df[df["stop_part"] + "_start"].T[0].values + 1 # this is wrong, just gives the locations for the first sequence
#   df[fw][fw_loc] = df["germline_alignment"][dna_loc]
# Compare the sequence to the germline and swap - need to change it in the "part" sequence as that's what gets taken forward
# Find the start and end of the region affected - fwrx_start and fwrx_end
# slice germline_alignment with these (remember to minus 1 from the start) - df.loc[seqid, "germline_alignment"].values[0][fwrx_start-1:fwrx_end]
# Translate germline and then only change the sequence where there is a difference



# seq1 = translate_sequence(sequence_part)
# seq2 = translate_sequence(germline_part)
# differences = compare_sequences(seq1, seq2)
# for diff in differences:
#     sequence_part[diff] = germline_part[diff]





def find_scfv_sequences(grouped_df, sequence_parts: list, linker_seq: str):
    """
    Finds the valid variable chain sequences in the grouped dataframe that contain all the parts defined in sequence-parts
    :param gdf: Dataframe grouped by sequence id, should be the output of process_igblast_results()
    :param sequence_parts: a list of all the parts that a valid sequence should contain, probably a minimum of ["cdr1", "cdr2", "cdr3"]
    :param linker_seq: dna sequence of the linker you want to separate heavy and light chains
    """
    sequence_ids = list(grouped_df.groups.keys())
    results_df = pd.DataFrame(index=sequence_ids, columns=["pre_seq", "post_seq", "ScFv_seq", "orientation"])
    interval = max(1, len(sequence_ids) // 50)
    results = []
    for i, (seqid, tdf) in enumerate(grouped_df):
        try:
            # tdf = grouped_df.get_group(seqid)
            pre_row = tdf[tdf["sequence_id"] == f"{seqid}_pre"].iloc[0]
            post_row = tdf[tdf["sequence_id"] == f"{seqid}_post"].iloc[0]

            pre_seq = "".join([pre_row[part] for part in sequence_parts])
            post_seq = "".join([post_row[part] for part in sequence_parts])

            # pre_seq = "".join(tdf.loc[tdf["sequence_id"] == f"{seqid}_pre", sequence_parts].values.flatten())

            pre_trans_frames = check_translation(pre_seq)
            if len(pre_trans_frames) == 1:
                pre_seq = pre_seq[pre_trans_frames[0]-1:]

            post_trans_frames = check_translation(post_seq)
            if len(post_trans_frames) == 1:
                post_seq = post_seq[post_trans_frames[0]-1:]

            # pre_seq_chain = tdf.loc[tdf["sequence_id"] == f"{seqid}_pre", "locus"].values[0]
            # post_seq = "".join(tdf.loc[tdf["sequence_id"] == f"{seqid}_post", sequence_parts].values.flatten())

            # post_seq_chain = tdf.loc[tdf["sequence_id"] == f"{seqid}_post", "locus"].values[0]
            scfv_seq = linker_seq.join([pre_seq, post_seq])
            # print(scfv_seq)
            orientation = f"{pre_row['locus']}-{post_row['locus']}"
            # results_df.loc[seqid] = [pre_seq, post_seq, scfv_seq, orientation]
            results.append({
              "sequence_id": seqid,
              "pre_seq": pre_seq,
              "post_seq": post_seq,
              "ScFv_seq": scfv_seq,
              "orientation": orientation
              })

            if i % interval == 0:
                print(f"{i}/{len(sequence_ids)} completed")
        except Exception as e:
            print(f"Exception with entry: {i}, error: {e}")
    results_df = pd.DataFrame(results).set_index("sequence_id")
    return results_df


def check_stops(df):
    """
    Post find check to make sure no stop codons have slipped through
    """
    df.loc[:, "translated_seq"] = df.apply(lambda x: translate_sequence(x["ScFv_seq"]))
    df.loc[:, "contains_stop"] = df["ScFv_seq"].str.contains("\*")
    print(f"{df['contains_stop'].sum()} number of final sequences with stop codons")


def main():
    # folder_path = get_folder_path()
    folder_path = "/Users/benjamindraper/Library/CloudStorage/Dropbox/01_UCL/Projects/PB00020_PHD_RS/Alk_ScFv_Library/Alk5-8"
    fastq = [x for x in os.listdir(folder_path) if ".fastq" in x]
    variables = {
        "min_length": 700, # minimum length of sequences you want to keep
        "max_length": 800, # maximum length of sequences you want to keep
        "quality": 15, # minimum quality of sequences you want to keep - fullcircle give you a graph showing quality, wouldn't recommend going lower than 10
        "fprimer": "GTCCCTGGCTCCACTGGA", # Sequence of the forward primer used to amplify the PCR fragment
        "rprimer": "GCGCTGGCGTCGTGGT", # Sequence of the reverse primer used to amplify the PCR fragment
        "fastq": f"{folder_path}/{fastq[0]}", # Path to the fastq file
        "reference": f"{data_path}/reference.fasta", # Path to the reference fasta file
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
    filtered_df = process_blast_df(
        min_match_length=40,
        min_pct_match=95,
        query_start=300,
        folder_path=folder_path
      )
    print("Filtering completed\nSplitting ScFv chains...")
    merged_df = split_df(filtered_df, folder_path)
    print("Splitting completed\nPreparing IGBLAST...")
    igblast = f"{data_path}/ncbi-igblast-1.22.0"
    if not check_igblast(igblast, data_path):
        print("IGBLAST successfully downloaded")
    print("Running IgBLAST...")
    env["IGDATA"] = igblast
    # run_igblast(igblast, folder_path)
    print("IgBLAST Successfully run")
    gdf = process_igblast_results(folder_path)
    print("IgBLAST Results filtered for complete variable sequences")

    sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    linker_seq = merged_df["ref_seq"].unique()[0]

    results_df = find_scfv_sequences(gdf, sequence_parts, linker_seq)
    check_stops(results_df)
    results_df.to_csv(f"{folder_path}/results_df.csv", index=False)
    print(results_df.shape[0])
    print(results_df["ScFv_seq"].unique().shape[0], "ScFv sequences found")




if __name__ == "__main__":

        main()

        print("All done!")











        # This is just for installing igblast and downloading the databases







        # Load the previous dataframe
        # merged_df = pd.read_csv(f"{folder_path}/merged_df.csv")
        # Filter by the cleaned sequence ids






        # # This goes through each sequence id and identifies the heavy and light chains
        # # It then concatenates the sequences to make a full length ScFv sequence
        # # TODO - make sure the nt sequences start in reading frame 1 - need to check vh and vl separately
        # # TODO - this can probably be cleaned up a bit but should work as is
        # sequence_parts = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
        # heavy_chain_id = "IGH"
        # light_chain_id = "IGK"
        # for i, seqid in enumerate(sequence_ids):
        #     tdf = gdf.get_group(seqid)

        #     pre_seq = "".join(tdf.loc[tdf["sequence_id"] == f"{seqid}_pre", sequence_parts].values.flatten())
        #     pre_trans_frames = check_translation(pre_seq)
        #     if len(pre_trans_frames) == 1:
        #         pre_seq = pre_seq[pre_trans_frames[0]-1:]
        #     pre_seq_chain = tdf.loc[tdf["sequence_id"] == f"{seqid}_pre", "locus"].values[0]


        #     post_seq = "".join(tdf.loc[tdf["sequence_id"] == f"{seqid}_post", sequence_parts].values.flatten())
        #     post_trans_frames = check_translation(post_seq)
        #     if len(post_trans_frames) == 1:
        #         post_seq = post_seq[post_trans_frames[0]-1:]
        #     post_seq_chain = tdf.loc[tdf["sequence_id"] == f"{seqid}_post", "locus"].values[0]

        #     scfv_seq = linker_seq.join([pre_seq, post_seq])

        #     # print(scfv_seq)

        #     orientation = "-".join([pre_seq_chain, post_seq_chain])

        #     merged_df.loc[merged_df["query_id"] == seqid, ["pre_seq", "post_seq", "ScFv_seq", "orientation"]] = [pre_seq, post_seq, scfv_seq, orientation]

        #     if i % 20 == 0:
        #         print(f"{i}/{len(sequence_ids)} completed")
        # # merged_df.head()

        # merged_df.to_csv(f"{folder_path}/merged_df.csv", index=False)

        # print(merged_df.shape[0])
        # print(merged_df["ScFv_seq"].unique().shape[0], "ScFv sequences found")
