import numpy as np
import os
import pandas as pd

def get_runs(reference_table):
    runs_file = f"{reference_table}.runs.txt"
    if os.path.exists(runs_file):
        df_replicates = pd.read_csv(runs_file, sep="\t")['0']
    else:
        df_replicates = pd.Series(pd.read_csv(reference_table, sep = "\t")["run"].unique())
        df_replicates.to_csv(runs_file, sep = "\t", index = False)
    return list(df_replicates)


def write_shortened_diann_file_w_channel_lib_pg_cutoff(diann_file, qval_cutoff):
    diann_df = pd.read_csv(diann_file, sep="\t")
    diann_df = diann_df[diann_df['Channel.Q.Value'] < qval_cutoff]
    diann_df = diann_df[diann_df['Lib.PG.Q.Value'] < 0.01]
    diann_df.to_csv(diann_file + f".filtered_lib_pg_ch_qval{qval_cutoff}.tsv", sep="\t", index=False)


def format_precursor_file_to_iq(precursor_file):
    df_input = pd.read_csv(precursor_file, sep = "\t")
    df_input = df_input.melt(id_vars=["protein", "ion"], value_name="quant", var_name="experiment")
    df_input = df_input[df_input["quant"]>1]
    outfile = f"{precursor_file}.formatted_for_iq.tsv"
    df_input.to_csv(outfile, sep="\t", index = None)

def format_list_of_precursor_files_to_iq(list_of_precursor_files):
    for precursor_file in list_of_precursor_files:
        format_precursor_file_to_iq(precursor_file)

def rescale_iq_files_from_log_to_linear(list_of_iq_files):
    iq_files = [f"{input}.formatted_for_iq.tsv.maxlfq_iq_protein_intensities.tsv" for input in list_of_iq_files]
    for processed_input in iq_files:
        df_input = pd.read_csv(processed_input, sep = "\t")
        df_input = (2**df_input.set_index("protein")).reset_index()
        df_input = df_input.replace(np.nan, 0)
        outfile = f"{processed_input}.linearized.tsv"
        df_input.to_csv(outfile, sep="\t", index = None)


def write_out_filtered_precursor_tables(list_of_precursor_files):
    for precursor_file in list_of_precursor_files:
        precursor_table = pd.read_csv(precursor_file, sep='\t')
        filtered_precursor_table = _filter_precursor_table_by_num_valid_precursors(precursor_table, 1000)
        filtered_precursor_table.to_csv(precursor_file + ".filtered.tsv", sep='\t', index=False)


def _filter_precursor_table_by_num_valid_precursors(precursor_table, min_num_valid_precursors):
    #get numeric columns
    numeric_columns = precursor_table.select_dtypes(include=['float64'])
    #count num non-null float values in each column
    num_valid_precursors = numeric_columns.count()
    #filter columns by num valid precursors
    columns_to_drop = num_valid_precursors[num_valid_precursors < min_num_valid_precursors].index
    filtered_precursor_table = precursor_table.drop(columns_to_drop, axis=1)
    return filtered_precursor_table