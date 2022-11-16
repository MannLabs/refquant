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
