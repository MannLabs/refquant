from ..reference_quantification.multi_run import multi_run_table_creation
from ..reference_quantification.multi_run import channel_bias_compensator
import directlfq.lfq_manager as lfq_manager
import pandas as pd


DIANN_FILE_NAME = None
PRECURSOR_REFQUANT_UNCORR = None
PRECURSOR_REFQUANT_CORR = None
PRECURSOR_MS1_UNCORR = None
PRECURSOR_MS1_CORR = None
PRECURSOR_DIANN_UNCORR = None
PRECURSOR_DIANN_CORR = None


def write_out_precursortable_in_multiple_variations(diann_file_name, single_labelled_precursors):
    
    init_file_locations(diann_file_name)
    
    table_creator_refquant = multi_run_table_creation.PrecursorQuantityTableCreator(single_labelled_precursors)
    table_creator_diann = multi_run_table_creation.PrecursorTableCreatorQuantityBySearchEngine(single_labelled_precursors)
    table_creator_ms1 = multi_run_table_creation.PrecursorTableCreatorQuantityByMS1(single_labelled_precursors)

    table_creator_refquant.precursorquantitytable.to_csv(PRECURSOR_REFQUANT_UNCORR, sep="\t", index = None)
    table_creator_ms1.precursorquantitytable.to_csv(PRECURSOR_MS1_UNCORR, sep="\t", index = None)
    table_creator_diann.precursorquantitytable.to_csv(PRECURSOR_DIANN_UNCORR, sep="\t", index = None)

    channel_bias_compensator.ChannelBiasCompensator(table_creator_refquant.precursorquantitytable).precursor_table_df.to_csv(PRECURSOR_REFQUANT_CORR, sep="\t", index = None)
    channel_bias_compensator.ChannelBiasCompensator(table_creator_ms1.precursorquantitytable).precursor_table_df.to_csv(PRECURSOR_MS1_CORR, sep="\t", index = None)
    channel_bias_compensator.ChannelBiasCompensator(table_creator_diann.precursorquantitytable).precursor_table_df.to_csv(PRECURSOR_DIANN_CORR, sep="\t", index = None)


def write_out_protein_tables():
    all_precursor_tables = [PRECURSOR_REFQUANT_UNCORR, PRECURSOR_REFQUANT_CORR, PRECURSOR_MS1_UNCORR, 
                            PRECURSOR_MS1_CORR, PRECURSOR_DIANN_UNCORR, PRECURSOR_DIANN_CORR]
    
    diann_additional_columns_df = pd.read_csv(DIANN_FILE_NAME, sep="\t", usecols = ["Protein.Group","Protein.Names", "Genes"]).drop_duplicates()
    
    for precursor_table in all_precursor_tables:
        lfq_manager.run_lfq(precursor_table, input_type_to_use = "directlfq")
        outfile_name = f"{precursor_table}.directlfq.protein_intensities.tsv"
        read_reformat_directlfq_output(outfile_name, diann_additional_columns_df)


    
def read_reformat_directlfq_output(directlfq_output, diann_additional_columns_df):
    df_directlfq = pd.read_csv(directlfq_output, sep="\t")
    lenght_before = len(df_directlfq.index)
    df_directlfq = df_directlfq.merge(diann_additional_columns_df, left_on="protein", right_on="Protein.Group", how="left")
    length_after = len(df_directlfq.index)
    
    assert lenght_before == length_after
    df_directlfq.to_csv(directlfq_output, sep="\t", index = None)

def init_file_locations(diann_file_name):
    global DIANN_FILE_NAME
    global PRECURSOR_REFQUANT_UNCORR
    global PRECURSOR_REFQUANT_CORR
    global PRECURSOR_MS1_UNCORR
    global PRECURSOR_MS1_CORR
    global PRECURSOR_DIANN_UNCORR
    global PRECURSOR_DIANN_CORR

    DIANN_FILE_NAME = diann_file_name

    PRECURSOR_REFQUANT_UNCORR = f"{diann_file_name}.precursors_uncorrected_refquant.tsv"
    PRECURSOR_REFQUANT_CORR = f"{diann_file_name}.precursors_corrected_refquant.tsv"
    
    PRECURSOR_MS1_UNCORR = f"{diann_file_name}.precursors_uncorrected_ms1.tsv"
    PRECURSOR_MS1_CORR = f"{diann_file_name}.precursors_corrected_ms1.tsv"

    PRECURSOR_DIANN_UNCORR = f"{diann_file_name}.precursors_uncorrected_diann.tsv"
    PRECURSOR_DIANN_CORR = f"{diann_file_name}.precursors_corrected_diann.tsv"


