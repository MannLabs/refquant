from ..reference_quantification.multi_run import multi_run_table_creation
from ..reference_quantification.multi_run import channel_bias_compensator

def write_out_precursortable_in_multiple_variations(diann_file_name, single_labelled_precursors):
    precursor_resfilename_uncorrected = f"{diann_file_name}.precursortable_uncorrected_refquant.tsv"
    precursor_resfilename_corrected = f"{diann_file_name}.precursortable_corrected_refquant.tsv"
    
    precursor_resfilename_ms1 = f"{diann_file_name}.precursortable_uncorrected_ms1.tsv"
    precursor_resfilename_ms1_corrected = f"{diann_file_name}.precursortable_corrected_ms1.tsv"

    precursor_resfilename_diann = f"{diann_file_name}.precursortable_uncorrected_diann.tsv"
    precursor_resfilename_diann_corrected = f"{diann_file_name}.precursortable_corrected_diann.tsv"
    
    table_creator_refquant = multi_run_table_creation.PrecursorQuantityTableCreator(single_labelled_precursors)
    table_creator_diann = multi_run_table_creation.PrecursorTableCreatorQuantityBySearchEngine(single_labelled_precursors)
    table_creator_ms1 = multi_run_table_creation.PrecursorTableCreatorQuantityByMS1(single_labelled_precursors)

    table_creator_refquant.precursorquantitytable.to_csv(precursor_resfilename_uncorrected, sep="\t", index = None)
    table_creator_ms1.precursorquantitytable.to_csv(precursor_resfilename_ms1, sep="\t", index = None)
    table_creator_diann.precursorquantitytable.to_csv(precursor_resfilename_diann, sep="\t", index = None)

    channel_bias_compensator.ChannelBiasCompensator(table_creator_refquant.precursorquantitytable).precursor_table_df.to_csv(precursor_resfilename_corrected, sep="\t", index = None)
    channel_bias_compensator.ChannelBiasCompensator(table_creator_ms1.precursorquantitytable).precursor_table_df.to_csv(precursor_resfilename_ms1_corrected, sep="\t", index = None)
    channel_bias_compensator.ChannelBiasCompensator(table_creator_diann.precursorquantitytable).precursor_table_df.to_csv(precursor_resfilename_diann_corrected, sep="\t", index = None)