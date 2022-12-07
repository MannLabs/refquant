import refquant.table_preparation.table_import as table_import
import refquant.reference_quantification.multi_precursor.loading_manager as loading_manager
import refquant.quant_table_writeout.writeout_manager as writeout_manager

def run_refquant(diann_file_qvalfiltered, use_multiprocessing=True):
    reformatter = table_import.TableReformatterDIANN(diann_file_qvalfiltered, quantitative_extraction_types=["diann_fragion_isotopes_mDIA_raw", "diann_precursors_mDIA"])
    refquant_reformatted_table = reformatter.filename_reformatted
    reference_quantified_precursors = loading_manager.get_all_single_labelled_precursors_in_dataset_diann(reference_table = refquant_reformatted_table, use_multiprocessing=use_multiprocessing)
    
    writeout_manager.write_out_precursortable_in_multiple_variations(diann_file_qvalfiltered, reference_quantified_precursors)

