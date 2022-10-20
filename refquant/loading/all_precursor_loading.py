
import multiprocessing
import refquant.refquant_classes as refquant_classes
import refquant.refquant_utils as utils
import pandas as pd

def get_all_single_labelled_precursors_in_dataset_diann(reference_table, diann_file, use_multiprocessing=False, number_of_cores = None):
    diann_qvaladder = utils.DIANNQvalueAdder(diann_file)
    run2df = get_run2df_dict(reference_table)
    print(f"processing {len(run2df.keys())} runs")
    if use_multiprocessing:
        return run_multiprocessing_diann(run2df, diann_qvaladder, number_of_cores)
    else:
        return run_linear_processing_diann(run2df, diann_qvaladder)


def run_multiprocessing_diann(run2df, diann_qvaladder, number_of_cores = None):
    multiprocessing.freeze_support()
    args = [(x, run2df, diann_qvaladder)   for x in run2df.keys()]
    if number_of_cores is None:
        number_of_cores = int(multiprocessing.cpu_count()/2)
    print(f"{number_of_cores} cores of {multiprocessing.cpu_count()} used")

    with multiprocessing.Pool(int(multiprocessing.cpu_count()/2)) as pool:
        single_labelled_precursors = pool.starmap(get_single_labelled_precursors_diann, args)
    #join list of lists
    single_labelled_precursors = [item for sublist in single_labelled_precursors for item in sublist]
    return single_labelled_precursors

def run_linear_processing_diann(run2df, diann_qvaladder):
    single_labelled_precursors = []
    for run in run2df.keys():
        single_labelled_precursors.extend(get_single_labelled_precursors_diann(run, run2df, diann_qvaladder))
    return single_labelled_precursors

def get_run2df_dict(reference_table):
    run2df_dict = {}
    reference_df = pd.read_csv(reference_table, sep = "\t").set_index("run")
    runs = reference_df.index.unique()
    for run in runs:
        run2df_dict[run] = reference_df.loc[run].reset_index()
    reference_df = None
    return run2df_dict

def get_single_labelled_precursors_diann(run, run2df, diann_qvaladder):
    print("run: ", run)
    single_labelled_precursors = []
    reference_df = run2df[run]
    precursor_loader = refquant_classes.PrecursorLoaderDIANNFromDf(reference_df, diann_qvaladder)
    label_combiner = refquant_classes.PrecursorCombinerToSingleReference(precursor_loader)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors


def get_single_labelled_precursors_spectronaut(run, run2df):
    print("run: ", run)
    single_labelled_precursors = []
    reference_df = run2df[run]
    precursor_loader = refquant_classes.PrecursorLoaderFromDf(reference_df)
    label_combiner = refquant_classes.PrecursorCombinerToSingleReference(precursor_loader)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors

