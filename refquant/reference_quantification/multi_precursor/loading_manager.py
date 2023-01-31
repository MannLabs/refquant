import multiprocessing
import pandas as pd

from . import precursor_combiner
from . import precursor_loader


def get_all_single_labelled_precursors_in_dataset_diann(reference_table, use_multiprocessing=False, number_of_cores = None):
    run2df = get_run2df_dict(reference_table)
    print(f"processing {len(run2df.keys())} runs")
    if use_multiprocessing:
        return run_multiprocessing(run2df, get_single_labelled_precursors_diann,number_of_cores)
    else:
        return run_linear_processing(run2df, get_single_labelled_precursors_diann)


def get_all_single_labelled_precursors_in_dataset_prm(reference_table, use_multiprocessing=False, number_of_cores = None):
    run2df = get_run2df_dict(reference_table)
    print(f"processing {len(run2df.keys())} runs")
    if use_multiprocessing:
        return run_multiprocessing(run2df, get_single_labelled_precursors_prm,number_of_cores)
    else:
        return run_linear_processing(run2df, get_single_labelled_precursors_prm)


def run_multiprocessing(run2df,  get_single_labelled_precursors_function,number_of_cores = None):
    args = [(x, run2df)   for x in run2df.keys()]
    pool = get_configured_multiprocessing_pool(number_of_cores)
    single_labelled_precursors = pool.starmap( get_single_labelled_precursors_function, args)
    pool.close()
    #join list of lists
    single_labelled_precursors = [item for sublist in single_labelled_precursors for item in sublist]
    return single_labelled_precursors


def get_configured_multiprocessing_pool(num_cores):
    multiprocessing.freeze_support()
    if num_cores is None:
        num_cores = int(multiprocessing.cpu_count()) if int(multiprocessing.cpu_count()/2) < 60 else 60 #windows upper thread limit
    pool = multiprocessing.Pool(num_cores)
    print(f"using {pool._processes} processes")
    return pool


def run_linear_processing(run2df, get_single_labelled_precursors_function):
    single_labelled_precursors = []
    for run in run2df.keys():
        single_labelled_precursors.extend(get_single_labelled_precursors_function(run, run2df))
    return single_labelled_precursors

def get_run2df_dict(reference_table):
    run2df_dict = {}
    reference_df = pd.read_csv(reference_table, sep = "\t").set_index("run")
    runs = reference_df.index.unique()
    for run in runs:
        run2df_dict[run] = reference_df.loc[run].reset_index()
    reference_df = None
    return run2df_dict

def get_single_labelled_precursors_diann(run, run2df):
    print("run: ", run)
    single_labelled_precursors = []
    reference_df = run2df[run]
    precursor_loader_initialized = precursor_loader.PrecursorLoaderDIANNFromDf(reference_df)
    label_combiner = precursor_combiner.PrecursorCombinerToSingleReference(precursor_loader_initialized)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors


def get_single_labelled_precursors_spectronaut(run, run2df):
    print("run: ", run)
    single_labelled_precursors = []
    reference_df = run2df[run]
    precursor_loader_initialized = precursor_loader.PrecursorLoaderFromDf(reference_df)
    label_combiner = precursor_combiner.PrecursorCombinerToSingleReference(precursor_loader_initialized)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors


def get_single_labelled_precursors_prm(run, run2df):
    print("run: ", run)
    single_labelled_precursors = []
    reference_df = run2df[run]
    precursor_loader_initialized = precursor_loader.PrecursorLoaderPRMFromDf(reference_df)
    label_combiner = precursor_combiner.PrecursorCombinerToPRMReference(precursor_loader_initialized)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors

