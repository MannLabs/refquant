import pandas as pd
import os
import pathlib


if "__file__" in globals():#only run in the translated python file, as __file__ is not defined with ipython
    INTABLE_CONFIG = os.path.join(pathlib.Path(__file__).parent.absolute(),"..", "configs", "intable_config.yaml") #the yaml config is located one directory below the python library files


# Cell
def get_condpairname(condpair):
    return f"{condpair[0]}_VS_{condpair[1]}"

# Cell

def get_quality_score_column(acquisition_info_df):
    if "FG.ShapeQualityScore" in acquisition_info_df.columns:
        param = "FG.ShapeQualityScore"
    elif "Quantity.Quality" in acquisition_info_df.columns:
        param = "Quantity.Quality"
    return param

# Cell
import os

def make_dir_w_existcheck(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

# Cell
import os
def get_results_plot_dir_condpair(results_dir, condpair):
    results_dir_plots = f"{results_dir}/{condpair}_plots"
    make_dir_w_existcheck(results_dir_plots)
    return results_dir_plots

# Cell
def get_middle_elem(sorted_list):
    nvals = len(sorted_list)
    if nvals==1:
        return sorted_list[0]
    middle_idx = nvals//2
    if nvals%2==1:
        return sorted_list[middle_idx]
    return 0.5* (sorted_list[middle_idx] + sorted_list[middle_idx-1])

# Cell
import numpy as np
def get_nonna_array(array_w_nas):
    res = []
    isnan_arr = np.isnan(array_w_nas)

    for idx in range(len(array_w_nas)):
        sub_res = []
        sub_array = array_w_nas[idx]
        na_array = isnan_arr[idx]
        for idx2 in range(len(sub_array)):
            if not na_array[idx2]:
               sub_res.append(sub_array[idx2])
        res.append(np.array(sub_res))
    return np.array(res)

# Cell
import numpy as np
def get_non_nas_from_pd_df(df):
    return {
        pep_name: sub_vals[~np.isnan(sub_vals)] for pep_name, sub_vals in
        zip( df.index.values, df.values)
    }

# Cell
import numpy as np
def get_ionints_from_pd_df(df):
    return {
        pep_name: sub_vals for pep_name, sub_vals in
        zip( df.index.values, df.values)
    }

# Cell
def invert_dictionary(my_map):
    inv_map = {}
    for k, v in my_map.items():
        inv_map[v] = inv_map.get(v, []) + [k]
    return inv_map

# Cell
import statistics

def get_z_from_p_empirical(p_emp,p2z):
    p_rounded = np.format_float_scientific(p_emp, 1)
    if p_rounded in p2z:
        return p2z.get(p_rounded)
    z = statistics.NormalDist().inv_cdf(float(p_rounded))
    p2z[p_rounded] = z
    return z

# Cell
def count_fraction_outliers_from_expected_fc(result_df, threshold, expected_log2fc):
    num_outliers = sum([abs(x-expected_log2fc)> threshold for x in result_df["log2fc"]])
    fraction_outliers = num_outliers/len(result_df["log2fc"])
    print(f"{round(fraction_outliers, 2)} outliers")
    return fraction_outliers

# Cell
import os
import shutil
def create_or_replace_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)

# Cell
def write_chunk_to_file(chunk, filepath ,write_header):
    """write chunk of pandas dataframe to a file"""
    chunk.to_csv(filepath, header=write_header, mode='a', sep = "\t", index = None)

# Cell
def index_and_log_transform_input_df(data_df):
    data_df = data_df.set_index(["protein", "ion"])
    return np.log2(data_df.replace(0, np.nan))

# Cell
def remove_allnan_rows_input_df(data_df):
    return data_df.dropna(axis = 0, how = 'all')




# Cell
import yaml
import itertools

def get_relevant_columns(protein_cols, ion_cols, sample_ID, quant_ID, filter_dict):
    filtcols = []
    for filtconf in filter_dict.values():
        filtcols.append(filtconf.get('param'))
    relevant_cols = protein_cols + ion_cols + [sample_ID] + [quant_ID] + filtcols
    relevant_cols = list(set(relevant_cols)) # to remove possible redudancies
    return relevant_cols


def get_relevant_columns_config_dict(config_typedict):
    filtcols = []
    dict_ioncols = []
    for filtconf in config_typedict.get('filters', {}).values():
        filtcols.append(filtconf.get('param'))

    if 'ion_hierarchy' in config_typedict.keys():
        for headr in config_typedict.get('ion_hierarchy').values():
            ioncols = list(itertools.chain.from_iterable(headr.get("mapping").values()))
            dict_ioncols.extend(ioncols)

    quant_ids = get_quant_ids_from_config_dict(config_typedict)
    sample_ids = get_sample_ids_from_config_dict(config_typedict)
    channel_ids = get_channel_ids_from_config_dict(config_typedict)
    relevant_cols = config_typedict.get("protein_cols") + config_typedict.get("ion_cols", []) + sample_ids + quant_ids + filtcols + dict_ioncols + channel_ids
    relevant_cols = list(set(relevant_cols)) # to remove possible redudancies
    return relevant_cols

def get_quant_ids_from_config_dict(config_typedict):
    quantID = config_typedict.get("quant_ID")
    if type(quantID) ==type("string"):
        return [config_typedict.get("quant_ID")]
    if quantID == None:
        return[]
    else:
        return list(config_typedict.get("quant_ID").values())

def get_sample_ids_from_config_dict(config_typedict):
    sampleID = config_typedict.get("sample_ID")
    if type(sampleID) ==type("string"):
        return [config_typedict.get("sample_ID")]
    if sampleID == None:
        return []
    else:
        return config_typedict.get("sample_ID")

def get_channel_ids_from_config_dict(config_typedict):
    return config_typedict.get("channel_ID", [])



def load_config(config_yaml):
    stream = open(config_yaml, 'r')
    config_all = yaml.safe_load(stream)
    return config_all

def get_type2relevant_cols(config_all):
    type2relcols = {}
    for type in config_all.keys():
        config_typedict = config_all.get(type)
        relevant_cols = get_relevant_columns_config_dict(config_typedict)
        type2relcols[type] = relevant_cols
    return type2relcols

# Cell

def filter_input(filter_dict, input):
    if filter_dict == None:
        return input
    for filtname,filterconf in filter_dict.items():
        param = filterconf.get('param')
        comparator = filterconf.get('comparator')
        value = filterconf.get('value')

        if comparator not in [">",">=", "<", "<=", "==", "!="]:
            raise TypeError(f"cannot identify the filter comparator of {filtname} given in the longtable config yaml!")

        if comparator=="==":
            input = input[input[param] ==value]
            continue
        try:
            input = input.astype({f"{param}" : "float"})
        except:
            pass

        if comparator==">":
            input = input[input[param].astype(type(value)) >value]

        if comparator==">=":
            input = input[input[param].astype(type(value)) >=value]

        if comparator=="<":
            input = input[input[param].astype(type(value)) <value]

        if comparator=="<=":
            input = input[input[param].astype(type(value)) <=value]

        if comparator=="!=":
            input = input[input[param].astype(type(value)) !=value]

    return input

# Cell
def merge_protein_and_ion_cols(input_df, config_dict):
    protein_cols =  config_dict.get("protein_cols")
    ion_cols = config_dict.get("ion_cols")
    input_df['protein'] = input_df.apply(lambda row : "_".join(row[protein_cols].astype('string')), axis = 1)
    input_df['ion'] = input_df.apply(lambda row : "_".join(row[ion_cols].astype('string')), axis = 1)
    input_df = input_df.rename(columns = {config_dict.get('quant_ID') : "quant_val"})
    return input_df

# Cell
import copy
def merge_protein_cols_and_ion_dict(input_df, config_dict):
    """[summary]

    Args:
        input_df ([pandas dataframe]): longtable containing peptide intensity data
        confid_dict ([dict[String[]]]): nested dict containing the parse information. derived from yaml file

    Returns:
        pandas dataframe: longtable with newly assigned "protein" and "ion" columns
    """
    protein_cols = config_dict.get("protein_cols")
    ion_hierarchy = config_dict.get("ion_hierarchy")
    splitcol2sep = config_dict.get('split_cols')
    quant_id_dict = config_dict.get('quant_ID')

    ion_dfs = []
    input_df['protein'] = input_df.apply(lambda row : "_".join(row[protein_cols].astype('string')), axis = 1)

    input_df = input_df.drop(columns = [x for x in protein_cols if x!='protein'])
    for hierarchy_type in ion_hierarchy.keys():
        df_subset = input_df.copy()
        ion_hierarchy_local = ion_hierarchy.get(hierarchy_type).get("order")
        ion_headers_merged, ion_headers_grouped = get_ionname_columns(ion_hierarchy.get(hierarchy_type).get("mapping"), ion_hierarchy_local) #ion headers merged is just a helper to select all relevant rows, ionheaders grouped contains the sets of ionstrings to be merged into a list eg [[SEQ, MOD], [CH]]
        quant_columns = get_quantitative_columns(df_subset, hierarchy_type, config_dict, ion_headers_merged)
        headers = list(set(ion_headers_merged + quant_columns + ['protein']))
        if "sample_ID" in config_dict.keys():
            headers+=[config_dict.get("sample_ID")]
        df_subset = df_subset[headers].drop_duplicates()

        if splitcol2sep is not None:
            if quant_columns[0] in splitcol2sep.keys(): #in the case that quantitative values are stored grouped in one column (e.g. msiso1,msiso2,msiso3, etc.), reformat accordingly
                df_subset = split_extend_df(df_subset, splitcol2sep)
            ion_headers_grouped = adapt_headers_on_extended_df(ion_headers_grouped, splitcol2sep)

        #df_subset = df_subset.set_index(quant_columns)

        df_subset = add_merged_ionnames(df_subset, ion_hierarchy_local, ion_headers_grouped, quant_id_dict, hierarchy_type)
        ion_dfs.append(df_subset)
    input_df = pd.concat(ion_dfs, ignore_index=True)
    return input_df


def get_quantitative_columns(input_df, hierarchy_type, config_dict, ion_headers_merged):
    naming_columns = ion_headers_merged + ['protein']
    if config_dict.get("format") == 'longtable':
        quantcol = config_dict.get("quant_ID").get(hierarchy_type)
        return [quantcol]

    if config_dict.get("format") == 'widetable':
        quantcolumn_candidates = [x for x in input_df.columns if x not in naming_columns]
        if "quant_prefix" in config_dict.keys():
            return [x for x in quantcolumn_candidates if x.startswith(config_dict.get("quant_prefix"))] # in the case that the quantitative columns have a prefix (like "Intensity " in MQ peptides.txt), only columns with the prefix are filtered
        else:
            return quantcolumn_candidates #in this case, we assume that all non-ionname/proteinname columns are quantitative columns


def get_ionname_columns(ion_dict, ion_hierarchy_local):
    ion_headers_merged = []
    ion_headers_grouped = []
    for lvl in ion_hierarchy_local:
        vals = ion_dict.get(lvl)
        ion_headers_merged.extend(vals)
        ion_headers_grouped.append(vals)
    return ion_headers_merged, ion_headers_grouped


def adapt_headers_on_extended_df(ion_headers_grouped, splitcol2sep):
    #in the case that one column has been split, we need to designate the "naming" column
    ion_headers_grouped_copy = copy.deepcopy(ion_headers_grouped)
    for vals in ion_headers_grouped_copy:
        if splitcol2sep is not None:
            for idx in range(len(vals)):
                if vals[idx] in splitcol2sep.keys():
                    vals[idx] = vals[idx] + "_idxs"
    return ion_headers_grouped_copy

def split_extend_df(input_df, splitcol2sep, value_threshold=10):
    """reformats data that is stored in a condensed way in a single column. For example isotope1_intensity;isotope2_intensity etc. in Spectronaut

    Args:
        input_df ([type]): [description]
        splitcol2sep ([type]): [description]
        value_threshold([type]): [description]

    Returns:
        Pandas Dataframe: Pandas dataframe with the condensed items expanded to long format
    """
    if splitcol2sep==None:
        return input_df

    for split_col, separator in splitcol2sep.items():
        idx_name = f"{split_col}_idxs"
        split_col_series = input_df[split_col].str.split(separator)
        input_df = input_df.drop(columns = [split_col])

        input_df[idx_name] = [list(range(len(x))) for x in split_col_series]
        exploded_input = input_df.explode(idx_name)
        exploded_split_col_series = split_col_series.explode()

        exploded_input[split_col] = exploded_split_col_series.replace('', 0) #the column with the intensities has to come after to column with the idxs

        exploded_input = exploded_input.astype({split_col: float})
        exploded_input = exploded_input[exploded_input[split_col]>value_threshold]
        #exploded_input = exploded_input.rename(columns = {'var1': split_col})
    return exploded_input



def add_merged_ionnames(df_subset, ion_hierarchy_local, ion_headers_grouped, quant_id_dict, hierarchy_type):
    """puts together the hierarchical ion names as a column in a given input dataframe"""
    all_ion_headers = list(itertools.chain.from_iterable(ion_headers_grouped))
    columns_to_index = [x for x in df_subset.columns if x not in all_ion_headers]
    df_subset = df_subset.set_index(columns_to_index)

    rows = df_subset[all_ion_headers].to_numpy()
    ions = []

    for row in rows: #iterate through dataframe
        count = 0
        ionstring = ""
        for lvl_idx in range(len(ion_hierarchy_local)):
            ionstring += f"{ion_hierarchy_local[lvl_idx]}"
            for sublvl in ion_headers_grouped[lvl_idx]:
                ionstring+= f"_{row[count]}_"
                count+=1
        ions.append(ionstring)
    df_subset['ion'] = ions
    df_subset = df_subset.reset_index()
    if quant_id_dict!= None:
        df_subset = df_subset.rename(columns = {quant_id_dict.get(hierarchy_type) : "quant_val"})
    return df_subset

# Cell
import os.path
def reformat_and_write_longtable_according_to_config(input_file, outfile_name, config_dict_for_type, sep = "\t",decimal = ".", enforce_largefile_processing = False, chunksize =1000_000):
    """Reshape a long format proteomics results table (e.g. Spectronaut or DIA-NN) to a wide format table.
    :param file input_file: long format proteomic results table
    :param string input_type: the configuration key stored in the config file (e.g. "diann_precursor")
    """
    filesize = os.path.getsize(input_file)/(1024**3) #size in gigabyte
    file_is_large = (filesize>20 and str(input_file).endswith(".zip")) or filesize>50 or enforce_largefile_processing

    if file_is_large:
        tmpfile_large = f"{input_file}.tmp.longformat.columnfilt.tsv" #only needed when file is large
        #remove potential leftovers from previous processings
        if os.path.exists(tmpfile_large):
            os.remove(tmpfile_large)
        if os.path.exists(outfile_name):
            os.remove(outfile_name)

    relevant_cols = get_relevant_columns_config_dict(config_dict_for_type)
    input_df_it = pd.read_csv(input_file, sep = sep, decimal=decimal, usecols = relevant_cols, encoding ='latin1', chunksize = chunksize)
    input_df_list = []
    header = True
    for input_df_subset in input_df_it:
        input_df_subset = adapt_subtable(input_df_subset, config_dict_for_type)
        if file_is_large:
            write_chunk_to_file(input_df_subset,tmpfile_large, header)
        else:
            input_df_list.append(input_df_subset)
        header = False

    if file_is_large:
        process_with_dask(tmpfile_columnfilt=tmpfile_large , outfile_name = outfile_name, config_dict_for_type=config_dict_for_type)
    else:
        input_df = pd.concat(input_df_list)
        input_reshaped = reshape_input_df(input_df, config_dict_for_type)
        input_reshaped.to_csv(outfile_name, sep = "\t", index = None)


def adapt_subtable(input_df_subset, config_dict):
    input_df_subset = filter_input(config_dict.get("filters", {}), input_df_subset)
    if "ion_hierarchy" in config_dict.keys():
        return merge_protein_cols_and_ion_dict(input_df_subset, config_dict)
    else:
        return merge_protein_and_ion_cols(input_df_subset, config_dict)


# Cell
import dask.dataframe as dd
import pandas as pd
import glob
import os
import shutil

def process_with_dask(*, tmpfile_columnfilt, outfile_name, config_dict_for_type):
    df = dd.read_csv(tmpfile_columnfilt, sep = "\t")
    allcols = df[config_dict_for_type.get("sample_ID")].drop_duplicates().compute() # the columns of the output table are the sample IDs
    allcols = extend_sample_allcolumns_for_mDIA_case(allcols_samples=allcols, config_dict_for_type=config_dict_for_type)
    allcols = ['protein', 'ion'] + sorted(allcols)
    df = df.set_index('protein')
    sorted_filedir = f"{tmpfile_columnfilt}_sorted"
    df.to_csv(sorted_filedir, sep = "\t")
    #now the files are sorted and can be pivoted chunkwise (multiindex pivoting at the moment not possible in dask)
    files_dask = glob.glob(f"{sorted_filedir}/*part")
    header = True
    for file in files_dask:
        input_df = pd.read_csv(file, sep = "\t")
        if len(input_df.index) <2:
            continue
        input_reshaped = reshape_input_df(input_df, config_dict_for_type)
        input_reshaped = sort_and_add_columns(input_reshaped, allcols)
        write_chunk_to_file(input_reshaped, outfile_name, header)
        header = False
    os.remove(tmpfile_columnfilt)
    shutil.rmtree(sorted_filedir)

def reshape_input_df(input_df, config_dict):
    input_df = input_df.astype({'quant_val': 'float'})
    input_df = adapt_input_df_columns_in_case_of_mDIA(input_df=input_df, config_dict_for_type=config_dict)
    input_reshaped = pd.pivot_table(input_df, index = ['protein', 'ion'], columns = config_dict.get("sample_ID"), values = 'quant_val', fill_value=0)

    input_reshaped = input_reshaped.reset_index()
    return input_reshaped


def sort_and_add_columns(input_reshaped, allcols):
    missing_cols = set(allcols) - set(input_reshaped.columns)
    input_reshaped[list(missing_cols)] = 0
    input_reshaped = input_reshaped[allcols]
    return input_reshaped


def extend_sample_allcolumns_for_mDIA_case(allcols_samples, config_dict_for_type):
    if is_mDIA_table(config_dict_for_type):
        new_allcols = []
        channels = ['Dimethyl-n-0', 'Dimethyl-n-4', 'Dimethyl-n-8']
        for channel in channels:
            for sample in allcols_samples:
                new_allcols.append(merge_channel_and_sample_string(sample, channel))
        return new_allcols
    else:
        return allcols_samples


# Cell
#mDIA case

def adapt_input_df_columns_in_case_of_mDIA(input_df,config_dict_for_type):
    if is_mDIA_table(config_dict_for_type):
        input_df = extend_sampleID_column_for_mDIA_case(input_df, config_dict_for_type)
        input_df = set_mtraq_reduced_ion_column_into_dataframe(input_df)
        return input_df
    else:
        return input_df


def extend_sampleID_column_for_mDIA_case(input_df,config_dict_for_type):
    channels_per_peptide = parse_channel_from_peptide_column(input_df)
    return merge_sample_id_and_channels(input_df, channels_per_peptide, config_dict_for_type)


def set_mtraq_reduced_ion_column_into_dataframe(input_df):
    new_ions = remove_mtraq_modifications_from_ion_ids(input_df['ion'])
    input_df['ion'] = new_ions
    return input_df

def remove_mtraq_modifications_from_ion_ids(ions):
    new_ions = []
    for ion in ions:
        new_ions.append( re.sub("\(Dimethyl-\w-\d\)","", ion))
    return new_ions


def is_mDIA_table(config_dict_for_type):
    return config_dict_for_type.get('channel_ID') == ['Channel.0', 'Channel.4']


import re
def parse_channel_from_peptide_column(input_df):
    channels = []
    for pep in input_df['Modified.Sequence']:
        pattern = "(.*)(\(Dimethyl-n-.\))(.*)"
        matched = re.match(pattern, pep)
        num_appearances = pep.count("Dimethyl-n-")
        if matched and num_appearances==1:
            channels.append(matched.group(2))
        else:
            channels.append("NA")
    return channels

def merge_sample_id_and_channels(input_df, channels, config_dict_for_type):
    sample_id = config_dict_for_type.get("sample_ID")
    sample_ids = list(input_df[sample_id])
    input_df[sample_id] = [merge_channel_and_sample_string(sample_ids[idx], channels[idx]) for idx in range(len(sample_ids))]
    return input_df

def merge_channel_and_sample_string(sample, channel):
    return f"{sample}_{channel}"


# Cell
def reformat_and_write_wideformat_table(peptides_tsv, outfile_name, config_dict):
    input_df = pd.read_csv(peptides_tsv,sep="\t", encoding ='latin1')
    filter_dict = config_dict.get("filters")
    protein_cols = config_dict.get("protein_cols")
    ion_cols = config_dict.get("ion_cols")
    input_df = filter_input(filter_dict, input_df)
    #input_df = merge_protein_and_ion_cols(input_df, config_dict)
    input_df = merge_protein_cols_and_ion_dict(input_df, config_dict)
    if 'quant_prefix' in config_dict.keys():
        quant_prefix = config_dict.get('quant_prefix')
        headers = ['protein', 'ion'] + list(filter(lambda x: x.startswith(quant_prefix), input_df.columns))
        input_df = input_df[headers]
        input_df = input_df.rename(columns = lambda x : x.replace(quant_prefix, ""))

    input_df = input_df.reset_index()

    input_df.to_csv(outfile_name, sep = '\t', index = None)

# Cell
import os
def check_for_processed_runs_in_results_folder(results_folder):
    contained_condpairs = []
    folder_files = os.listdir(results_folder)
    result_files = list(filter(lambda x: "results.tsv" in x ,folder_files))
    for result_file in result_files:
        res_name = result_file.replace(".results.tsv", "")
        if ((f"{res_name}.normed.tsv" in folder_files) & (f"{res_name}.results.ions.tsv" in folder_files)):
            contained_condpairs.append(res_name)
    return contained_condpairs



def import_data(input_file, input_type_to_use = None, samples_subset = None, results_dir = None):
    """
    Function to import peptide level data. Depending on available columns in the provided file,
    the function identifies the type of input used (e.g. Spectronaut, MaxQuant, DIA-NN), reformats if necessary
    and returns a generic wide-format dataframe
    :param file input_file: quantified peptide/ion -level data
    :param file results_folder: the folder where the directlfq outputs are stored
    """

    samples_subset = add_ion_protein_headers_if_applicable(samples_subset)
    if "aq_reformat" in input_file:
        file_to_read = input_file
    else:
        file_to_read = reformat_and_save_input_file(input_file=input_file, input_type_to_use=input_type_to_use)

    input_reshaped = pd.read_csv(file_to_read, sep = "\t", encoding = 'latin1', usecols=samples_subset)
    input_reshaped = input_reshaped.drop_duplicates(subset='ion')
    return input_reshaped


def reformat_and_save_input_file(input_file, input_type_to_use = None):

    input_type, config_dict_for_type, sep = get_input_type_and_config_dict(input_file, input_type_to_use)
    print(f"using input type {input_type}")
    format = config_dict_for_type.get('format')
    outfile_name = f"{input_file}.{input_type}.aq_reformat.tsv"

    if format == "longtable":
        reformat_and_write_longtable_according_to_config(input_file, outfile_name,config_dict_for_type, sep = sep)
    elif format == "widetable":
        reformat_and_write_wideformat_table(input_file, outfile_name, config_dict_for_type)
    else:
        raise Exception('Format not recognized!')
    return outfile_name




def add_ion_protein_headers_if_applicable(samples_subset):
    if samples_subset is not None:
        return samples_subset + ["ion", "protein"]
    else:
        return None






# Cell
import pandas as pd
import os.path
import pathlib

def get_input_type_and_config_dict(input_file, input_type_to_use = None):
    #parse the type of input (e.g. Spectronaut Fragion+MS1Iso) out of the input file


    config_dict = load_config(INTABLE_CONFIG)
    type2relevant_columns = get_type2relevant_cols(config_dict)

    if "aq_reformat.tsv" in input_file:
        input_file = get_original_file_from_aq_reformat(input_file)

    filename = str(input_file)
    if '.csv' in filename:
        sep=','
    if '.tsv' in filename:
        sep='\t'
    if '.txt' in filename:
        sep='\t'

    if 'sep' not in locals():
        raise TypeError(f"neither of the file extensions (.tsv, .csv, .txt) detected for file {input_file}! Your filename has to contain one of these extensions. Please modify your file name accordingly.")



    uploaded_data_columns = set(pd.read_csv(input_file, sep=sep, nrows=1, encoding ='latin1').columns)

    for input_type in type2relevant_columns.keys():
        if (input_type_to_use is not None) and (input_type!=input_type_to_use):
            continue
        relevant_columns = type2relevant_columns.get(input_type)
        relevant_columns = [x for x in relevant_columns if x] #filter None values
        if set(relevant_columns).issubset(uploaded_data_columns):
            config_dict_type =  config_dict.get(input_type)
            return input_type, config_dict_type, sep
    raise TypeError("format not specified in intable_config.yaml!")

import re
def get_original_file_from_aq_reformat(input_file):
    matched = re.match("(.*)(\..*\.)(aq_reformat\.tsv)",input_file)
    return matched.group(1)

# Cell
def import_config_dict():
    config_dict = load_config(INTABLE_CONFIG)
    return config_dict

# Cell

import pandas as pd

def load_samplemap(samplemap_file):
    file_ext = os.path.splitext(samplemap_file)[-1]
    if file_ext=='.csv':
        sep=','
    if (file_ext=='.tsv') | (file_ext=='.txt'):
        sep='\t'

    if 'sep' not in locals():
        print(f"neither of the file extensions (.tsv, .csv, .txt) detected for file {samplemap_file}! Trying with tab separation. In the case that it fails, please add the appropriate extension to your file name.")
        sep = "\t"

    return pd.read_csv(samplemap_file, sep = sep, encoding ='latin1', dtype='str')

# Cell
def prepare_loaded_tables(data_df, samplemap_df):
    """
    Integrates information from the peptide/ion data and the samplemap, selects the relevant columns and log2 transforms intensities.
    """
    samplemap_df = samplemap_df[samplemap_df["condition"]!=""] #remove rows that have no condition entry
    filtvec_not_in_data = [(x in data_df.columns) for x in samplemap_df["sample"]] #remove samples that are not in the dataframe
    samplemap_df = samplemap_df[filtvec_not_in_data]
    headers = ['protein'] + samplemap_df["sample"].to_list()
    data_df = data_df.set_index("ion")
    for sample in samplemap_df["sample"]:
        data_df[sample] = np.log2(data_df[sample].replace(0, np.nan))
    return data_df[headers], samplemap_df

# Cell


#export
class LongTableReformater():
    """Generic class to reformat tabular files in chunks. For the specific cases you can inherit the class and specify reformat and iterate function
    """
    def __init__(self, input_file):
        self._input_file = input_file
        self._reformatting_function = None
        self._iterator_function = self.__initialize_df_iterator__
        self._concat_list = []

    def reformat_and_load_acquisition_data_frame(self):

        input_df_it = self._iterator_function()

        input_df_list = []
        for input_df_subset in input_df_it:
            input_df_subset = self._reformatting_function(input_df_subset)
            input_df_list.append(input_df_subset)
        input_df = pd.concat(input_df_list)

        return input_df

    def reformat_and_save_acquisition_data_frame(self, output_file):

        input_df_it = self._iterator_function()
        write_header = True

        for input_df_subset in input_df_it:
            input_df_subset = self._reformatting_function(input_df_subset)
            self.__write_reformatted_df_to_file__(input_df_subset, output_file, write_header)
            write_header = False

    def __initialize_df_iterator__(self):
        return pd.read_csv(self._input_file, sep = "\t", encoding ='latin1', chunksize=1000000)

    @staticmethod
    def __write_reformatted_df_to_file__(reformatted_df, filepath ,write_header):
        reformatted_df.to_csv(filepath, header=write_header, mode='a', sep = "\t", index = None)


# Cell

import os
import re

class AcquisitionTableHandler():
    def __init__(self, results_dir, samples):
        self._table_infos = AcquisitionTableInfo(results_dir=results_dir)
        self._header_infos = AcquisitionTableHeaders(self._table_infos)
        self._samples = self.__reformat_samples_if_necessary(samples)

    def get_acquisition_info_df(self):
        return self.__get_reformated_df__()

    def save_dataframe_as_new_acquisition_dataframe(self):
        self._output_paths = AcquisitionTableOutputPaths(self._table_infos)
        self.__remove_possible_pre_existing_ml_table__(self._output_paths.output_file_name)
        df_reformater = AcquisitionTableReformater(table_infos = self._table_infos, header_infos=self._header_infos, samples = self._samples, dataframe_already_preformated=False)
        df_reformater.reformat_and_save_acquisition_data_frame(self._output_paths.output_file_name)

    def update_ml_file_location_in_method_parameters_yaml(self):
        method_params = load_method_parameters(self._table_infos._results_dir)
        if self._output_paths == None:
            raise Exception("output paths not initialized! This could be because no dataframe was saved before")
        method_params[self._output_paths.ml_file_accession_in_yaml] = self._output_paths.output_file_name
        save_dict_as_yaml(method_params, self._output_paths.method_parameters_yaml_path)

    def __get_reformated_df__(self):
        df_reformater = AcquisitionTableReformater(table_infos = self._table_infos, header_infos=self._header_infos, samples = self._samples, dataframe_already_preformated=True)
        df = df_reformater.reformat_and_load_acquisition_data_frame()
        return df.convert_dtypes()

    def __reformat_samples_if_necessary(self, samples):
        if "mDIA" in  self._table_infos._input_type:
            return self.__get_mDIA_samplenames__(samples)
        else:
            return samples

    def __get_mDIA_samplenames__(self, samples):
        new_samples = []
        for sample in samples:
            new_samples.append(self.__get_samplename_without_mtraq_tag__(sample))
        return new_samples

    @staticmethod
    def __get_samplename_without_mtraq_tag__(samplename):
        pattern = "(.*)(_\(Dimethyl-n-.\))"
        matched = re.match(pattern, samplename)
        return matched.group(1)

    @staticmethod
    def __remove_possible_pre_existing_ml_table__(output_file_name):
        if os.path.exists(output_file_name):
            os.remove(output_file_name)
            print(f"removed pre existing {output_file_name}")


class AcquisitionTableInfo():
    def __init__(self, results_dir, sep = "\t", decimal = "."):
        self._results_dir = results_dir
        self._sep = sep
        self._decimal = decimal
        self._method_params_dict = load_method_parameters(results_dir)
        self._input_file = self.__get_input_file__()
        self._file_ending_of_formatted_table = ".ml_info_table.tsv"
        self.already_formatted =  self.__check_if_input_file_is_already_formatted__()
        self._input_type, self._config_dict = self.__get_input_type_and_config_dict__()
        self._sample_column = self.__get_sample_column__()
        self.last_ion_level_to_use = self.__get_last_ion_level_to_use__()

    def __get_input_file__(self):
        if self._method_params_dict.get('ml_input_file') is None:
            return self.__get_location_of_original_file__()
        else:
            return self._method_params_dict.get('ml_input_file')

    def __check_if_input_file_is_already_formatted__(self):
        if self._file_ending_of_formatted_table in self._input_file:
            return True
        else:
            return False

    def __get_input_type_and_config_dict__(self):
        if self.already_formatted:
            original_file = self.__get_location_of_original_file__()
        else:
            original_file = self._input_file
        input_type, config_dict, _ = get_input_type_and_config_dict(original_file)
        return input_type, config_dict

    def __get_location_of_original_file__(self):
        input_file = self._method_params_dict.get('input_file')
        return self.__get_original_filename_from_input_file__(input_file)

    @staticmethod
    def __get_original_filename_from_input_file__(input_file):
        pattern = "(.*\.tsv|.*\.csv|.*\.txt)(\..*)(.aq_reformat.tsv)"
        m = re.match(pattern=pattern, string=input_file)
        if m:
            return m.group(1)
        else:
            return input_file


    def __get_sample_column__(self):
        return self._config_dict.get("sample_ID")

    def __get_last_ion_level_to_use__(self):
        return self._config_dict["ml_level"]





class AcquisitionTableHeaders():
    def __init__(self, acquisition_table_info):

        self._table_info = acquisition_table_info

        self._ion_hierarchy = self.__get_ordered_ion_hierarchy__()
        self._included_levelnames = self.__get_included_levelnames__()
        self._ion_headers_grouped = self.__get_ion_headers_grouped__()
        self._ion_headers = self.__get_ion_headers__()
        self._numeric_headers = self.__get_numeric_headers__()
        self._relevant_headers = self.__get_relevant_headers__()

    def __get_ordered_ion_hierarchy__(self):
        ion_hierarchy = self._table_info._config_dict.get("ion_hierarchy")
        hier_key = 'fragion' if 'fragion' in ion_hierarchy.keys() else list(ion_hierarchy.keys())[0]
        ion_hierarchy_on_chosen_key = ion_hierarchy.get(hier_key)
        return ion_hierarchy_on_chosen_key

    def __get_included_levelnames__(self):
        levelnames = self.__get_all_levelnames__(self._ion_hierarchy)
        last_ionlevel_idx = levelnames.index(self._table_info.last_ion_level_to_use)
        return levelnames[:last_ionlevel_idx+1]

    @staticmethod
    def __get_all_levelnames__(ion_hierarchy):
        return  ion_hierarchy.get('order')

    def __get_ion_headers_grouped__(self):
        mapping_dict = self.__get_levelname_mapping_dict(self._ion_hierarchy)
        return [mapping_dict.get(x) for x in self._included_levelnames]#on each level there can be multiple names, so it is a list of lists

    @staticmethod
    def __get_levelname_mapping_dict(ion_hierarchy):
        return ion_hierarchy.get('mapping')

    def __get_ion_headers__(self):
        return list(itertools.chain(*self._ion_headers_grouped))


    def __get_relevant_headers__(self):
        relevant_headers = self._numeric_headers+self._ion_headers + [self._table_info._sample_column]
        return self.__remove_possible_none_values_from_list__(relevant_headers)

    @staticmethod
    def __remove_possible_none_values_from_list__(list):
        return [x for x in list if x is not None]

    def __get_numeric_headers__(self):
        df_sample = pd.read_csv(self._table_info._input_file, sep = self._table_info._sep, decimal = self._table_info._decimal, encoding='latin1', nrows=3000) #sample 3000 rows from the df to assess the types of each row
        df_sample = df_sample.replace({False: 0, True: 1})
        numeric_headers =  list(df_sample.select_dtypes(include=np.number).columns)
        numeric_headers = AcquisitionTableHeaderFilter().filter_numeric_headers_if_specified(input_type = self._table_info._input_type, numeric_headers = numeric_headers)
        return numeric_headers


class AcquisitionTableOutputPaths():
    def __init__(self, table_info):
        self._table_info = table_info
        self.output_file_name = self.__get_output_file_name__()
        self.method_parameters_yaml_path = self.__get_method_parameters_yaml_path__()
        self.ml_file_accession_in_yaml = "ml_input_file"

    def __get_output_file_name__(self):
        old_file_name = self._table_info._input_file
        new_file_name = old_file_name+self._table_info._file_ending_of_formatted_table
        return new_file_name

    def __get_method_parameters_yaml_path__(self):
        return f"{self._table_info._results_dir}/aq_parameters.yaml"


class AcquisitionTableReformater(LongTableReformater):
    def __init__(self, table_infos, header_infos, samples, dataframe_already_preformated = False):

        LongTableReformater.__init__(self, table_infos._input_file)
        self._table_infos = table_infos
        self._header_infos = header_infos
        self._samples = samples
        self._dataframe_already_preformated = dataframe_already_preformated

        #set the two functions that specify the explicit reformatting
        self._reformatting_function = self.__reformatting_function__
        self._iterator_function = self.__initialize_iterator_with_specified_columns__

    def __reformatting_function__(self, input_df_subset):
        input_df_subset = input_df_subset.drop_duplicates()
        input_df_subset = self.__filter_reformated_df_if_necessary__(input_df_subset)
        if not self._dataframe_already_preformated:
            input_df_subset = add_merged_ionnames(input_df_subset, self._header_infos._included_levelnames, self._header_infos._ion_headers_grouped, None, None)
        return input_df_subset

    def __filter_reformated_df_if_necessary__(self, reformatted_df):
        if 'spectronaut' in self._table_infos._input_type or 'diann' in self._table_infos._input_type:
            return self.__filter_reformatted_dataframe_to_relevant_samples__(reformatted_df)
        else:
            return reformatted_df

    def __filter_reformatted_dataframe_to_relevant_samples__(self, input_df_subset):
        return input_df_subset[[x in self._samples for x in input_df_subset[self._table_infos._sample_column]]]

    def __initialize_iterator_with_specified_columns__(self):
        cols_to_use = self.__get_cols_to_use__()
        return pd.read_csv(self._table_infos._input_file, sep = self._table_infos._sep, decimal=self._table_infos._decimal, usecols = cols_to_use, encoding ='latin1', chunksize=1000000)

    def __get_cols_to_use__(self):
        cols_to_use = self._header_infos._relevant_headers
        if self._dataframe_already_preformated:
            return cols_to_use+['ion']
        else:
            return cols_to_use




class AcquisitionTableHeaderFilter():
    def __init__(self):
        self._spectronaut_header_filter = lambda x : (("EG." in x) | ("FG." in x)) and ("Global" not in x)
        self._maxquant_header_filter = lambda x : ("Intensity" not in x) and ("Experiment" not in x)

    def filter_numeric_headers_if_specified(self, input_type, numeric_headers):
        if 'spectronaut' in input_type:
            return [x for x in numeric_headers if self._spectronaut_header_filter(x)]
        elif 'maxquant' in input_type:
            return [x for x in numeric_headers if self._maxquant_header_filter(x)]
        else:
            return numeric_headers





# Cell

def merge_acquisition_df_parameter_df(acquisition_df, parameter_df, groupby_merge_type = 'mean'):
    """acquisition df contains details on the acquisition, parameter df are the parameters derived from the tree
    """
    merged_df = parameter_df.merge(acquisition_df, how = 'left', on = 'ion')
    if groupby_merge_type == 'mean':
        merged_df = merged_df.groupby('ion').mean().reset_index()
    if groupby_merge_type == 'min':
        merged_df = merged_df.groupby('ion').min().reset_index()
    if groupby_merge_type == 'max':
        merged_df = merged_df.groupby('ion').max().reset_index()
    merged_df = merged_df.dropna(axis=1, how='all')
    return merged_df

# Cell

import matplotlib.pyplot as plt
from scipy import stats
import itertools

def plot_withincond_fcs(normed_intensity_df, cut_extremes = True):
    """takes a normalized intensity dataframe and plots the fold change distribution between all samples. Column = sample, row = ion"""

    samplecombs = list(itertools.combinations(normed_intensity_df.columns, 2))

    for spair in samplecombs:#compare all pairs of samples
        s1 = spair[0]
        s2 = spair[1]
        diff_fcs = normed_intensity_df[s1].to_numpy() - normed_intensity_df[s2].to_numpy() #calculate fold changes by subtracting log2 intensities of both samples

        if cut_extremes:
            cutoff = max(abs(np.nanquantile(diff_fcs,0.025)), abs(np.nanquantile(diff_fcs, 0.975))) #determine 2.5% - 97.5% interval, i.e. remove extremes
            range = (-cutoff, cutoff)
        else:
            range = None
        plt.hist(diff_fcs,80,density=True, histtype='step',range=range) #set the cutoffs to focus the visualization
        plt.xlabel("log2 peptide fcs")

    plt.show()