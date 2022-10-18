import pandas as pd
import os
import re
import table_utils
from abc import ABC, abstractmethod

class TableReformatter(ABC):
    def __init__(self, input_file : str, quantitative_extraction_types : list):
        self._input_file = input_file
        self._quantitative_extraction_types = quantitative_extraction_types
        
        self._tmpdir_refquant = None
        self._per_channel_subfiles = None
        self._merged_aq_reformatted_df = None


        self.spectronaut_table_reformatted = None
        self.outfile_name = None

        self._create_refquant_tmp_dir()
    
    def _create_refquant_tmp_dir(self):
        #get the directory of the input file
        input_file_dir = os.path.dirname(os.path.abspath(self._input_file))
        print("input file")
        print(input_file_dir)
        print(self._input_file)
        self._tmpdir_refquant = f"{input_file_dir}/tmpdir_refquant"
        if not os.path.exists(self._tmpdir_refquant):
            os.makedirs(self._tmpdir_refquant)


class TableReformatterSpectronaut(TableReformatter):
    def __init__(self, input_file : str, quantitative_extraction_types : list = ["spectronaut_precursor_v3", "spectronaut_fragion_isotopes"]):
        self._input_file = input_file
        self._quantitative_extraction_types = quantitative_extraction_types
        
        self._tmpdir_refquant = None
        self._per_channel_subfiles = None
        self._merged_aq_reformatted_df = None


        self.spectronaut_table_reformatted = None
        self.outfile_name = None

        super()._create_refquant_tmp_dir()
        self._define_per_channel_files()
        self._define_merged_aq_reformatted_df()
        self._define_spectronaut_table_reformatted()
        self._save_reformatted_table_and_define_outfile_name()

    def _create_refquant_tmp_dir(self):
        #get the directory of the input file
        input_file_dir = os.path.dirname(os.path.abspath(self._input_file))
        self._tmpdir_refquant = f"{input_file_dir}/tmpdir_refquant"
        if not os.path.exists(self._tmpdir_refquant):
            os.makedirs(self._tmpdir_refquant)

    def _define_per_channel_files(self):
        self._per_channel_subfiles = SpectronautPerChannelSubfileWriter(self._input_file, self._tmpdir_refquant).per_channel_subfiles #writes out the channel files and stores them in a list

    def _define_merged_aq_reformatted_df(self):
        self._merged_aq_reformatted_df = TableReformatterAndMerger(self._per_channel_subfiles, self._quantitative_extraction_types).merged_aq_reformatted_df

    def _define_spectronaut_table_reformatted(self):
        self.spectronaut_table_reformatted = AQTableReformatterSpectronaut(self._merged_aq_reformatted_df).reformatted_df
    
    def _save_reformatted_table_and_define_outfile_name(self):
        self.outfile_name = f"{self._input_file}.reformatted_for_refquant.tsv"
        self.spectronaut_table_reformatted.to_csv(self.outfile_name, sep="\t", index=False)


import table_utils
class TableReformatterDIANN(TableReformatter):
    def __init__(self, input_file : str, quantitative_extraction_types : list = ["diann_fragion_isotopes_mDIA_raw", "diann_precursors_mDIA"]):
        self._input_file = input_file
        self._quantitative_extraction_types = quantitative_extraction_types
        
        self._tmpdir_refquant = None
        self._per_channel_subfiles = None
        self._aq_reformatted_df = None


        self.diann_table_reformatted = None
        self.outfile_name = None

        super()._create_refquant_tmp_dir()
        self._define_merged_aq_reformatted_df()
        self._define_diann_table_reformatted()
        self._save_reformatted_table_and_define_outfile_name()

    def _define_merged_aq_reformatted_df(self):
        reformatted_dfs = []
        for quantitative_extraction_type in self._quantitative_extraction_types:
            aq_reformatted_df = table_utils.import_data(self._input_file, results_dir=self._tmpdir_refquant, input_type_to_use=quantitative_extraction_type)
            reformatted_dfs.append(aq_reformatted_df)
        self._aq_reformatted_df = pd.concat(reformatted_dfs, ignore_index=True)
        self._aq_reformatted_df = DIANNfragionIDAdder(diann_original_input_file=self._input_file, aq_reformatted_df=self._aq_reformatted_df).aq_reformatted_df

    def _define_diann_table_reformatted(self):
        self.diann_table_reformatted = AQTableReformatterDIANN(self._aq_reformatted_df).reformatted_df
    
    def _save_reformatted_table_and_define_outfile_name(self):
        self.outfile_name = f"{self._input_file}.reformatted_for_refquant.tsv"
        self.diann_table_reformatted.to_csv(self.outfile_name, sep="\t", index=False)


class TableReformatterDIANNCh8Ref(TableReformatterDIANN):
    def _define_diann_table_reformatted(self):
        self.diann_table_reformatted = AQTableReformatterDIANNCh8Ref(self._aq_reformatted_df).reformatted_df

    


class SpectronautPerChannelSubfileWriter():
    channels = [0, 4, 8]

    def __init__(self, input_file, tmpdir_refquant):
        self._input_file = input_file
        self._tmpdir_refquant = tmpdir_refquant

        self._input_filename = None
        self._spectronaut_df = None

        self.per_channel_subfiles = []

        self._define_input_filename()
        self._load_spectronaut_df()
        self._write_subfiles_per_channel()

    def _define_input_filename(self):
        self._input_filename = os.path.basename(self._input_file)

    def _load_spectronaut_df(self):
        self._spectronaut_df = pd.read_csv(self._input_file, sep = "\t")

    def _write_subfiles_per_channel(self):
        for channel in self.channels:
            channel_file = f'{self._tmpdir_refquant}/{self._input_filename}.channel{channel}.tsv.zip'
            self._subset_spectronaut_df_to_channel_and_save(channel, channel_file)
            self.per_channel_subfiles.append(channel_file)

    def _subset_spectronaut_df_to_channel_and_save(self, channel_number, channel_file):
        is_correct_channel_number = [self._get_channel_number_from_fg_id(fg_id) == channel_number for fg_id in self._spectronaut_df['FG.Id']]
        df_spectronaut_channel = self._spectronaut_df[is_correct_channel_number]
        df_spectronaut_channel.to_csv(channel_file, sep='\t',compression='zip')

    @staticmethod
    def _get_channel_number_from_fg_id(fg_id):
        #search for regex pattern in string
        match = re.search('(_\[DimethNter)(.)(\].*)', fg_id)
        return int(match.group(2))

    def get_table_file(folder_containing_table):
        for file in os.listdir(folder_containing_table):
            if file.endswith("Report.tsv"):
                return os.path.join(folder_containing_table, file)


from functools import reduce
class TableReformatterAndMerger():
    def __init__(self, per_channel_subfiles :list, quantitative_extraction_types : list):
        self._per_channel_subfiles : list = per_channel_subfiles

        self._quantitative_extraction_types : list = quantitative_extraction_types

        self._channel_table_infos : list(ChannelTableInfo) = []

        self.merged_aq_reformatted_df = None

        self._define_channel_table_infos()
        self._define_merged_df()

    def _define_channel_table_infos(self):
        for per_channel_subfile in self._per_channel_subfiles:
            channel_table_info = ChannelTableInfoCollector(per_channel_subfile, self._quantitative_extraction_types).channel_table_info
            self._channel_table_infos.append(channel_table_info)
    
    def _define_merged_df(self):
        self.merged_aq_reformatted_df = pd.concat([channel_table_info.merged_reformatted_df for channel_table_info in self._channel_table_infos], axis=1, sort=False)
        self.merged_aq_reformatted_df = self.merged_aq_reformatted_df.reset_index()


class ChannelTableInfoCollector():
    def __init__(self,  channel_subfile : str,quantitative_extraction_types : list):
        self._channel_subfile : str = channel_subfile
        self._quantitative_extraction_types : list = quantitative_extraction_types

        self.channel_table_info : ChannelTableInfo = None

        self._define_channel_table_info()

    def _define_channel_table_info(self):
        reformatted_tables = self._get_reformatted_tables(self._channel_subfile)
        channel_name = self._get_channel_name_from_subfile(self._channel_subfile)
        self.channel_table_info = ChannelTableInfo(reformatted_dfs_based_on_quantinfos= reformatted_tables, channel_name= channel_name)

          
    def _get_reformatted_tables(self, subfile):
        reformatted_tables = []
        for quantitative_extraction_type in self._quantitative_extraction_types:
            reformatted_table = self._reformat_table_to_specified_format(subfile, quantitative_extraction_type)
            reformatted_tables.append(reformatted_table)
        return reformatted_tables

    def _get_channel_name_from_subfile(self, subfile):
        return re.search('(.*)\.(channel.)\.tsv\.zip', subfile).group(2)

    def _reformat_table_to_specified_format(self, subfile, format_to_use_for_reformatting):
        reformatted_table = table_utils.import_data(subfile, input_type_to_use=format_to_use_for_reformatting)
        reformatted_table = reformatted_table.set_index(["protein", "ion"])
        return reformatted_table


class ChannelTableInfo():
    def __init__(self, reformatted_dfs_based_on_quantinfos : list, channel_name : str):
        
        self._reformatted_dfs_based_on_quantinfos : list = reformatted_dfs_based_on_quantinfos
        self.channel_name :str = channel_name

        self.merged_reformatted_df = reformatted_dfs_based_on_quantinfos

        self._define_merged_reformatted_df()
        self._add_channel_suffixes()

    def _define_merged_reformatted_df(self):
        self.merged_reformatted_df = pd.concat(self._reformatted_dfs_based_on_quantinfos)
    
    def _add_channel_suffixes(self):
        self.merged_reformatted_df.columns = [col + '_' + self.channel_name for col in self.merged_reformatted_df.columns]

import numpy as np
class AQTableReformatter(ABC):
    ordered_headers =  ["run", "protein", "precursor","ion", "reference", "target4", "target8"]

    def __init__(self, aq_format_df : pd.DataFrame):
        self._aq_format_df : pd.DataFrame = aq_format_df

        self.reformatted_df : pd.DataFrame = None

        self._filter_and_annotate_aq_format_df()
        self._define_reformatted_df()

    def _filter_and_annotate_aq_format_df(self):
        self._filter_aq_df_to_meaningful_ions()
        self._add_precursor_column_to_aq_format_df()
        self._change_zero_to_nan_in_aq_format_df()
        self._rename_columns_in_aq_format_df()
    
    def _filter_aq_df_to_meaningful_ions(self):
        self._aq_format_df  = self._aq_format_df[[self._check_if_ion_is_meaningful(ion) for ion in  self._aq_format_df['ion']]]
    
    def _add_precursor_column_to_aq_format_df(self):
        self._aq_format_df["precursor"] = [self._parse_precursor_from_ion(ion) for ion in self._aq_format_df["ion"]]

    def _change_zero_to_nan_in_aq_format_df(self):
        self._aq_format_df = self._aq_format_df.replace(0, np.nan)

    def _rename_columns_in_aq_format_df(self):
        renamed_columns = []
        for column in self._aq_format_df.columns:
            column_renamed = column
            for channel, name in self.channel2name.items():
                column_renamed = column_renamed.replace(channel, name)
            renamed_columns.append(column_renamed)
        self._aq_format_df.columns = renamed_columns

    def _define_reformatted_df(self):
        self.reformatted_df = self._melt_aq_reformatted_df()
        self._annotate_channel_to_reformatted_df()
        self._annotate_run_to_reformatted_df()
        self._unstack_channels()
        self._reduce_to_ordered_headers()
        self._move_channels_into_same_row()
        self._drop_ions_w_no_data_in_target_channels()
        
    
    def _melt_aq_reformatted_df(self):
        return pd.melt(self._aq_format_df, id_vars=["protein", "ion", "precursor"], var_name="run_channel", value_name="intensity").dropna().drop_duplicates(subset=["ion", "run_channel"])

    def _annotate_channel_to_reformatted_df(self):
        self.reformatted_df["channel"] = [self._get_channel_name(x) for x in self.reformatted_df["run_channel"]]

    def _annotate_run_to_reformatted_df(self):
        self.reformatted_df["run"] = [self._get_run_name(x) for x in self.reformatted_df["run_channel"]]

    def _drop_run_channel_column(self):
        self.reformatted_df =  self.reformatted_df.drop(columns=["run_channel"])

    def _unstack_channels(self):
        self.reformatted_df = pd.pivot(self.reformatted_df, index=["protein", "ion", "precursor", "run"], columns="channel", values="intensity").reset_index()
        
    def _reduce_to_ordered_headers(self):
        self.reformatted_df = self.reformatted_df[self.ordered_headers]
        self.reformatted_df = self.reformatted_df.drop_duplicates()

    def _move_channels_into_same_row(self):
        self.reformatted_df =self.reformatted_df.groupby(['run', 'protein', 'precursor', 'ion']).sum().reset_index()

    def _drop_ions_w_no_data_in_target_channels(self):
        target_channels = [x for x in self.ordered_headers if "target" in x]
        if len(target_channels) >2:
            raise ValueError("More than 2 target channels")
        at_least_one_not_na = np.logical_or(self.reformatted_df[target_channels[0]].replace(0, np.nan).notna().values, self.reformatted_df[target_channels[1]].replace(0, np.nan).notna().values)
        self.reformatted_df = self.reformatted_df[at_least_one_not_na]

    @staticmethod
    def _parse_precursor_from_ion(ion):
        if "FRGION" in ion:
            return ion.split("FRGION")[0]
        elif "MS1ISOTOPES" in ion:
            return ion.split("MS1ISOTOPES")[0]
        else:
            return ion

    @staticmethod
    def _check_if_ion_is_meaningful(ion):
        if "FRGION" not in ion and "MS1ISOTOPES" not in ion:
            return True
        if "MS1ISOTOPES" in ion:
            return True
        elif "K[DimethLys0]_" in ion: #if the end is labelled, all fragions are meaningful
            return True
        elif "FRGION_b" in ion:
            return True
        else:
            return False
        
    @staticmethod
    def _get_channel_name(run_channel):
        channelname =  run_channel.split("_")[-1]
        if (not "target" in channelname) and (not "reference" in channelname):
            raise ValueError("Channel name not recognized")
        return channelname

    @staticmethod
    def _get_run_name(run_channel):
        channelname =  run_channel.split("_")[-1]
        run_name = run_channel.replace(f"_{channelname}", "")
        return run_name


class AQTableReformatterSpectronaut(AQTableReformatter):
    channels = ["channel0", "channel4", "channel8"]
    channel_names = [ "reference", "target4", "target8"]
    channel2name = dict(zip(channels, channel_names))


class AQTableReformatterDIANN(AQTableReformatter):
    channels = ["(Dimethyl-n-0)", "(Dimethyl-n-4)", "(Dimethyl-n-8)"]
    channel_names = [ "reference", "target4", "target8"]
    channel2name = dict(zip(channels, channel_names))

    def __init__(self, aq_format_df : pd.DataFrame):
        self._aq_format_df : pd.DataFrame = aq_format_df

        self.reformatted_df : pd.DataFrame = None

        self._remove_decoy_columns_if_present()
        super()._filter_and_annotate_aq_format_df()
        super()._define_reformatted_df()
    
    def _remove_decoy_columns_if_present(self):
        decoy_columns = [x for x in self._aq_format_df.columns if x.endswith("_NA")]
        self._aq_format_df = self._aq_format_df.drop(columns=decoy_columns)

class AQTableReformatterDIANNCh8Ref(AQTableReformatterDIANN):
    ordered_headers =  ["run", "protein", "precursor","ion", "reference", "target0", "target4"]
    channels = ["(Dimethyl-n-0)", "(Dimethyl-n-4)", "(Dimethyl-n-8)"]
    channel_names = [ "target0", "target4", "reference"]
    channel2name = dict(zip(channels, channel_names))


class DIANNfragionIDAdder():
    def __init__(self, diann_original_input_file, aq_reformatted_df):
        self._diann_original_input_file = diann_original_input_file
        self.aq_reformatted_df = aq_reformatted_df

        self._diann_df_relevant_cols = None
        self._ionid2fragionid = {}

        self._load_diann_df_relevant_cols()
        self._define_ionid2fragion()
        self._add_fragion_id_to_aq_reformatted_df()
    
    def _load_diann_df_relevant_cols(self):
        self._diann_df_relevant_cols = pd.read_csv(self._diann_original_input_file, sep="\t", usecols=['Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Fragment.Info']).drop_duplicates()

    def _define_ionid2fragion(self):
        #iterate over rows and generate idstrings
        for idx, row in self._diann_df_relevant_cols.iterrows():
            self._ionid2fragionid.update(self._generate_idstrings_from_row(row))

    def _add_fragion_id_to_aq_reformatted_df(self):
        self.aq_reformatted_df["ion"] = [self._ionid2fragionid.get(x, x) for x in self.aq_reformatted_df["ion"]]


    def _generate_idstrings_from_row(self, row):
        stripped_seq = row['Stripped.Sequence']
        mod_seq = row['Modified.Sequence']
        charge = row['Precursor.Charge']
        ion = row['Fragment.Info']
        stripped_seq = self._remove_mtraq_modifications_from_ion_id(stripped_seq)
        mod_seq = self._remove_mtraq_modifications_from_ion_id(mod_seq)

        idstring_stump = f"SEQ_{stripped_seq}_MOD_{mod_seq}_CHARGE_{charge}_FRGION_"
        idstring_w_number2idstring_w_fragion = self._go_through_fragment_infos_and_generate_id_mappings(ion, idstring_stump)
        return idstring_w_number2idstring_w_fragion

    @staticmethod
    def _remove_mtraq_modifications_from_ion_id(ion):
        all_mtraq_tags = ["(Dimethyl-K-0)", "(Dimethyl-K-4)", "(Dimethyl-K-8)", "(Dimethyl-n-0)", "(Dimethyl-n-4)", "(Dimethyl-n-8)"]
        for tag in all_mtraq_tags:
            ion = ion.replace(tag, "")
        return ion

    @staticmethod
    def _go_through_fragment_infos_and_generate_id_mappings(fragment_info_string, idstring_stump):
        idstring_w_number2idstring_w_fragion = {}
        fragment_info_string_split = fragment_info_string.split(";")
        for idx in range(len(fragment_info_string_split)):
            fragment_info = fragment_info_string_split[idx]
            fragion = fragment_info.split("^")[0]
            idstring_w_number = f"{idstring_stump}{idx}_"
            idstring_w_fragion = f"{idstring_stump}{fragion}"
            idstring_w_number2idstring_w_fragion[idstring_w_number] = idstring_w_fragion
        return idstring_w_number2idstring_w_fragion

