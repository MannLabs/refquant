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





import pandas as pd
import os
class DIANNQvalueAdder():
    channel2name = {"(Dimethyl-n-0)": "reference", "(Dimethyl-n-4)": "target4", "(Dimethyl-n-8)": "target8"}
    def __init__(self, diann_input_table):
        self._diann_input_table = diann_input_table
        self._preprocessed_file = f"{self._diann_input_table}.qvalshortened.tsv"

        self._preprocessed_file_exists = None

        self._input_df = None
        
        self._precursor2run2channel2qvalue = {}

        self._check_if_preprocessed_file_exists()
        self._load_input_df()
        self._save_precprocessed_file_if_it_does_not_exist()
        self._define_precuror2run2qvalue()

    def add_qvalue_infos(self, single_labelled_precursor):
        key = (single_labelled_precursor.name, single_labelled_precursor.replicate_name, single_labelled_precursor.channel_name)
        is_in  = key in self._precursor2run2channel2qvalue.keys()
        if is_in:
            single_labelled_precursor.channel_qvalue = self._precursor2run2channel2qvalue[key]
    
    def _check_if_preprocessed_file_exists(self):
        self._preprocessed_file_exists = os.path.exists(self._preprocessed_file)

    def _load_input_df(self):
        if self._preprocessed_file_exists:
            self._input_df = pd.read_csv(self._preprocessed_file, sep="\t")
        else:
            self._input_df = pd.read_csv(self._diann_input_table, sep="\t", usecols=['Run','Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Channel.Q.Value']).drop_duplicates()
            self._add_precursor_column_to_input_df()
            self._add_channel_column_to_input_df()
    
    def _add_precursor_column_to_input_df(self):
        self._input_df['precursor'] = self._input_df.apply(self._generate_precursor_from_columnvals, axis=1)

    def _add_channel_column_to_input_df(self):
        self._input_df['channel'] = [self._get_channel_column(x) for x in self._input_df['Modified.Sequence']]
    
    def _get_channel_column(self, mod_seq):
        if "(Dimethyl-n-0)" in mod_seq:
            return self.channel2name["(Dimethyl-n-0)"]
        elif "(Dimethyl-n-4)" in mod_seq:
            return self.channel2name["(Dimethyl-n-4)"]
        elif "(Dimethyl-n-8)" in mod_seq:
            return self.channel2name["(Dimethyl-n-8)"]
        
    def _save_precprocessed_file_if_it_does_not_exist(self):
        self._input_df.to_csv(self._preprocessed_file, sep="\t", index=False)

    def _define_precuror2run2qvalue(self):
        qval_df = self._input_df[['precursor', 'Run', 'channel','Channel.Q.Value']]
        qval_df = qval_df.groupby(['precursor', 'Run', 'channel']).max()

        self._precursor2run2channel2qvalue = qval_df.to_dict()['Channel.Q.Value']


    def _generate_precursor_from_columnvals(self, row: dict):
        stripped_seq = row['Stripped.Sequence']
        mod_seq = row['Modified.Sequence']
        charge = row['Precursor.Charge']
        stripped_seq = self._remove_mtraq_modifications_from_ion_id(stripped_seq)
        mod_seq = self._remove_mtraq_modifications_from_ion_id(mod_seq)

        return f"SEQ_{stripped_seq}_MOD_{mod_seq}_CHARGE_{charge}_"
        

    @staticmethod
    def _remove_mtraq_modifications_from_ion_id(ion):
        all_mtraq_tags = ["(Dimethyl-K-0)", "(Dimethyl-K-4)", "(Dimethyl-K-8)", "(Dimethyl-n-0)", "(Dimethyl-n-4)", "(Dimethyl-n-8)"]
        for tag in all_mtraq_tags:
            ion = ion.replace(tag, "")
        return ion

class DIANNQvalueAdderRef8(DIANNQvalueAdder):
    channel2name = {"(Dimethyl-n-0)": "target0", "(Dimethyl-n-4)": "target4", "(Dimethyl-n-8)": "reference"}