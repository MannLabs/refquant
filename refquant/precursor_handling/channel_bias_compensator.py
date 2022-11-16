import numpy as np

class ChannelBiasCompensator():
    def __init__(self, precursor_table_df):
        self.precursor_table_df = precursor_table_df.copy()
        
        self._channel4_columns = None
        self._channel8_columns = None
        self._median_ratios_compensator_channel = None

        self._adapt_precursor_table_df()
        self._define_channel_columns()
        self._define_median_ratios_between_each_channel()
        self._compensate_channel()
        self._rescale_precursor_table_df()
    
    def _adapt_precursor_table_df(self):
        self.precursor_table_df = self.precursor_table_df.replace(0, np.nan)
        self.precursor_table_df = np.log2(self.precursor_table_df.set_index(["protein", "ion"])).reset_index()


    def _define_channel_columns(self):
        self._channel4_columns = self.precursor_table_df.columns[self.precursor_table_df.columns.str.contains("target4")]
        self._channel8_columns = self.precursor_table_df.columns[self.precursor_table_df.columns.str.contains("target8")]
    
    def _define_median_ratios_between_each_channel(self):
        median_channel4 = self.precursor_table_df[self._channel4_columns].median(axis = 1).to_numpy()
        median_channel8 = self.precursor_table_df[self._channel8_columns].median(axis = 1).to_numpy()
        self._median_ratios_compensator_channel = median_channel8 - median_channel4


    def _compensate_channel(self):
        self.precursor_table_df[self._channel4_columns] = self.precursor_table_df[self._channel4_columns].to_numpy() + self._median_ratios_compensator_channel.reshape(-1, 1)

    def _rescale_precursor_table_df(self):
        self.precursor_table_df = (2**self.precursor_table_df.set_index(["protein", "ion"])).reset_index()
        self.precursor_table_df = self.precursor_table_df.fillna(0)



class MildChannelBiasCompensator(ChannelBiasCompensator):
    def __init__(self, precursor_table_df):
        super().__init__(precursor_table_df)

    def _compensate_channel(self):
        self.precursor_table_df[self._channel4_columns] = self.precursor_table_df[self._channel4_columns].to_numpy() + np.nanmedian(self._median_ratios_compensator_channel)
