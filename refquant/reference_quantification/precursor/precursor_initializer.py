import numpy as np
import pandas as pd
from . import precursor_classes



class PrecursorsFromDataframeInititalizer():
    channel_names = ["reference", "target4", "target8"]
    def __init__(self, precursor_name, precursor_df):
        self._precursor_name = precursor_name
        self._precursor_df = precursor_df

        self.single_labelled_precursors = []

        self._initialize_single_labelled_precursors()

    def _initialize_single_labelled_precursors(self):
        for channel_name in self.channel_names:
            single_labelled_precursor = SingleLabelledPrecursorIntitializer(self._precursor_name, self._precursor_df, channel_name).single_labelled_precursor
            self.single_labelled_precursors.append(single_labelled_precursor)

class PrecursorsFromDataframeInititalizerSpikeInRuns(PrecursorsFromDataframeInititalizer):
    channel_names = ["reference4", "reference8", "target4", "target8"]

class PrecursorsFromDataframeInititalizerReference8(PrecursorsFromDataframeInititalizer):
    channel_names = ["target0", "target4", "reference"]

class PrecursorsFromDataframeInitializerPRM(PrecursorsFromDataframeInititalizer):
    channel_names = ["target", "reference"]

class SingleLabelledPrecursorIntitializer():
    def __init__(self, precursor_name, precursor_df, channel_name):
        self._precursor_name = precursor_name
        self._precursor_df = precursor_df
        self._channel_name = channel_name

        self.single_labelled_precursor = precursor_classes.SingleLabelledPrecursor()

        self._initialize_single_labelled_precursor()
    
    def _initialize_single_labelled_precursor(self):
        self.single_labelled_precursor.channel_name = self._channel_name
        self.single_labelled_precursor.name = self._precursor_name
        self.single_labelled_precursor.protein_name = self._get_protein_name()
        self.single_labelled_precursor.replicate_name = self._get_replicate_name()
        self.single_labelled_precursor.fragion2quantity = self._get_fragment_ion_quantities_for_channel_name()
        self.single_labelled_precursor.search_engine_derived_quantity = self._get_quantity_estimated_by_diann()

    def _get_protein_name(self):
        return self._precursor_df["protein"].values[0]

    def _get_replicate_name(self):
        return self._precursor_df["run"].values[0]

    def _get_fragment_ion_quantities_for_channel_name(self):
        fragion2quantity = dict(zip(self._precursor_df["ion"], self._precursor_df[self._channel_name]))
        fragion2quantity = {key:value for key, value in fragion2quantity.items() if self._is_valid_dict_entry(key, value)}
        return fragion2quantity
    
    def _is_valid_dict_entry(self, key, value):
        return not np.isnan(value) #and ("FRGION" in key or "MS1ISO" in key)

    def _get_quantity_estimated_by_diann(self):
        ion_vals = self._precursor_df["ion"].values
        quant_vals = self._precursor_df[self._channel_name].values
        #get index of first true value in numpy array
        index_of_precursor = np.where(ion_vals == self._precursor_name)[0]
        if len(index_of_precursor) == 0:
            self.single_labelled_precursor.search_engine_derived_quantity_not_provided = True
            return np.nan# self._estimate_quantity_from_available_ions()
        return quant_vals[index_of_precursor[0]]

        

    def _estimate_quantity_from_available_ions(self):
        return self._precursor_df[self._channel_name].max()-10


