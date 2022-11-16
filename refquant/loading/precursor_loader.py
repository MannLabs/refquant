import numpy as np
import pandas as pd

from ..precursor_handling import precursor_objects
from ..precursor_handling import precursor_initializer

class PrecursorLoader():
    def __init__(self, formatted_input_file, replicate_name):
        self._formatted_input_file = formatted_input_file
        self._replicate_name = replicate_name
        
        self._merged_precusor_and_fragion_df = None
        self._precursor2df = {}

        self.precursors_w_all_labels = []

        self._define_merged_precusor_and_fragion_df()
        self._define_precursor2df()
        self._go_through_precursors_and_initialize_precursors_w_all_labels()


    def _define_merged_precusor_and_fragion_df(self):
        self._merged_precusor_and_fragion_df = pd.read_csv(self._formatted_input_file, sep="\t")
        self._subset_merged_precusor_and_fragion_df_to_replicate()
        self._adapt_formatting_of_dataframe()
        #self._merged_precusor_and_fragion_df = self._merged_precusor_and_fragion_df.set_index("precursor")
    
    def _adapt_formatting_of_dataframe(self):
        self._merged_precusor_and_fragion_df["run"] = self._merged_precusor_and_fragion_df["run"].astype(str)
        self._merged_precusor_and_fragion_df =  self._merged_precusor_and_fragion_df.replace(0, np.nan)
        numeric_cols = self._merged_precusor_and_fragion_df.select_dtypes(include=[np.number]).columns
        self._merged_precusor_and_fragion_df[numeric_cols] = np.log2(self._merged_precusor_and_fragion_df[numeric_cols])
    
    def _subset_merged_precusor_and_fragion_df_to_replicate(self):
        #subset datframe to rows containing the replicate name
        self._merged_precusor_and_fragion_df = self._merged_precusor_and_fragion_df[self._merged_precusor_and_fragion_df["run"] == self._replicate_name]

    def _define_precursor2df(self):
        self._precursor2df = {k:v for k, v in self._merged_precusor_and_fragion_df.groupby('precursor')}
        self._clear_merged_precusor_and_fragion_df()
    
    def _go_through_precursors_and_initialize_precursors_w_all_labels(self):
        num_precursors = len(self._precursor2df.keys())
        prec_counter = 0
        for precursor in self._precursor2df.keys():
            if prec_counter%1000==0:
                print("{}/{}".format(prec_counter, num_precursors))
            prec_counter += 1
            precursor_df = self._precursor2df[precursor]
            list_of_single_labelled_precursors = self._get_list_of_single_labelled_precursors(precursor, precursor_df)
            precursor_w_all_labels = precursor_objects.PrecursorWithAllLabels(list_of_single_labelled_precursors)
            self.precursors_w_all_labels.append(precursor_w_all_labels)
        self._clear_precursor2df()
    
    def _get_list_of_single_labelled_precursors(self, precursor, precursor_df):
        return precursor_initializer.PrecursorsFromDataframeInititalizer(precursor, precursor_df).single_labelled_precursors
    
    def _clear_merged_precusor_and_fragion_df(self):
        self._merged_precusor_and_fragion_df = None
    
    def _clear_precursor2df(self):
        self._precursor2df = None

class PrecursorLoaderFromDf(PrecursorLoader):
    def __init__(self, merged_precusor_and_fragion_df):

        self._merged_precusor_and_fragion_df = merged_precusor_and_fragion_df
        self._precursor2df = {}

        self.precursors_w_all_labels = []

        self._define_merged_precusor_and_fragion_df()
        self._define_precursor2df()
        self._go_through_precursors_and_initialize_precursors_w_all_labels()
    

class PrecursorLoaderDIANN(PrecursorLoader):
    def __init__(self, formatted_input_file, replicate_name):
        super().__init__(formatted_input_file, replicate_name)
   

class PrecursorLoaderDIANNRef8(PrecursorLoader):
    def __init__(self, formatted_input_file, replicate_name):
        super().__init__(formatted_input_file, replicate_name)

    def _get_list_of_single_labelled_precursors(self, precursor, precursor_df):
        return precursor_initializer.PrecursorsFromDataframeInititalizerReference8(precursor, precursor_df).single_labelled_precursors



class PrecursorLoaderDIANNFromDf(PrecursorLoaderDIANN):
    def __init__(self, merged_precusor_and_fragion_df):

        self._merged_precusor_and_fragion_df = merged_precusor_and_fragion_df

        self._precursor2df = {}
        
        self.precursors_w_all_labels = []
        self._adapt_formatting_of_dataframe()
        self._define_precursor2df()
        self._go_through_precursors_and_initialize_precursors_w_all_labels()
    

