# %% [markdown]
# ## Define labelled precursor objects

# %%
import re
from timeit import repeat
import numpy as np
import refquant.geometric_ratio as geometric_ratio
import refquant.refquant_utils as utils

class SingleLabelledPrecursor():
    def __init__(self):
        self.name = None
        self.protein_name = None
        self.channel_name = None
        self.replicate_name = None

        self.fragion2quantity = {}
        
        self.comparison_derived_quantity = np.nan
        
        self.ratio_to_reference = np.nan
        self.median_ratio_to_reference = np.nan
        self.min_ratio_to_reference = np.nan
        self.ratio_to_reference_intensity_based = np.nan
        self.ratio_of_most_abundant_fragion_to_reference = np.nan
        self.geometric_ratio = np.nan
        self.search_engine_derived_quantity = np.nan

        self.number_of_ratios_used = np.nan
        self.number_of_fragment_ions_used = np.nan
        self.cosine_similarity = np.nan
        self.decoy_cosine_similarity = np.nan
        self.search_engine_derived_quantity_not_provided = np.nan

        self.search_engine_derived_reference_quantity = np.nan
        self.search_engine_derived_reference_quantity_normed = np.nan
        self.derived_reference_quantity_static = np.nan


class PrecursorWithAllLabels():
    def __init__(self, list_of_single_labelled_precursors):
        self.list_of_single_labelled_precursors = list_of_single_labelled_precursors

class PrecursorWithMatchedLabels():
    def __init__(self,  reference_precursor, list_of_target_precursors):     
        self.reference_precursor = reference_precursor
        self.list_of_target_precursors = list_of_target_precursors

        self.annotate_precursors()

    def annotate_precursors(self):
        for target_precursor in self.list_of_target_precursors:
            TargetPrecursorAnnotator(self.reference_precursor, target_precursor)


class TargetPrecursorAnnotator():
    def __init__(self, reference_precursor, target_precursor):
        self.reference_precursor = reference_precursor
        self.target_precursor = target_precursor

        self._intensities_target = None
        self._intensities_reference = None
        self._list_of_intersection_ions = None
        self._ratios_to_reference = None

        self._define_intersecting_fragment_ions()
        self._define_intensities_of_target_and_reference()
        self._define_ratios_to_reference()

        self._annotate_number_of_ratios_used_to_precursor()
        self._annotate_intensity_based_reference_ratio()
        self._annotate_search_engine_derived_reference_quantity()
        self._annotate_ms1_reference_quantity()
        self._annotate_summed_top5_reference_quantity()

        self._annotate_comparison_derived_quantity_to_precursor()
        self._annotate_ms1_ratio_and_intensity()
        self._annotate_derived_ratio()
        
    def _define_intersecting_fragment_ions(self):
        self._list_of_intersection_ions = list(set(self.reference_precursor.fragion2quantity.keys()).intersection(set(self.target_precursor.fragion2quantity.keys())))
    
    def _define_intensities_of_target_and_reference(self):
        self._intensities_target = np.array([self.target_precursor.fragion2quantity[fragion] for fragion in self._list_of_intersection_ions])
        self._intensities_reference = np.array([self.reference_precursor.fragion2quantity[fragion] for fragion in self._list_of_intersection_ions])
        if np.nan in self._intensities_target or np.nan in self._intensities_reference:
            raise Exception("Nans not filtered as expected!") 
    
    def _define_ratios_to_reference(self):
        self._ratios_to_reference = self._intensities_target - self._intensities_reference #the intensities need be be log2 transformed



    def _annotate_number_of_ratios_used_to_precursor(self):
        self.target_precursor.number_of_ratios_used = len(self._list_of_intersection_ions)
    
    def _annotate_derived_ratio(self):
        if self.target_precursor.number_of_ratios_used == 0:
            self.target_precursor.derived_ratio = np.nan
            return
        self.target_precursor.median_ratio_to_reference = np.median(self._ratios_to_reference)
        sorted_ratios = np.sort(self._ratios_to_reference)
        idx_quantile_min = self._get_index_of_quantile(0.1)
        idx_quantile = self._get_index_of_quantile(0.25)
        self.target_precursor.min_ratio_to_reference = sorted_ratios[idx_quantile_min]
        self.target_precursor.ratio_to_reference = sorted_ratios[idx_quantile]


    def _get_index_of_quantile(self,quantile):
        return int(quantile * len(self._ratios_to_reference))

    def _annotate_ms1_ratio_and_intensity(self):
        is_ms1 = ["MS1" in x for x in self._list_of_intersection_ions]
        if sum(is_ms1)==1:
            ms1_ratio = self._ratios_to_reference[is_ms1][0]
            ms1_intensity = self._intensities_target[is_ms1][0]
        elif sum(is_ms1) == 0:
            ms1_ratio = np.nan
            ms1_intensity = np.nan
        else:
            raise ValueError("More than one MS1 ion in intersection")

        self.target_precursor.ms1_ratio_to_reference = ms1_ratio
        self.target_precursor.ms1_intensity = ms1_intensity
    
    def _annotate_search_engine_derived_reference_quantity(self):
        self.target_precursor.search_engine_derived_quantity_reference = self.reference_precursor.search_engine_derived_quantity

    def _annotate_ms1_reference_quantity(self):
        is_ms1 = ["MS1" in x for x in self._list_of_intersection_ions]
        if sum(is_ms1)==1:
            self.target_precursor.ms1_quantity_reference = self._intensities_reference[is_ms1][0]
        else:
            self.target_precursor.ms1_quantity_reference = np.nan

    def _annotate_summed_top5_reference_quantity(self):
        sorted_intensities_descending = np.sort(self._intensities_reference)[::-1]
        self.target_precursor.summed_quantity_reference = np.log2(np.sum(2**sorted_intensities_descending[:5]))

    def _annotate_intensity_based_reference_ratio(self):
        if self.target_precursor.search_engine_derived_quantity is not None and self.reference_precursor.search_engine_derived_quantity is not None:
            self.target_precursor.ratio_to_reference_intensity_based = self.target_precursor.search_engine_derived_quantity - self.reference_precursor.search_engine_derived_quantity
    
    def _annotate_comparison_derived_quantity_to_precursor(self):
        self.target_precursor.comparison_derived_quantity = self.target_precursor.ratio_to_reference + self.reference_precursor.search_engine_derived_quantity
        
    def _annotate_ratio_of_most_abundant_fragion_to_reference(self):
        is_fragion = ["FRGION" in x for x in self._list_of_intersection_ions]
        if sum(is_fragion)>0:
            ratios_to_reference_just_fragions = self._ratios_to_reference[is_fragion]
            intensities_target_just_fragions = self._intensities_target[is_fragion]
            self.target_precursor.ratio_of_most_abundant_fragion_to_reference = ratios_to_reference_just_fragions[np.argmax(intensities_target_just_fragions)]
        
    def _annotate_cosine_similarity(self):
        self.target_precursor.cosine_similarity = CosineSimilarityCalculator(self._intensities_target, self._intensities_reference).get_cosine_similarity()
        self.target_precursor.decoy_cosine_similarity = CosineSimilarityCalculatorShuffled(self._intensities_target, self._intensities_reference).get_cosine_similarity()

    def _annotate_number_of_fragment_ions_available(self):
        self.target_precursor.number_of_fragment_ions_used = len(set(filter(lambda x : "FRGION" in x, self._list_of_intersection_ions)))

    def _annotate_geometric_ratio(self):
        if len(self._list_of_intersection_ions)>2:
            calculator = geometric_ratio.GeometricRatioCalculator(self._list_of_intersection_ions, self.target_precursor.fragion2quantity, self.reference_precursor.fragion2quantity)
            self.target_precursor.geometric_ratio = calculator.ratio_to_reference



class CosineSimilarityCalculator():
    def __init__(self, intensities_target_logged, intensities_reference_logged):
        self._intensities_target_logged = intensities_target_logged
        self._intensities_reference_logged = intensities_reference_logged
        
        self._intensities_target_linscaled = None
        self._intensities_reference_linscaled = None
        self._define_intensities_target_and_reference_linscaled()

    
    def _define_intensities_target_and_reference_linscaled(self):
        self._intensities_target_linscaled = 2**self._intensities_target_logged
        self._intensities_reference_linscaled = 2**self._intensities_reference_logged

    def get_cosine_similarity(self):
        dotscore = np.dot(self._intensities_target_linscaled, self._intensities_reference_linscaled)
        normalization = np.linalg.norm(self._intensities_target_linscaled) * np.linalg.norm(self._intensities_reference_linscaled)
        cosine_similarity = dotscore / normalization
        return cosine_similarity


class CosineSimilarityCalculatorShuffled(CosineSimilarityCalculator):
    def __init__(self, intensities_target_logged, intensities_reference_logged):
        super().__init__(intensities_target_logged, intensities_reference_logged)
        self._shuffle_target_intensities()
    
    def _shuffle_target_intensities(self):
        self._intensities_target_linscaled = np.random.permutation(self._intensities_target_linscaled)


import numpy as np
import pandas as pd

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
            precursor_w_all_labels = PrecursorWithAllLabels(list_of_single_labelled_precursors)
            self.precursors_w_all_labels.append(precursor_w_all_labels)
        self._clear_precursor2df()
    
    def _get_list_of_single_labelled_precursors(self, precursor, precursor_df):
        return PrecursorsFromDataframeInititalizer(precursor, precursor_df).single_labelled_precursors
    
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
        return PrecursorsFromDataframeInititalizerReference8(precursor, precursor_df).single_labelled_precursors



class PrecursorLoaderDIANNFromDf(PrecursorLoaderDIANN):
    def __init__(self, merged_precusor_and_fragion_df):

        self._merged_precusor_and_fragion_df = merged_precusor_and_fragion_df

        self._precursor2df = {}
        
        self.precursors_w_all_labels = []
        self._adapt_formatting_of_dataframe()
        self._define_precursor2df()
        self._go_through_precursors_and_initialize_precursors_w_all_labels()
    




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

class SingleLabelledPrecursorIntitializer():
    def __init__(self, precursor_name, precursor_df, channel_name):
        self._precursor_name = precursor_name
        self._precursor_df = precursor_df
        self._channel_name = channel_name

        self.single_labelled_precursor = SingleLabelledPrecursor()

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




# %% [markdown]
# ## combine precursor objects

# %%
from abc import ABC, abstractmethod
class PrecursorCombiner(ABC):
    def __init__(self, precursorloader : PrecursorLoader):
        self._precursorloader = precursorloader

        self.list_of_precursors_with_matched_labels = []

        self._define_list_of_precursors_with_matched_labels()
    
    def _define_list_of_precursors_with_matched_labels(self):
        for precursor_w_all_labels in self._precursorloader.precursors_w_all_labels:
            self._extend_list_of_precursors_w_matched_labels(precursor_w_all_labels)
    
    @abstractmethod
    def _extend_list_of_precursors_w_matched_labels(self):
        pass

class PrecursorCombinerToSingleReference(PrecursorCombiner):
    def _extend_list_of_precursors_w_matched_labels(self, precursor_w_all_labels):
        matched_precursor = self._define_precursor_w_matched_labels_using_target_and_reference(precursor_w_all_labels)
        self.list_of_precursors_with_matched_labels.append(matched_precursor)


    def _define_precursor_w_matched_labels_using_target_and_reference(self, precursor_w_all_labels):
        list_of_target_precursors = list(filter(lambda x : "target" in x.channel_name , precursor_w_all_labels.list_of_single_labelled_precursors))
        reference_precursor = list(filter(lambda x : "reference" in x.channel_name , precursor_w_all_labels.list_of_single_labelled_precursors))[0]
        precursor_with_matched_labels = PrecursorWithMatchedLabels(reference_precursor=reference_precursor, list_of_target_precursors=list_of_target_precursors)
        return precursor_with_matched_labels

class PrecursorCombinerPairWiseToReference(PrecursorCombiner):
    def _extend_list_of_precursors_w_matched_labels(self, precursor_w_all_labels):
        matched_precursor_sn_run1 = self._define_precursor_w_matched_labels_using_target_and_reference(precursor_w_all_labels, target_column="target4", reference_column="reference4")
        matched_precursor_sn_run2 = self._define_precursor_w_matched_labels_using_target_and_reference(precursor_w_all_labels, target_column="target8", reference_column="reference8")
        self.list_of_precursors_with_matched_labels.append(matched_precursor_sn_run1)
        self.list_of_precursors_with_matched_labels.append(matched_precursor_sn_run2)


    def _define_precursor_w_matched_labels_using_target_and_reference(self, precursor_w_all_labels, target_column, reference_column):
        target_precursor = list(filter(lambda x : x.channel_name == target_column, precursor_w_all_labels.list_of_single_labelled_precursors))[0]
        reference_precursor = list(filter(lambda x : x.channel_name == reference_column, precursor_w_all_labels.list_of_single_labelled_precursors))[0]
        precursor_with_matched_labels = PrecursorWithMatchedLabels(reference_precursor=reference_precursor, list_of_target_precursors=[target_precursor])
        return precursor_with_matched_labels


import multiprocess
def get_all_single_labelled_precursors_in_dataset(reference_table, use_multiprocessing=False):
    single_labelled_precursors = []
    runs = utils.get_runs(reference_table)
    print(f"processing {len(runs)} runs: {runs}")
    if use_multiprocessing:
        multiprocess.freeze_support()
        args = [(x, reference_table)   for x in runs]
        with multiprocess.Pool(int(multiprocess.cpu_count()/2)) as pool:
            single_labelled_precursors = pool.starmap(get_single_labelled_precursors, args)
        #join list of lists
        single_labelled_precursors = [item for sublist in single_labelled_precursors for item in sublist]
    else:
        for run in runs:
            single_labelled_precursors.extend(get_single_labelled_precursors(run, reference_table))
    return single_labelled_precursors

def get_single_labelled_precursors(run, reference_table):
    print("run: ", run)
    single_labelled_precursors = []
    precursor_loader = PrecursorLoader(reference_table, run)
    label_combiner = PrecursorCombinerToSingleReference(precursor_loader)
    for matched_prec in label_combiner.list_of_precursors_with_matched_labels:
        for target_precursor in matched_prec.list_of_target_precursors:
            single_labelled_precursors.append(target_precursor)
    return single_labelled_precursors


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
