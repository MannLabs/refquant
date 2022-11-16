import numpy as np

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
