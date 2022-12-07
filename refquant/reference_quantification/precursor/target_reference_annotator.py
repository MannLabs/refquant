import numpy as np
from . import precursor_classes

class PrecursorWMatchedLabelsAnnotator():
    def __init__(self, precursor_with_matched_labels : precursor_classes.PrecursorWithMatchedLabels):
        self._reference_precursor = precursor_with_matched_labels.reference_precursor
        self._list_of_target_precursors = precursor_with_matched_labels.list_of_target_precursors

    def annotate_precursor(self):
        for target_precursor in self._list_of_target_precursors:
            TargetReferenceAnnotator(self._reference_precursor, target_precursor)


class TargetReferenceAnnotator():
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
        idx_quantile = self._get_index_of_quantile(0.4)
        self.target_precursor.min_ratio_to_reference = sorted_ratios[idx_quantile_min]
        self.target_precursor.ratio_to_reference = np.mean(sorted_ratios[:idx_quantile])


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
        
    def _annotate_number_of_fragment_ions_available(self):
        self.target_precursor.number_of_fragment_ions_used = len(set(filter(lambda x : "FRGION" in x, self._list_of_intersection_ions)))


