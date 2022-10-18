import numpy as np

class GeometricRatioCalculator():

    def __init__(self, list_of_intersection_ions, fragion2quantity_target,fragion2quantity_reference, min_number_ions = 3):
        self._list_of_intersected_ions = list_of_intersection_ions
        self._fragion2quantity_target = fragion2quantity_target
        self._fragion2quantity_reference = fragion2quantity_reference
        self._min_number_ions = min_number_ions
        self._transform_fragion2quantities_from_log2_to_linear()

        self._intensity_sorted_ions = None
        self._reference_slope = None
        self._min_slopediff = np.inf
        self._best_cutoff_idx = min_number_ions

        self.ratio_to_reference = None
    
        self._define_intensity_sorted_ions()
        self._define_reference_slope()
        self._define_best_cutoff_idx()
        self._define_ratio_to_reference()
    
    def _transform_fragion2quantities_from_log2_to_linear(self):
        self._fragion2quantity_reference = {k:2**v for k,v in self._fragion2quantity_reference.items()}
        self._fragion2quantity_target = {k:2**v for k,v in self._fragion2quantity_target.items()}

    def _define_intensity_sorted_ions(self):
        self._intensity_sorted_ions = list(sorted(self._list_of_intersected_ions, key=lambda x : self._fragion2quantity_reference.get(x), reverse=True))
    
    def _define_reference_slope(self):
        self._reference_slope = self._calculate_slope_for_n_first_ions(fragion2quantity_dict=self._fragion2quantity_reference, n=len(self._intensity_sorted_ions))       

    def _define_best_cutoff_idx(self):
        for idx in range(self._min_number_ions, len(self._intensity_sorted_ions)):
            target_slope = self._calculate_slope_for_n_first_ions(fragion2quantity_dict=self._fragion2quantity_target, n=idx)
            slopediff = self._calculate_difference_to_reference_slope(target_slope)
            self._upate_best_cutoff_and_min_slopediff_if_applicable(slopediff, idx)
    
    def _calculate_slope_for_n_first_ions(self, fragion2quantity_dict, n):
        n_first_ions = self._intensity_sorted_ions[:n]
        sorted_intensities = np.array([fragion2quantity_dict[fragion] for fragion in n_first_ions])
        idxs_sorted_intensities = np.array(range(len(sorted_intensities)))
        slope, _ = np.polyfit(idxs_sorted_intensities,sorted_intensities,1)
        return slope
    
    def _calculate_difference_to_reference_slope(self, target_slope):
        return np.abs(1-target_slope/self._reference_slope) 

    def _upate_best_cutoff_and_min_slopediff_if_applicable(self, slopediff, idx):
        if slopediff < self._min_slopediff:
            self._min_slopediff = slopediff
            self._best_cutoff_idx = idx
    
    def _define_ratio_to_reference(self):
        selected_intensity_ions = self._intensity_sorted_ions[:self._best_cutoff_idx]
        intensities_target = np.array([self._fragion2quantity_target[x] for x in selected_intensity_ions])
        intensities_reference = np.array([self._fragion2quantity_reference[x] for x in selected_intensity_ions])
        ratios_to_references = np.log2(intensities_target/intensities_reference)
        self.ratio_to_reference = np.mean(ratios_to_references)
