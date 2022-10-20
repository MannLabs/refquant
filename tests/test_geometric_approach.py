import refquant.refquant_classes as rqc
import unittest
import pandas as pd
import numpy as np

class TestGeometricCalculation(unittest.TestCase):

    def test_shortening_of_df_dict(self):
        for idx in range(10):
            spectrum_generator = TargetReferenceSpectrumGenerator()

            ratio_target_reference = rqc.GeometricRatioCalculator(spectrum_generator.list_of_intersecting_ions, spectrum_generator.target_dict, spectrum_generator.reference_dict).ratio_to_reference
            ratio_target_reference_alternative_calculation = alternative_implementation_of_ratio_calculation(spectrum_generator.reference_dict, spectrum_generator.target_dict, minimum_number_required_fragions=3)

            # print(ratio_target_reference)
            # print(ratio_target_reference_alternative_calculation)
            # print("iteration through")
            self.assertAlmostEqual(ratio_target_reference, ratio_target_reference_alternative_calculation)
            print(f"ratio {ratio_target_reference}\nratio alternative{ratio_target_reference_alternative_calculation}")



class TargetReferenceSpectrumGenerator():
    def __init__(self):
        self.target_dict = None
        self.reference_dict = None
        self.list_of_intersecting_ions = None
        self._generate_target_and_reference_spectrum()
        self._generate_list_of_intersecting_ions()
    
    def _generate_target_and_reference_spectrum(self):
        num_peaks =  np.random.randint(3, 13)
        self.target_dict = {f"{x}" : np.log2(abs(np.random.random()+0.2)*100.000) for x in range(num_peaks)}
        self.reference_dict = {f"{x}" : np.log2(abs(np.random.random()) * 100.000) for x in range(num_peaks)}
    
    def _generate_list_of_intersecting_ions(self):
        self.list_of_intersecting_ions = list(set(self.target_dict.keys()) & set(self.reference_dict.keys()))
    



def alternative_implementation_of_ratio_calculation(fragion2intensity_dict_spec_lib, fragion2intensity_dict_target, minimum_number_required_fragions):
    fragion2intensity_dict_spec_lib = {k : 2**v for k,v in fragion2intensity_dict_spec_lib.items()}
    fragion2intensity_dict_target = {k : 2**v for k,v in fragion2intensity_dict_target.items()}

    fragion2intensity_dict_spec_lib_sorted = {k: v for k, v in sorted(fragion2intensity_dict_spec_lib.items(), key=lambda item: item[1])[::-1]}
    frag_lookup_idx = {k:idx for idx, k in enumerate(fragion2intensity_dict_spec_lib_sorted.keys())}
    shared_frags = [_ for _ in fragion2intensity_dict_spec_lib_sorted.keys() if _ in fragion2intensity_dict_target.keys()]
    
    y = np.fromiter(fragion2intensity_dict_spec_lib_sorted.values(), dtype=float)
    x = np.arange(len(y))

    y_ = np.array([fragion2intensity_dict_target[_] for _ in shared_frags])
    x_ = np.array([frag_lookup_idx[_] for _ in shared_frags])

    slope, intercept = np.polyfit(x,y,1)

    best_idx = minimum_number_required_fragions
    best_diff = np.inf

    for i in range(minimum_number_required_fragions, len(x_)):
        slope_, intercept_ = np.polyfit(x_[:i], y_[:i], 1)
        
        diff = np.abs(1-slope_/slope)
        if diff < best_diff:
            best_idx = i
            best_diff = diff

    print(f"alternative best idx {best_idx}")
    shortened_fragion2intensity_dict_target = {k: fragion2intensity_dict_target[k] for k in shared_frags[:best_idx]}
    ratio_to_reference = -np.mean(np.log2(([fragion2intensity_dict_spec_lib[k]/v for k,v in shortened_fragion2intensity_dict_target.items()])))


    return ratio_to_reference


