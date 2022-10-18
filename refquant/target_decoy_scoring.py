from re import T
import numpy as np

class TargetDecoyScorer():
    def __init__(self, list_of_single_labelled_precursor):
        self._list_of_single_labelled_precursors = list_of_single_labelled_precursor
        self._list_of_ranked_precursors = []
        self._precursorname2fdr = {}


        self._filter_list_of_single_labelled_precursors_to_valid_precursors()
        self._define_list_of_ranked_precursors()
        self._define_precursor2fdr()
        self._add_fdr_to_precursors()

    
    def _filter_list_of_single_labelled_precursors_to_valid_precursors(self):
        self._list_of_single_labelled_precursors = [x for x in self._list_of_single_labelled_precursors if x.cosine_similarity is not None]

    def _define_list_of_ranked_precursors(self):
        for precursor in self._list_of_single_labelled_precursors:
            ranked_prec_target = RankedPrecursor(name = precursor.name, score = precursor.cosine_similarity, is_decoy = False)
            ranked_prec_decoy = RankedPrecursor(name = precursor.name, score = precursor.decoy_cosine_similarity, is_decoy = True)
            self._list_of_ranked_precursors.append(ranked_prec_target)
            self._list_of_ranked_precursors.append(ranked_prec_decoy)
        self._list_of_ranked_precursors = sorted(self._list_of_ranked_precursors, key=lambda x: x.score, reverse=True)
        self._add_ranking_to_rankedprecursors()

    
    def _add_ranking_to_rankedprecursors(self):
        num_decoys = 0
        for idx in range(len(self._list_of_ranked_precursors)):
            ranked_precursor = self._list_of_ranked_precursors[idx]
            ranked_precursor.list_position = idx+1
            ranked_precursor.num_decoys = num_decoys
            ranked_precursor.calculate_empirical_FDR()
            if ranked_precursor.is_decoy:
                num_decoys += 1

    def _define_precursor2fdr(self):
        for rankedprec in self._list_of_ranked_precursors:
            if not rankedprec.is_decoy:
                self._precursorname2fdr[rankedprec.name] = rankedprec.empirical_fdr

    def _add_fdr_to_precursors(self):
        for precursor in self._list_of_single_labelled_precursors:
            precursor.empirical_fdr = self._precursorname2fdr[precursor.name]

class RankedPrecursor():
    def __init__(self, name, score, is_decoy):
        self.name = name
        self.is_decoy = is_decoy
        self.score = score
        self.list_position = None
        self.num_decoys = None

        self.empirical_fdr = None

    def calculate_empirical_FDR(self):
        self.empirical_fdr = 2*self.num_decoys / self.list_position