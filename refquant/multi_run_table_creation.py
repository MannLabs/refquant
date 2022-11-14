import re
import pandas as pd
import os

class PrecursorQuantityTableCreator():
    def __init__(self, list_of_single_labelled_precursors):
        self._list_of_single_labelled_precursors = list_of_single_labelled_precursors

        self._precursor2protein = {}
        self._replicate_and_channel_to_precursors = {}
        self._list_of_precursors_series = []

        self.precursorquantitytable = None

        self._add_static_cross_table_quantities_to_single_labelled_precursors()
        self._define_precursor2protein()
        self._define_run_and_channel_to_precursors()
        self._reformat_run_and_channel_precursors_to_precursor_indexed_series()
        self._define_precursor2quantitytable()

    def _add_static_cross_table_quantities_to_single_labelled_precursors(self):
        StaticReferenceChannelQuantityAdder(self._list_of_single_labelled_precursors)

    def _define_precursor2protein(self):
        for precursor in self._list_of_single_labelled_precursors:
            self._precursor2protein[precursor.name] = precursor.protein_name

    def _define_run_and_channel_to_precursors(self):
        for precursor in self._list_of_single_labelled_precursors:
            replicate_and_channel = self._get_replicate_and_channel_name(precursor)
            if replicate_and_channel not in self._replicate_and_channel_to_precursors:
                self._replicate_and_channel_to_precursors[replicate_and_channel] = []
            self._replicate_and_channel_to_precursors[replicate_and_channel].append(precursor)
    
    def _reformat_run_and_channel_precursors_to_precursor_indexed_series(self):
        for replicate_and_channel, list_of_precursors in self._replicate_and_channel_to_precursors.items():
            series_values = self._get_quantitative_values(list_of_precursors)
            series_index = [precursor.name for precursor in list_of_precursors]
            series_quantities = pd.Series(series_values, index = series_index)
            series_quantities.name = replicate_and_channel
            self._list_of_precursors_series.append(series_quantities)

    def _get_quantitative_values(self, list_of_precursors):
        return [precursor.derived_reference_quantity_static + precursor.ratio_to_reference for precursor in list_of_precursors]

    def _define_precursor2quantitytable(self):
        self.precursorquantitytable = pd.concat(self._list_of_precursors_series, axis = 1)
        self.precursorquantitytable.index.name = "precursor"
        self.precursorquantitytable["protein"] = [self._precursor2protein[precursor] for precursor in self.precursorquantitytable.index]
        self.precursorquantitytable = self.precursorquantitytable.sort_values(by = "protein")
        self.precursorquantitytable = self.precursorquantitytable.reset_index().set_index(["protein", "precursor"])
        self.precursorquantitytable = 2**self.precursorquantitytable
        self.precursorquantitytable = self.precursorquantitytable.fillna(0)
        self.precursorquantitytable = self.precursorquantitytable.reset_index()
        self.precursorquantitytable = self.precursorquantitytable.rename(columns = {"precursor": "ion"})
    
    @staticmethod
    def _get_replicate_and_channel_name(precursor):
        return f"{precursor.replicate_name}_{precursor.channel_name}"


class PrecursorTableCreatorQuantityBySearchEngine(PrecursorQuantityTableCreator):
    def _get_quantitative_values(self, list_of_precursors):
        return [precursor.search_engine_derived_quantity for precursor in list_of_precursors]


class PrecursorTableCreatorQuantityByMS1(PrecursorQuantityTableCreator):
    def _get_quantitative_values(self, list_of_precursors):
        return [precursor.ms1_intensity for precursor in list_of_precursors]

class PrecursorTableCreatorMS1RatioToReference(PrecursorQuantityTableCreator):
    def _get_quantitative_values(self, list_of_precursors):
        return [precursor.derived_reference_quantity_static + precursor.ms1_ratio_to_reference for precursor in list_of_precursors]

class PrecursorTableCreatorNoStaticReference(PrecursorQuantityTableCreator):
    def _get_quantitative_values(self, list_of_precursors):
        return [precursor.search_engine_derived_quantity_reference + precursor.ratio_to_reference for precursor in list_of_precursors]




class StaticReferenceChannelQuantityAdder():
    def __init__(self, list_of_single_labelled_precursors):

        self._quantities_to_use = ["ms1_quantity_reference", "search_engine_derived_quantity_reference", "summed_quantity_reference"]
        self._list_of_single_labelled_precursors = list_of_single_labelled_precursors

        self._iteratively_add_reference_quantities()

    def _iteratively_add_reference_quantities(self):
        unassigned_precursors = [x for x in self._list_of_single_labelled_precursors]
        for quantity in self._quantities_to_use:
            unassigned_precursors = self._get_unassigned_precursors(unassigned_precursors)
            unassigned_precursors_w_quantity_values =  self._filter_for_precursors_w_quantity_values(unassigned_precursors, quantity)
            print(f"Trying to assign {quantity} to {len(unassigned_precursors_w_quantity_values)} precursors")
            if len(unassigned_precursors_w_quantity_values) >100:
                StaticReferenceChannelQuantityAdderForSelectedQuantity(unassigned_precursors, quantity)

    def _get_unassigned_precursors(self, single_labelled_precursors):
        precursors_w_no_static_quantity = [precursor for precursor in single_labelled_precursors if np.isnan(precursor.derived_reference_quantity_static)]
        return precursors_w_no_static_quantity

    def _filter_for_precursors_w_quantity_values(self, unassigned_precursors, quantity):
        return [precursor for precursor in unassigned_precursors if not np.isnan(getattr(precursor, quantity))]

class StaticReferenceChannelQuantityAdderForSelectedQuantity():
    def __init__(self, list_of_single_labelled_precursors, selected_quantity):
        
        self._list_of_single_labelled_precursors = list_of_single_labelled_precursors
        self._selected_quantity = selected_quantity
        self._run2singlelabelledPrecursors = {}
        self._list_of_annotated_series = []

        self._reference_intensity_dataframe : pd.DataFrame = None

        self._define_run2singlelabelledPrecursors()
        self._define_list_of_annotated_series()
        if len(self._list_of_annotated_series) > 0:
            self._merge_list_of_annotated_series_to_reference_intensity_dataframe()
            self._normalize_reference_intensity_and_annotate_precursors()
            self._add_static_quantities()

    def _define_run2singlelabelledPrecursors(self):
        for precursor in self._list_of_single_labelled_precursors:
            if precursor.replicate_name not in self._run2singlelabelledPrecursors.keys():
                self._run2singlelabelledPrecursors[precursor.replicate_name] = []
            self._run2singlelabelledPrecursors[precursor.replicate_name].append(precursor)

    def _define_list_of_annotated_series(self):
        for run, list_of_precursors in self._run2singlelabelledPrecursors.items():
            list_of_precursors = [x for x in list_of_precursors if hasattr(x, self._selected_quantity)]
            list_of_precursors = [x for x in list_of_precursors if np.isfinite(getattr(x, self._selected_quantity))]
            series_values = [getattr(precursor, self._selected_quantity) for precursor in list_of_precursors]
            series_index = [precursor.name for precursor in list_of_precursors]
            series_quantities = pd.Series(series_values, index = series_index, dtype='float')
            duplicate_indices = series_quantities.index[series_quantities.index.duplicated()]
            series_quantities = series_quantities.drop(duplicate_indices)
            series_quantities.name = run
            if len(series_quantities.index)>0:
                self._list_of_annotated_series.append(series_quantities)

    def _merge_list_of_annotated_series_to_reference_intensity_dataframe(self):
        self._reference_intensity_dataframe = pd.concat(self._list_of_annotated_series, axis = 1, join="outer")

    def _normalize_reference_intensity_and_annotate_precursors(self):
        self._reference_intensity_dataframe = ReferenceChannelNormalizer(self._reference_intensity_dataframe, self._run2singlelabelledPrecursors).reference_intensity_dataframe

    def _add_static_quantities(self):
        ReferenceChannelQuantityDeriver(self._reference_intensity_dataframe, self._list_of_single_labelled_precursors).add_static_quantities()


import numpy as np
class ReferenceChannelNormalizer():
    def __init__(self, reference_intensity_dataframe : pd.DataFrame, run2singlelabelledPrecursors : dict):
        self.reference_intensity_dataframe = reference_intensity_dataframe
        self._run2singlelabelledPrecursors = run2singlelabelledPrecursors

        self._run_id_to_norm_towards = None
        self._run2shift = {}
        self._define_run_id_to_norm_towards()
        self._define_run2shift()
        self._apply_run2shift()
    
    def _define_run_id_to_norm_towards(self):
        self._run_id_to_norm_towards = self.reference_intensity_dataframe.median(axis = 0).idxmax()

    def _define_run2shift(self):
        runs_to_shift = [run for run in self.reference_intensity_dataframe.columns if run != self._run_id_to_norm_towards]
        for run in runs_to_shift:
            fcs_between_runs = self.reference_intensity_dataframe[run].values - self.reference_intensity_dataframe[self._run_id_to_norm_towards].values
            self._run2shift[run] = np.nanmedian(fcs_between_runs)

    def _apply_run2shift(self):
        for run, shift in self._run2shift.items():
            self._normalize_df_column(run, shift)
            self._add_normalized_quantity_to_precursor(run, shift)

    def _normalize_df_column(self, run, shift):
        self.reference_intensity_dataframe[run] = self.reference_intensity_dataframe[run] + shift
    
    def _add_normalized_quantity_to_precursor(self, run, shift):
        for precursor in self._run2singlelabelledPrecursors[run]:
            precursor.search_engine_derived_reference_quantity_normed = precursor.search_engine_derived_reference_quantity + shift

import refquant.refquant_classes as refquant_classes
class ReferenceChannelQuantityDeriver():
    def __init__(self, reference_intensity_dataframe : pd.DataFrame, list_of_single_labelled_precursors : list([refquant_classes.SingleLabelledPrecursor])):
        self._reference_intensity_dataframe = reference_intensity_dataframe
        self._list_of_single_labelled_precursors = list_of_single_labelled_precursors

        self._precursorname2staticquantity = {}
        self._define_precursorname2staticquantity()

    def add_static_quantities(self):
        for singlelabelledprecursor in self._list_of_single_labelled_precursors:
            self._annotate_static_reference_quantity(singlelabelledprecursor)
            self._annotate_target_quantity_based_on_comparison_to_static(singlelabelledprecursor)
        
    def _annotate_static_reference_quantity(self, singlelabelledprecursor):
        if singlelabelledprecursor.name in self._precursorname2staticquantity.keys():
            singlelabelledprecursor.derived_reference_quantity_static = self._precursorname2staticquantity[singlelabelledprecursor.name]


    def _annotate_target_quantity_based_on_comparison_to_static(self, singlelabelledprecursor):
        if singlelabelledprecursor.ratio_to_reference is not None and hasattr(singlelabelledprecursor, "derived_reference_quantity_static"):
            singlelabelledprecursor.comparison_derived_quantity_static = self._calculate_comparison_derived_quantity_static(singlelabelledprecursor)
        
    def _calculate_comparison_derived_quantity_static(self, precursor):
        return precursor.derived_reference_quantity_static + precursor.ratio_to_reference

    def _define_precursorname2staticquantity(self):
        self._reference_intensity_dataframe["median"] = self._reference_intensity_dataframe.median(axis = 1)
        self._precursorname2staticquantity = dict(zip(self._reference_intensity_dataframe.index, self._reference_intensity_dataframe["median"]))
    
