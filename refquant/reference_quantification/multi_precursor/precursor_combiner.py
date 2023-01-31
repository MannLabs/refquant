from abc import ABC, abstractmethod
from  . import precursor_loader
from ..precursor import precursor_classes
from ..precursor import target_reference_annotator

class PrecursorCombiner(ABC):
    def __init__(self, precursorloader : precursor_loader.PrecursorLoader):
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
        precursor_with_matched_labels = precursor_classes.PrecursorWithMatchedLabels(reference_precursor=reference_precursor, list_of_target_precursors=list_of_target_precursors)
        target_reference_annotator.PrecursorWMatchedLabelsAnnotator(precursor_with_matched_labels).annotate_precursor()
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
        precursor_with_matched_labels = precursor_classes.PrecursorWithMatchedLabels(reference_precursor=reference_precursor, list_of_target_precursors=[target_precursor])
        target_reference_annotator.PrecursorWMatchedLabelsAnnotator(precursor_with_matched_labels).annotate_precursor()
        return precursor_with_matched_labels

class PrecursorCombinerToPRMReference(PrecursorCombinerToSingleReference):
    def _define_precursor_w_matched_labels_using_target_and_reference(self, precursor_w_all_labels):
        list_of_target_precursors = list(filter(lambda x : "target" in x.channel_name , precursor_w_all_labels.list_of_single_labelled_precursors))
        reference_precursor = list(filter(lambda x : "reference" in x.channel_name , precursor_w_all_labels.list_of_single_labelled_precursors))[0]
        precursor_with_matched_labels = precursor_classes.PrecursorWithMatchedLabels(reference_precursor=reference_precursor, list_of_target_precursors=list_of_target_precursors)
        target_reference_annotator.PrecursorWMatchedLabelsAnnotator(precursor_with_matched_labels).annotate_precursor_prm()
        return precursor_with_matched_labels

