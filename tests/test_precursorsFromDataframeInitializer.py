
import refquant.refquant_classes as rqc
import unittest
import pandas as pd

class TestDataframeInit(unittest.TestCase):


    def test_PrecursorsFromDataframeInititalizer(self):
        test_file =  "Spectronaut/test_df_w_one_precursor.tsv"
        test_df  = pd.read_csv(test_file, sep="\t")
        precursor_name = test_df["precursor"][0]
        test_df = test_df.set_index("precursor")
        prec_df_init = rqc.PrecursorsFromDataframeInititalizer(precursor_name, test_df)
        single_labelled_precursors = prec_df_init.single_labelled_precursors
        self.assert_that_all_expected_precursors_are_loaded( precursor_name,single_labelled_precursors)
        self.assert_that_quant_values_are_loaded_correctly(precursor_name, single_labelled_precursors)

    def assert_that_all_expected_precursors_are_loaded(self, precursor_name, single_labelled_precursors):
        assert len(single_labelled_precursors) == 4
        channel_names = ["reference4", "reference8", "target4", "target8"]
        channel_names_single_labelled_precursors = [single_labelled_precursor.channel_name for single_labelled_precursor in single_labelled_precursors]
        assert channel_names == channel_names_single_labelled_precursors
        assert single_labelled_precursors[0].name == precursor_name

    def assert_that_quant_values_are_loaded_correctly(self, precursor_name, single_labelled_precursors):
        ions = ["SEQ_NLFEEYR_MOD__[DimethNter0]NLFEEYR__CHARGE_2_FRGION_b6__noloss__1_","SEQ_NLFEEYR_MOD__[DimethNter0]NLFEEYR__CHARGE_2_MS1ISOTOPES_0_", "SEQ_NLFEEYR_MOD__[DimethNter0]NLFEEYR__CHARGE_2_MS1ISOTOPES_1_",
        "SEQ_NLFEEYR_MOD__[DimethNter0]NLFEEYR__CHARGE_2_MS1ISOTOPES_2_", "SEQ_NLFEEYR_MOD__[DimethNter0]NLFEEYR__CHARGE_2_MS1ISOTOPES_3_"]
        quantities_target4 = [9.188793312,12.70238916,13.41560974,13.20747196,13.03118426]
        quantities_target8 = [7.674903863,12.58871464,13.35493805,12.97727992,14.643067]
        quantities_reference4 = [7.663307381,13.91653263,13.6614441,12.36714169,13.15180936]
        quantities_reference8 = [7.341600259,13.91653263,13.6614441,12.36714169,12.67131991]
        list_of_quantities = [quantities_reference4, quantities_reference8, quantities_target4, quantities_target8]
        quant_dicts = self.compare_quant_dicts(single_labelled_precursors ,ions, list_of_quantities)

    def compare_quant_dicts(self,single_labelled_precursors , ions, list_of_quantities):
        for idx in range(len(single_labelled_precursors)):
            single_labelled_precursor = single_labelled_precursors[idx]
            quant_dict_reference = self.get_quant_dict(ions, list_of_quantities[idx])
            self.assert_that_dicts_are_close(quant_dict_reference, single_labelled_precursor.fragion2quantity)
            print("dicts are close")

    @staticmethod
    def get_quant_dict(ions, quantities):
        return dict(zip(ions, quantities))

    #check that two dicts are close enough
    @staticmethod
    def assert_that_dicts_are_close(dict1, dict2, tol=1e-6):
        for key in dict1:
            assert abs(dict1[key] - dict2[key]) < tol


class TestPrecursorLoading(unittest.TestCase):

    def test_PrecursorLoader(self):


        replictate_no = 3970
        input_file = "Spectronaut/test_df_w_many_precursors.tsv"
        prec_loader = rqc.PrecursorLoader(input_file, replictate_no)
        assert len(prec_loader.precursors_w_all_labels) == len(pd.read_csv(input_file, sep="\t")["precursor"].unique())
        assert len(prec_loader.precursors_w_all_labels[0].list_of_single_labelled_precursors) == 4
        assert prec_loader.precursors_w_all_labels[0].list_of_single_labelled_precursors[0].channel_name == "reference4"
        assert prec_loader.precursors_w_all_labels[0].list_of_single_labelled_precursors[0].name == "SEQ_ALALLEDEER_MOD__[DimethNter0]ALALLEDEER__CHARGE_2_"


    
