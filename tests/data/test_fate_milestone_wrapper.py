import pytest
import cfe

import pandas as pd
from ..test_util import compare_dataframes


def setup_method_data():
    """ Create data for testing, convinient for other test file reuse """
    # id_list = ["W", "X", "Y", "Z", "A"]
    milestone_network = pd.DataFrame(
        columns=["from", "to", "length", "directed"],
        data=[
            ["W", "X", 1.0, True],
            ["X", "Y", 1.0, True],
            ["X", "Z", 1.0, True],
            ["Z", "A", 2.0, True]
        ]
    )
    divergence_regions = pd.DataFrame(
        columns=["divergence_id", "milestone_id", "is_start"],
        data=[
            ["XYZ", "X", True],
            ["XYZ", "Y", False],
            ["XYZ", "Z", False]
        ]
    )
    milestone_percentages = pd.DataFrame(
        columns=["cell_id", "milestone_id", "percentage"],
        data=[
            ["a", "W", 1.0],
            ["b", "W", 0.2],
            ["b", "X", 0.8],
            ["c", "X", 0.8],
            ["c", "Z", 0.2],
            ["d", "Z", 1.0],
            ["e", "X", 0.3],
            ["e", "Y", 0.2],
            ["e", "Z", 0.5],
            ["f", "Z", 0.8],
            ["f", "A", 0.2],
        ]
    )
    milestone_wrapper = cfe.data.MilestoneWrapper(
        milestone_network=milestone_network,
        divergence_regions=divergence_regions,
        milestone_percentages=milestone_percentages,
    )
    return milestone_wrapper


class TestMilestoneWrapper:
    def setup_method(self):
        self.milestone_wrapper = setup_method_data()

    def test_magic_method(self):
        """test __***__ methods"""
        mw = self.milestone_wrapper

        # test __contains__
        assert "id" in mw, "id should in mw"

        # test __getitem__
        assert mw["id"] == mw.id, "mw['id'] should be the same as mw.id"

        # test keys
        mw_dict = dict(mw)
        attribute_name_list = ["id"]
        assert set(attribute_name_list).issubset(set(mw_dict.keys())), f"{attribute_name_list} should be the keys of the dict: {mw_dict}"

    def test_milestone_network(self):
        mw = self.milestone_wrapper
        id_list = mw.id_list
        milestone_network = mw.milestone_network
        assert (set(milestone_network["from"].unique()) | set(milestone_network["to"].unique())) == set(id_list), \
            "every id should show in 'from' or 'to' column in  milestone_network dataframe"

    def test_convert_milestone_percentages_to_progressions(self):
        mw = self.milestone_wrapper
        # static method can be called with class or instance
        progression = mw.convert_milestone_percentages_to_progressions(mw.milestone_network, mw.milestone_percentages)

        expected_progression = pd.DataFrame(
            columns=["cell_id", "from", "to", "percentage"],
            data=[
                ["a", "W", "W", 1],
                ["b", "W", "X", 0.8],
                ["c", "X", "Z", 0.2],
                ["d", "Z", "Z", 1],
                ["e", "X", "Y", 0.2],
                ["e", "X", "Z", 0.5],
                ["f", "Z", "A", 0.2],
            ]
        )
        assert isinstance(progression, pd.DataFrame), "progression should be a dataframe"
        assert compare_dataframes(progression, expected_progression, on_columns=["cell_id", "from", "to"])

    # this test case can't execute the method
    # def test_convert_progressions_to_milestone_percentages(self):
    #     pass

    def test_pipeline(self):
        mw = self.milestone_wrapper

        mw.pipeline()

        assert mw.progressions is not None
        assert mw.milestone_percentages is not None


if __name__ == "__main__":
    pytest.main(["-v", __file__])
