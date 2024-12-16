import pytest
import cfe

import pandas as pd
import anndata as ad

from ..test_util import compare_dataframes


class TestFateAnnData:
    def setup_method(self):
        self.fadata = cfe.data.FateAnnData()

    def test_parent_class(self):
        assert isinstance(self.fadata, ad.AnnData)

    def test_add_prior_information(self):
        self.fadata.add_prior_information()

    def test_add_trajectory(self):
        from .test_fate_milestone_wrapper import setup_method_data
        milestone_wrapper = setup_method_data()
        self.fadata.add_trajectory(milestone_wrapper)
        assert self.fadata.is_wrapped_with_trajectory

    def test_add_waypoints(self):
        from .test_fate_milestone_wrapper import setup_method_data
        milestone_wrapper = setup_method_data()
        self.fadata.add_waypoints(milestone_wrapper)
        assert self.fadata.is_wrapped_with_waypoints

    def test_add_branch_trajectory(self):
        branch_network = pd.DataFrame(
            columns=["from", "to"],
            data=[
                ["A", "B"],
                ["A", "C"],
                ["B", "D"],
            ],
        )
        branch_progressions = pd.DataFrame(
            columns=["cell_id", "branch_id", "percentage"],
            data=[
                ["a", "A", 0.0],
                ["b", "A", 0.8],
                ["c", "B", 0.2],
                ["d", "B", 1.0],
                # compared to "test_add_trajectory" test case, cell "e" is moved to branch "C" from divergence region
                ["e", "C", 0.2],
                ["f", "D", 0.2],
            ]
        )
        branches = pd.DataFrame(
            columns=["branch_id", "length", "directed"],
            data=[
                ["A", 1.0, True],
                ["B", 1.0, True],
                ["C", 1.0, True],
                ["D", 2.0, True],
            ]
        )

        self.fadata.add_branch_trajectory(branch_network, branch_progressions, branches)

        # 预期构造结果
        expected_milestone_network = pd.DataFrame(
            columns=["from", "to", "length", "directed"],
            data=[
                ["1", "2", 1.0, True],
                ["2", "3", 1.0, True],
                ["2", "4", 1.0, True],
                ["3", "5", 2.0, True],
            ]
        )
        expected_progressions = pd.DataFrame(
            columns=["cell_id", "from", "to", "percentage"],
            data=[
                ["a", "1", "2", 0.0],
                ["b", "1", "2", 0.8],
                ["c", "2", "3", 0.2],
                ["d", "2", "3", 1.0],
                ["e", "2", "4", 0.2],
                ["f", "3", "5", 0.2],
            ]
        )

        assert compare_dataframes(self.fadata.cfe_dict["milestone_wrapper"]["milestone_network"], expected_milestone_network, on_columns=["from", "to"])
        assert compare_dataframes(self.fadata.cfe_dict["milestone_wrapper"]["progressions"], expected_progressions, on_columns=["cell_id", "from", "to"])


if __name__ == "__main__":
    pytest.main(["-v", __file__])
