import pytest
import cfe

import os
from scipy.sparse import csc_matrix
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from ..test_util import compare_dataframes, compare_dataframes_closely

# test case for add_trajectory, add_waypoints


def setup_method_data():
    counts = np.array([
        [0, 10],
        [8, 10],
        [12, 12],
        [20, 20],
        [15, 16],
        [22, 20],
    ])

    counts = csc_matrix(counts)

    fadata = cfe.data.FateAnnData(X=counts)
    fadata.obs.index = ["a", "b", "c", "d", "e", "f"]
    fadata.obs["clusters"] = [1, 1, 2, 2, 2, 3]
    fadata.var.index = ["g1", "g2"]
    fadata.layers["counts"] = counts
    fadata.layers["expression"] = counts.copy()
    fadata.obsm["X_emb"] = counts.toarray().copy()

    return fadata


class TestFateAnnData:
    def setup_method(self):
        self.fadata = setup_method_data()

    def test_init(self):
        assert isinstance(self.fadata, ad.AnnData)
        assert self.fadata.shape == (6, 2)
        assert "cfe" in self.fadata.uns.keys()

    def test_from_anndata(self):
        # data source: https://github.com/theislab/cellrank_reproducibility/blob/master/data/dyngen_simulated_data/bifurcating.h5ad
        adata = sc.read_h5ad(f"{os.path.dirname(__file__)}/bifurcating.h5ad")
        fadata = cfe.data.FateAnnData.from_anndata(adata)
        assert fadata.id is not None

    @pytest.mark.skipif(not cfe.settings.r_available, reason="R is not available")
    def test_read_dynverse_simulation_data(self):
        fadata = cfe.data.FateAnnData.read_dynverse_simulation_data()
        assert fadata.is_wrapped_with_trajectory

    def test_add_model_name(self):
        # test in test_read_dynverse_simulation_data
        pass

    def test_get_all_model_name(self):
        self.test_add_trajectory()
        fadata = self.fadata
        # first model
        fadata.add_model_name("first model")
        # second model
        milestone_wrapper = fadata.milestone_wrapper
        fadata.add_trajectory(
            milestone_network=milestone_wrapper.milestone_network,
            divergence_regions=milestone_wrapper.divergence_regions,
            milestone_percentages=milestone_wrapper.milestone_percentages
        )
        from cfe.util import random_time_string
        fadata.add_model_name(random_time_string("second model")) # radom_time_string for parsing

        model_name_list = self.fadata.get_all_model_name()
        assert sorted(model_name_list) == sorted(["first model", "second model"])

    def test_get_item(self):
        pass

    def test_add_prior_information(self):
        self.fadata.add_prior_information(start_id="a", group_id=self.fadata.obs["clusters"].tolist())
        self.fadata.add_prior_information(end_id="f")
        assert set(["start_id", "group_id", "end_id"]) <= set(self.fadata.prior_information.keys())

    def test_add_trajectory(self):
        from .test_fate_milestone_wrapper import setup_method_data
        milestone_wrapper = setup_method_data()
        self.fadata.add_trajectory(
            milestone_network=milestone_wrapper.milestone_network,
            divergence_regions=milestone_wrapper.divergence_regions,
            milestone_percentages=milestone_wrapper.milestone_percentages,
            # progressions=milestone_wrapper.progressions
        )
        assert self.fadata.is_wrapped_with_trajectory
        # test write_h5ad
        self.fadata.write_h5ad("test_fate_anndata.h5ad")
        fadata = cfe.data.read_h5ad("test_fate_anndata.h5ad")
        assert fadata.milestone_wrapper is not None

    def test_add_waypoints(self):
        # from .test_fate_milestone_wrapper import setup_method_data
        # milestone_wrapper = setup_method_data()
        self.test_add_trajectory()
        self.fadata.add_waypoints()
        assert self.fadata.is_wrapped_with_waypoints
        # test write_h5ad
        self.fadata.write_h5ad("test_fate_anndata.h5ad")
        fadata = cfe.data.read_h5ad("test_fate_anndata.h5ad")
        assert fadata.waypoint_wrapper is not None

    def test_add_trajectory_branch(self):
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

        self.fadata.add_trajectory_branch(branch_network, branch_progressions, branches)

        # expected results
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

    def test_add_trajectory_linear(self):
        # new test case: pseudotime and FateAnnData
        name = "test_add_trajectory_linear"
        cell_ids = ["a", "b", "c", "d", "e", "f"]
        pseudotime = [0.0, 0.1, 0.4, 0.5, 0.8, 1.0]

        expression = np.tile(pseudotime, (2, 1)).T
        fadata = cfe.data.FateAnnData(X=expression, name=name)
        fadata.obs.index = cell_ids
        fadata.layers["expression"] = expression.copy()

        fadata.add_trajectory_linear(pseudotime)

        expected_milestone_ids = ["milestone_begin", "milestone_end"]
        expected_milestone_network = pd.DataFrame({
            "from": "milestone_begin",
            "to": "milestone_end",
            "length": 1,
            "directed": False,
        }, index=[0])
        expected_progressions = pd.DataFrame({
            "cell_id": cell_ids,
            "from": "milestone_begin",
            "to": "milestone_end",
            "percentage": pseudotime,
        })

        # 构造的milestone_network和progressions与预期对比
        assert fadata.milestone_wrapper["id_list"] == expected_milestone_ids
        assert fadata.milestone_wrapper["milestone_network"].equals(expected_milestone_network)
        assert fadata.milestone_wrapper["progressions"].equals(expected_progressions)

    def test_add_trajectory_velocity(self):
        # TODO
        pass

    def test_group_onto_trajectory_edges(self):
        self.test_add_trajectory()  # reuse test case from test_add_trajectory
        fadata = self.fadata

        cluster_key = "group"
        fadata.group_onto_trajectory_edges(cluster_key=cluster_key)
        excepted_group = ["W->W", "W->X", "X->Z", "Z->Z", "X->Z", "Z->A"]

        cluster_key = "group"
        assert cluster_key in self.fadata.obs.columns
        assert excepted_group == self.fadata.obs[cluster_key].tolist()

    def test_group_onto_nearest_milestones(self):
        self.test_add_trajectory()  # reuse test case from test_add_trajectory
        fadata = self.fadata

        cluster_key = "group"
        fadata.group_onto_nearest_milestones(cluster_key=cluster_key)
        excepted_group = ["W", "X", "X", "Z", "Z", "Z"]

        assert cluster_key in self.fadata.obs.columns
        assert excepted_group == self.fadata.obs[cluster_key].tolist()

    def test_simplify_trajectory(self):
        # create FateAnnData object
        id = "directed_linear"
        cell_ids = ["a", "b", "c", "d", "e"]
        milestone_network = pd.DataFrame(
            data=[
                ["A", "B", 1, True],
                ["B", "C", 1, True],
                ["C", "D", 1, True]
            ],
            columns=["from", "to", "length", "directed"],
        )
        progressions = pd.DataFrame(
            data=[
                ["a", "A", "B", 0.3],
                ["b", "A", "B", 0.6],
                ["c", "B", "C", 0.2],
                ["d", "B", "C", 0.8],
                ["e", "C", "D", 0.4],
            ],
            columns=["cell_id", "from", "to", "percentage"]
        )
        fadata = cfe.data.FateAnnData(name=id, X=np.zeros((len(cell_ids), 2)))
        fadata.add_trajectory(
            milestone_network=milestone_network,
            progressions=progressions
        )

        # simpify trajectory
        simplified_milestone_wrapper = fadata.simplify_trajectory()

        # assert results
        expected_milestone_network = pd.DataFrame(
            data=[["A", "D", 3, True]],
            columns=["from", "to", "length", "directed"],
        )
        expected_progressions = pd.DataFrame(
            data=[
                ["a", "A", "D", 0.1],
                ["b", "A", "D", 0.2],
                ["c", "A", "D", 0.4],
                ["d", "A", "D", 0.6],
                ["e", "A", "D", 0.8],
            ],
            columns=["cell_id", "from", "to", "percentage"]
        )

        assert simplified_milestone_wrapper.milestone_network.equals(expected_milestone_network)
        assert compare_dataframes_closely(simplified_milestone_wrapper.progressions, expected_progressions, on_columns="percentage")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
