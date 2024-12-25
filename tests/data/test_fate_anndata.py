import pytest
import cfe

import os
from scipy.sparse import csc_matrix
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from ..test_util import compare_dataframes

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
        from .test_fate_milestone_wrapper import setup_method_data
        milestone_wrapper = setup_method_data()
        self.fadata.add_waypoints(milestone_wrapper)
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
        # TODO: linear
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


if __name__ == "__main__":
    pytest.main(["-v", __file__])
