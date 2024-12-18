import pytest
import cfe

import os
import scanpy as sc


class TestCFPAGA():
    def setup_method(self):
        adata = sc.read_h5ad(f"{os.path.dirname(__file__)}/../../data/bifurcating.h5ad")
        self.fadata = cfe.data.FateAnnData.from_anndata(adata)
        self.fadata.obs.index = self.fadata.obs["cell_id"]

    def test_paga(self):
        # add priority and parameeters
        prior_information = {
            "start_id": "cell1",
            "groups_id": self.fadata.obs["lineage"].tolist()
        }
        parameters = {"connectivity_cutoff": 0.8, "filter_features": False}
        self.fadata.add_prior_information(**prior_information)  # add prior information to fadata
        fadata = cfe.method.cf_paga(self.fadata, parameters)  # add parameters when inferring trajectory
        assert fadata.is_wrapped_with_trajectory


if __name__ == "__main__":
    pytest.main(["-v", __file__])
