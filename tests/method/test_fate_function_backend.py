import pytest
import cfe

import os.path
import scanpy as sc

function_name = "cf_paga"


class TestFunctionBackend:

    def setup_method(self):
        self.function_backend = cfe.method.FunctionBackend(function_name)

    def test_load_backend(self):
        assert self.function_backend.function_name == function_name

    def test_run(self):
        adata = sc.read(f"{os.path.dirname(__file__)}/../data/bifurcating.h5ad")
        fadata = cfe.data.FateAnnData.from_anndata(adata)
        fadata.layers["counts"] = fadata.X.copy()
        fadata.layers["expression"] = fadata.X.copy()
        cluster_key = "lineage"
        fadata.obs.index = fadata.obs["cell_id"]
        # prior_information,  parameters
        prior_information = {
            "start_id": "cell1",
            "groups_id": fadata.obs[cluster_key].tolist()
        }
        parameters = {"filter_features": False}
        fadata.add_prior_information(**prior_information)  # add prior information to fadata

        self.function_backend.run(fadata, parameters)

        assert fadata.is_wrapped_with_trajectory

    def test_load_definition(self):
        # _load_definition is called in test_load_backend, which is called in __init__
        definition = self.function_backend.definition
        assert isinstance(definition, cfe.method.Definition)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
