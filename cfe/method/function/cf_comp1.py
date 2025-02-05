import pickle
import anndata as ad
import scanpy as sc


def cf_comp1(
    adata: ad.AnnData,
    prior_information: dict = {},
    parameters: dict = {}
):
    sc.pp.pca(adata, n_comps=parameters["ndim"])

    pseudotime = adata.obsm["X_pca"][:, parameters["component"] - 1]
    trajectory_dict = {
        "pseudotime": pseudotime,
        # "trajectory_type": "linear",
    }

    return trajectory_dict


if __name__ == "__main__":

    from parse_args import parse_args

    adata, prior_information, parameters, output_filename = parse_args()

    trajectory_dict = cf_comp1(adata, prior_information, parameters)

    with open(output_filename, "wb") as f:
        pickle.dump(trajectory_dict, f)
    print("Comp1 Finish!")
