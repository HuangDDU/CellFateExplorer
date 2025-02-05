import pytest
import cfe

import os
import matplotlib.pyplot as plt


def get_plot_fadata():
    # add trajectory
    from ..data.test_fate_milestone_wrapper import setup_method_data as get_milestone_wrapper
    from ..data.test_fate_anndata import setup_method_data as get_fadata
    fadata = get_fadata()
    fadata.obsm["X_umap"] = fadata.obsm["X_emb"]
    milestone_wrapper = get_milestone_wrapper()
    fadata.add_trajectory(
        milestone_network=milestone_wrapper.milestone_network,
        divergence_regions=milestone_wrapper.divergence_regions,
        milestone_percentages=milestone_wrapper.milestone_percentages,
    )
    return fadata


def test_plot_trajectory():
    fadata = get_plot_fadata()
    cfe.plot.plot_trajectory(fadata, color="clusters", basis="umap")
    plt.savefig(f"{os.path.dirname(__file__)}/img/test_plot_trajectory.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
