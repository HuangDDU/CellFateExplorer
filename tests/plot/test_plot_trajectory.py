import pytest
import cfe


def test_plot_trajectory():
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

    # plot
    cfe.plot.plot_trajectory(fadata, color="clusters", basis="umap", save="test_plot_trajectory.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
