import pytest
import cfe

import os
import matplotlib.pyplot as plt


def test_plot_graph():
    fadata = cfe.data.FateAnnData.read_dynverse_simulation_data()
    cfe.plot.plot_graph(fadata)
    plt.savefig(f"{os.path.dirname(__file__)}/img/test_plot_graph.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
