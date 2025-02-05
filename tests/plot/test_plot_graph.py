import pytest
import cfe

import os
import matplotlib.pyplot as plt

from .test_plot_trajectory import get_plot_fadata

def test_plot_graph():
    fadata = get_plot_fadata()
    cfe.plot.plot_graph(fadata)
    plt.savefig(f"{os.path.dirname(__file__)}/img/test_plot_graph.png")


if __name__ == "__main__":
    pytest.main(["-v", __file__])
