import pytest
import cfe


import pandas as pd

import cfe.metric


def test_metric_isomorphic():
    net1 = pd.DataFrame(
        data=[
            ["A", "B", 1, True,],
            ["B", "C", 1, True,],
            ["C", "D", 1, True,],
        ],
        columns=["from", "to", "length", "direction"]
    )
    net2 = pd.DataFrame(
        data=[
            ["A", "B", 1, True,],
            ["B", "C", 1, True,],
            ["B", "D", 1, True,],
        ],
        columns=["from", "to", "length", "direction"]
    )
    assert cfe.metric.calc_isomorphic(net1, net1) == 1
    assert cfe.metric.calc_isomorphic(net1, net2) == 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])
