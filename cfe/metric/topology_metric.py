import pandas as pd
import networkx as nx


def calc_isomorphic(net1: pd.DataFrame, net2: pd.DataFrame):
    """Judge if two milestone network are  isomorphic

    Args:
        net1 (pd.DataFrame): reference milestone network
        net2 (pd.DataFrame): predict milestone network

    Returns:
        int: 0 is not isomorphic, 1 is isomorphic
    """
    # 图同构
    graph1 = nx.from_pandas_edgelist(net1, source="from", target="to")
    graph2 = nx.from_pandas_edgelist(net2, source="from", target="to")
    if nx.is_isomorphic(graph1, graph2):
        return 1
    else:
        return 0


def calculate_edge_flip(net1, net2):
    # TODO: 边反转，最大公共边子图
    return 0


def calculate_him(net1, net2):
    # TODO: HIM
    return 0
