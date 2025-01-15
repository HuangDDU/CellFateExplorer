import networkx as nx


def calc_isomorphic(net1, net2):
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
