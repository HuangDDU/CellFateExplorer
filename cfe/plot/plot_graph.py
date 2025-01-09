import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import networkx as nx
import scanpy as sc
import seaborn as sns

from ..data import FateAnnData


def plot_graph(
    fadata: FateAnnData,
    save: str = None
):
    """Plot DAG base on milestone network

    Args:
        fadata (FateAnnData): _description_
        save (str, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """

    # extract milestone network
    milestone_wrapper = fadata.milestone_wrapper
    milestone_id_list = milestone_wrapper.id_list
    milestone_network = milestone_wrapper.milestone_network
    milestone_percentages = milestone_wrapper.milestone_percentages
    is_directed = milestone_network["directed"].any()

    # color and position fo milestone
    milestone_color_list = sns.color_palette("Set3")[:len(milestone_id_list)]  # color rgb
    milestone_color_dict = dict(zip(milestone_id_list, milestone_color_list))
    G = nx.from_pandas_edgelist(
        milestone_network,
        source="from",
        target="to",
        edge_attr=True,
        create_using=nx.DiGraph if is_directed else nx.Graph
    )
    milestone_emb_dict = nx.nx_agraph.graphviz_layout(G, prog="dot")  # position
    ax = plt.subplots(1, 1)[1]
    nx.draw(G,
            milestone_emb_dict,
            with_labels=True,
            node_color=[milestone_color_dict[node] for node in G.nodes],
            ax=ax
            )

    # color and position fo cell
    milestone_color_df = pd.DataFrame(milestone_color_dict, index=["r", "g", "b"]).T

    def mix_color(mpg):
        mpg_color = milestone_color_df.loc[mpg["milestone_id"]]
        mix_color_array = mpg_color.apply(lambda rgb_channel: (rgb_channel.array * mpg["percentage"].array).sum())
        return mix_color_array
    cell_color_df = milestone_percentages.groupby("cell_id").apply(lambda mpg: mix_color(mpg))

    milestone_emb_df = pd.DataFrame(milestone_emb_dict).T

    def mix_emb(mpg):
        mpg_emb = milestone_emb_df.loc[mpg["milestone_id"]]
        return mpg_emb.apply(lambda emb_dim: (emb_dim.array * mpg["percentage"].array)).sum()
    cell_emb_df = milestone_percentages.groupby("cell_id").apply(lambda mpg: mix_emb(mpg))

    basis = "milestone_network_emb"
    cell_color_key = "cell_color"
    fadata.obsm[basis] = cell_emb_df.loc[fadata.obs.index].values
    fadata.obs[cell_color_key] = pd.Categorical(fadata.obs.index, categories=fadata.obs.index.tolist())
    fadata.uns[f"{cell_color_key}_colors"] = rgb2hex(cell_color_df.loc[fadata.obs.index].values)

    sc.pl.embedding(fadata, basis=basis, color=cell_color_key, ax=ax, title="", legend_loc=None, save=save)


def rgb2hex(palette):
    return [mcolors.to_hex(color) for color in palette]


def hex2rgb(palette):
    return [mcolors.to_rgb(color) for color in palette]
