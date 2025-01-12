import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import networkx as nx
import scanpy as sc
import seaborn as sns

from ..data import FateAnnData
from .add_color import add_milestone_color, add_milestone_cell_color


def plot_graph(
    fadata: FateAnnData,
    color: str | list = "milestone",
    save: str = None
):
    """Plot DAG base on milestone network

    Args:
        fadata (FateAnnData): FateAnnData
        save (str, optional): img path.

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
    milestone_color_list = add_milestone_color(len(milestone_id_list))  # color rgb
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
            width=5,  # TODO: adjusted by cell size
            edge_color="gray",
            arrowstyle="simple",
            arrowsize=30,   # TODO: adjusted by cell size
            ax=ax,
            )

    # position fo cell
    milestone_emb_df = pd.DataFrame(milestone_emb_dict).T

    def mix_emb(mpg):
        # mix related milestone emb to get position for a cell
        mpg_emb = milestone_emb_df.loc[mpg["milestone_id"]]
        return mpg_emb.apply(lambda emb_dim: (emb_dim.array * mpg["percentage"].array)).sum()
    basis = "milestone_network_emb"
    cell_emb_df = milestone_percentages.groupby("cell_id").apply(lambda mpg: mix_emb(mpg))
    fadata.obsm[basis] = cell_emb_df.loc[fadata.obs.index].values

    if color == "milestone":
        # color of cells
        cell_color_key = "milestone"
        cell_color_df = add_milestone_cell_color(milestone_color_dict, milestone_percentages)
        fadata.obs[cell_color_key] = pd.Categorical(fadata.obs.index, categories=fadata.obs.index.tolist())
        fadata.uns[f"{cell_color_key}_colors"] = cell_color_df.loc[fadata.obs.index].values

    sc.pl.embedding(
        fadata,
        basis=basis,
        color=color,
        ax=ax,
        title="",
        legend_loc=None,
        save=save
    )
