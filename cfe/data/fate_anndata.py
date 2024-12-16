import networkx as nx
import pandas as pd
import anndata as ad
from .._logging import logger
from ..util import random_time_string

from .fate_milestone_wrapper import MilestoneWrapper
from .fate_waypoint_wrapper import WaypointWrapper


class FateAnnData(ad.AnnData):
    """
    AnnData object for CellFateExplorer, related data are stored in the object.uns[cfe] attribute.
    """

    def __init__(self, name="FateAnnData", *args, **kwargs):
        logger.debug("FateAnnData __init__")
        self.id = random_time_string(name)
        super().__init__()

        # adata.uns["cfe"] is synchronized with cfe_dict
        self.cfe_dict = {}
        self.uns["cfe"] = self.cfe_dict
        # self.cfe_dict["trajectory"] = {}
        # self.cfe_dict["waypoint"] = {}

        self.is_wrapped_with_trajectory = False
        self.is_wrapped_with_waypoints = False

    def add_prior_information(self):
        pass

    def add_trajectory(self, milestone_wrapper: MilestoneWrapper):
        logger.debug("FateAnnData add_trajectory")
        # trajectory = self.cfe_dict.get("milestone_wrapper")

        # if not trajectory is None:
        #     return

        milestone_wrapper.pipeline()
        self.cfe_dict["milestone_wrapper"] = milestone_wrapper
        self.is_wrapped_with_trajectory = True

    def add_embedding(self):
        pass

    def add_waypoints(self, waypoint_wrapper: WaypointWrapper):
        logger.debug("FateAnnData add_trajectory")
        # waypoint = self.cfe_dict.get("waypoint_wrapper")

        # if not waypoint is None:
        #     return

        waypoint_wrapper.pipeline()
        self.cfe_dict["waypoint_wrapper"] = waypoint_wrapper
        self.is_wrapped_with_waypoints = True

    def add_branch_trajectory(
            self,
            branch_network: pd.DataFrame,
            branch_progressions: pd.DataFrame,
            branches: pd.DataFrame
    ):
        cell_id_list = branch_progressions["cell_id"].unique().tolist()
        branch_id_list = branches["branch_id"]

        milestone_network = pd.DataFrame({
            "from": map(lambda x: f"{x}_from", branch_id_list),
            "to": map(lambda x: f"{x}_to", branch_id_list),
            "branch_id": branch_id_list
        })
        milestone_mapper_network = pd.concat(
            [
            # single from node
            pd.DataFrame({
                "from": map(lambda x: f"{x}_from", branch_id_list),
                "to": map(lambda x: f"{x}_from", branch_id_list),
            }),
            # connected node, if "A->B" in branch_network , then "A_to->B_from" in here,
            pd.DataFrame({
                "from": map(lambda x: f"{x}_to", branch_network["from"]),
                "to": map(lambda x: f"{x}_from", branch_network["to"]),
            }),
            # single to node
            pd.DataFrame({
                "from": map(lambda x: f"{x}_to", branch_id_list),
                "to": map(lambda x: f"{x}_to", branch_id_list),
            }),
        ])
        # transform node name to connected component id
        mapper = {}
        graph = nx.from_pandas_edgelist(milestone_mapper_network, source="from", target="to")
        connected_components = nx.connected_components(graph)
        for component_index, component in enumerate(connected_components):
            for node in component:
                mapper[node] = str(component_index+1)  # milestone序号从1开始
        milestone_network["from"] = milestone_network["from"].apply(lambda x: mapper[x])
        milestone_network["to"] = milestone_network["to"].apply(lambda x: mapper[x])
        milestone_network = pd.merge(milestone_network, branches, on="branch_id")

        
        progressions = pd.merge(branch_progressions, milestone_network, on="branch_id")[["cell_id", "from", "to", "percentage"]]
        
        milestone_network = milestone_network[["from", "to", "length", "directed"]]

        milestone_wrapper = MilestoneWrapper(milestone_network=milestone_network, progressions=progressions)
        self.add_trajectory(milestone_wrapper)
