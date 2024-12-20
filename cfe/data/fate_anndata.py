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

    def __init__(
            self,
            name: str = "FateAnnData",
            *args,
            **kwargs
    ):
        # logger.debug("FateAnnData __init__")
        self.id = random_time_string(name)
        super().__init__(*args, **kwargs)

        self.cfe_dict = {}
        self.uns["cfe"] = self.cfe_dict
        self.prior_information = {}

        self.is_wrapped_with_trajectory = False
        self.is_wrapped_with_waypoints = False

    @classmethod
    def from_anndata(cls, adata: ad.AnnData):
        """
        Create a FateAnnData object from an existing AnnData object.
        """
        logger.debug("Creating FateAnnData from existing AnnData")

        fate_adata = cls(name=adata.name if hasattr(adata, "name") else "FateAnnData",
                         X=adata.X,
                         obs=adata.obs,
                         var=adata.var,
                         uns=adata.uns,
                         obsm=adata.obsm,
                         varm=adata.varm,
                         layers=adata.layers)

        return fate_adata

    def add_prior_information(self, **kwargs):
        """
        ref: pydynverse/wrap/wrap_add_prior_information add_prior_information
        """
        self.prior_information.update(kwargs)

    def add_trajectory(
        self,
        milestone_network: pd.DataFrame,
        divergence_regions: pd.DataFrame = None,
        milestone_percentages: pd.DataFrame = None,
        progressions: pd.DataFrame = None,
    ):
        """
        create MilestoneWrapper object
        """
        logger.debug("FateAnnData add_trajectory")

        milestone_wrapper = MilestoneWrapper(
            milestone_network=milestone_network,
            divergence_regions=divergence_regions,
            milestone_percentages=milestone_percentages,
            progressions=progressions
        )
        milestone_wrapper.pipeline()
        self.milestone_wrapper = milestone_wrapper
        self.cfe_dict["milestone_wrapper"] = milestone_wrapper
        self.is_wrapped_with_trajectory = True

    def add_waypoints(self, milestone_wrapper: MilestoneWrapper):
        """
        create WaypointWrapper object
        """
        logger.debug("FateAnnData add_waypoints")

        waypoint_wrapper = WaypointWrapper(milestone_wrapper)
        waypoint_wrapper.pipeline()
        self.waypoint_wrapper = waypoint_wrapper
        self.cfe_dict["waypoint_wrapper"] = waypoint_wrapper
        self.is_wrapped_with_waypoints = True

    def add_branch_trajectory(
            self,
            branch_network: pd.DataFrame,
            branch_progressions: pd.DataFrame,
            branches: pd.DataFrame
    ):
        """
        ref: PyDynverse/pydynverse/wrap/wrap_add_branch_trajectory.add_branch_trajectory
        """
        logger.debug("FateAnnData add_branch_trajectory")

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
                mapper[node] = str(component_index+1)  # milestone id starts from 1
        milestone_network["from"] = milestone_network["from"].apply(lambda x: mapper[x])
        milestone_network["to"] = milestone_network["to"].apply(lambda x: mapper[x])
        milestone_network = pd.merge(milestone_network, branches, on="branch_id")

        progressions = pd.merge(branch_progressions, milestone_network, on="branch_id")[["cell_id", "from", "to", "percentage"]]

        milestone_network = milestone_network[["from", "to", "length", "directed"]]

        self.add_trajectory(milestone_network=milestone_network, progressions=progressions)

    def __getitem__(self, key):
        sub_adata = super().__getitem__(key)
        sub_fadata = self.from_anndata(sub_adata)
        # TODO: add sub operation for all other attributes, such as prior_information, milestone_wrapper, wayppoint_wrapper, etc.
        return sub_fadata
