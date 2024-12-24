import networkx as nx
import pandas as pd
import anndata as ad
import scanpy as sc

from .._logging import logger
from ..util import random_time_string

from .fate_milestone_wrapper import MilestoneWrapper
from .fate_waypoint_wrapper import WaypointWrapper


class FateAnnData(ad.AnnData):
    """AnnData object for CellFateExplorer, related data are stored in the object.uns["cfe"] attribute.
    """

    def __init__(
            self,
            name: str = "FateAnnData",
            *args,
            **kwargs
    ):
        """Initialize the FateAnnData class.

        Args:
            name (str, optional): name of the FateAnnData object.
        """
        # logger.debug("FateAnnData __init__")
        self.id = random_time_string(name)
        super().__init__(*args, **kwargs)

        cfe_dict = self.uns.get("cfe", {})  # try to get the stored FateAnnData information

        self.prior_information = cfe_dict.get("prior_information", {})
        cfe_dict["prior_information"] = self.prior_information

        self.milestone_wrapper = cfe_dict.get("milestone_wrapper", None)
        self.waypoint_wrapper = cfe_dict.get("waypoint_wrapper", None)

        # NOTE: Other attributes will be added later.
        self.is_wrapped_with_trajectory = False
        self.is_wrapped_with_waypoints = False

        self.cfe_dict = cfe_dict
        self.uns["cfe"] = self.cfe_dict

    @classmethod
    def from_anndata(cls, adata: ad.AnnData) -> "FateAnnData":
        """Create a FateAnnData object from an existing AnnData object.

        Args:
            adata (ad.AnnData): existing AnnData object

        Returns:
            fadata (cfe.data.FateAnnData): generated FateAnnData object
        """

        logger.debug("Create a FateAnnData object from an existing AnnData object.")

        fadata = cls(name=adata.name if hasattr(adata, "name") else "FateAnnData",
                     X=adata.X,
                     obs=adata.obs,
                     var=adata.var,
                     uns=adata.uns,
                     obsm=adata.obsm,
                     varm=adata.varm,
                     layers=adata.layers)

        return fadata

    def add_prior_information(self, **kwargs) -> None:
        """Add prior information to the FateAnnData object.

        ref: pydynverse/wrap/wrap_add_prior_information add_prior_information
        """
        self.prior_information.update(kwargs)

    def add_trajectory(
        self,
        milestone_network: pd.DataFrame,
        divergence_regions: pd.DataFrame = None,
        milestone_percentages: pd.DataFrame = None,
        progressions: pd.DataFrame = None,
    ) -> None:
        """Create MilestoneWrapper object as trajectory

        Args:
            milestone_network (pd.DataFrame): milestone network with column list: ["from", "to", "length", "directed"]
            divergence_regions (pd.DataFrame, optional): divergence regions with column list: ["divergence_id", "milestone_id", "is_start"].
            milestone_percentages (pd.DataFrame, optional): milestone percentage with column list: ["cell_id", "milestone_id", "percentage"].
            progressions (pd.DataFrame, optional): progressions with column list: ["cell_id", "from", "to", "percentage"].
        """

        logger.debug("FateAnnData add_trajectory")

        milestone_wrapper = MilestoneWrapper(
            milestone_network=milestone_network,
            divergence_regions=divergence_regions,
            milestone_percentages=milestone_percentages,
            progressions=progressions
        )
        self.milestone_wrapper = milestone_wrapper
        self.cfe_dict["milestone_wrapper"] = milestone_wrapper
        self.is_wrapped_with_trajectory = True

    def add_trajectory_by_type(self, trajectory_dict: dict) -> None:
        """Call the trajectory addition method based on specific trajectory types

        Args:
            trajectory_dict (dict): trajectory dict result based on specific trajectory types
        """

        # TODO : 暂时只是paga
        trajectory_type = "paga"
        self.add_trajectory_branch(
            branch_network=trajectory_dict["branch_network"],
            branches=trajectory_dict["branches"],
            branch_progressions=trajectory_dict["branch_progressions"]
        )

    def add_waypoints(self, milestone_wrapper: MilestoneWrapper) -> None:
        """Create WaypointWrapper object

        Args:
            milestone_wrapper (MilestoneWrapper): trajectory wrapper object
        """
        logger.debug("FateAnnData add_waypoints")

        waypoint_wrapper = WaypointWrapper(milestone_wrapper)
        self.waypoint_wrapper = waypoint_wrapper
        self.cfe_dict["waypoint_wrapper"] = waypoint_wrapper
        self.is_wrapped_with_waypoints = True

    def add_trajectory_branch(
            self,
            branch_network: pd.DataFrame,
            branch_progressions: pd.DataFrame,
            branches: pd.DataFrame
    ) -> None:
        """Add branch trajectory,such as paga

        ref: PyDynverse/pydynverse/wrap/wrap_add_branch_trajectory.add_branch_trajectory

        Args:
            branch_network (pd.DataFrame): branch network with column list: ["from", "to"]
            branch_progressions (pd.DataFrame): branch progressions with column list: ["cell_id", "branch_id", "percentage"
            branches (pd.DataFrame): branches with column list: ["branch_id", "length", "directed"]
        """
        logger.debug("FateAnnData add_trajectory_branch")

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
                mapper[node] = str(component_index + 1)  # milestone id starts from 1
        milestone_network["from"] = milestone_network["from"].apply(lambda x: mapper[x])
        milestone_network["to"] = milestone_network["to"].apply(lambda x: mapper[x])
        milestone_network = pd.merge(milestone_network, branches, on="branch_id")

        progressions = pd.merge(branch_progressions, milestone_network, on="branch_id")[["cell_id", "from", "to", "percentage"]]

        milestone_network = milestone_network[["from", "to", "length", "directed"]]

        self.add_trajectory(milestone_network=milestone_network, progressions=progressions)

    def write_h5ad(self, *args, **kwargs):
        if self.cfe_dict.get("milestone_wrapper", None) is not None:
            self.cfe_dict["milestone_wrapper"] = dict(self.cfe_dict["milestone_wrapper"])
        if self.cfe_dict.get("waypoint_wrapper", None) is not None:
            self.cfe_dict["waypoint_wrapper"] = dict(self.cfe_dict["waypoint_wrapper"])
            self.cfe_dict["waypoint_wrapper"]["milestone_wrapper"] = None  # milestone_wrapper is redundent
            waypoints = self.cfe_dict["waypoint_wrapper"]["waypoints"]
            self.cfe_dict["waypoint_wrapper"]["waypoints"] = waypoints.fillna("")  # "" replace None
        return super().write_h5ad(*args, **kwargs)

    def __getitem__(self, key):
        sub_adata = super().__getitem__(key)
        sub_fadata = self.from_anndata(sub_adata)
        # TODO: add sub operation for all other attributes, such as prior_information, milestone_wrapper, wayppoint_wrapper, etc.
        return sub_fadata


def read_h5ad(*args, **kwargs):
    adata = sc.read_h5ad(*args, **kwargs)
    fadata = FateAnnData.from_anndata(adata)
    return fadata
