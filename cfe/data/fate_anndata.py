import numpy as np
import pandas as pd
import networkx as nx
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

        # TODO: parse milestone_wrapper, waypoint_wrapper from dict, move to read_h5ad
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

        fadata = cls(
            name=adata.name if hasattr(adata, "name") else "FateAnnData",
            X=adata.X,
            obs=adata.obs,
            var=adata.var,
            uns=adata.uns,
            obsm=adata.obsm,
            varm=adata.varm,
            layers=adata.layers
        )

        return fadata

    @classmethod
    def read_dynverse_simulation_data(cls, data_filename="synthetic/dyntoy/bifurcating_1.rds", data_dir="/home/huang/RCode/scrna_tools/dynbenchmark/data"):
        # read dynverse simulation data and create FateAnnData object
        import rpy2.robjects as ro
        from ..util import rpy2_read  # rpy2 data structure transfer automatically

        r_script = f"""
        dataset <- readRDS("{data_dir}/{data_filename}")
        dataset
        """
        dataset = ro.r(r_script)

        # crreate FateAnnData object base expression and count matrix
        layers = {}
        if "expression" in dataset:
            X = dataset["expression"]
            layers["expression"] = dataset["expression"]
        if "count" in dataset:
            X = dataset["count"]
            layers["count"] = dataset["count"]
        fadata = cls(name=dataset["id"], X=X)
        fadata.layers = layers

        # other Anndata attributes
        # if dataset.has_key("cell_info"):
        #     fadata.obs = dataset["cell_info"]
        fadata.obs = dataset.get("cell_info", fadata.obs)  # equal to above
        fadata.obs.index = dataset["cell_ids"]
        fadata.var = dataset.get("feature_info", fadata.obs)
        fadata.var.index = dataset.get("feature_ids", fadata.var.index)

        # call FateAnnData object method
        if "prior_information" in dataset:
            fadata.add_prior_information(**dataset["prior_information"])
        if "milestone_network" in dataset:
            milestone_network = dataset["milestone_network"]
            milestone_percentages = dataset["milestone_percentages"]
            divergence_regions = dataset["divergence_regions"]
            progressions = dataset["progressions"]
            fadata.add_trajectory(
                milestone_network=milestone_network,
                divergence_regions=divergence_regions,
                milestone_percentages=milestone_percentages,
                # progressions=progressions # may cover milestone_percentages
            )
        # TODO: waypoint添加

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
        # TODO: save multiple trajectory in cfe_dict
        self.is_wrapped_with_trajectory = True

    def add_trajectory_by_type(self, trajectory_dict: dict) -> None:
        """Call the trajectory addition method based on specific trajectory types

        Args:
            trajectory_dict (dict): trajectory dict result based on specific trajectory types
        """

        # TODO: for more wrapper
        if "pseudotime" in trajectory_dict.keys():
            self.add_trajectory_linear(trajectory_dict["pseudotime"])
        elif "branch_network" in trajectory_dict.keys():
            self.add_trajectory_branch(
                branch_network=trajectory_dict["branch_network"],
                branches=trajectory_dict["branches"],
                branch_progressions=trajectory_dict["branch_progressions"]
            )
        else:
            # defult direct output
            self.add_trajectory(**trajectory_dict)

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
        """Add branch trajectory,such as PAGA

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

    def add_trajectory_linear(
        self,
        pseudotime: list,
        directed: bool = False,
        do_scale_minmax: bool = True,
    ) -> None:
        """add linear trajectory, such as Palantir, Cytotrace.

        Args:
            pseudotime (list): pseudotime sequence.
        """
        pseudotime = np.array(pseudotime)

        # min-max scale pseudotime to [0, 1]
        if do_scale_minmax:
            pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
        else:
            assert (pseudotime >= 0).all() and (pseudotime <= 1).all()
        milestone_ids = ["milestone_begin", "milestone_end"]
        # milestone_network datframe construction, length=1
        milestone_network = pd.DataFrame({
            "from": milestone_ids[0],
            "to": milestone_ids[1],
            "length": 1,
            "directed": directed,
        }, index=[0])  # all scalar, need "index" to show sample num
        # progressions datafram construction， percentage=pseudotime
        progressions = pd.DataFrame({
            "cell_id": self.obs.index,
            "from": milestone_ids[0],
            "to": milestone_ids[1],
            "percentage": pseudotime,
        })
        self.add_trajectory(
            milestone_network=milestone_network,
            divergence_regions=None,
            progressions=progressions
        )

    def add_trajectory_velocity():
        # TODO: add velocity trajectory, such as scVelo, VeloAE
        pass

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
    """_summary_

    Returns:
        _type_: _description_
    """
    adata = sc.read_h5ad(*args, **kwargs)
    fadata = FateAnnData.from_anndata(adata)
    return fadata
