import pandas as pd

from .._logging import logger
from ..util import random_time_string
from .fate_wrapper import FateWrapper


class MilestoneWrapper(FateWrapper):
    def __init__(
        self,
        milestone_network: pd.DataFrame,
        divergence_regions: pd.DataFrame = None,
        milestone_percentages: pd.DataFrame = None,
        progressions: pd.DataFrame = None
    ):
        self.id = random_time_string("MilestoneWrapper")
        self.milestone_network = milestone_network
        self.id_list = milestone_network[["from", "to"]].stack().unique().tolist()

        if divergence_regions is None:
            self.divergence_regions = pd.DataFrame(columns=["divergence_id", "milestone_id", "is_start"])
        else:
            self.divergence_regions = divergence_regions

        # choose milestone_percentages or progressions
        if (milestone_percentages is None) == (progressions is None):
            if milestone_percentages is not None:
                logger.warning("Both milestone_percentages and progressions are given, will only use progressions")
                milestone_percentages = None
            else:
                raise ValueError("Exactly one of milestone_percentages or progressions, must be defined, the other must be None")
        if milestone_percentages is not None:
            self.cell_id_list = milestone_percentages["cell_id"].unique().tolist()
        else:
            self.cell_id_list = progressions["cell_id"].unique().tolist()
        self.milestone_percentages = milestone_percentages
        self.progressions = progressions

        self.milestone_network_class = "N"
        self.directed = False

    @staticmethod
    def convert_milestone_percentages_to_progressions(
        milestone_network: pd.DataFrame,
        milestone_percentages: pd.DataFrame
    ):
        """
        milestone_percentages -> progressions, "add_trajectory" test case
        ref: pydynverse/wrap/convert_milestone_percentages_to_progressions.convert_milestone_percentages_to_progressions
        """
        # part1: for cells that have 2 or more milestones
        # first merge based on "to" key result in many invalid cell_id-form relationship
        df1 = pd.merge(milestone_network, milestone_percentages, left_on="to", right_on="milestone_id")
        # second merge based on "to" key
        df2 = pd.merge(df1, milestone_percentages[["cell_id", "milestone_id"]], left_on=["from", "cell_id"], right_on=["milestone_id", "cell_id"])
        # TODO: if the two step merge can be done simutaneously?
        progr_part1 = df2[["cell_id", "from", "to", "percentage"]]

        # for cells that have just 1 milestone
        # TODO: only simple reserve cells with one milestone
        progr_part2 = milestone_percentages.groupby("cell_id").filter(lambda x: len(x) == 1)
        progr_part2["from"] = progr_part2["milestone_id"]
        progr_part2["to"] = progr_part2["milestone_id"]
        progr_part2 = progr_part2[["cell_id", "from", "to", "percentage"]]

        # progressions = pd.concat([progr_part1], ignore_index=True)
        progressions = pd.concat([progr_part1, progr_part2], ignore_index=True).reset_index(drop=True)

        return progressions

    @staticmethod
    def convert_progressions_to_milestone_percentages(
        milestone_network: pd.DataFrame,
        progressions: pd.DataFrame
    ):
        """
        progressions -> milestone_percentages, "add_branch_trajectory" test case.
        resuse for "wrap_add_waypoint", written as the static method of MilestoneWrapper Class

        ref: pydynverse/wrap/convert_progressions_to_milestone_percentages.convert_progressions_to_milestone_percentages
        """
        # self loops
        selfs = progressions[progressions["from"] == progressions["to"]]
        selfs = selfs[["cell_id", "from"]].copy().rename(columns={"from": "milestone_id"})
        selfs["percentage"] = 1

        # not self loops
        progressions = progressions[~(progressions["from"] == progressions["to"])]
        # percentage for "from milestone"
        froms = progressions[["cell_id", "from", "percentage"]].copy().rename(columns={"from": "milestone_id"})
        froms["percentage"] = 1 - froms["percentage"]
        # percentage for "to milestone"
        tos = progressions[["cell_id", "to", "percentage"]].copy().rename(columns={"to": "milestone_id"})

        milestone_percentages = pd.concat([selfs, froms, tos]).reset_index(drop=True)

        return milestone_percentages

    def pipeline(self):
        # NOTE: no redundant check
        if (self.milestone_percentages is None) == (self.progressions is None):
            if self.milestone_percentages is not None:
                logger.warning("Both milestone_percentages and progressions are given, will only use progressions")
                self.milestone_percentages = None
            else:
                raise ValueError("Exactly one of milestone_percentages or progressions, must be defined, the other must be None")

        if self.progressions is None:
            # milestone_percentages -> progressions, 'add_trajectory' test case
            self.progressions = MilestoneWrapper.convert_milestone_percentages_to_progressions(self.milestone_network, self.milestone_percentages)
        else:
            # progressions -> milestone_percentages, 'add_branch_trajectory' test case
            self.milestone_percentages = MilestoneWrapper.convert_progressions_to_milestone_percentages(self.milestone_network, self.progressions)

        self.classify_milestone_network()

    def classify_milestone_network(self):
        """
        milestone network classification

        ref: pydynverse/wrap/wrap_add_trajectory.changed_topology
        """
        # TODO: PyDynverse and CFE implementation
        self.milestone_network_class = "N"
        self.directed = False

    def gather_cells_at_milestones(self):
        """
        move cells to their nearest milestone

        ref: pydynverse/wrap/wrap_gather_cells_at_milestones.gather_cells_at_milestones
        """
        # gather all cells to their nearest milestone
        pass
