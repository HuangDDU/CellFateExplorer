from .._logging import logger
from ..data import FateAnnData
from ..method import FateMethod


def plot_trajectory(
        fdata: FateAnnData = None,
        fmethod: FateMethod = None,
        ax=None,
):
    logger.debug("plot_trajectory")
    # NOTE: a fdata, a method
    # TODO: a fdata, many methods
