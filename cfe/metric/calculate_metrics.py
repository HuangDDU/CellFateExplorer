from .topology_metric import calc_isomorphic, calculate_edge_flip, calculate_him
from .cluster_metric import calculate_mapping_branches, calculate_mapping_milestones


def calculate_metrics(
        fadata,
        now_model=None,
        ref_model=None,
        metrics=["isomorphic"]
):
    summary_dict = {}

    ref_simplified_milestone_wrapper = fadata.simplify_trajectory(ref_model)
    simplified_milestone_wrapper = fadata.simplify_trajectory(now_model)

    net1 = ref_simplified_milestone_wrapper.milestone_network  # ref model
    net2 = simplified_milestone_wrapper.milestone_network

    # topology metric
    if "isomorphic" in metrics:
        summary_dict["isomorphic"] = calc_isomorphic(net1, net2)
    if "edge_flip" in metrics:
        summary_dict["edge_flip"] = calculate_edge_flip(net1, net2)
    if "him" in metrics:
        summary_dict["him"] = calculate_him(net1, net2)

    # cluster metric
    if "F1_branch" in metrics:
        summary_dict["F1_branch"] = calculate_mapping_branches()
    if "F1_milestone" in metrics:
        summary_dict["F1_milestone"] = calculate_mapping_milestones()

    # TODO: other metric
    return summary_dict
