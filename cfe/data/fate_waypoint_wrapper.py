import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics.pairwise import pairwise_distances

from ..util import random_time_string
from .fate_wrapper import FateWrapper
from .fate_milestone_wrapper import MilestoneWrapper


class WaypointWrapper(FateWrapper):
    def __init__(
        self,
        milestone_wrapper: MilestoneWrapper,
        name="WaypointWrapper"
    ):
        self.id = random_time_string(name)
        self.milestone_wrapper = milestone_wrapper

    def pipeline(self):
        self._select_waypoints()

    def _select_waypoints(
            self,
            n_waypoints=200,
            transform=lambda x: x,  # edge length transform function
            resolution=None):
        """
        select waypoints base milestone network edge length and resolution parameter
        ref: pydynverse/wrap/wrap_add_waypoints.select_waypoints
        """
        mr = self.milestone_wrapper

        if resolution is None:
            # compute resolution automaticall based on the sum of milestone network length after transformation
            resolution = mr.milestone_network["length"].apply(lambda x: transform(x)).sum() / n_waypoints

        # percentage list construction and explode
        def waypoint_id_from_progressions_row(row):
            # get waypoint_id by considering a row comprehensively
            match row["percentage"]:
                case 0:
                    return f"MILESTONE_BEGIN_W{row['from']}_{row['to']}"
                case 1:
                    return f"MILESTONE_END_W{row['from']}_{row['to']}"
                case _:
                    return f"W{row.name+1}"  # waypoint id start from 1
        waypoint_progressions = mr.milestone_network.copy()
        waypoint_progressions["percentage"] = waypoint_progressions["length"].apply(lambda x: [i / x for i in np.arange(0, x, resolution)] + [1])
        waypoint_progressions = waypoint_progressions[["from", "to", "percentage"]]
        waypoint_progressions = waypoint_progressions.explode("percentage").reset_index(drop=True)
        waypoint_progressions["percentage"] = waypoint_progressions["percentage"].astype("float")
        waypoint_progressions["waypoint_id"] = waypoint_progressions.apply(waypoint_id_from_progressions_row, axis=1)
        self.waypoint_progressions = waypoint_progressions

        self.id_list = waypoint_progressions["waypoint_id"].unique().tolist()

        # progressions -> percentages
        waypoint_progressions_tmp = waypoint_progressions.copy()
        waypoint_progressions_tmp = waypoint_progressions_tmp.rename(columns={"waypoint_id": "cell_id"})  # reuse pre column name
        # tmp "cell_id" column name for MilestoneWrapper.reuse convert_progressions_to_milestone_percentages
        waypoint_milestone_percentages = MilestoneWrapper.convert_progressions_to_milestone_percentages(
            milestone_network=mr.milestone_network,
            progressions=waypoint_progressions_tmp
        ).rename(columns={"cell_id": "waypoint_id"})
        self.waypoint_milestone_percentages = waypoint_milestone_percentages

        self.waypoint_geodesic_distances = self._calculate_geodesic_distances().loc[waypoint_progressions["waypoint_id"]]

        waypoint_network = waypoint_progressions\
            .sort_values(by=["from", "to", "percentage"])\
            .groupby(["from", "to"])\
            .apply(lambda group: group.assign(
                from_waypoint=group["waypoint_id"],
                to_waypoint=group["waypoint_id"].shift(-1),
            ))\
            .dropna()\
            .reset_index(drop=True)  # 分组后, 组内按照percentage排序, lead函数向当前看下一行的元素, 如果是组内最后一个元素则获得NULL
        waypoint_network = waypoint_network[["from_waypoint", "to_waypoint", "from", "to"]]
        waypoint_network.columns = ["from", "to", "from_milestone_id", "to_milestone_id"]
        self.waypoint_network = waypoint_network

        waypoints = waypoint_milestone_percentages.iloc[waypoint_milestone_percentages.groupby("waypoint_id")["percentage"].idxmax()].reset_index(drop=True)
        waypoints["milestone_id"] = waypoints.apply(lambda x: x["milestone_id"] if x["percentage"] == 1 else None, axis=1)  # 不在里程碑上的waypoint的milestone_id为None
        waypoints = waypoints[["waypoint_id", "milestone_id"]]
        self.waypoints = waypoints

    def _calculate_geodesic_distances(self):
        """
        calculate geodesic distances between cells and waypoints/milestones
        overall idea:
            1. calculate the full path of the target point within each divergent region separately
            2. merge and calculate the distance on the overall graph

        ref: PyDynverse/pydynverse/wrap/calculate_geodesic_distances.py
        """
        # attribute in the MilestoneWrapper
        cell_id_list = self.milestone_wrapper.cell_id_list
        milestone_id_list = self.milestone_wrapper.id_list
        milestone_network = self.milestone_wrapper.milestone_network
        milestone_percentages = self.milestone_wrapper.milestone_percentages
        divergence_regions = self.milestone_wrapper.divergence_regions
        directed = self.milestone_wrapper.directed

        waypoint_id_list = self.id_list
        waypoint_milestone_percentages = self.waypoint_milestone_percentages

        milestone_percentages = pd.concat([
            milestone_percentages,
            waypoint_milestone_percentages.rename(columns={"waypoint_id": "cell_id"})
        ])

        # remae all milestone ids to MILESTONE_ID
        def milestone_trafo_fun(x): return f"MILESTONE_{x}"
        milestone_network["from"] = milestone_network["from"].apply(milestone_trafo_fun)
        milestone_network["to"] = milestone_network["to"].apply(milestone_trafo_fun)
        milestone_id_list = list(map(milestone_trafo_fun, milestone_id_list))
        milestone_percentages["milestone_id"] = milestone_percentages["milestone_id"].apply(milestone_trafo_fun)
        divergence_regions["milestone_id"] = divergence_regions["milestone_id"].apply(milestone_trafo_fun)

        # 添加extra发散区域, 正常的边也被当作发散区域
        extra_divergences = milestone_network.copy()
        extra_divergences = extra_divergences[~(extra_divergences["from"] == extra_divergences["to"])]
        # extra_divergences = extra_divergences.query("not from == to") # query更加优雅
        # in_divergence判断当前边是否在已有的发散区域内，标准的延迟承诺区域
        divergence_regions_set_list = divergence_regions.groupby("divergence_id")["milestone_id"].apply(set).tolist()

        def is_milestone_in_divergence(milestone_set, divergence_regions_set_list):
            for divergence_regions_set in divergence_regions_set_list:
                if milestone_set.issubset(divergence_regions_set):
                    return True
            return False
        extra_divergences["in_divergence"] = extra_divergences.apply(lambda x: is_milestone_in_divergence({x["from"], x["to"]}, divergence_regions_set_list), axis=1)
        extra_divergences = extra_divergences[~extra_divergences["in_divergence"]]  # 只保留新的发散区域
        extra_divergences["divergence_id"] = extra_divergences.apply(lambda x: f"{x['from']}__{x['to']}", axis=1)
        extra_divergences = pd.concat([
            # 添加新的milestone_id和is_start列
            extra_divergences.assign(milestone_id=extra_divergences["from"], is_start=True),
            extra_divergences.assign(milestone_id=extra_divergences["to"], is_start=False)
        ])[["divergence_id", "milestone_id", "is_start"]]

        # 合并发散区域
        divergence_regions = pd.concat([divergence_regions, extra_divergences]).reset_index(drop=True)
        divergence_ids = divergence_regions["divergence_id"].unique()

        # 准备使用NetworkX, 构造相关数据, 从DataFrame开始构造
        milestone_graph = nx.from_pandas_edgelist(milestone_network, source="from", target="to", edge_attr="length")

        # NOTE: 1. 分别计算
        # 计算发散区域内部细胞间距离
        def calc_divergence_inner_distance_df(did):
            dir = divergence_regions[divergence_regions["divergence_id"] == did]
            mid = dir[dir["is_start"]]["milestone_id"].tolist()  # 该区域起点milestone_id
            tent = dir["milestone_id"].tolist()  # 该区域所有milestone_id
            tent_distances = pd.DataFrame(index=mid, columns=tent, data=np.zeros((len(mid), len(tent))))  # 区域内起点到所有milestone的距离
            # 从图中提取对应的边
            for i in mid:
                for j in tent:
                    if i == j:
                        tent_distances.loc[i, j] = 0
                    else:
                        tent_distances.loc[i, j] = milestone_graph.edges[(i, j)]["length"]
            # 此处复用来找相关的点的cell_id
            relevant_pct_cell_id_list = milestone_percentages.groupby("cell_id")["milestone_id"].apply(lambda x: is_milestone_in_divergence(set(x), [set(tent)]))
            relevant_pct_cell_id_list = relevant_pct_cell_id_list[relevant_pct_cell_id_list].index.to_list()
            relevant_pct = milestone_percentages[milestone_percentages["cell_id"].apply(lambda x: x in relevant_pct_cell_id_list)]
            if relevant_pct.shape[0] <= 1:
                return None

            scaled_dists = relevant_pct.copy()
            scaled_dists["dist"] = scaled_dists.apply(lambda x: x["percentage"]*tent_distances.loc[mid, x["milestone_id"]], axis=1)
            tent_distances_long = tent_distances.melt(var_name="from", value_name="length")  # 宽数据转化为长数据
            tent_distances_long["to"] = tent_distances_long["from"]

            pct_mat = pd.concat([
                scaled_dists[["cell_id", "milestone_id", "dist"]].rename(columns={"cell_id": "from", "milestone_id": "to", "dist": "length"}),
                tent_distances_long
            ]).pivot(index="from", columns="to", values="length").fillna(0)  # (n_cell+n_milestone+n_waypoint)*n_milestone, 长数据转宽数据, 索引名为from

            wp_cells = list(set(pct_mat.index) & set(waypoint_id_list))

            if directed:
                # TODO: 暂时不管有向图
                pass

            distances = pairwise_distances(pct_mat, pct_mat.loc[wp_cells+tent], metric="manhattan")
            distances = pd.DataFrame(index=pct_mat.index, columns=wp_cells+tent, data=distances)
            distances = distances.reset_index().melt(id_vars="from", var_name="to", value_name="length")  # 宽数据转化为长数据
            distances = distances[~(distances["from"] == distances["to"])]
            return distances

        cell_in_tent_distances = pd.concat([calc_divergence_inner_distance_df(did) for did in divergence_ids])

        if directed:
            # TODO: 暂时不管有向图
            pass

        # NOTE: 2. 合并计算
        # 合并两个图到一张图上
        graph = pd.concat([milestone_network, cell_in_tent_distances]).groupby(["from", "to"]).agg({"length": "min"}).reset_index()  # 合并后提取最短边，目前没什么效果，可能对环图有用
        graph = nx.from_pandas_edgelist(graph, source="from", target="to", edge_attr="length")
        # 选择后续有向图的最短距离模式
        if directed or directed == "forward":
            mode = "out"
        elif directed == "reverse":
            mode = "in"
        else:
            mode = "all"
        mode

        # 调用dijkstra计算, 计算waypoint和cell之间的
        out = pd.DataFrame(np.zeros((len(waypoint_id_list), len(cell_id_list))), index=waypoint_id_list, columns=cell_id_list)
        for source in waypoint_id_list:
            length_dict = nx.single_source_dijkstra_path_length(graph, source=source, weight="length")
            for target in cell_id_list:
                if target in length_dict:
                    out.loc[source, target] = length_dict[target]
                else:
                    out.loc[source, target] = np.inf

        # TODO: 过滤一些细胞
        cell_ids_filtered_list = []
        if len(cell_ids_filtered_list) > 0:
            pass

        return out.loc[waypoint_id_list, cell_id_list]
