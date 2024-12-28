# Welcome to Cell Fate Explorer🔍 Document

## Introduction

**Cell Fate Explorer(cfe)** is a integration platform for *inferring*, *visualizing* and *benchmarking* cell fate trajectory for single-cell RNA-seq data.

## Key Concept

**FateAnnData**: Unified trajectory inference data structure base AnnData[^1], includes two types of key external data:

- **Milestone**: The key nodes in the trajectory of cell fate are milestones, and the milestone network is a simplification of the cell trajectory.
- **Waypoint**: Sampling points for cell fate and trajectories, can be used to simplify calculations and visualize trajectories in embedding space.

**FateMethod**: Unified trajectory inference method interface, includes three types of backend ( Users can select one of the calls during the execution process ), we highly recommend **Python Function** or **CFE Docker** backend:

|Backend|Description|Advantage|Disadvantage|
| ---- | ---- | ---- | ---- |
|**Python Function**|The function developed by this project incorporates the latest trajectory inference methods in recent years, making it particularly well-suited to the project's framework.|1. New Methods in recenty years. <br>| 1. Different trajectory inference package versions in the same Python environment may conflict. |
|**Dynverse Docker**|Docker image for trajectory inference refers to dynverse[^2].|1. The ease of use of Docker |1.Methods on R language not be compatible. <br> 2.Methods are old relatively. <br> 3. Docker environment is need.
|**CFE Docker**|Docker image for trajectory inference are developed by this project.|1. New Methods in recenty years. <br>2. The ease of use of Docker. |1. Docker environment is need. |

## Toc

- [Tutorial](./tutorial.md)
- [API document](./api.md)
- [Change log](./change_log.md)

[^1]: AnnData: https://anndata.readthedocs.io/en/stable/

[^2]: Saelens, Wouter, et al. "A comparison of single-cell trajectory inference methods." Nature biotechnology 37.5 (2019): 547-554.
