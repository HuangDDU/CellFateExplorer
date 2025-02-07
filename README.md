# CellFateExplorer: An integrated platform for exploring cell fate
[![test](https://github.com/HuangDDU/CellFateExplorer/actions/workflows/dev_branch_test.yml/badge.svg)](https://github.com/HuangDDU/CellFateExplorer/actions/workflows/dev_branch_test.yml)
[![document](https://readthedocs.org/projects/cellfateexplorer/badge/?version=latest)](https://cellfateexplorer.readthedocs.io/en/latest/)

**Cell Fate Explorer(cfe)** is a integration platform for *inferring*, *visualizing* and *benchmarking* cell fate trajectory for single-cell RNA-seq data.

## Framework

![CellFateExplorer Framework](./docs/img/framework.png)

## Installation

1. clone the repository and enter the directory.
2. Create a conda environment and install the dependencies by running the following commands.

    ```bash
    conda create -n cfe python=3.10
    pip install -r requirements.txt
    ```

3. Pypi package will be released soon.

## Quick Start

You can run the [quickstart.ipynb](https://cellfateexplorer.readthedocs.io/en/latest/tutorial/quickstart/) using jupyter noboker to learn the basic function of tools quickly.

## Document

1. links: For [`User`](https://cellfateexplorer-cellfateexplorer.readthedocs-hosted.com/en/latest/api/), for [`Developer`](https://cellfateexplorer-cellfateexplorer.readthedocs-hosted.com/en/latest/api/)

2. if you want to build the docs locally, run the following command in now conda environment.

    ```bash
    pip install -r docs/requirements.txt
    mkdir -p docs/_build/html
    ```
