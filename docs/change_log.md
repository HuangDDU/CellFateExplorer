# Change Log

## 2024.12.11: initial repository

## 2024.12.24: 0.1.0 version, project framework finish

### File structure

1. `cfe`: **source code**
   - `data`:
     - base: FateAnnData.
     - two trajectory related wrapper: MilestoneWrapper, WaypointWrapper.
   - `method`:
     - base: FateMethod.
     - three method backend: Function, CFE Docker, Dynverse Docker(based R).
   - `plot`:
     - plot_trajectory: plot embdding and trajectory
   - `util`: common utils
   - `_logging.py`: logging module
   - `_settings.py`: settings module
2. `docs`: **document** based on **mkdocs**, link to [CellFateExplorer](https://cellfateexplorer-cellfateexplorer.readthedocs-hosted.com/en/latest//)
3. `notebooks`: **jupyter notebook** for quickstart using Slingshot and PAGA
4. `test`: **unit test** with case

### Workflow for new requirement

1. create new *xxx_feature* branch based on *dev* branch(Administrator).
2. write source(`cfe`) and test`test` code locally on PC(Developer).
3. update document(`docs`)(Developer).
4. commit and push to remote(Developer).
5. merge to *dev* branch, check Github Action and readthedocs(Administrator)
6. merge to *main* branch, tag it, check Github Action and readthedocs (Administrator)

### Version Management

|Version|Aim and task|
| ---- | ---- |
| 0.1.0 | Project framework |
| 0.2.0 | TI methods and visualization |
| 0.3.0 | Metric construction |
| 0.4.0 | Simulation data and benchmark |
| 0.5.0 | Web UI for Explorer |
| 0.6.0 | Downstream analysis demo|
| 1.0.0 | The First paper |
| N.0.0 | The Second paper |
