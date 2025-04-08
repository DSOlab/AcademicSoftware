# AcademicSoftware

This repository is meant to accompany Geodesy and Satellite Geodesy courses offered by DSO. It consists of:
 - a *core* python module named `dsoclasses`, on top of which
 - a list of [Jupyter Notebooks](https://jupyter.org/) are available, showcasing trivial space geodesy data analysis tasks.

To use the Notebooks, you need to install the core python module.

## Installation of `dsoclasses` module

This guide walks you through installing the [AcademicSoftware](https://github.com/DSOlab/AcademicSoftware) Python module from source using [Hatch](https://hatch.pypa.io/).


0. Prerequisites

Make sure your system has the following installed:

- Python 3.8 or higher
- `git`, `pip`, and optionally `venv`

 1. Clone the repository

```bash
git clone https://github.com/DSOlab/AcademicSoftware.git
cd AcademicSoftware
```

2. Create a Virtual Environment (Optional)

```bash
python3 -m venv .venv
source .venv/bin/activate
```

3. Install  [Hatch](https://hatch.pypa.io/)

```bash
pip install hatch
```

4. Build and Install the package

```bash
hatch build
pip install dist/*.whl
```

### Install in Editable Mode (Optional)

Editable mode allows you to make changes to your project and have them reflected 
immediately without having to reinstall it every time. This is great for development.

To install in editable mode, run the following command from the root of the 
project (where the `pyproject.toml` is located): `pip install -e .`

## Jupyter Notebooks

The notebooks are placed under the `JupiterLab` folder. Hence, assuming jupyterlab is available on your system (if not, `pip install jupyterlab` would do it) the following command should 
launch a local web server and open JupyterLab in your browser `jupyter lab --notebook-dir=JupiterLab/` (from the top-level directory).

## A note on data

The notebooks use data to showcase different analysis tasks. The data though are not distributed with the project, and you'll have to download them youself. You can conviniently place them at a `data/` top-level directory, or anywhere else you see fit. The actual filenames that appear on the notebooks are not (supposed to be) binding; its up to the user selection.
