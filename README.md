# AcademicSoftware

This repository is meant to accompany Geodesy and Satellite Geodesy courses offered by DSO. It consists of:
 - a *core* python module named `dsoclasses`, on top of which
 - a list of [Jupyter Notebooks](https://jupyter.org/) are available, showcasing trivial space geodesy data analysis tasks.

To use the Notebooks, you need to install the core python module.

## Installation of `dsoclasses` module

This guide walks you through installing the `dsoclasses` Python module from source.

0. Prerequisites

Make sure your system has the following installed:

- Python 3.8 or higher
- `git`, `pip`, and optionally `venv`

 1. Clone the repository

```bash
git clone https://github.com/DSOlab/AcademicSoftware.git
cd AcademicSoftware
```

2. Create a Virtual Environment

**Note that depending on your OS and setup, you may need to replace `python` with `python3`.**

```bash
python -m venv .venv
source .venv/bin/activate
```

For **Windows** users the above should be replaced with:
```bash
python -m venv .venv
.venv\Scripts\activate
```

3. Upgrade build tools

**Note that depending on your OS and setup, you may need to replace `pip` with `pip3`.**

```bash
python -m pip install --upgrade pip
```

4. Build and Install the package

```bash
pip install -e .
```

### Updating

To fetch latest changes and/or additions to the online repository, you will need to 
run the following command from the root of the project (where the `pyproject.toml` is 
located): `git pull origin`. No other step should be needed.


## Jupyter Notebooks

The notebooks are placed under the `JupyterLab` folder. Hence, assuming jupyterlab 
is available on your system (if not, `pip install jupyterlab` would do it) the following 
command should launch a local web server and open JupyterLab in your browser 
`jupyter lab --notebook-dir=JupyterLab/` (from the top-level directory).

## A note on data

This repository comes with a few data files that are needed to run the examples presented in 
the notebooks. Alternate or updated data should be seeked at the dedicated web repositories.

## Bugs & Troubleshooting

If you encounter any problems, please contact:
* Xanthos Papanikolaou xanthos@mail.ntua.gr,
* Prof. Dimitris Anastasiou danastasiou@mail.ntua.gr
