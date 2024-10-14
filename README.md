# AcademicSoftware

## Installation

Building the project is handled by [pip](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/) and  
[hatch](https://hatch.pypa.io/latest/). 

Steps:

*  make sure you have Hatch installed; if not, you can install it via `pip install hatch`
*  go to the top-level directory (i.e. where the `pyproject.toml` file is located) and run: `hatch build`. This will create a `dist/` directory with the build artifacts (such as a `.whl` file). 
*  then, you can install the project using: `pip install dist/dsoacademicsoftware-0.1.2-py3-none-any.whl`; note that the exact name of the file may differ depending on version.
*  once installed, you can test if the executable works by running it from the command line, i.e. `dso_groundtrack --help`.

### Install in Editable Mode

Editable mode allows you to make changes to your project and have them reflected 
immediately without having to reinstall it every time. This is great for development.

To install in editable mode, run the following command from the root of the 
project (where the `pyproject.toml` is located): `pip install -e`.

## List of programs

* dso_groundtrack
* dso_skyplot
