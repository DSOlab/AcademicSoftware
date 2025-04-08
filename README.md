# AcademicSoftware

## Installation of `dsoclasses` module

This guide walks you through installing the [AcademicSoftware](https://github.com/DSOlab/AcademicSoftware) Python module from source using [Hatch](https://hatch.pypa.io/).


0. Prerequisites
---

Make sure your system has the following installed:

- Python 3.8 or higher
- `git`, `pip`, and optionally `venv`

1. Clone the repository
---

```bash
git clone https://github.com/DSOlab/AcademicSoftware.git
cd AcademicSoftware
```

2. Create a Virtual Environment (Optional)
---

```bash
python3 -m venv .venv
source .venv/bin/activate
```

3. Install  [Hatch](https://hatch.pypa.io/)
---

```bash
pip install hatch
```

4. Build and Install the package
---

```bash
hatch build
pip install dist/*.whl
```

### Install in Editable Mode (Optional)

Editable mode allows you to make changes to your project and have them reflected 
immediately without having to reinstall it every time. This is great for development.

To install in editable mode, run the following command from the root of the 
project (where the `pyproject.toml` is located): `pip install -e .`


