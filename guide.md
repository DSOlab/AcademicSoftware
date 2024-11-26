## Install Git

### Download Git:

    - Visit the official Git website: https://git-scm.com/.
    - Download the latest version of Git for Windows.

### Run the Installer:

    - Double-click the downloaded .exe file to start the installation.
    - During installation:
        - Choose "Use Git from the Command Prompt" to make Git available system-wide.
        - Select your preferred editor (e.g., Vim, Notepad++, VS Code).
        - Leave the default options for "Adjusting your PATH environment."
        - Choose "Git Credential Manager" for managing repository credentials.

### Verify Git Installation:

    - Open Command Prompt or PowerShell and run: `git --version`

## Install Python and Pip

### Download Python:

    - Visit the official Python website: https://www.python.org/.
    - Download the latest version for Windows.

### Run the Installer:

    - Double-click the downloaded .exe file to start the installation.
    - Ensure you check the box "Add Python to PATH" at the bottom of the installer.
    - Choose "Customize Installation" and enable these optional features:
        - pip
        - tcl/tk and IDLE
        - Python test suite
    - Complete the installation process.

### Verify Python and pip Installation:

    - Open Command Prompt or PowerShell and run:
    ```
    python --version
    pip --version
    ```

## Install Jupyter Lab

    - Install Jupyter Notebook Using pip:
        - Open Command Prompt or PowerShell and run: `pip install jupyterlab`

### Verify Installation:

    - Run the following to verify Jupyter Notebook is installed: `jupyter lab --version`. This should output the Jupyter version.

## Clone the Repository

    - Clone the Repository:
        Navigate to the directory where you want to store the project. Then run: `git clone https://github.com/DSOlab/AcademicSoftware.git`

    - Navigate to the Project Directory and Change into the cloned repository's directory: `cd AcademicSoftware`

## Build and Install the Python Module

    - Build and Install Using pip: `pip install .`

    - Verify Installation: After installation, test if the module is accessible: `python -c "from dsoclasses.rinex.gnss.rinex import GnssRinex"`

## Run Jupyter Notebooks Locally

    - Start Jupyter Notebook: From the cloned repository's directory, start Jupyter Notebook: `jupyter lab --notebook-dir=JupiterLab/`



# Installation with Anaconda

```
conda create -n myenv python=3.10
conda activate myenv
conda install git
conda install -c conda-forge jupyterlab
pip install .
conda install -c conda-forge build-essential
git clone https://github.com/DSOlab/AcademicSoftware.git
cd AcademicSoftware
pip install .
jupyter lab --notebook-dir=JupiterLab/
```
