# 1. Use a Python base image
FROM python:3.10

# 2. Set the working directory
WORKDIR /app

# 3. Install build tools for building the Python module
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 4. Copy the pyproject.toml and other build files
COPY pyproject.toml ./

# 5. Install dependencies using pip or poetry
# Assuming pyproject.toml follows PEP 517/518 and works with pip:
RUN pip install --no-cache-dir build \
    && python -m build --wheel \
    && pip install --no-cache-dir dist/*.whl

# 6. Copy the entire project to the container
COPY . .

# 7. Install Jupyter Notebook
RUN pip install --no-cache-dir notebook

# 8. Expose Jupyter Notebook default port
EXPOSE 8888

# 9. Set the default command to start Jupyter Notebook
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
