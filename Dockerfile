FROM ghcr.io/prefix-dev/pixi:0.60.0-bullseye-slim AS build

# Use bash as shell
SHELL ["/bin/bash", "-c"]

WORKDIR /app

COPY pixi.toml pixi.lock README.md pyproject.toml /app/
COPY src /app

RUN pixi build


FROM condaforge/mambaforge:24.9.2-0 AS install

# Use bash as shell
SHELL ["/bin/bash", "-c"]

WORKDIR /analysis

COPY --from=build /app/*.conda /app/coidb.conda

RUN conda install -c bioconda -c conda-forge snk-cli python && \
    conda install /app/coidb.conda && \
    conda clean -ay && \
    rm /app/coidb.conda

ENTRYPOINT [ "coidb" ]