FROM ghcr.io/prefix-dev/pixi:0.60.0-bullseye-slim AS build

# Use bash as shell
SHELL ["/bin/bash", "-c"]

WORKDIR /app

ENV PATH=$PATH:/opt/conda/bin:/opt/conda/condabin

RUN pixi global install pixi-install-to-prefix pixi-inject

COPY src /app
COPY pixi.lock pyproject.toml README.md /app/

RUN pixi-install-to-prefix -l pixi.lock /opt/conda && \
    pixi build --locked --path pyproject.toml && \
    pixi inject --prefix /opt/conda --package coidb*.conda


FROM condaforge/mambaforge:24.9.2-0 AS production

SHELL ["/bin/bash", "-c"]

WORKDIR /app

COPY --from=build /opt/conda /opt/conda

EXPOSE 8000

ENTRYPOINT [ "coidb" ]