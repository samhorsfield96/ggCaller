FROM mambaorg/micromamba:latest

USER root

# create a project directory inside user home
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# create a project directory inside user home
# (this isn't used with a clone running snakemake)
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
# copy the code in
COPY . $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build conda env
ENV ENV_PREFIX $PROJECT_DIR/env
COPY --chown=$user:$user docker/environment_docker.yml /tmp/environment_docker.yml

COPY --chown=$user:$user docker/entrypoint.sh /usr/local/bin/
RUN chmod u+x /usr/local/bin/entrypoint.sh

RUN micromamba install -y -n base -f /tmp/environment_docker.yml && \
    micromamba clean --all --yes && \
    python -m pip install --no-deps --ignore-installed . && \
    mkdir ggc_db && \
    ggcaller --balrog-db ggc_db && \
    PATH=$PATH:/opt/conda/bin

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]
