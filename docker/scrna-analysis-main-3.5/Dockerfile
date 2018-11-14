## Build stage 1: Get private git repositories

FROM rocker/rstudio:3.5.0 as intermediate

# docker build --no-cache -t scrna-analysis-main-3.5:latest -t scrna-analysis-main-3.5:v1.2 .
# docker tag scrna-analysis-main-3.5:v1.2 alzhang/scrna-analysis-main-3.5:v1.2
# docker tag scrna-analysis-main-3.5:latest alzhang/scrna-analysis-main-3.5:latest
# docker push alzhang/scrna-analysis-main-3.5:v1.2
# docker push alzhang/scrna-analysis-main-3.5:latest

# Install git

RUN apt-get update
RUN apt-get install -y git ssh

# Git credentials
ARG GITHUB_SSH_KEY=keys/github_private_key
RUN mkdir /root/.ssh
ADD ${GITHUB_SSH_KEY} /root/.ssh/id_rsa

RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

RUN mkdir /root/repositories && \
    cd /root/repositories && \
    git clone git@github.com:Irrationone/cellassign.utils.git && \
    git clone git@github.com:Irrationone/scrnaseq.git && \
    git clone git@github.com:Irrationone/scrnatools.git && \
    git clone git@github.com:Irrationone/cellassign.git && \
    git clone git@github.com:Irrationone/vdj.utils.git && \
    git clone git@github.com:Irrationone/scvis.git && \
    git clone git@github.com:Irrationone/splatter.git && \
    git clone git@github.com:Irrationone/seurat.git

## Build stage 2: Actually build the Docker image

FROM alzhang/scrna-analysis-base-3.5:latest

# Copy over private repositories from previous image
COPY --from=intermediate /root/repositories/scrnaseq /scrnaseq
COPY --from=intermediate /root/repositories/cellassign.utils /cellassign.utils
COPY --from=intermediate /root/repositories/scrnatools /scrnatools
COPY --from=intermediate /root/repositories/cellassign /cellassign
COPY --from=intermediate /root/repositories/vdj.utils /vdj.utils
COPY --from=intermediate /root/repositories/scvis /scvis
COPY --from=intermediate /root/repositories/splatter /splatter
COPY --from=intermediate /root/repositories/seurat /seurat

RUN /bin/bash -c "PATH=$HOME/miniconda/bin:$PATH && source activate r-tensorflow && pip install /scrnatools"

USER rstudio

RUN /bin/bash -c "PATH=$HOME/miniconda/bin:$PATH && source activate r-tensorflow && pip install /scrnatools"

USER root

RUN Rscript -e "devtools::install('/scrnaseq/packages/R_packages/scrna.sceutils')" && \
    Rscript -e "devtools::install('/scrnaseq/packages/R_packages/scrna.utils')" && \
    Rscript -e "devtools::install('/cellassign.utils')" && \
    Rscript -e "devtools::install('/cellassign')" && \
    Rscript -e "devtools::install('/vdj.utils')" && \
    Rscript -e "devtools::install('/scvis')" && \
    Rscript -e "devtools::install('/splatter')" && \
    Rscript -e "devtools::install('/seurat')"

RUN rm -rf "/scrnaseq" && \
    rm -rf "/cellassign.utils" && \
    rm -rf "/scrnatools" && \
    rm -rf "/cellassign" && \
    rm -rf "/vdj.utils" && \
    rm -rf "/scvis" && \
    rm -rf "/splatter" && \
    rm -rf "/seurat"

RUN Rscript -e "install.packages('lmerTest',dependencies=TRUE)" && \
    Rscript -e "install.packages('argparse',dependencies=TRUE)" && \
    Rscript -e "install.packages('dunn.test',dependencies=TRUE)" && \
    Rscript -e "install.packages('nestedRanksTest',dependencies=TRUE)"

RUN echo 'BiocManager::install(c("ReactomePA"), update = FALSE)' >> /opt/packages.r && \
    Rscript /opt/packages.r
