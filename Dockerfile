FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="Image with tools used in demuxlink" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/demuxlink/bin:$PATH
