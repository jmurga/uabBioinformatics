# Galaxy uabGenomics. Based on galaxy-ngs-preprocessing and galaxy-exome-seq Galaxy Docker Images (https://github.com/bgruening/docker-galaxy-ngs-preprocessing, https://github.com/bgruening/docker-galaxy-exome-seq)

FROM bgruening/galaxy-stable:19.01

MAINTAINER Jesus Murga Moreno, jesus.murga@uab.cat

ENV GALAXY_CONFIG_BRAND uabGenomics

# Install tools
ADD ngsPreprocessing.yml $GALAXY_ROOT/ngsTools.yaml
RUN install-tools $GALAXY_ROOT/ngsTools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf
