# Galaxy uabGenomics. Based on galaxy-ngs-preprocessing and galaxy-exome-seq 
## Galaxy Docker Images (https://github.com/bgruening/docker-galaxy-ngs-preprocessing,
## https://github.com/bgruening/docker-galaxy-exome-seq)

FROM bgruening/galaxy-stable

MAINTAINER Jesus Murga Moreno, jesus.murga@uab.cat

ENV GALAXY_CONFIG_BRAND UAB-Genomics

# Install tools

ADD ngsPreprocessing.yaml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf
