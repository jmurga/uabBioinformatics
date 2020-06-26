# Galaxy uabBioinformatics. Based on docker-galaxy-stable and galaxy-rna-workbench
## https://github.com/bgruening/docker-galaxy-stable
## https://github.com/bgruening/galaxy-rna-workbench

FROM bgruening/galaxy-stable

MAINTAINER Jesus Murga Moreno, jesus.murga@uab.cat

ENV GALAXY_CONFIG_BRAND UAB-Bioinformatics
ENV GALAXY_CONFIG_TOOL_PATH=/galaxy-central/tools/ \ 
    GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# WORKDIR /galaxy-central

# Install tools

#RUN install-biojs msa

ADD galaxy/dataManagers.yaml $GALAXY_ROOT/dataManagers.yaml
RUN install-tools $GALAXY_ROOT/dataManagers.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf

ADD galaxy/bioinformatics.yaml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf

# Container Style.
# Same structure of galaxy-rna-workbench. Bootstrap framework maintained customizing welcome page
ADD galaxy/assets/img/bioinfoCloud.png $GALAXY_CONFIG_DIR/web/welcome_image.png
ADD galaxy/welcome.html $GALAXY_CONFIG_DIR/web/welcome.html

# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Add workflows to the Docker image
# ADD ./rna-workbench-workflow/* $GALAXY_HOME/workflows/

# ENV GALAXY_CONFIG_TOOL_PATH=/galaxy-central/tools/

# Download training data and populate the data library
# RUN startup_lite && \ galaxy-wait && \ #workflow-install --workflow_path $GALAXY_HOME/workflows/ -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \  #setup-data-libraries -i $GALAXY_ROOT/library_data.yaml -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
