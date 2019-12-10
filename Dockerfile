# Galaxy uabGenomics. Based on galaxy-ngs-preprocessing and galaxy-exome-seq 
## Galaxy Docker Images (https://github.com/bgruening/docker-galaxy-ngs-preprocessing,
## https://github.com/bgruening/docker-galaxy-exome-seq)

FROM bgruening/galaxy-stable

MAINTAINER Jesus Murga Moreno, jesus.murga@uab.cat

ENV GALAXY_CONFIG_BRAND UAB-Genomics

WORKDIR /galaxy-central


# Install tools

ADD dataManagers.yaml $GALAXY_ROOT/managers.yaml
RUN install-tools $GALAXY_ROOT/managers.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf

ADD genomics.yaml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml && \
    /tool_deps/_conda/bin/conda clean --tarballs && \
    rm /export/galaxy-central/ -rf



# Add workflows to the Docker image
# ADD ./rna-workbench-workflow/* $GALAXY_HOME/workflows/

# ENV GALAXY_CONFIG_TOOL_PATH=/galaxy-central/tools/

# Download training data and populate the data library
# RUN startup_lite && \
#    galaxy-wait && \
    #workflow-install --workflow_path $GALAXY_HOME/workflows/ -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \
    #setup-data-libraries -i $GALAXY_ROOT/library_data.yaml -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    #  
# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]