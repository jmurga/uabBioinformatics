FROM debian:latest

RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

ENV LANG en_US.utf8

ENV PATH /opt/conda/bin:$PATH
ARG condaEnv=bioinformatics

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git nano vim libxtst6 libxrender1 libxext6 libexpat1 fonts-dejavu-core fontconfig-config libfontconfig1 libfreetype6 libpng16-16 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

RUN conda update -n base -c defaults conda

ADD debian/bioinformatics.yml /tmp/environment1.yml
ADD debian/genomics.yml /tmp/environment2.yml
ADD debian/artemis.yml /tmp/environment3.yml

RUN conda env create -n bioinformatics -f /tmp/environment1.yml 
ENV PATH /opt/conda/envs/bioinformatics/bin:$PATH
RUN conda env create -n genomics -f /tmp/environment2.yml 
ENV PATH /opt/conda/envs/genomics/bin:$PATH
RUN conda env create -n artemis -f /tmp/environment3.yml 
ENV PATH /opt/conda/envs/artemis/bin:$PATH

RUN conda init bash 
# RUN /opt/conda/bin/conda activate bioinformatics
# RUN conda config --set auto_activate_base false
ENV PATH /opt/conda/envs/bioinformatics/bin:$PATH

RUN mkdir -p /data
WORKDIR /data

RUN wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz && \
    tar -xvf exonerate-2.2.0-x86_64.tar.gz && \
    cp -r exonerate-2.2.0-x86_64/bin* /bin/ && \
    rm exonerate-2.2.0-x86_64.tar.gz

RUN echo "LC_ALL=C" >> ~/.bashrc
RUN echo "LANG=C" >> ~/.bashrc
RUN echo "export LC_ALL LANG" >> ~/.bashrc

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate bioinformatics" >> ~/.bashrc

RUN echo "source activate $(head -1 /tmp/environment1.yml | cut -d' ' -f2)" >> ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH



# CMD [ "/bin/bash" ]

