FROM ubuntu:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
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

ADD conda/test.yml /tmp/environment1.yml

RUN conda env create -n bioinformatics -f /tmp/environment1.yml 
ENV PATH /opt/conda/envs/bioinformatics/bin:$PATH

RUN groupadd -g 2000 student \
&& useradd -m -u 2001 -g student student

USER student
WORKDIR /home/student

RUN conda init bash 
RUN mkdir -p session8 session9
ADD session8/ /home/student/session8/
ADD session9/ /home/student/session9/

CMD [ "/bin/bash" ]