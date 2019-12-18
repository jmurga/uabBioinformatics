## Bioinformatics - Universitat Autònoma de Barcelona

Repository including material developed to complete the [bioinformatics course](https://stepik.org/course/1171/) from Genetics bachelor's degree at Universitat Autonoma de Barcelona. We recommend use Galaxy instance or Debian image included at the repository, both dockerized. They are prepared to include software related with pre-processing NGS data, sequence statistics, pairwise alignment, multiple alignment sequence and variant calling. Please check [Galaxy-Training](https://training.galaxyproject.org/training-material/) when to understand pipelines in Galaxy!

### Galaxy instance

Galaxy Docker Image is based on [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable), other [Galaxy flavours](https://github.com/bgruening/docker-galaxy-stable#List-of-Galaxy-flavours) and [Galaxy Training](https://galaxyproject.github.io/training-material/). To include Galaxy Interactive Enviroment you need to pull [Docker-Jupyter-Notebook](https://hub.docker.com/r/bgruening/docker-jupyter-notebook). To run other Docker images from Galaxy Instance, you need to give  execution permission from *docker* group on docker.sock and run Docker Image properly.

Some connections to docker container are still in development. For that reason some scripts installed from [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu) are not working. You may take this into account if you want to expand the Galaxy flavour.

```bash 
# Download galaxy 
docker pull jmurga/uab-bioinformatics
# To include docker jupyter and rstudio
docker pull bgruening/docker-jupyter-notebook
docker pull erasche/docker-rstudio-notebook
# To execute docker containers from galaxy
chmod g+rx /var/run/docker.sock
# To run docker image
docker run -d -t -p 8080:80 -p 8021:21 -p 8800:8800 \
 --privileged=True \
 -e GALAXY_DOCKER_ENABLED=True \
 -v /var/run/docker.sock:/var/run/docker.sock \
 -v /home/${USER}/galaxy_storage/:/export/ \
 jmurga/uab-bioinformatics
```

### Debian image
Debian image conda enviroment is included on folder *debian/*. Based on [miniconda image](https://hub.docker.com/r/continuumio/miniconda3). Enviroment *bioinformatics* and user *student* are activated by default when running the image. *student* user has no root permission to avoid problems with specific *root* paths and commands. */home/student* folder is the main working directory where user *student* has full permisions. In addition *student* cannot create or modificate conda enviroments.

Include two images 

To build or pull the image run the following commands
```bash 
docker pull jmurga/uab-debian-bioinformatics
# or
docker build -t uab/debian-bioinformatics -f debian/Dockerfile .
```

To run the images with jupyter notebook on [localhost:8888](http://localhost:8888)

```bash
# Run docker bash interactive session
docker run -i -t -p 8888:8888 \
	-v /tmp/.X11-unix:/tmp/.X11-unix \
    -e DISPLAY=unix$DISPLAY \
    ubuntu/bioinfo /bin/bash
# Run only jupyter notebook from docker image
docker run -i -t -p 8888:8888 uab/bioinfo /bin/bash -c "jupyter notebook --ip='*' --port=8888 --no-browser"
```