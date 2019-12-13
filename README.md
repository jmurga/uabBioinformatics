# Galaxy Docker Image - Bioinformatics Universitat Autònoma de Barcelona
Repository including a Galaxy Docker Image (Galaxy flavor). Galaxy instance is prepared to include software related with pre-processing NGS data, sequence statistics, pairwise alignment, multiple alignment sequence and variant calling. We deposited all the required software to complete the [bioinformatics course](https://stepik.org/course/1171/) from Genetics bachelor's degree at Universitat Autonoma de Barcelona. Please check [Galaxy-Training](https://training.galaxyproject.org/training-material/) when to understand pipelines in Galaxy!

Galaxy Docker Image is based on [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable), other [Galaxy flavours](https://github.com/bgruening/docker-galaxy-stable#List-of-Galaxy-flavours) and [Galaxy Training](https://galaxyproject.github.io/training-material/).

To include Galaxy Interactive Enviroment you need to pull [Docker-Jupyter-Notebook](https://hub.docker.com/r/bgruening/docker-jupyter-notebook). To run other Docker images from Galaxy Instance, you need to give  execution permission from *docker* group on docker.sock and run Docker Image properly.

Some connections to docker container are still in development. For that reason there are some scripts installed from [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu) are not working.


```bash 
# To include docker jupyter and rstudio
docker pull bgruening/docker-jupyter-notebook
docker pull quay.io/erasche/docker-rstudio-notebook
# To execute docker containers from galaxy
sudo chmod g+rx /var/run/docker.sock
# To run docker image
docker run -i -t -p 8080:80 -p 8021:21 -p 8800:8800 \
 --privileged=True \
 -e GALAXY_DOCKER_ENABLED=True \
 -v /var/run/docker.sock:/var/run/docker.sock \
 -v /home/${USER}/galaxy_storage/:/export/ \
 uab/bioinformatics
```