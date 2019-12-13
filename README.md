# uabBioinformatics
Repository including a Galaxy Docker Image (Galaxy flavor). Galaxy instance is prepared to include software related with pre-processing NGS data, sequence statistics estimation, pairwise aligment, multiple alignment sequence and variant calling.

Galaxy Docker Image is based on [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable), other [Galaxy flavours](https://github.com/bgruening/docker-galaxy-stable#List-of-Galaxy-flavours) and (Galaxy Training)[https://galaxyproject.github.io/training-material/].

To include Galaxy Interactive Enviroment we include [Jupyter-Notebook](https://hub.docker.com/r/bgruening/docker-jupyter-notebook). To run other Docker images from Galaxy Instance, you need to give  execution permission from *docker* group on docker.sock and run Docker Image properly.

```bash 
sudo chmod g+rx /var/run/docker.sock
# To run docker image
docker run -i -t -p 8080:80 -p 8021:21 -p 8800:8800 \
 --privileged=True \
 -e GALAXY_DOCKER_ENABLED=True \
 -v /var/run/docker.sock:/var/run/docker.sock \
 -v /home/jmurga/galaxy_storage/:/export/ \
 uab/bioinformatics

```


Some connections to docker container are still in development. For that reason there are some scripts installed from [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu) are not working.
