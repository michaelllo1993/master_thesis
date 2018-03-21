# Single amino acid repeats (SAAR) evolutionary analysis pipeline

This repository consists of a suite of R and Perl scripts performing the evolutionary SAAR analysis. The workflow is automatized and maintained by [Snakemake](https://snakemake.readthedocs.io/en/stable/) software.

For the convenience of the portability and ease of use it is also distributed as a [Docker](https://www.docker.com/) container, which can be accessed at [DockerHub repository](https://hub.docker.com/r/mstolarczyk/saarpipeline/).

## Usage
**1. Install the Docker software**
  
  The Docker installation guide can be found [here](https://docs.docker.com/install/)

**2. Downolad the container with the pipeline setup and all dependencies installed**

`(sudo) docker pull mstolarczyk/saarpipeline`

**3. Make sure the image is available**

`(sudo) docker images`

The line above should output similar to the one presented below

```
REPOSITORY                         TAG                 IMAGE ID            CREATED             SIZE
mstolarczyk/saarpipeline           latest              1142d796e7ad        3 minutes ago       6.27GB
```

**4. Run the image**

`(sudo) docker run -it mstolarczyk/saarpipeline `

**5. Perform the analysis**

The commands that will run the analysis are `Snakemake` commands:

`snakemake <target file>`

To test the part of workflow to be executed run

`snakemake -np <target file>`

The line above will show the execution plan instead of actually perform the steps. The `-p` flag instructs Snakemake to also print the resulting shell command for illustration

_More details regarding the usage will be provided soon_

After completion of the analysis exit the container

`exit`

**6. Save the docker container state and copy the results**

In order to pick up when you have left save the state of the container

`(sudo) docker commit <container ID> <tag>`

e.g. `sudo docker commit 1bf3b5f36b94 mstolarczyk/saarpipeline_changed`
`(sudo)`

To copy the results of the analysis use the Docker `docker cp` command, for more details see the [website](https://docs.docker.com/engine/reference/commandline/cp/#description)

```
docker cp [OPTIONS] CONTAINER:SRC_PATH DEST_PATH|-
docker cp [OPTIONS] SRC_PATH|- CONTAINER:DEST_PATH
```

