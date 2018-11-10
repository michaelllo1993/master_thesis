FROM ubuntu:xenial as OS

MAINTAINER Michal Stolarczyk (stolarczyk.michal93@gmail.com)

# Add cran R backport
RUN apt-get -y update
RUN apt-get -y install apt-transport-https
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN echo "deb https://mirrors.ebi.ac.uk/CRAN/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list

# Update & upgrade sources
RUN apt-get -y update
RUN apt-get -y dist-upgrade

#Install system requirements for the software packages below
RUN apt-get install -y build-essential libssl-dev libffi-dev python-dev libcurl4-openssl-dev libxml2-dev wget libx11-dev libcairo2-dev cpanminus vim libcairo2-dev libjpeg-dev libgif-dev
RUN apt -y install python3-pip

#Install snakemake 
RUN pip3 install snakemake

#Install R and R packages
RUN apt install -y r-base
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('seqinr')"
RUN Rscript -e "install.packages('RCurl')"
RUN Rscript -e "install.packages('httr')"
RUN Rscript -e "install.packages('XML')"
RUN Rscript -e "install.packages('seqinr')"
RUN Rscript -e "install.packages('svglite')"
RUN Rscript -e "install.packages('ape')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('svglite')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('Biostrings')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('biomaRt')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('msa')"
#Install required Perl modules
RUN cpan List::MoreUtils
#Install git
RUN apt-get install -y git
#Clone the git repository
RUN git clone https://github.com/michalstolarczyk/SAARpipeline.git

#copy the example data
#COPY ./Pipeline/data SAARpipeline/Pipeline/data

#download and install EMBOSS software
WORKDIR SAARpipeline/Pipeline
RUN wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
RUN gunzip EMBOSS-6.6.0.tar.gz 
RUN tar xvf EMBOSS-6.6.0.tar
RUN rm EMBOSS-6.6.0.tar
RUN mv EMBOSS-6.6.0 software
WORKDIR software/EMBOSS-6.6.0
RUN ./configure
RUN make
WORKDIR ../..
RUN mkdir tmp

