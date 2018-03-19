FROM ubuntu:latest as OS

MAINTAINER Michal Stolarczyk (mjs5kd@virginia.edu)

# Add cran R backport
RUN apt-get -y update
RUN apt-get -y install apt-transport-https
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN echo "deb https://mirrors.ebi.ac.uk/CRAN/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list

# Update & upgrade sources
RUN apt-get -y update
RUN apt-get -y dist-upgrade









#Install other R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('markdown')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('shiny')"
RUN Rscript -e "install.packages('shinythemes')"
RUN Rscript -e "install.packages('shinyBS')"
RUN Rscript -e "install.packages('intergraph')"
RUN Rscript -e "install.packages('visNetwork')"
RUN Rscript -e "install.packages('xtable')"
RUN Rscript -e "install.packages('sna')"
RUN Rscript -e "install.packages('rPython')"



#Clone the git repository
RUN git clone https://github.com/michalstolarczyk/master_thesis.git
