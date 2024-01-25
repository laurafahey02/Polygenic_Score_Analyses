FROM ubuntu:18.04
RUN apt-get update && apt-get install -y wget unzip
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=America/New_York
RUN apt update -y && apt install -y software-properties-common && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' && \
    apt-get autoclean && apt-get autoremove && rm -rf /var/lib/apt/lists/* && \
    wget http://se.archive.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb && \
    dpkg -i libpng12-0_1.2.54-1ubuntu1_amd64.deb
RUN apt-get update -y && apt install -y --allow-unauthenticated libpng12-0 libreadline7 r-base-core
RUN wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip && unzip PRSice_linux.zip -d /opt/ && rm PRSice_linux.zip && chmod +x /opt/
RUN R -e "install.packages(c('ggplot2', 'optparse', 'car', 'data.table','methods', 'tools',  'grDevices', 'RColorBrewer', 'Rtsne', 'FactoMineR'))"
ENV PATH="/opt/:${PATH}"
