FROM ubuntu 

LABEL MAINTAINER "Lourdes Peña-Castillo <lourdes@mun.ca>" \
      VERSION "1.0" \
      DESCRIPTION "sRNA-Characterization-Pipeline-nf"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y gnupg2
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9


#install neeeded tools
#ADD CRAN to the list of repositories
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

#bedtools v 2.27, R 3.4.4, boost 1.58
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    bedtools \
    libboost-all-dev \
    curl  \
    wget \
    gawk \
    unzip \
    build-essential \
    g++ \
    r-base \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#install ViennaRNA
# -f fail silently
RUN curl -f https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz | tar xz && \
   cd ViennaRNA-2.4.18 && \
   ./configure && \
   make && \
   make install

#install CentroidFold v 0.0.15
# -L follows redirection, -o save into a file
RUN curl -sL https://github.com/satoken/centroid-rna-package/archive/master.zip -o master.zip && \
    unzip master.zip && \
    rm master.zip && \
    cd centroid-rna-package-master && \
    ./configure && make && make install

#install transtermHP 2.09
RUN curl -fsL http://transterm.cbcb.umd.edu/transterm_hp_v2.09.zip -o transterm_hp_v2.09.zip && \
    unzip transterm_hp_v2.09.zip && \
    rm transterm_hp_v2.09.zip && \
    cd transterm_hp_v2.09 
    # make clean transterm

# bprom requires 32-bit libraries
RUN dpkg --add-architecture i386 && apt-get update \
    && apt-get install -y libc6:i386

# install Python3 and its libraries
RUN apt-get update --fix-missing -qq && apt-get -y install -y -q \
    software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y python3 && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-skbio

ENV PATH $PATH:/transterm_hp_v2.09:/lin/

ENV TERM_DATA="/transterm_hp_v2.09/expterm.dat"
ENV TSS_DATA="/lin/data"
