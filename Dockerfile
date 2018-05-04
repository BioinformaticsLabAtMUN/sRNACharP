FROM ubuntu 

LABEL MAINTAINER "Lourdes Peña-Castillo <lourdes@mun.ca>" \
      VERSION "1.0" \
      DESCRIPTION "sRNA-Characterization-Pipeline-nf"

#install neeeded tools
#ADD CRAN to the list of repositories
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

#bedtools v 2.26, R 3.4.4, boost 1.58
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
RUN curl -f https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.5.tar.gz | tar xz && \
   cd ViennaRNA-2.4.5 && \
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
    cd transterm_hp_v2.09 && \
    make clean transterm 

# bprom requires 32-bit libraries
RUN dpkg --add-architecture i386 && apt-get update \
    && apt-get install -y libc6:i386

#install bprom
#RUN curl -fsSL http://linux5.softberry.com/downloads/ea7abecd8f86195f4772082e6bf17ead/bprom.tar.bz2 | tar xvj 
COPY bprom.tar.bz2 /bprom.tar.bz2
RUN tar -xvjf /bprom.tar.bz2 && chmod +x /lin/bprom && rm /bprom.tar.bz2
    
ENV PATH $PATH:/transterm_hp_v2.09:/lin/

ENV TERM_DATA="/transterm_hp_v2.09/expterm.dat"
ENV TSS_DATA="/lin/data"






    






