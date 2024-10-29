FROM ubuntu:24.04
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y autoconf automake bzip2 cmake g++ libboost-all-dev libbz2-dev libcurl4-openssl-dev liblzma-dev make python3 wget zlib1g-dev && \

    # install Google Sparsehash (needed for BioBloom)
    wget -qO- "https://github.com/sparsehash/sparsehash/archive/refs/tags/sparsehash-2.0.4.tar.gz" | tar -zx && \
    cd sparsehash-* && ./configure && make && make install && cd .. && rm -rf sparsehash-* && \

    # install sdsl-lite (needed for BioBloom)
    wget -qO- "https://github.com/simongog/sdsl-lite/releases/download/v2.1.1/sdsl-lite-2.1.1.tar.gz.offline.install.gz" | tar -zx && \
    cd sdsl-lite-* && ./install.sh /usr/local/ && cd .. && rm -rf sdsl-lite-* && \

    # install BioBloom
    wget -qO- "https://github.com/bcgsc/biobloom/releases/download/2.3.5/biobloomtools-2.3.5.tar.gz" | tar -zx && \
    cd biobloomtools-* && sed -i 's/c++11/c++14/g' configure.ac && ./configure && make && make install && cd .. && rm -rf biobloomtools-* && \

    # clean up
    apt-get clean && rm -rf /root/.cache /tmp/*
