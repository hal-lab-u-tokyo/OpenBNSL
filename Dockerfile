ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE}

# ARG statement must be defined after FROM statement since FROM resets ARG
# https://docs.docker.jp/engine/reference/builder.html#understand-how-arg-and-from-interact
ARG INSTALL_R=false
ARG INSTALL_GUROBI=false

# Set timezone and locale for R
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo
RUN apt-get update && \
    apt-get install -y tzdata && \
    ln -sf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y locales \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# Install packages
RUN apt update && apt install -y \
    build-essential \
    cmake \
    ccache \
    wget \
    git \
    clang-format \
    doxygen \
    graphviz \
    libgraphviz-dev \
    python3.10 \
    python3.10-dev \
    python3-pip \
    libomp-dev \
    libboost-all-dev \
    libgtest-dev \
    && rm -rf /var/lib/apt/lists/*

# Set Python3.10 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1 && \
    update-alternatives --set python3 /usr/bin/python3.10

# Build Google Test
RUN cd /usr/src/googletest && \
    cmake CMakeLists.txt && \
    make && \
    cp lib/*.a /usr/lib

# Install Python packages 
# pgmpy depends on networkx, numpy, pandas
RUN pip3 install --upgrade pip && pip3 install \
    pgmpy==1.0.0 \
    matplotlib==3.9.2 \
    graphviz==0.21 \
    pygraphviz==1.14 \
    notebook==7.2.2 \
    pulp==2.9.0 \
    pytest==8.3.3 \
    black==24.10.0 \
    sphinx==8.1.3 \
    breathe==4.35.0

# Install R and packages if needed
RUN if [ "$INSTALL_R" = "true" ] ; then \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /etc/apt/trusted.gpg.d/cran.gpg && \
    echo "deb [signed-by=/etc/apt/trusted.gpg.d/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" > /etc/apt/sources.list.d/cran.list && \
    apt-get update && apt-get install -y r-base && \
    R -e "install.packages('bnlearn', repos='https://cloud.r-project.org/')" && \
    R -e "install.packages('IRkernel', repos='https://cloud.r-project.org/')" && \
    R -e "IRkernel::installspec()" ; \
    fi
    
# Install Gurobi if needed
RUN if [ "$INSTALL_GUROBI" = "true" ] ; then \
    wget https://packages.gurobi.com/11.0/gurobi11.0.0_linux64.tar.gz && \
    tar -xvf gurobi11.0.0_linux64.tar.gz -C /opt && \
    rm gurobi11.0.0_linux64.tar.gz && \
    pip3 install gurobipy==11.0.0 && \
    echo "export GUROBI_HOME=/opt/gurobi1100/linux64" >> ~/.bashrc && \
    echo 'export PATH=$PATH:$GUROBI_HOME/bin' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$GUROBI_HOME/lib' >> ~/.bashrc ; \
    fi

WORKDIR /workspace
