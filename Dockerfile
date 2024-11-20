FROM nvidia/cuda:12.6.2-devel-ubuntu22.04

ARG INSTALL_GUROBI=false

RUN apt update && apt install -y \
    build-essential \
    cmake \
    git \
    libgtest-dev \
    python3.10 \
    python3.10-dev \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set Python3.10 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1

# Build Google Test
RUN cd /usr/src/googletest && \
    cmake CMakeLists.txt && \
    make && \
    cp lib/*.a /usr/lib

# Install Python packages
RUN pip3 install --upgrade pip
RUN pip3 install \
    pandas==2.2.3 \
    networkx==3.4.2 \
    matplotlib==3.9.2 \
    notebook==7.2.2 \
    pulp==2.9.0 \
    pytest==8.3.3

# Install Gurobi if needed
RUN if [ "$INSTALL_GUROBI" = "true" ] ; then \
    apt update && apt install -y wget && \
    wget https://packages.gurobi.com/11.0/gurobi11.0.0_linux64.tar.gz && \
    tar -xvf gurobi11.0.0_linux64.tar.gz -C /opt && \
    rm gurobi11.0.0_linux64.tar.gz && \
    pip3 install gurobipy==11.0.0 && \
    echo "export GUROBI_HOME=/opt/gurobi1100/linux64" >> ~/.bashrc && \
    echo 'export PATH=$PATH:$GUROBI_HOME/bin' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$GUROBI_HOME/lib' >> ~/.bashrc ; \
    fi

WORKDIR /workspace

# docker build --build-arg INSTALL_GUROBI=true -t openbn-image .
# docker run -it --rm --memory=8g --shm-size=4g --cpus=16 --gpus='device=0' -v ${PWD}:/workspace --name openbn-container openbn-image
