ARG BASE_IMAGE=ubuntu:24.04
FROM ${BASE_IMAGE}

# ARG statement must be defined after FROM statement since FROM resets ARG
# https://docs.docker.jp/engine/reference/builder.html#understand-how-arg-and-from-interact
ARG INSTALL_R=false
ARG INSTALL_GUROBI=false

# Set timezone and locale for R
# ENV DEBIAN_FRONTEND=noninteractive
# ENV TZ=Asia/Tokyo
# RUN apt-get update && \
#     apt-get install -y tzdata && \
#     ln -sf /usr/share/zoneinfo/$TZ /etc/localtime && \
#     echo $TZ > /etc/timezone

# RUN apt-get update && apt-get install -y locales \
#     && locale-gen en_US.UTF-8 \
#     && update-locale LANG=en_US.UTF-8
# ENV LANG=en_US.UTF-8 \
#     LANGUAGE=en_US:en \
#     LC_ALL=en_US.UTF-8

# Install packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ccache \
    curl \
    ca-certificates \
    clang-format \
    doxygen \
    file \
    git \
    graphviz \
    libgraphviz-dev \
    libomp-dev \
    libboost-all-dev \
    libgtest-dev \
    python3 \
    python3-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN curl -fsSL https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Create uv virtual environment
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python packages 
# # pgmpy depends on networkx, numpy, pandas
RUN uv pip install \
    pgmpy==1.0.0 \
    matplotlib \
    graphviz \
    pygraphviz \
    notebook \
    pulp \
    pytest \
    black \
    sphinx \
    breathe
    
# Future work: Install R and packages if needed
# RUN if [ "$INSTALL_R" = "true" ] ; then \
#     wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /etc/apt/trusted.gpg.d/cran.gpg && \
#     echo "deb [signed-by=/etc/apt/trusted.gpg.d/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" > /etc/apt/sources.list.d/cran.list && \
#     apt-get update && apt-get install -y r-base && \
#     R -e "install.packages('bnlearn', repos='https://cloud.r-project.org/')" && \
#     R -e "install.packages('IRkernel', repos='https://cloud.r-project.org/')" && \
#     R -e "IRkernel::installspec()" ; \
#     fi
    
# Future work: Install Gurobi if needed
# RUN if [ "$INSTALL_GUROBI" = "true" ] ; then \
#     wget https://packages.gurobi.com/11.0/gurobi11.0.0_linux64.tar.gz && \
#     tar -xvf gurobi11.0.0_linux64.tar.gz -C /opt && \
#     rm gurobi11.0.0_linux64.tar.gz && \
#     uv pip install gurobipy==11.0.0 && \
#     echo "export GUROBI_HOME=/opt/gurobi1100/linux64" >> ~/.bashrc && \
#     echo 'export PATH=$PATH:$GUROBI_HOME/bin' >> ~/.bashrc && \
#     echo 'export LD_LIBRARY_PATH=$GUROBI_HOME/lib' >> ~/.bashrc ; \
#     fi

WORKDIR /workspace
