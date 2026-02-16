FROM centos:7
LABEL maintainer="ccdl@alexslemonade.org"

# Install spaceranger
WORKDIR /opt
COPY spaceranger-1.3.1.tar.gz spaceranger-1.3.1.tar.gz
RUN tar xzf spaceranger-1.3.1.tar.gz && rm spaceranger-1.3.1.tar.gz
ENV PATH="/opt/spaceranger-1.3.1:${PATH}"

WORKDIR /home
