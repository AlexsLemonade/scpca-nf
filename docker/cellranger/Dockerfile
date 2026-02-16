FROM almalinux:9
LABEL maintainer="ccdl@alexslemonade.org"

# Install dependencies to use with Nextflow
RUN dnf install -y procps-ng && dnf clean all

# Install bcl2fastq
WORKDIR /tmp
COPY bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm bcl2fastq2.rpm
RUN rpm -i bcl2fastq2.rpm

# Install cellranger
COPY cellranger-9.0.1 /opt/cellranger
ENV PATH="/opt/cellranger:${PATH}"

WORKDIR /home
