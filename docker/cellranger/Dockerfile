FROM almalinux:9
ARG CELLRANGER_VERSION=9.0.1

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="scpca-nf/cellranger"
LABEL org.opencontainers.image.description="Docker image to run Cell Ranger in scpca-nf"
LABEL org.opencontainers.image.authors="Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation <scpca@ccdatalab.org>"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/scpca-nf/tree/main/docker/cellranger"

# Install dependencies to use with Nextflow
RUN dnf install -y procps-ng && dnf clean all

# Install Cell Ranger
COPY cellranger-${CELLRANGER_VERSION} /opt/cellranger-${CELLRANGER_VERSION}
ENV PATH="/opt/cellranger-${CELLRANGER_VERSION}:${PATH}"

RUN cellranger telemetry disable

WORKDIR /home
