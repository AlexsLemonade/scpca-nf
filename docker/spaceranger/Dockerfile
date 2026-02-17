FROM almalinux:9
ARG SPACERANGER_VERSION=4.0.1

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="scpca-nf/spaceranger"
LABEL org.opencontainers.image.description="Docker image to run Space Ranger in scpca-nf"
LABEL org.opencontainers.image.authors="Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation <scpca@ccdatalab.org>"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/scpca-nf/tree/main/docker/spaceranger"

# Install dependencies to use with Nextflow
RUN dnf install -y procps-ng && dnf clean all

# Install Space Ranger
COPY spaceranger-${SPACERANGER_VERSION} /opt/spaceranger-${SPACERANGER_VERSION}
ENV PATH="/opt/spaceranger-${SPACERANGER_VERSION}:${PATH}"

RUN spaceranger telemetry disable

WORKDIR /home
