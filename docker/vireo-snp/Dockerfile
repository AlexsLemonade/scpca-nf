FROM python:3.10.19-slim

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="scpca-nf/vireo-snp"
LABEL org.opencontainers.image.description="Docker image to run the Vireo SNP in scpca-nf"
LABEL org.opencontainers.image.authors="Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation <scpca@ccdatalab.org>"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/scpca-nf/tree/main/docker/vireo-snp"

## Install ps for nextflow
RUN apt-get update \
    && apt-get install -y procps \
    && rm -rf /var/lib/apt/lists/*

# Install vireo using pip
RUN pip install vireoSNP==0.5.7 \
    && pip cache purge

WORKDIR /home
CMD ["bash"]
