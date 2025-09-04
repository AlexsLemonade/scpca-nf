FROM python:3.10.7-slim
LABEL maintainer = "ccdl@alexslemonade.org"
LABEL org.opencontainers.image.source = "https://github.com/AlexsLemonade/alsf-scpca"

## Install ps for nextflow
RUN apt-get update \
    && apt-get install -y procps \
    && rm -rf /var/lib/apt/lists/*

# Install vireo using pip
RUN pip install vireoSNP==0.5.7 \
    && pip cache purge

WORKDIR /home
