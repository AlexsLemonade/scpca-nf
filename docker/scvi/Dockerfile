FROM ghcr.io/prefix-dev/pixi:0.83.0 AS build

# change to CPU-only at build time with --build-arg PIXI_ENV=cpu
ARG PIXI_ENV=gpu
# the `gpu` environment requires a CUDA virtual package to solve; this lets
# you build the gpu image on a host with no GPU (e.g. CI runners). The CUDA
# runtime libs are bundled by conda regardless of host, so the resulting env
# imports and runs fine on CPU-only machines too.
ARG CONDA_OVERRIDE_CUDA=12.9

COPY pixi.* /app/
WORKDIR /app
RUN pixi install --locked --environment ${PIXI_ENV}
RUN pixi shell-hook --environment ${PIXI_ENV} > /shell-hook.sh
RUN echo 'exec "$@"' >> /shell-hook.sh


FROM ubuntu:24.04 AS production

ARG PIXI_ENV=gpu

LABEL org.opencontainers.image.title="scpca-nf/scvi"
LABEL org.opencontainers.image.description="scvi-tools Docker image for scpca-nf"
LABEL org.opencontainers.image.authors="Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation <scpca@ccdatalab.org>"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/scpca-nf/tree/main/docker/scvi"

WORKDIR /app
COPY --from=build /app/.pixi/envs/${PIXI_ENV} /app/.pixi/envs/${PIXI_ENV}
COPY --from=build /shell-hook.sh /shell-hook.sh

ENV PATH="/app/.pixi/envs/${PIXI_ENV}/bin:${PATH}"

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["bash"]
