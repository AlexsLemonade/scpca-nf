FROM ghcr.io/prefix-dev/pixi:0.67.0 AS build

ARG PIXI_ENV=cpu

COPY pixi.* /app/
WORKDIR /app
RUN pixi install --locked --environment ${PIXI_ENV}
RUN pixi shell-hook --environment ${PIXI_ENV} > /shell-hook.sh
RUN echo 'exec "$@"' >> /shell-hook.sh


FROM ubuntu:24.04 AS production

# change to GPU at runtime with --build-arg PIXI_ENV=gpu
ARG PIXI_ENV=cpu

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
