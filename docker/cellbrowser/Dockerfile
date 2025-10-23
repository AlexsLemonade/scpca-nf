FROM ghcr.io/prefix-dev/pixi:0.57.0 AS build

COPY pixi.* /app/
WORKDIR /app
RUN pixi install --locked
# build shell hook for activation without pixi installed
RUN pixi shell-hook > /shell-hook.sh
RUN echo 'exec "$@"' >> /shell-hook.sh

# use multi-stage build to remove pixi, cache, etc. in final image
FROM ubuntu:24.04 AS production

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="scpca-nf/cellbrowser"
LABEL org.opencontainers.image.description="Docker image to run the UCSC Cell Browser in scpca-nf"
LABEL org.opencontainers.image.authors="Childhood Cancer Data Lab, Alex's Lemonade Stand Foundation <scpca@ccdatalab.org>"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/scpca-nf/tree/main/docker/cellbrowser"

WORKDIR /app
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /shell-hook.sh /shell-hook.sh

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["bash"]
