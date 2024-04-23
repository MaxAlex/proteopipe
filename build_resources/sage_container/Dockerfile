# The Sage container has stripped out its unzip tool and AWS only provides
# their CLI utility as a zip file, so we instead grab it from another container.
FROM java_components AS aux
FROM ghcr.io/lazear/sage:latest

COPY --from=aux aws aws
RUN ./aws/install

