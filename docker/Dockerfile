FROM rocker/shiny:4.1.0

MAINTAINER keith.hughitt@nih.gov

RUN mkdir -p /scripts/

COPY ./scripts/install_tidyverse.sh /scripts/install_tidyverse.sh
COPY ./scripts/install_custom.sh /scripts/install_custom.sh

RUN /scripts/install_tidyverse.sh
RUN /scripts/install_custom.sh
