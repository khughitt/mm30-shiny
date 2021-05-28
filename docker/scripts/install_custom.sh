
#!/bin/bash
#
# install2.r       - cran
# installBioc.r    - bioconductor
# installGithub.r  - github
#

## build ARGs
NCPUS=${NCPUS:-1}

set -e

# install apache arrow
# https://arrow.apache.org/install/
# sudo apt update
# sudo apt install -y -V ca-certificates lsb-release wget
# wget https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
# sudo apt install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
# sudo apt update
# sudo apt install -y -V libarrow-dev

# install r packages
install2.r --error --skipinstalled -r $CRAN -n $NCPUS \
    DT \
    ggdark \
    ggforce \
    ggrepel \
    gridExtra \
    heatmaply \
    Hmisc \
    plotly \
    shinythemes \
    shinycssloaders \
    survival \
    survminer \
    uwot \
    yaml

installGithub.r stephenturner/annotables

rm -rf /tmp/downloaded_packages

