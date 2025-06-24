FROM bioconductor/bioconductor_docker:RELEASE_3_21

MAINTAINER emanuel.sonder@hest.ethz.ch

RUN R -e "if (!requireNamespace('devtools', quietly=TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"

RUN R -e "devtools::install_github('https://github.com/ETHZ-INS/TFBlearner')"
