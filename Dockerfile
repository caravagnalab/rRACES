FROM r-base:4.4.0

RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y build-essential cmake git \
       ssh g++ r-cran-dplyr r-cran-crayon r-cran-rcpp \
       r-cran-devtools r-cran-ggplot2 r-cran-ggraph \
       r-cran-knitr r-cran-r.utils r-cran-tidygraph \
       r-cran-hexbin r-cran-rcolorbrewer r-cran-cli \
    && apt-get clean

RUN R -e 'install.packages(c("ggmuller"))'

RUN R -e 'devtools::install_github("caravagnalab/rRACES")'