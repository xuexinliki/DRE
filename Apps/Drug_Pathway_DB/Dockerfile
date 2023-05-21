FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("BiocManager","shiny","plotly","tidyverse","shinycssloaders","DT","shinydlplot","shinyjs","ggplot2","RColorBrewer","reshape2","pheatmap","stringr","GGally","ggparallel","ggthemes","markdown","devtools"))'
RUN Rscript -e 'BiocManager::install(c("fgsea"),ask = F)'
RUN Rscript -e 'devtools::install_github(repo = "davidsjoberg/ggsankey",force = T)'

COPY /DB/ /srv/shiny-server/drugpathwaydb/DB 
COPY global.R /srv/shiny-server/drugpathwaydb/
COPY server.R /srv/shiny-server/drugpathwaydb/
COPY ui.R /srv/shiny-server/drugpathwaydb/

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/drugpathwaydb

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
