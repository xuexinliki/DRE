FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("BiocManager","shiny","shinyscreenshot","ggplot2","shinydlplot","shinyjs","shinycssloaders","plotly","DT","dplyr","randomcoloR","data.table","RColorBrewer"))'
RUN Rscript -e 'BiocManager::install(c("qusage","GeneOverlap","GSEABase"),ask = F)'

COPY /DB/ /srv/shiny-server/genesetsimilaritytool/DB 
COPY global.R /srv/shiny-server/genesetsimilaritytool/
COPY server.R /srv/shiny-server/genesetsimilaritytool/
COPY ui.R /srv/shiny-server/genesetsimilaritytool/

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/genesetsimilaritytool

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
