FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("BiocManager","shiny","shinyscreenshot","ggplot2","shinydlplot","shinyjs","fgsea","shinycssloaders","plotly","DT","dplyr","devtools","RColorBrewer","reshape2","stringr","markdown"))'

COPY /DB/ /srv/shiny-server/immunedb/DB 
COPY global.R /srv/shiny-server/immunedb/
COPY server.R /srv/shiny-server/immunedb/
COPY ui.R /srv/shiny-server/immunedb/

RUN cd /srv/shiny-server/ && \
    sudo chown -R shiny:shiny /srv/shiny-server/immunedb

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

ENV SHINY_LOG_STDERR=1

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
