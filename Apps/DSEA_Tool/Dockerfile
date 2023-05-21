FROM rocker/shiny:latest

RUN apt-get update && \
apt-get upgrade -y && \
apt-get clean && \
apt-get install -y git libxml2-dev libmagick++-dev && \
rm -rf /var/lib/apt/lists/*
  
RUN Rscript -e 'options(install.packages.compile.from.source = "always");install.packages(c("BiocManager","shiny","shinyscreenshot","ggplot2","shinydlplot","shinyjs","shinycssloaders","plotly","DT","dplyr"))'
RUN Rscript -e 'BiocManager::install(c("fgsea"),ask = F)'

COPY /DB/ /srv/shiny-server/dseatool/DB 
COPY global.R /srv/shiny-server/dseatool/
COPY server.R /srv/shiny-server/dseatool/
COPY ui.R /srv/shiny-server/dseatool/
  
RUN cd /srv/shiny-server/ && \
sudo chown -R shiny:shiny /srv/shiny-server/dseatool

COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
