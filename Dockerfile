FROM rocker/tidyverse:4.0.1 AS builder

RUN apt-get update && \
    apt-get install --yes python3-pip

COPY requirements.txt .

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install -r requirements.txt

RUN R -e 'chooseCRANmirror(ind=52); install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db"))'

RUN apt-get install --yes libglpk-dev
RUN R -e 'devtools::install_github("eddelbuettel/rcppspdlog")'
RUN R -e 'devtools::install_github("kharchenkolab/pagoda2")'
RUN R -e 'devtools::install_github("kharchenkolab/conos")'

# ----

FROM builder AS prod
WORKDIR /src
CMD ["/data-ingest/src/docker-entrypoint.sh"]

# ----

FROM builder AS dev

WORKDIR /data-ingest

RUN pip3 install -U jedi radian
RUN R -e 'install.packages("languageserver", repos="http://cran.r-project.org")'

COPY data-ingest.code-workspace .
CMD ["tail", "-f", "/dev/null"]