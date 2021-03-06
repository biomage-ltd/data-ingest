FROM rocker/tidyverse:4.0.1 AS builder

RUN apt-get update && \
    apt-get install --yes python3-pip

COPY requirements.txt .

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install -r requirements.txt

# DropletUtils needs v57+ to work, but needs to be done manually:
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/matrixStats_0.58.0.tar.gz", repos=NULL, type="source")'
RUN R -e 'chooseCRANmirror(ind=52); install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "DropletUtils"))'

RUN apt-get install --yes libglpk-dev
RUN R -e 'devtools::install_github("eddelbuettel/rcppspdlog")'

RUN apt-get install --yes libxt-dev
RUN R -e 'devtools::install_github("kharchenkolab/conos@v1.3.1")'

RUN R -e 'install.packages("gprofiler2")'
RUN R -e 'install.packages("RJSONIO")'
RUN R -e 'install.packages("Seurat")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("MASS")'
RUN R -e 'install.packages("R.utils")'
RUN R -e 'install.packages("testthat")'
RUN R -e 'BiocManager::install("plger/scDblFinder")'

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
