FROM dynverse/dynwrapr:v0.1.0

ARG GITHUB_PAT

RUN R -e 'devtools::install_cran("destiny")'

COPY definition.yml run.R example.R /code/

ENTRYPOINT ["/code/run.R"]
