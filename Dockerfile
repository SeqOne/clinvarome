FROM registry.gitlab.seq.one/devops/dockerfiles/java8-python3-7-htslib:3.7-hts1.10.2
LABEL maintainer='denis.bertrand@seqone.com'
ARG SEQONE_PYPI_REPOSITORY
ARG SEQONE_PYPI_USERNAME
ARG SEQONE_PYPI_PASSWORD
RUN pip install poetry
WORKDIR /root
COPY . .
RUN poetry config repositories.seqone $SEQONE_PYPI_REPOSITORY
RUN poetry config http-basic.seqone $SEQONE_PYPI_USERNAME $SEQONE_PYPI_PASSWORD
RUN poetry install
ENTRYPOINT ["poetry", "run", "./clinvarome/clinvarome_annotation.py"]