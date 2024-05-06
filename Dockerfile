FROM python:3.9.19-bookworm
workdir app
COPY . .
RUN ./download-conda.sh
RUN ./create-conda-env.sh
