#!/bin/bash

mkdir -p db/minikraken-DB

wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz -P db/minikraken-DB
tar xzf db/minikraken-DB/minikraken_8GB_202003.tgz -C db/minikraken-DB
rm minikraken-DB/minikraken_8GB_202003.tgz