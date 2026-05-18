#!/bin/bash
set -e

pip3 install --user sphinx sphinx_rtd_theme --break-system-packages

mkdir -p source/_static source/_templates build/html

python3 generate_doc.py

cp conf.py source/

make html
