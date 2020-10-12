#! /bin/bash

for f in $(ls *.pdf); do pdfcrop "$f" "$f"; done
