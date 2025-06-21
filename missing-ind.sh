#!/bin/bash

awk 'NR > 1 && $6 > 0.03 {print $1, $2}' onco.imiss > fail-missing-snp.txt


