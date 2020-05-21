#!/usr/bin/env bash

# Loop through each gam file
Graphname="graph"
CoverageDir="coverage"

echo "Using ${Graphname}.xg and gam files in ${Filepath}"

while read Filepath;
do
    # paramter expansion http://mywiki.wooledge.org/BashFAQ/073
    Filename=${Filepath##*/}
    FileStub=${Filename%.*}

    echo "packing ${CoverageDir}/${FileStub}.pack.table"
    vg pack \
       -d \
       -x ${Graphname}.xg \
       -g ${Filepath} \
       > ${CoverageDir}/${FileStub}.pack.table
done
