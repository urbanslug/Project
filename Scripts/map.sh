#!/usr/bin/env bash

# Loop through each file in the sequence list and call vg to map it to the graph
Graphname="graph"
GAMDir="gams"

while read Filepath;
do
    # paramter expansion http://mywiki.wooledge.org/BashFAQ/073
    Filename=${Filepath##*/}
    FileStub=${Filename%_*}

    echo "Mapping $FileStub"
    vg map \
       -f ${Filepath} \
       -x ${Graphname}.xg \
       -g ${Graphname}.gcsa \
       > ${GAMDir}/${FileStub}.gam
done

echo "Done"
