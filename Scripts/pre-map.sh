#!/usr/bin/env bash

# premap.sh

InputGFA=$1
BuildVG="build-odgi.vg"
ChoppedVG="chopped-odgi.vg"
SortedVG="sorted-odgi.vg"
ViewGFA="view-odgi.gfa"

echo "Using $InputGFA"

odgi build \
     -s \
     -g ${InputGFA} \
     -o ${BuildVG}

echo "Chopping to size 1024"

odgi chop \
     -i build-odgi.vg \
     -c 1024 \
     -o ${ChoppedVG}

echo "Sorting"

odgi sort \
     -i chopped-odgi.vg \
     -o ${SortedVG}

echo "Generating GFA $ViewGFA"

odgi view \
     -i sorted-odgi.vg \
     -g \
     > ${ViewGFA}
