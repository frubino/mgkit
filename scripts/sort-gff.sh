#!/bin/sh

if [ -z '${ARGS}' ]; then
    ARGS="";
fi

sort -s -k 1,1 -k 7,7 $ARGS $1
