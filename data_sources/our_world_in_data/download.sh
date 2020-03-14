#!/bin/bash
wget cowid.netlify.com/data/full_data.csv
DIFF=`diff full_data.csv owid_ts.csv`
if [ "$DIFF" != "" ] 
then
    mv full_data.csv owid_ts.csv
    echo "data updated."
else
    rm full_data.csv
    echo "data not yet updated."
fi
