#!/bin/bash
wget -q -O temp.csv "https://docs.google.com/spreadsheets/d/1avGWWl1J19O_Zm0NGTGy2E-fOG05i4ljRfjl87P7FiA/gviz/tq?tqx=out:csv&sheet=COVID-19"
DIFF=`diff temp.csv tableau_ts.csv`
if [ "$DIFF" != "" ] 
then
    mv temp.csv tableau_ts.csv
    echo "data updated."
else
    rm temp.csv
    echo "data not yet updated."
fi
