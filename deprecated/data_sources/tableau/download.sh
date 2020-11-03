#!/bin/bash
CURR_DATA="tableau_ts.csv"
NEW_DATA="temp.csv"

wget -q -O $NEW_DATA "https://docs.google.com/spreadsheets/d/1avGWWl1J19O_Zm0NGTGy2E-fOG05i4ljRfjl87P7FiA/gviz/tq?tqx=out:csv&sheet=COVID-19"

if [ ! -f "$CURR_DATA" ]; then
    if [ -f "$NEW_DATA" ]; then
	mv $NEW_DATA $CURR_DATA 
	echo "data downloaded."
    else
	echo "data download error."
    fi 
else
    DIFF=`diff $NEW_DATA $CURR_DATA`
    if [ "$DIFF" != "" ] 
    then
	mv $NEW_DATA $CURR_DATA
	echo "data updated."
    else
	rm $NEW_DATA
	echo "data not yet updated."
    fi
fi
