#!/bin/bash
CURR_DATA="owid_ts.csv"
NEW_DATA="full_data.csv"

wget -q -O $NEW_DATA https://covid.ourworldindata.org/data/ecdc/full_data.csv

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

