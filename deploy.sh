#if [ -f "dataNineNewRegionAdded.json" ]; 
#	then mv dataNineNewRegionAdded.json output/cntry_stat_owid.json; 
#fi
cd scripts
python prepare_owid_viz_data.py
python nyt_state_data.py
python prepare_state_level_viz_data.py 
cd ..
git add .
msg="updating data on $(date)" 
git commit -m "$msg"
git push


