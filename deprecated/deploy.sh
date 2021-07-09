#if [ -f "dataNineNewRegionAdded.json" ]; 
#	then mv dataNineNewRegionAdded.json output/cntry_stat_owid.json; 
#fi
julia scripts/prepare_owid_viz_data.jl

julia scripts/prepare_state_level_viz_data.jl

git add .
msg="updating data on $(date)" 
git commit -m "$msg"
git push


