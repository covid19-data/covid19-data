mv dataNineNewRegionAdded.json output/cntry_stat_owid.json
bash scripts/dataUpdate.sh
git add .
msg="updating data on $(date)" 
git commit -m "$msg"
git push
