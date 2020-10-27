mv dataNineNewRegionAdded.json output/cntry_stat_owid.json
bash scripts/update.sh
git add .
msg="updating data on $(date)" 
git commit -m "$msg"
git push
