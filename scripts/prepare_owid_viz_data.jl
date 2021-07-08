using DataFrames, CSV, HTTP, Dates, JSON, DataStructures

OUT_FNAME = string(dirname(pwd()), "/output/cntry_stat_owid.json")

# Thanks to https://stackoverflow.com/a/24841669/13716814 & https://stackoverflow.com/a/59473612/13716814

# import OWID raw data
OWID_DATA = "https://covid.ourworldindata.org/data/owid-covid-data.csv"
df_raw = CSV.read(HTTP.get(OWID_DATA).body, 
    DataFrame)[:,[
        "date", "iso_code", "continent", "location", 
        "total_cases", "total_deaths", "population"
        ]]
df = rename(df_raw, :iso_code => :country_code, :location => :country_name)

# import metadata which contains each country/region's corresponding world_region
REGION_DATA = "https://raw.githubusercontent.com/hongtaoh/covid19-data/master/data_sources/metadata/worldbank/country_metadata.csv"
metadata_raw = CSV.read(HTTP.get(REGION_DATA).body, DataFrame)[
    :, ["Country Code", "Region"]
]
metadata = rename(metadata_raw, "Country Code" => :country_code, "Region" => :world_region)

# The following function does the following:
# 1. Normalize each country/region by the same data range
# 2. For each country/region, if the data for total_cases/total_deaths on 2020-01-01 was missing,
#       Change it to 0. This is a necessary procedure for the next step. 
# 3. Change missing values in total_cases & total_deaths to previous available data. 
#       This procedure is called "forward filling".
# 4. Change forward filled data to "missing" for left behind country/region. For example, On 2021-07-02, 
#       Anguilla's most updated data in on 2021-06-25. In the previous stap, forward filling 
#       has filled cases & datahs data for [2021-06-26, 2021-07-02] to be that on 2021-06-25.
#       To be more accurate, we need to change them back to "missing."

function normalize_fill(df) # input should be df
	dfs = [] # initiate the array of dfs
	const_cols = ["country_code", "continent", "country_name", "population"]
	date_range = minimum(df.date):Day(1):maximum(df.date) #get all dates
	# convert all dates into a data frame for leftjoin
	date_range_df = DataFrame(:date => date_range) 
	 # forward filling algorithm
	ffill(a) = a[accumulate(max, [i*!ismissing(a[i]) for i in 1:length(a)], init = 1)]
	for group in groupby(df, :country_code)
		# get the index of the last available date for each country/region
		last_avail_date_index = findfirst(isequal(maximum(group.date)), date_range)
		# normalize each country with the same date rnage, then sort the data frame
		# the sorting will sort by date automatically
		cntry_df = sort(leftjoin(date_range_df, group, on = "date"))
		# if the first case/death is missing, convert it to 0. This is necessary for forward filling
		cntry_df.total_cases[begin] = coalesce.(cntry_df.total_cases[begin], 0)
		cntry_df.total_deaths[begin] = coalesce.(cntry_df.total_deaths[begin], 0)
		# Fill missing values in first row caused by normalization
		# Fill with group[!, i][end] rather than cntry_df[!, i][end]; the latter might be missing itself
		for i in 1:length(const_cols)
			cntry_df[!, i][begin] = coalesce.(cntry_df[!, i][begin], group[!, i][end])
		end
		# Fill missing values with forward filling
		cntry_df = DataFrame([ffill(cntry_df[!, c]) for c in names(cntry_df)], names(cntry_df))
		# Convert unavailable data back to missing
		if last_avail_date_index != length(date_range)
			cntry_df[!, [:total_cases, :total_deaths]][last_avail_date_index+1 : end, :] .= missing
		end 
		push!(dfs, cntry_df)
	end
	return dfs # output is dfs
end

# The following snippet does the following:
#   1. To vertically concatenate dfs, and then merge the concat_df with metadata
#   2. To fill missing values in world_region

function concat_fill(df) #input shouuld be dfs
	concat_df = vcat(df..., cols=:union)
	# wb means world bank where world_region info is obtained
	concat_df_wb = leftjoin(concat_df, metadata, on = "country_code")
	code_region_dict = Dict(
    "AIA" => "Latin America & Caribbean", 
    "BES" => "Latin America & Caribbean", 
    "COK" => "East Asia & Pacific", 
    "FLK" => "Latin America & Caribbean", 
    "GGY" => "Europe & Central Asia", 
    "JEY" => "Europe & Central Asia",
    "OWID_KOS" => "Europe & Central Asia", 
    "MSR" => "Latin America & Caribbean", 
    "OWID_CYN" => "Europe & Central Asia", 
    "PCN" => "East Asia & Pacific",
    "SHN" => "Sub-Saharan Africa", 
    "TWN" => "East Asia & Pacific", 
    "VAT" => "Europe & Central Asia", 
    "WLF" => "East Asia & Pacific", 
    "OWID_WRL" => "World")
	for i in eachrow(concat_df_wb)
		if i.country_code in collect(keys(code_region_dict))
			i.world_region = coalesce.(i.world_region, code_region_dict[i.country_code])
		end
	end
	return concat_df_wb # outout is concat_df_wb
end

# Make each country's data as a dictionary to be used for json conversion
function prepare_data_structure(df) # input is concat_df_wb
	data = []
	for group in groupby(df, :country_code)
		country_data = OrderedDict(
			"country_code" => group.country_code[end],
			"country_name" => group.country_name[end],
			"population" => group.population[end],
			"region" => group.world_region[end],
			"confirmed" => collect(zip(group.date, group.total_cases)),
			"deaths" => collect(zip(group.date, group.total_deaths))
			)
		push!(data, country_data)
	end
	return data
end 

dfs = normalize_fill(df)

concat_df_wb = concat_fill(dfs)

data = prepare_data_structure(concat_df_wb)

open(OUT_FNAME, "w") do f 
	JSON.print(f, data)
end

# alternatively:

# json_string = JSON.json(data)
# open(OUT_FNAME, "w") do f 
# 	write(f, json_string)
# end

