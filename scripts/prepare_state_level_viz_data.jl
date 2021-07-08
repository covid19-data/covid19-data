using HTTP, CSV, Dates, DataFrames, DataStructures, JSON

NYT_STATE_DATA = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"
OUT_FNAME = string(dirname(pwd()), "/output/us_state_nyt.json")

us_state_abbrev = Dict(
    "Alabama"=> "AL",
    "Alaska"=> "AK",
    "American Samoa"=> "AS",
    "Arizona"=> "AZ",
    "Arkansas"=> "AR",
    "California"=> "CA",
    "Colorado"=> "CO",
    "Connecticut"=> "CT",
    "Delaware"=> "DE",
    "District of Columbia"=> "DC",
    "Florida"=> "FL",
    "Georgia"=> "GA",
    "Guam"=> "GU",
    "Hawaii"=> "HI",
    "Idaho"=> "ID",
    "Illinois"=> "IL",
    "Indiana"=> "IN",
    "Iowa"=> "IA",
    "Kansas"=> "KS",
    "Kentucky"=> "KY",
    "Louisiana"=> "LA",
    "Maine"=> "ME",
    "Maryland"=> "MD",
    "Massachusetts"=> "MA",
    "Michigan"=> "MI",
    "Minnesota"=> "MN",
    "Mississippi"=> "MS",
    "Missouri"=> "MO",
    "Montana"=> "MT",
    "Nebraska"=> "NE",
    "Nevada"=> "NV",
    "New Hampshire"=> "NH",
    "New Jersey"=> "NJ",
    "New Mexico"=> "NM",
    "New York"=> "NY",
    "North Carolina"=> "NC",
    "North Dakota"=> "ND",
    "Northern Mariana Islands"=>"MP",
    "Ohio"=> "OH",
    "Oklahoma"=> "OK",
    "Oregon"=> "OR",
    "Pennsylvania"=> "PA",
    "Puerto Rico"=> "PR",
    "Rhode Island"=> "RI",
    "South Carolina"=> "SC",
    "South Dakota"=> "SD",
    "Tennessee"=> "TN",
    "Texas"=> "TX",
    "Utah"=> "UT",
    "Vermont"=> "VT",
    "Virgin Islands"=> "VI",
    "Virginia"=> "VA",
    "Washington"=> "WA",
    "West Virginia"=> "WV",
    "Wisconsin"=> "WI",
    "Wyoming"=> "WY"
)

name_to_code = x -> us_state_abbrev[x]

function get_df(url) # input is NYT_STATE_DATA
	df_raw = CSV.read(HTTP.get(url).body, DataFrame)
	df = rename(df_raw, 
    :cases => :total_cases, 
    :deaths => :total_deaths,
    :state => :state_name
    )
    df[!, :state_code] = name_to_code.(
    	df.state_name[i] for i in 1:length(df.state_name))
    return df # output is df
end

function normalize_fill_concat(df) # input is df
    #= Normalize each state data by using the same date range. =#
    dfs = []
    date_range = Date(2019, 12, 31):Day(1):maximum(df.date)
    date_range_df = DataFrame(:date => date_range)
    for g in groupby(df, :state_code)
        state_df = sort(leftjoin(date_range_df, g, on = "date"))
        state_df.total_cases = coalesce.(state_df.total_cases, 0)
        state_df.total_deaths = coalesce.(state_df.total_deaths, 0)
        state_df.state_code .= g.state_code[end]
        state_df.state_name .= g.state_name[end]
        state_df.fips .= g.fips[end]
        push!(dfs, state_df[:, 
                [:date, :state_code, :fips, :state_name, :total_cases, :total_deaths]
                ])
    end
    return sort(vcat(dfs..., cols=:union), :state_code)
end

function prepare_data_structure(df, gby = :state_code) # input is concat_df
    data = []
    for g in groupby(df, gby)
        state_data = OrderedDict(
            :code => g.state_code[end],
            :name => g.state_name[end],
            :confirmed => collect(zip(g.date, g.total_cases)),
            :deaths => collect(zip(g.date, g.total_deaths))
            )
        push!(data, state_data)
    end 
    return data 
end

df = get_df(NYT_STATE_DATA)

concat_df = normalize_fill_concat(df)

data = prepare_data_structure(concat_df)

# CSV.write(OUT_FNAME, concat_df)

open(OUT_FNAME, "w") do f
    JSON.print(f, data)
end 




