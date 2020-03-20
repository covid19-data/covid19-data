import json
import pandas as pd
import requests

URL_TEMPLATE = (
    "https://en.wikipedia.org/w/api.php?action=query&format=json"
    "&prop=revisions&titles={}&formatversion=2&rvprop=content&rvslots=*"
)

def get_content_from_url(url):
    r = requests.get(url)
    return json.loads(r.text)["query"]["pages"][0]["revisions"][0]["slots"]["main"][
        "content"
    ]


def get_case_and_death_dict(content, cidx_total_case, cidx_total_death):
    data = {}
    for line in content.splitlines():
        line = line.lower()
        if not line.startswith("{{medical cases chart/row|") and not line.startswith("{{bar stacked|"):
            continue
        temp  = [x.strip() for x in line.split("|")]
        if temp[1]:
            if temp[cidx_total_case].startswith("{{#expr:") or temp[cidx_total_death].startswith("{{#expr:"):
                temp[cidx_total_case] = temp[cidx_total_case].split("/")[0].split("&nbsp;")[0].lstrip("{{#expr:").strip() 
                temp[cidx_total_death] = temp[cidx_total_death].split("/")[0].split("&nbsp;")[0].lstrip("{{#expr:").strip()               
            data[temp[1]] = (
                int(eval(temp[cidx_total_case])) if temp[cidx_total_case] else 0,
                int(eval(temp[cidx_total_death])) if temp[cidx_total_death] else 0,
            )                
    return data

def get_confirmed_and_deaths(content, cidx_total_case, cidx_total_death):
    data = get_case_and_death_dict(content, cidx_total_case, cidx_total_death)
    df = pd.DataFrame.from_dict(
        data, orient="index", columns=["total_cases", "total_deaths"]
    )
    df.index = pd.DatetimeIndex(df.index, dayfirst=True)
    return (
        df.reindex(pd.date_range(df.index[0], df.index[-1]), method="pad")
        .reset_index()
        .rename(columns={"index": "date"})
    )

for idx, row in enumerate(pd.read_csv(snakemake.input[0]).itertuples()):
    url = URL_TEMPLATE.format(row.page_name)
    cidx_total_case = int(row.cidx_total_case)
    cidx_total_death = int(row.cidx_total_death)
    content = get_content_from_url(url)
    df = get_confirmed_and_deaths(content, cidx_total_case, cidx_total_death)
    df["country_code"] = row.country_code
    df["country_name"] = row.country_name
    df[["date", "country_code", "country_name", "total_cases", "total_deaths"]].to_csv(
        snakemake.output[idx], index=False
    )
