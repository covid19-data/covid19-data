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


def get_case_and_death_dict(content):
    data = {}
    for line in content.splitlines():
        line = line.lower()
        if not line.startswith("{{medical cases chart/row|") and not line.startswith("{{bar stacked|"):
            continue
        temp = line.split("|")
        temp = [tempitem.strip() for tempitem in temp]
        if temp[1]:
            if temp[4].startswith("{{#expr:"):
                temp[4] = temp[4].split("/")[0].split("&nbsp;")[0].lstrip("{{#expr:").strip()
                temp[2] = temp[2].split("/")[0].split("&nbsp;")[0].lstrip("{{#expr:").strip() 
                data[temp[1]] = (
                    int(temp[2]) if temp[2] else 0,
                    int(temp[4]) if temp[4] else 0,
                )
            else:
                data[temp[1]] = (
		    int(eval(temp[4])) if temp[4] else 0,
                    int(eval(temp[2])) if temp[2] else 0,
                )
    return data


def get_confirmed_and_deaths(content):
    data = get_case_and_death_dict(content)
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
    print(url)
    content = get_content_from_url(url)
    df = get_confirmed_and_deaths(content)
    df["country_code"] = row.country_code
    df["country_name"] = row.country_name
    df[["date", "country_code", "country_name", "total_cases", "total_deaths"]].to_csv(
        snakemake.output[idx], index=False
    )
    print("Done")

