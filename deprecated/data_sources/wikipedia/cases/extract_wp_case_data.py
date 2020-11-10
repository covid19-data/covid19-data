import json
import logging
import sys

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO)


class WikiCaseParser:
    def __init__(self, code, page_name, cidx_total_case, cidx_total_death, cidx_total_recovery):
        self.url = (
            "https://en.wikipedia.org/w/api.php?action=query&format=json"
            "&prop=revisions&"
            "titles={}"
            "&formatversion=2&rvprop=content&rvslots=*"
        ).format(page_name)
        if isinstance(code, list):
            self.country_code = code[0]
            self.state_code = code[1]
        else:
            self.country_code = code
        self.cidx_total_case = int(cidx_total_case)
        self.cidx_total_death = int(cidx_total_death)
        self.cidx_total_recovery = int(cidx_total_recovery)

    def download_content(self):
        r = requests.get(self.url)
        return json.loads(r.text)["query"]["pages"][0]["revisions"][0]["slots"]["main"][
            "content"
        ]

    def parse_entry(self, line):
        """Parse one line entry that contains the date, total cases and deaths."""

        def parse_case_str(items, idx):
            s = items[idx].split("/")[0].split("&nbsp;")[0].lstrip("{{#expr:").strip()
            return int(eval(s)) if s else 0

        items = [x.strip() for x in line.split("|")]

        date = items[1]
        total_cases = parse_case_str(items, self.cidx_total_case)
        total_deaths = parse_case_str(items, self.cidx_total_death)
        if self.cidx_total_recovery >= 0:
            total_recoveries = parse_case_str(items, self.cidx_total_recovery)
        else:
            total_recoveries = None

        return date, (total_cases, total_deaths, total_recoveries)

    def get_case_and_death_dict(self, content):
        def is_valid(line):
            line = line.lower()
            if (
                line.startswith("{{medical cases chart/row|")
                or line.startswith("{{bar stacked|")
            ) and line.split("|")[1].strip():
                return True

        return dict(
            self.parse_entry(line) for line in content.splitlines() if is_valid(line)
        )

    def get_confirmed_and_deaths(self):
        data = self.get_case_and_death_dict(self.download_content())
        df = pd.DataFrame.from_dict(
            data, orient="index", columns=["total_cases", "total_deaths", "total_recoveries"]
        )

        df.index = (
            pd.DatetimeIndex(df.index, dayfirst=True)
            if self.country_code == "DNK"
            else pd.DatetimeIndex(df.index, dayfirst=False)
        )

        return (
            df.sort_index()
            .reindex(pd.date_range(df.index[0], df.index[-1]), method="pad")
            .reset_index()
            .rename(columns={"index": "date"})
        )


for idx, row in enumerate(pd.read_csv(snakemake.input[0]).itertuples()):
    df = WikiCaseParser(
        row.country_code, row.page_name, row.cidx_total_case, row.cidx_total_death, row.cidx_total_recovery
    ).get_confirmed_and_deaths()

    columns = ["date", "country_code"]
    df["country_code"] = row.country_code

    if hasattr(row, 'state_code'):
        columns.append("state_code")
        columns.append("state_name")
        df["state_code"] = row.state_code
        df["state_name"] = row.state_name
    else:
        columns.append("country_name")
        df["country_name"] = row.country_name

    columns.append("total_cases")
    columns.append("total_deaths")
    columns.append("total_recoveries")

    df[columns].to_csv(
        snakemake.output[idx], index=False
    )

    logging.info("%s data downloaded & parsed.", row.page_name)
