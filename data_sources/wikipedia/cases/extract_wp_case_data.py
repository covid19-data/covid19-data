import json
import logging

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO)


class WikiCaseParser:
    def __init__(self, country_code, page_name, cidx_total_case, cidx_total_death):
        self.url = (
            "https://en.wikipedia.org/w/api.php?action=query&format=json"
            "&prop=revisions&"
            "titles={}"
            "&formatversion=2&rvprop=content&rvslots=*"
        ).format(page_name)
        self.country_code = country_code
        self.cidx_total_case = int(cidx_total_case)
        self.cidx_total_death = int(cidx_total_death)

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

        return date, (total_cases, total_deaths)

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
            data, orient="index", columns=["total_cases", "total_deaths"]
        )

        df.index = (
            pd.DatetimeIndex(df.index, dayfirst=True)
            if self.country_code == "DNK"
            else pd.DatetimeIndex(df.index, dayfirst=False)
        )

        return (
            df.reindex(pd.date_range(df.index[0], df.index[-1]), method="pad")
            .reset_index()
            .rename(columns={"index": "date"})
        )


for idx, row in enumerate(pd.read_csv(snakemake.input[0]).itertuples()):
    df = WikiCaseParser(
        row.country_code, row.page_name, row.cidx_total_case, row.cidx_total_death
    ).get_confirmed_and_deaths()

    df["country_code"] = row.country_code
    df["country_name"] = row.country_name
    df[["date", "country_code", "country_name", "total_cases", "total_deaths"]].to_csv(
        snakemake.output[idx], index=False
    )
    logging.info("%s data downloaded & parsed.", row.country_code)
