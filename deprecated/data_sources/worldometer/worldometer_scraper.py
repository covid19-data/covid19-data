"""
get country data from worldometers
live page from domain and archive from webarchive
"""

import json
from datetime import datetime

import urllib3

from bs4 import BeautifulSoup

live_url = "https://www.worldometers.info/coronavirus/#countries"

# webarchive
url2 = "http://web.archive.org/web/20200319012758/https://www.worldometers.info/coronavirus/"
archive_url = "http://web.archive.org/web/20200318022844/https://www.worldometers.info/coronavirus/"


def get_page(url, sdate, table_name):
    http = urllib3.PoolManager()
    r = http.request("GET", url)
    soup = BeautifulSoup(r.data, "html.parser")
    daily_table = soup.find(id=table_name)
    rows = daily_table.find_all("tr")
    country_data = {}
    country_data_csv = []

    for row in rows[:-1]:
        try:
            cols = row.find_all("td")
            country = cols[0].text.strip()
            total_cases = cols[1].text.strip()
            total_cases = total_cases.replace(",", "")
            new_cases = cols[2].text.strip()
            new_cases = new_cases.replace(",", "")
            new_cases = new_cases.replace("+", "")
            td = cols[3].text.strip().replace(",", "")
            nd = cols[4].text.strip().replace(",", "")
            nd = nd.replace("+", "")

            country_data[country] = total_cases
            s = "\t".join([country, total_cases, new_cases, td, nd])
            country_data_csv.append(s)

        # TODO: except should catch specific exception
        except:
            continue

    # with open('worldometer_country.json','w') as f:
    #    f.write(json.dumps(country_data))

    # TODO: it's probably better to use pandas or csv module.
    with open("worldometer_country_" + sdate + ".csv", "w") as f:
        h = "\t".join(
            ["country", "total cases", "new cases", "total deaths", "new deaths"]
        )
        f.write(h + "\n")
        for r in country_data_csv:
            f.write(r + "\n")


# TODO: Ideally we want to have a way to automatically get historical data.
def get_live():
    sdate = datetime.today().strftime("%Y-%m-%d")
    table_name = "main_table_countries_yesterday"
    get_page(live_url, sdate, table_name)


def get_archive():
    sdate = "20200318022844"
    # table_name = "main_table_countries_yesterday"
    # table_name = 'main_table_countries_today'
    table_name = "main_table_countries"
    get_page(archive_url, sdate, table_name)


get_live()
