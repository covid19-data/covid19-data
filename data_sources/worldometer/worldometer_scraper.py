"""
get country data from worldometers
"""

import urllib3
import json
from bs4 import BeautifulSoup


url = "https://www.worldometers.info/coronavirus/#countries" 

def get_page():
    http = urllib3.PoolManager()
    r = http.request('GET', url)
    soup = BeautifulSoup(r.data, 'html.parser')
    daily_table = soup.find(id='main_table_countries_today')
    rows = daily_table.find_all('tr')
    country_data = {}
    country_data_csv = []
    
    for row in rows[:-1]:
        try:
            cols = row.find_all('td')
            #print (cols)
            country = cols[0].text.strip()
            total_cases = cols[1].text.strip()
            total_cases = total_cases.replace(",","")
            new_cases = cols[2].text.strip()
            new_cases = new_cases.replace(",","")
            new_cases = new_cases.replace("+","")
            td = cols[3].text.strip().replace(",","")
            nd = cols[4].text.strip().replace(",","")
            nd = nd.replace("+","")

            country_data[country] = total_cases
            s = "\t".join([country, total_cases, new_cases, td, nd])
            country_data_csv.append(s)
        except:
            continue

    #with open('worldometer_country.json','w') as f:
    #    f.write(json.dumps(country_data))

    with open('worldometer_country.csv','w') as f:
        h = "\t".join(["country", "total cases", "new cases", "total deaths", "new deaths"])
        f.write(h +'\n')
        for r in country_data_csv:
            f.write(r + '\n')


get_page()

