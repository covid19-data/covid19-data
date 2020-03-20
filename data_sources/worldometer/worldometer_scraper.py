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
    for row in rows[:-1]:
        try:
            cols = row.find_all('td')
            #print (cols)
            country = cols[0].text.strip()
            total_cases = cols[1].text.strip()
            total_cases = total_cases.replace(",","")
            print (country,total_cases)
            country_data[country] = total_cases
        except:
            continue

    with open('worldometer_country.json','w') as f:
        f.write(json.dumps(country_data))

get_page()
