import requests

response = requests.get(
    "https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}".format(
        "1avGWWl1J19O_Zm0NGTGy2E-fOG05i4ljRfjl87P7FiA", "COVID-19"
    )
)
assert response.status_code == 200, "Wrong status code"
open("tableau_ts.csv", "w").write(response.content.decode("ascii"))
