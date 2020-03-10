import os

import atoma
import requests


def notify(title, text):
    os.system(
        """osascript -e 'display notification "{}" with title "{}"''""".format(
            text, title
        )
    )


timestamp_lastcommit = atoma.parse_atom_bytes(
    requests.get(
        "https://github.com/CSSEGISandData/COVID-19/commits/master.atom"
    ).content
).updated.timestamp()

if os.path.getmtime("cntry_stat.json") < timestamp_lastcommit:
    notify(
        title="COVID-19 data updated",
        text="There is a new commit in the COVID-19 data repository!",
    )
