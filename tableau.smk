rule download_tableau:
    output: TABLEAU_TS
    shell: "wget -O {output} https://docs.google.com/spreadsheets/d/1avGWWl1J19O_Zm0NGTGy2E-fOG05i4ljRfjl87P7FiA/gviz/tq?tqx=out:csv&sheet=COVID-19"
