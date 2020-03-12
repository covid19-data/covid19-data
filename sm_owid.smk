rule download_owid:
    output: OWID_TS
    shell: "wget -O {output} cowid.netlify.com/data/full_data.csv"
