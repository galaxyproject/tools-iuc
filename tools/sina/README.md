**Add reference libraries**

[Reference databases](https://www.arb-silva.de/download/arb-files/) provided by the [SILVA project](https://www.arb-silva.de/) can be used within this wrapper by adding the corresponding files to the tool-data directory and editing` tool-data/sina_references.loc.sample` as follows: 

    LSU_Parc	${__HERE__}/SILVA_132_LSUParc_12_12_17_opt.arb
    LSU_Ref	${__HERE__}/SILVA_132_LSURef_07_12_17_opt.arb
    SSU_Ref	${__HERE__}/SILVA_132_SSURef_12_12_17_opt.arb
    Ref_NR_99	${__HERE__}/SILVA_132_SSURef_NR99_13_12_17_opt.arb