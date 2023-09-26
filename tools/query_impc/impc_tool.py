import sys

import mygene
import pandas as pd
import requests


impc_api_url = "https://www.ebi.ac.uk/mi/impc/bulkdata-api"
impc_api_search_url = f"{impc_api_url}/genes"
impc_api_gene_bundle_url = f"{impc_api_url}/geneBundles"


def stop_err(msg):
    sys.exit(msg)


def main():
    inp = str(sys.argv[1])
    query = str(sys.argv[3])

    try:
        if query == "7":
            g_out = str(sys.argv[5])
            full_gene_table(g_out)
            sys.exit(0)

        if str(sys.argv[5]) == "txt":
            s = str(sys.argv[6])
            if s == "t":
                sep = "\t"
            elif s == "s":
                sep = " "
            elif s in ",;.":
                sep = s
            else:
                sys.exit("Separator not valid, please change it.")
            inp = pd.read_csv(inp, header=None, delimiter=sep)
            if len(inp.columns) == 1:
                inp = inp.to_csv(header=None,
                                 index=False).strip("\n").split("\n")
                inp = ",".join(inp)
            else:
                inp = inp.to_csv(header=None,
                                 index=False).strip(sep).split(sep)
                inp = ",".join(inp)

        if query == "8":
            if str(sys.argv[5]) == "txt":
                g_out = str(sys.argv[7])
            else:
                g_out = str(sys.argv[6])
            genes_in_pipeline(inp, g_out)
            sys.exit(0)
        elif query == "10":
            par_pip_ma(inp)
            sys.exit(0)
        elif query == "11":
            par_gen(inp)
            sys.exit(0)
        elif query == "2" or query == "4":
            final_list = pheno_mapping(inp)
        else:
            final_list = gene_mapping(inp)
        inp = ",".join(final_list)

        if query == "1":
            get_pheno(inp)
            sys.exit(0)
        elif query == "2":
            if str(sys.argv[5]) == "txt":
                g_out = str(sys.argv[7])
            else:
                g_out = str(sys.argv[6])
            get_genes(inp, g_out)
            sys.exit(0)
        elif query == "3":
            gene_set(inp)
            sys.exit(0)
        elif query == "4":
            extr_img(inp)
            sys.exit(0)
        elif query == "5":
            parameters(inp)
            sys.exit(0)
        elif query == "6":
            sign_par(inp)
            sys.exit(0)
        elif query == "9":
            if str(sys.argv[5]) == "txt":
                g_out = str(sys.argv[7])
            else:
                g_out = str(sys.argv[6])
            sign_mp(inp, g_out)
            sys.exit(0)
        else:
            stop_err("Error, non-implemented query selected: " + query)
    except Exception as ex:
        stop_err("Error running impc_tool.py:\n" + str(ex))


# 1-Given a gene id, retrieve all the phenotypes related to it (id and name)
def get_pheno(inp):
    head = sys.argv[4]
    mgi_accession_id = inp

    gene_url = f"{impc_api_search_url}/{mgi_accession_id}"
    gene_data = requests.get(gene_url).json()

    p_list = []
    id_list = []

    if gene_data["significantMpTerms"] is None:
        stop_err("No significant MP terms found for this gene")
    else:
        for x in gene_data["significantMpTerms"]:
            p_list.append(x["mpTermId"])
            id_list.append(x["mpTermName"])

    df = pd.DataFrame()
    df["MP term name"] = p_list
    df["MP term id"] = id_list

    if head == "True":
        df.to_csv(sys.argv[2], header=True, index=False,
                  sep="\t", index_label=False)
    else:
        df.to_csv(sys.argv[2], header=False, index=False,
                  sep="\t", index_label=False)


# 3-Extract all genes having a particular phenotype or a set of phenotypes
# (e.g. relevant to a disease)
def get_genes(inp, g_out):
    head = sys.argv[4]
    target_mp_terms = inp

# All the data is paginated using the page and size parameters,
# by default the endpoint returns the first 20 hits
    gene_by_phenotypes_query = f"{impc_api_search_url}" \
                               f"/search/findAllBySignificantMpTermIdsContains" \
                               f"?mpTermIds={target_mp_terms}&page=0&size=20"
    genes_with_clinical_chemistry_phen = \
        requests.get(gene_by_phenotypes_query).json()
    print(f"Genes with {target_mp_terms}: "
          f"{genes_with_clinical_chemistry_phen['page']['totalElements']}")
    acc = []
    name = []
    url = []

    for gene in genes_with_clinical_chemistry_phen["_embedded"]["genes"]:
        acc.append(gene["mgiAccessionId"])
        name.append(gene["markerName"])
        url.append(gene["_links"]["geneBundle"]["href"])

    if g_out == "sym":
        list_of_genes = pd.DataFrame(columns=["Gene symbol id", "Gene name",
                                              "Gene bundle url"])
        list_of_genes["Gene symbol id"] = mgi_sym_map(acc)
    else:
        list_of_genes = pd.DataFrame(columns=["Gene accession id",
                                              "Gene name", "Gene bundle url"])
        list_of_genes["Gene accession id"] = acc
    list_of_genes["Gene name"] = name
    list_of_genes["Gene bundle url"] = url

    if head == "True":
        list_of_genes.to_csv(sys.argv[2], header=True, index=False,
                             sep="\t", index_label=False)
    else:
        list_of_genes.to_csv(sys.argv[2], header=False, index=False,
                             sep="\t", index_label=False)


# 4. Extract all phenotypes which are present in a particular gene set
# (e.g. genes together in a pathway)
def gene_set(inp):
    head = sys.argv[4]
    target_genes = inp

    genes_in_gene_list_query = f"{impc_api_search_url}/search/" \
                               f"findAllByMgiAccessionIdIn?" \
                               f"mgiAccessionIds={target_genes}"

    genes_in_gene_list = requests.get(genes_in_gene_list_query).json()
    mp_terms_vs_gene_idx = {}

    for gene in genes_in_gene_list["_embedded"]["genes"]:
        mp_terms = gene["significantMpTerms"]
        gene_acc_id = gene["mgiAccessionId"]
        if mp_terms is None:
            continue
        for mp_term_name in mp_terms:
            if mp_term_name["mpTermId"] not in mp_terms_vs_gene_idx:
                mp_terms_vs_gene_idx[mp_term_name["mpTermId"]] = \
                    {"mp_term": mp_term_name["mpTermId"],
                     "mp_name": mp_term_name["mpTermName"], "genes": []}
            mp_terms_vs_gene_idx[mp_term_name["mpTermId"]]["genes"].\
                append(gene_acc_id)
    genes_by_mp_term = list(mp_terms_vs_gene_idx.values())

    df = pd.DataFrame()
    terms = []
    names = []
    genes = []
    for i in genes_by_mp_term:
        terms.append(i["mp_term"])
        names.append(i["mp_name"])
        genes.append(",".join(i["genes"]))

    df["mp_term"] = terms
    df["mp_name"] = names
    df["genes"] = genes

    if head == "True":
        df.to_csv(sys.argv[2], header=True, index=False,
                  sep="\t", index_label=False)
    else:
        df.to_csv(sys.argv[2], header=False, index=False,
                  sep="\t", index_label=False)


# 7. Extract images with a particular phenotype or a set of phenotypes
def extr_img(inp):
    head = sys.argv[4]
    target_mp_terms = inp  # ["MP:0002110", "MP:0000559"]

# All the data is paginated using the page and size parameters,
# by default the endpoint returns the first 20 hits
    gene_by_phenotypes_query = f"{impc_api_search_url}/search/" \
                               f"findAllBySignificantMpTermIdsContains?" \
                               f"mpTermIds={target_mp_terms}&page=0&size=20"
    genes_with_morph_mps = requests.get(gene_by_phenotypes_query).json()
    list_of_gene_bundle_urls = [
        gene["_links"]["geneBundle"]["href"] for gene in
            genes_with_morph_mps["_embedded"]["genes"]
        ]

    gene_bundles = []
    for gene_bundle_url in list_of_gene_bundle_urls:
        gene_bundle = requests.get(gene_bundle_url).json()
        gene_bundles.append(gene_bundle)

    images_with_morphology_mps = []

    # Doing just the first 20 and filtering out fields on the images
    display_fields = ["geneSymbol", "parameterName", "biologicalSampleGroup",
                      "colonyId", "zygosity", "sex", "downloadUrl",
                      "externalSampleId", "thumbnailUrl"]

    for gene_bundle in gene_bundles[:20]:
        if len(gene_bundle) == 4:
            continue
        if gene_bundle["geneImages"] is not None:
            images = gene_bundle["geneImages"]
            for image in images:
                display_image = {k: v for k, v in image.items()
                                 if k in display_fields}
                images_with_morphology_mps.append(display_image)

    images_table = []
    print(f"Images related to phenotype {target_mp_terms}: "
          f"{len(images_with_morphology_mps)}")
    # Displaying just the first 20 images
    for i in images_with_morphology_mps[:20]:
        row = [f"<img src='{i['thumbnailUrl']}' />"] + list(i.values())
        images_table.append(row)

    df = pd.DataFrame()
    externalSampleId = []
    geneSymbol = []
    biologicalSampleGroup = []
    sex = []
    colonyId = []
    zygosity = []
    parameterName = []
    downloadUrl = []
    thumbnailUrl = []

    for i in images_table:
        externalSampleId.append(i[1])
        geneSymbol.append(i[2])
        biologicalSampleGroup.append(i[3])
        sex.append(i[4])
        colonyId.append(i[5])
        zygosity.append(i[6])
        parameterName.append(i[7])
        downloadUrl.append(i[8])
        thumbnailUrl.append(i[9])

    df["externalSampleId"] = externalSampleId
    df["geneSymbol"] = geneSymbol
    df["biologicalSampleGroup"] = biologicalSampleGroup
    df["sex"] = sex
    df["colonyId"] = colonyId
    df["zygosity"] = zygosity
    df["parameterName"] = parameterName
    df["downloadUrl"] = downloadUrl
    df["thumbnailUrl"] = thumbnailUrl

    if head == "True":
        df.to_csv(sys.argv[2], header=True, index=False,
                  sep="\t", index_label=False)
    else:
        df.to_csv(sys.argv[2], header=False, index=False,
                  sep="\t", index_label=False)


# 11- Which parameters have been measured for a particular knockout
def parameters(inp):
    head = sys.argv[4]
    knockout = inp  # "MGI:104636"
    gene_info = requests.get(impc_api_search_url + "/" + knockout).json()

    if gene_info["phenotypingDataAvailable"]:
        geneBundle = requests.get(gene_info["_links"]["geneBundle"]["href"])\
            .json()
        gen_imgs = geneBundle["geneImages"]
        par_list = []
        lis = {}
        for i in gen_imgs:
            lis = {"Parameter Name": i["parameterName"]}
            if lis not in par_list:
                par_list.append(lis)
        df = pd.DataFrame()
        li = []

        for i in par_list:
            li.append(i["Parameter Name"])

        df["Parameter"] = li
        if head == "True":
            df.to_csv(sys.argv[2], header=True, index=False,
                      sep="\t", index_label=False)
        else:
            df.to_csv(sys.argv[2], header=False, index=False,
                      sep="\t", index_label=False)

    else:
        stop_err("No parameters available for this knockout gene")


# 12- Which parameters identified a significant finding for a particular
# knockout line (colony)
def sign_par(inp):
    head = sys.argv[4]
    knockout = inp  # "MGI:104636"

    gene_info = requests.get(f"{impc_api_url}statisticalResults/search/"
                             f"findAllByMarkerAccessionIdIsAndSignificantTrue?"
                             f"mgiAccessionId=" + knockout).json()
    gene_stats = gene_info["_embedded"]["statisticalResults"]

    if len(gene_stats) == 0:
        stop_err("No statistically relevant parameters found "
                 "for this knockout gene")
    else:
        df = pd.DataFrame()
        n = []
        p = []
        for g in gene_stats:
            n.append(g["parameterName"])
            p.append(g["pvalue"])

        df["Parameter name"] = n
        df["p-value"] = p
        if head == "True":
            df.to_csv(sys.argv[2], header=True, index=False,
                      sep="\t", index_label=False)
        else:
            df.to_csv(sys.argv[2], header=False, index=False,
                      sep="\t", index_label=False)


# 13- List of genes names and ID measured in a pipeline
def genes_in_pipeline(inp, g_out):
    head = sys.argv[4]
    pip = inp

    g_in_p_query = f"{impc_api_search_url}/search/" \
                   f"findAllByTestedPipelineId?pipelineId={pip}&" \
                   f"page=0&size=1000"
    genes_in_pip = requests.get(g_in_p_query).json()
    pages = genes_in_pip["page"]["totalPages"]
    max_elem = genes_in_pip["page"]["totalElements"]

    print(f"Genes with {pip}: {genes_in_pip['page']['totalElements']}")
    list_d = []
    acc = []
    name = []

    if max_elem > 1000:
        g_in_p_query = genes_in_pip["_embedded"]["genes"]
        for i in range(1, pages):
            gl = requests.get(f"{impc_api_search_url}/search/"
                              f"findAllByTestedPipelineId?pipelineId={pip}&"
                              f"page={i}&"
                              f"size=1000").json()["_embedded"]["genes"]
            g_in_p_query += gl
    else:
        g_in_p_query = genes_in_pip["_embedded"]["genes"]

    for g in g_in_p_query:
        d = {"Gene Accession ID": g["mgiAccessionId"],
             "Gene Name": g["markerName"]}
        list_d.append(d)

    for i in list_d:
        acc.append(i["Gene Accession ID"])
        name.append(i["Gene Name"])
    if g_out == "sym":
        list_of_genes = pd.DataFrame(columns=["Gene symbol", "Gene name"])
        list_of_genes["Gene symbol"] = mgi_sym_map(acc)
    else:
        list_of_genes = pd.DataFrame(columns=["Gene accession id",
                                              "Gene name"])
        list_of_genes["Gene accession id"] = acc
    list_of_genes["Gene name"] = name

    if head == "True":
        list_of_genes.to_csv(sys.argv[2], header=True, index=False,
                             sep="\t", index_label=False)
    else:
        list_of_genes.to_csv(sys.argv[2], header=False, index=False,
                             sep="\t", index_label=False)


# 14- Extract all genes and corresponding phenotypes related to a
# particular organ system (eg: significatMPTerm)
def sign_mp(inp, g_out):
    head = sys.argv[4]
    mp_term = inp  # ["MP:0005391"]

    gene_by_mpterm_query = f"{impc_api_search_url}/search/" \
                           f"findAllBySignificantMpTermIdsContains?" \
                           f"mpTermIds={mp_term}&size=1000"
    genes_with_mpterm = requests.get(gene_by_mpterm_query).json()

    pages = genes_with_mpterm["page"]["totalPages"]
    genes_info = genes_with_mpterm["_embedded"]["genes"]

    for pn in range(1, pages):
        pq = f"{impc_api_search_url}/search/" \
             f"findAllBySignificantMpTermIdsContains?" \
             f"mpTermIds={mp_term}&page={pn}&size=1000"
        g = requests.get(pq).json()["_embedded"]["genes"]
        genes_info += g

    list_d = []
    d = {}
    for g in genes_info:
        names = []
        ids = []
        for s in g["significantMpTerms"]:
            names.append(s["mpTermName"])
            ids.append(s["mpTermId"])
        d = {"Gene": g["mgiAccessionId"], "mpTermId": ids, "mpTermName": names}
        list_d.append(d)

    g = []
    ids = []
    names = []
    for i in list_d:
        g.append(i["Gene"])
        ids.append(i["mpTermId"])
        names.append(i["mpTermName"])

    df = pd.DataFrame()
    if g_out == "sym":
        df["Gene symbol"] = mgi_sym_map(g)
    else:
        df["Gene Id"] = g
    df["Significant MP terms Ids"] = ids
    df["Significant MP terms Names"] = names

    if head == "True":
        df.to_csv(sys.argv[2], header=True, index=False,
                  sep="\t", index_label=False)
    else:
        df.to_csv(sys.argv[2], header=False, index=False,
                  sep="\t", index_label=False)


# 16- Full table of genes and all identified phenotypes
def full_gene_table(g_out):
    head = sys.argv[4]
    gene_list = requests.get(impc_api_search_url + "?page=0&size=1000").json()
    pages = gene_list["page"]["totalPages"]
    genes_info = gene_list["_embedded"]["genes"]

    for pn in range(1, pages):
        gp = requests.get(impc_api_search_url
                          + f"?page={pn}&"
                            f"size=1000").json()["_embedded"]["genes"]
        genes_info += gp

    d = {}
    list_d = []

    for i in genes_info:
        if i["significantMpTerms"] is None:
            d = {"Gene": i["mgiAccessionId"], "Identified phenotypes": "None"}
        else:
            d = {"Gene": i["mgiAccessionId"],
                 "Identified phenotypes": [
                     sub["mpTermId"] for sub in i["significantMpTerms"]
            ]}
        list_d.append(d)

    df = pd.DataFrame()
    g = []
    p = []
    for i in list_d:
        g.append(i["Gene"])
        p.append(i["Identified phenotypes"])

    if g_out == "sym":
        df["Gene symbol"] = mgi_sym_map(g)
    else:
        df["MGI id"] = g
    df["MP term list"] = p

    for i in range(0, len(df)):
        if df["MP term list"][i] != "None":
            df["MP term list"][i] = str(
                df["MP term list"][i]
            )[1:-1].replace("'", "")

    if str(sys.argv[1]) == "True":
        if head == "True":
            df.to_csv(sys.argv[2], header=True, index=False,
                      sep="\t", index_label=False)
        else:
            df.to_csv(sys.argv[2], header=False, index=False,
                      sep="\t", index_label=False)
    else:
        df = df[df["MP term list"] != "None"]
        df.reset_index(drop=True, inplace=True)
        if head == "True":
            df.to_csv(sys.argv[2], header=True, index=False,
                      sep="\t", index_label=False)
        else:
            df.to_csv(sys.argv[2], header=False, index=False,
                      sep="\t", index_label=False)


# 18- Extract measurements and analysis for a parameter or pipeline
def par_pip_ma(inp):
    head = sys.argv[4]
    id = inp

    if id[0:4] == "IMPC":
        par = True
        ma_query = f"{impc_api_search_url}/search/" \
                   f"findAllByTestedParameterId?" \
                   f"parameterId={id}&page=0&size=1000"
    else:
        ma_query = f"{impc_api_search_url}/search/" \
                   f"findAllByTestedPipelineId?" \
                   f"pipelineId={id}&page=0&size=1000"
        par = False

    ma_in_pip = requests.get(ma_query).json()
    pages = ma_in_pip["page"]["totalPages"]
    max_elem = ma_in_pip["page"]["totalElements"]

    print(f"Genes with {id}: {ma_in_pip['page']['totalElements']}")
    list_d = []
    list_of_genes = pd.DataFrame(columns=["Measurements", "Analysis"])
    mes = []
    an = []

    if max_elem > 1000:

        ma_in_pip = ma_in_pip["_embedded"]["genes"]
        for pn in range(1, pages):
            if par:
                pip = requests.get(f"{impc_api_search_url}/search/"
                                   f"findAllByTestedParameterId?"
                                   f"parameterId={id}&page={pn}&"
                                   f"size=1000").json()["_embedded"]["genes"]
            else:
                pip = requests.get(f"{impc_api_search_url}/search/"
                                   f"findAllByTestedPipelineId?"
                                   f"pipelineId={id}&page={pn}&"
                                   f"size=1000").json()["_embedded"]["genes"]
            ma_in_pip += pip

    else:
        ma_in_pip = ma_in_pip["_embedded"]["genes"]

    for g in ma_in_pip:
        d = {"Measurements": g[""], "Analysis": g[""]}
        list_d.append(d)

    for i in list_d:
        mes.append(i[""])
        an.append(i[""])

    list_of_genes["Analysis"] = an
    list_of_genes["Measurements"] = mes

    if head == "True":
        list_of_genes.to_csv(sys.argv[2], header=True, index=False,
                             sep="\t", index_label=False)
    else:
        list_of_genes.to_csv(sys.argv[2], header=False, index=False,
                             sep="\t", index_label=False)


# 19- Get all genes and measured values for a particular parameter
def par_gen(inp, g_out):
    head = sys.argv[4]
    id = inp

    pa_query = f"{impc_api_search_url}/search/" \
               f"findAllByTestedParameterId?parameterId={id}&page=0&size=1000"

    gm_par = requests.get(pa_query).json()
    pages = gm_par["page"]["totalPages"]
    max_elem = gm_par["page"]["totalElements"]

    print(f"Genes with {id}: {gm_par['page']['totalElements']}")
    list_d = []
    gen = []
    mes = []

    if max_elem > 1000:

        gm_par = gm_par["_embedded"]["genes"]

        for pn in range(1, pages):
            pip = requests.get(f"{impc_api_search_url}/search/"
                               f"findAllByTestedParameterId?"
                               f"parameterId={id}&page={pn}&"
                               f"size=1000").json()["_embedded"]["genes"]
            gm_par += pip

    else:
        gm_par = gm_par["_embedded"]["genes"]

    for g in gm_par:
        d = {"Genes": g["mgiAccessionId"], "Measured Values": g[""]}
        list_d.append(d)

    for i in list_d:
        gen.append(i["Genes"])
        mes.append(i["Measured Values"])

    if g_out == "sym":
        list_of_genes = pd.DataFrame(columns=["Gene symbol",
                                              "Measured Values"])
        list_of_genes["Gene symbol"] = mgi_sym_map(gen)
    else:
        list_of_genes = pd.DataFrame(columns=["Gene accession id",
                                              "Measured Values"])
        list_of_genes["Gene accession id"] = gen
    list_of_genes["Measured Values"] = mes

    if head == "True":
        list_of_genes.to_csv(sys.argv[2], header=True, index=False,
                             sep="\t", index_label=False)
    else:
        list_of_genes.to_csv(sys.argv[2], header=False, index=False,
                             sep="\t", index_label=False)


# Function to map gene symbol to MGI ids
def gene_mapping(inp):
    tmp = inp.split(",")
    final_list = []
    sym_list = []
    for i in tmp:
        if "MGI:" in i:
            final_list.append(i)
        else:
            sym_list.append(i)
    del i

    # symbol for symbols, mgi for MGI :
    # https://docs.mygene.info/en/latest/doc/query_service.html#available-fields
    if len(sym_list) != 0:
        mg = mygene.MyGeneInfo()
        ginfo = mg.querymany(sym_list, scopes="symbol", fields="symbol,MGI",
                             species="mouse")
        empty = True
        discarded = []
        for i in ginfo:
            try:
                final_list.append(i["MGI"])
                empty = False
            except KeyError:
                discarded.append(i["query"])
        if empty and len(final_list) == 0:
            stop_err("Error: it was not possible to map the input.")
        elif empty:
            print("Warning: it was not possible to map any of the symbol ids. "
                  "Only MGI ids will be used.")
        elif len(discarded) != 0:
            print("Warning: it was not possible to map these elements: "
                  "" + ",".join(discarded) + "\n")

    return final_list


# Function to map phenotypes ids to names
def pheno_mapping(inp):
    tmp = inp.split(",")
    final_list = []
    sym_list = []
    for i in tmp:
        if "MP:" in i:
            final_list.append(i)
        else:
            sym_list.append(i)
    del i
    if len(sym_list) != 0:
        url = "https://raw.githubusercontent.com/AndreaFurlani/" \
              "hp_mp_mapping_test/main/hp_mp_mapping.csv"
        mapper = pd.read_csv(url, header=0, index_col=2)
        empty = True
        discarded = []
        for i in sym_list:
            try:
                final_list.append(mapper.loc[i]["mpId"])
                empty = False
            except KeyError:
                discarded.append(i)
                continue
        if empty and len(final_list) == 0:
            stop_err("Error: it was not possible to map the input.")
        elif empty:
            print("Warning: it was not possible to map any of the "
                  "HP term entries. Only MP entries will be used.")
        elif len(discarded) != 0:
            print("Warning: it was not possible to "
                  "map these elements: " + ",".join(discarded) + "\n")
    return final_list


# Function to map MGI ids to Gene Symbols
def mgi_sym_map(mgi_list):
    sym_list = []
    mg = mygene.MyGeneInfo()
    ginfo = mg.querymany(mgi_list, scopes="MGI", fields="symbol,MGI",
                         species="mouse")
    discarded = []
    for i in ginfo:
        try:
            sym_list.append(i["symbol"])
        except KeyError:
            sym_list.append(i["query"])
            discarded.append(i["query"])
    if len(discarded) != 0:
        print("It was not possible to map these genes: " + ",".join(discarded))
    return sym_list


if __name__ == "__main__":
    main()
