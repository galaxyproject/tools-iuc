import argparse
import os

import requests


# doi
def get_openalex_id_from_doi(doi):
    url = f'https://api.openalex.org/works/https://doi.org/{doi}'
    response = requests.get(url)
    response.raise_for_status()
    return response.json()['id'].split('/')[-1]


# title
def get_openalex_id_from_title(title):
    url = f'https://api.openalex.org/works?search={title}'
    response = requests.get(url)
    response.raise_for_status()
    results = response.json().get('results', [])
    if not results:
        raise ValueError("No paper found with the given title.")
    return results[0]['id'].split('/')[-1]


# fetch papers
def fetch_citing_papers(openalex_id, max_citations=None):
    all_citing_papers = []
    per_page = 200
    page = 1

    work_url = f'https://api.openalex.org/works/{openalex_id}'
    response = requests.get(work_url)
    response.raise_for_status()
    work_data = response.json()

    cited_by_url = work_data.get('cited_by_api_url')
    if not cited_by_url:
        raise ValueError("This work has no citing papers.")

    while True:
        paged_url = f"{cited_by_url}&per_page={per_page}&page={page}"
        response = requests.get(paged_url)
        response.raise_for_status()
        data = response.json()

        results = data.get('results', [])
        if not results:
            break

        all_citing_papers.extend(results)

        if max_citations and len(all_citing_papers) >= max_citations:
            all_citing_papers = all_citing_papers[:max_citations]
            break

        if len(results) < per_page:
            break

        page += 1

    return all_citing_papers


def download_pdf(url, title, folder_name):
    try:
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        response = requests.get(url)
        if response.status_code == 200:
            safe_title = "".join(x for x in title if x.isalnum() or x in " _-").rstrip()
            file_path = os.path.join(folder_name, f"{safe_title}.pdf")
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"[✓] Downloaded: {file_path}")
        else:
            print(f"[x] Failed to download: {url}")
    except Exception as e:
        print(f"[!] Error downloading {url}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Fetch citing papers from OpenAlex")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--id', help='OpenAlex ID of the paper (e.g., W2088676066)')
    group.add_argument('--doi', help='DOI of the paper')
    group.add_argument('--title', help='Title of the paper')

    parser.add_argument('--download', action='store_true', help='Download available OA PDFs')
    parser.add_argument('--max-citations', type=str, default="50", dest='max_citations', help="Max citing papers to fetch or 'all'")
    parser.add_argument('--output-dir', default='.', help='Directory to save output files')
    args = parser.parse_args()

    output_dir = args.output_dir
    summary_path = os.path.join(output_dir, "summary.txt")
    tsv_path = os.path.join(output_dir, "citing_papers.tsv")
    download_dir = os.path.join(output_dir, "downloads")

    if args.max_citations.lower() == "all":
        max_citations = None
    else:
        max_citations = int(args.max_citations)

    try:
        if args.title:
            openalex_id = get_openalex_id_from_title(args.title)
        elif args.doi:
            openalex_id = get_openalex_id_from_doi(args.doi)
        else:
            openalex_id = args.id

        citing_papers = fetch_citing_papers(openalex_id, max_citations=max_citations)

        is_oa = 0
        is_not_oa = 0

        for paper in citing_papers:
            if not paper['locations']:
                continue
            location = paper['locations'][0]
            is_open = location.get('is_oa', False)
            landing_url = location.get('landing_page_url', 'No URL')

            if is_open:
                is_oa += 1
                print("[OA]", landing_url)
                if args.download:
                    pdf_url = location.get('pdf_url')
                    if pdf_url:
                        download_pdf(pdf_url, paper['title'], download_dir)
                    else:
                        print(f"[!] No direct PDF URL for: {paper['title']}")

            else:
                is_not_oa += 1
                print("[Closed]", landing_url)

        print("\nSummary:")
        print("Total citing papers:", len(citing_papers))
        print("Open Access papers:", is_oa)
        print("Closed Access papers:", is_not_oa)

        # save summary
        with open(summary_path, "w") as f:
            f.write(f"Total citing papers: {len(citing_papers)}\n")
            f.write(f"Open Access papers: {is_oa}\n")
            f.write(f"Closed Access papers: {is_not_oa}\n")

        # save  citing papers to a TSV file
        with open(tsv_path, "w", encoding="utf-8") as f:
            f.write("Title\tDOI\tIs_OA\n")
            for paper in citing_papers:
                raw_title = paper.get("title") or "N/A"
                title = raw_title.replace("\t", " ")
                doi = paper.get("doi", "N/A")
                location = paper['locations'][0] if paper['locations'] else {}
                is_oa = location.get("is_oa", False)
                # landing_url = location.get("landing_page_url", "N/A")
                # pdf_url = location.get("pdf_url", "N/A")

                f.write(f"{title}\t{doi}\t{is_oa}\n")

    except Exception as e:
        print(f"[!] Error: {e}")


if __name__ == '__main__':
    main()
