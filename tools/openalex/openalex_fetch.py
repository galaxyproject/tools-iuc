import argparse
import os
import sys
import time

import requests


# OpenAlex API base URL
OPENALEX_API_BASE = "https://api.openalex.org"

# Retry settings
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds


def make_request(url, params=None, email=None):
    """Make a request to the OpenAlex API with retry logic and optional polite pool."""
    headers = {}
    request_params = dict(params or {})
    if email:
        request_params["mailto"] = email

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, params=request_params, headers=headers, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                wait_time = RETRY_DELAY * (attempt + 1) * 2
                print(f"[!] Rate limited, waiting {wait_time}s before retry...")
                time.sleep(wait_time)
            elif response.status_code == 404:
                raise ValueError(f"Work not found: {url}")
            elif response.status_code >= 500:
                wait_time = RETRY_DELAY * (attempt + 1)
                print(f"[!] Server error {response.status_code}, retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise
        except requests.exceptions.ConnectionError:
            wait_time = RETRY_DELAY * (attempt + 1)
            print(f"[!] Connection error, retrying in {wait_time}s...")
            time.sleep(wait_time)
        except requests.exceptions.Timeout:
            wait_time = RETRY_DELAY * (attempt + 1)
            print(f"[!] Request timed out, retrying in {wait_time}s...")
            time.sleep(wait_time)

    raise RuntimeError(f"Failed after {MAX_RETRIES} retries for URL: {url}")


def resolve_openalex_id(identifier, id_type, email=None):
    """Resolve an identifier to an OpenAlex work ID."""
    if id_type == "id":
        return identifier
    elif id_type == "doi":
        url = f"{OPENALEX_API_BASE}/works/https://doi.org/{identifier}"
        data = make_request(url, email=email)
        return data['id'].split('/')[-1]
    elif id_type == "title":
        url = f"{OPENALEX_API_BASE}/works"
        params = {"search": identifier}
        data = make_request(url, params=params, email=email)
        results = data.get('results', [])
        if not results:
            raise ValueError("No paper found with the given title.")
        return results[0]['id'].split('/')[-1]
    else:
        raise ValueError(f"Unknown identifier type: {id_type}")


def fetch_work_details(openalex_id, email=None):
    """Fetch full details for a single OpenAlex work."""
    url = f"{OPENALEX_API_BASE}/works/{openalex_id}"
    return make_request(url, email=email)


def build_work_filters(base_filter, args):
    """Build OpenAlex works filter clauses from wrapper options."""
    filters = [base_filter]
    if args.open_access != "any":
        filters.append(f"open_access.is_oa:{args.open_access}")
    if args.publication_year_from:
        filters.append(f"from_publication_date:{args.publication_year_from}-01-01")
    if args.publication_year_to:
        filters.append(f"to_publication_date:{args.publication_year_to}-12-31")
    if args.work_type != "any":
        filters.append(f"type:{args.work_type}")
    return ",".join(filters)


def filter_work(work, args):
    """Apply the same work filters client-side for dereferenced works."""
    if args.open_access != "any":
        is_oa = work.get("open_access", {}).get("is_oa", False)
        if is_oa != (args.open_access == "true"):
            return False
    year = work.get("publication_year")
    if args.publication_year_from and (year is None or year < args.publication_year_from):
        return False
    if args.publication_year_to and (year is None or year > args.publication_year_to):
        return False
    if args.work_type != "any" and work.get("type") != args.work_type:
        return False
    return True


def sort_works(works, sort_by):
    """Sort client-side work lists with the same choices exposed in Galaxy."""
    if sort_by == "none":
        return works
    reverse = sort_by.endswith(":desc")
    field = sort_by.split(":", 1)[0]
    if field == "relevance_score":
        return works
    return sorted(works, key=lambda work: work.get(field) or 0, reverse=reverse)


def fetch_citing_papers(openalex_id, max_citations=None, email=None, args=None):
    """Fetch citing papers for a given OpenAlex work ID."""
    all_citing_papers = []
    per_page = 100
    page = 1

    work_data = fetch_work_details(openalex_id, email=email)
    cited_by_count = work_data.get('cited_by_count', 0)
    if cited_by_count == 0:
        raise ValueError("This work has no citing papers.")

    params = {
        "filter": build_work_filters(f"cites:{openalex_id}", args),
        "per_page": per_page,
        "page": page,
        "select": "id,title,doi,publication_year,publication_date,cited_by_count,type,open_access,locations",
    }
    if args.sort != "none":
        params["sort"] = args.sort
    if args.include_xpac:
        params["include_xpac"] = "true"

    while True:
        params["page"] = page
        data = make_request(f"{OPENALEX_API_BASE}/works", params=params, email=email)

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


def fetch_referenced_works(openalex_id, max_works=None, email=None, args=None):
    """Fetch works referenced by a given OpenAlex work (its bibliography)."""
    work_data = fetch_work_details(openalex_id, email=email)
    referenced_works = work_data.get('referenced_works', [])

    if not referenced_works:
        raise ValueError("This work has no referenced works.")

    details = []
    for ref_id in referenced_works:
        ref_openalex_id = ref_id.split('/')[-1]
        try:
            ref_data = fetch_work_details(ref_openalex_id, email=email)
            if filter_work(ref_data, args):
                details.append(ref_data)
                if max_works and len(details) >= max_works:
                    break
        except ValueError:
            continue

    return sort_works(details, args.sort)


def fetch_related_works(openalex_id, max_works=None, email=None, args=None):
    """Fetch works related to a given OpenAlex work."""
    url = f"{OPENALEX_API_BASE}/works"
    params = {
        "filter": build_work_filters(f"related_to:{openalex_id}", args),
        "per_page": min(max_works or 50, 100),
        "select": "id,title,doi,publication_year,publication_date,cited_by_count,type,open_access,locations",
    }
    if args.sort != "none":
        params["sort"] = args.sort
    if args.include_xpac:
        params["include_xpac"] = "true"
    data = make_request(url, params=params, email=email)
    results = data.get('results', [])

    if not results:
        raise ValueError("No related works found.")

    if max_works:
        results = results[:max_works]

    return results


def fetch_author_details(openalex_id, email=None):
    """Fetch author details for a given OpenAlex work."""
    work_data = fetch_work_details(openalex_id, email=email)
    authorships = work_data.get('authorships', [])

    if not authorships:
        raise ValueError("No author information found for this work.")

    authors = []
    for authorship in authorships:
        author = authorship.get('author', {})
        author_info = {
            'name': author.get('display_name', 'N/A'),
            'orcid': author.get('orcid', 'N/A'),
            'openalex_id': author.get('id', 'N/A'),
            'affiliation': 'N/A',
            'country': 'N/A'
        }
        institutions = authorship.get('institutions', [])
        if institutions:
            inst = institutions[0]
            author_info['affiliation'] = inst.get('display_name', 'N/A')
            author_info['country'] = inst.get('country_code', 'N/A')
        authors.append(author_info)

    return authors


def write_summary(work_data, citing_papers, output_path):
    """Write a summary file."""
    is_oa = 0
    is_not_oa = 0
    for paper in citing_papers:
        locations = paper.get('locations') or []
        if not locations:
            is_not_oa += 1
        else:
            location = locations[0]
            if location.get('is_oa', False):
                is_oa += 1
            else:
                is_not_oa += 1

    with open(output_path, "w") as f:
        f.write(f"Total citing papers: {len(citing_papers)}\n")
        f.write(f"Open Access papers: {is_oa}\n")
        f.write(f"Closed Access papers: {is_not_oa}\n")


def write_citing_papers_tsv(citing_papers, output_path):
    """Write citing papers to a TSV file."""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Title\tDOI\tIs_OA\n")
        for paper in citing_papers:
            raw_title = paper.get("title") or "N/A"
            title = raw_title.replace("\t", " ")
            doi = paper.get("doi", "N/A")
            locations = paper.get('locations') or []
            location = locations[0] if locations else {}
            is_oa_val = location.get("is_oa", False)
            f.write(f"{title}\t{doi}\t{is_oa_val}\n")


def write_referenced_works_tsv(referenced_works, output_path):
    """Write referenced works (bibliography) to a TSV file."""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Title\tDOI\tYear\tIs_OA\n")
        for work in referenced_works:
            raw_title = work.get("title") or "N/A"
            title = raw_title.replace("\t", " ")
            doi = work.get("doi", "N/A")
            year = work.get("publication_year", "N/A")
            is_oa = work.get("open_access", {}).get("is_oa", False)
            f.write(f"{title}\t{doi}\t{year}\t{is_oa}\n")


def write_related_works_tsv(related_works, output_path):
    """Write related works to a TSV file."""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Title\tDOI\tYear\tCited_By_Count\tIs_OA\n")
        for work in related_works:
            raw_title = work.get("title") or "N/A"
            title = raw_title.replace("\t", " ")
            doi = work.get("doi", "N/A")
            year = work.get("publication_year", "N/A")
            cited_by = work.get("cited_by_count", 0)
            is_oa = work.get("open_access", {}).get("is_oa", False)
            f.write(f"{title}\t{doi}\t{year}\t{cited_by}\t{is_oa}\n")


def write_authors_tsv(authors, output_path):
    """Write author details to a TSV file."""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Name\tORCID\tAffiliation\tCountry\n")
        for author in authors:
            name = author['name'].replace("\t", " ")
            orcid = author.get('orcid', 'N/A')
            affiliation = author.get('affiliation', 'N/A').replace("\t", " ")
            country = author.get('country', 'N/A')
            f.write(f"{name}\t{orcid}\t{affiliation}\t{country}\n")


def write_concepts_tsv(work_data, output_path):
    """Write concepts/topics to a TSV file."""
    concepts = work_data.get('concepts', [])
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Concept\tScore\tWikidata_ID\n")
        for concept in concepts:
            name = concept.get('display_name', 'N/A').replace("\t", " ")
            score = f"{concept.get('score', 0):.4f}"
            wikidata = concept.get('wikidata', 'N/A')
            f.write(f"{name}\t{score}\t{wikidata}\n")


def download_pdf(url, title, folder_name):
    """Download a PDF from a URL and save it to the specified folder."""
    try:
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            safe_title = "".join(x for x in title if x.isalnum() or x in " _-").rstrip()
            file_path = os.path.join(folder_name, f"{safe_title}.pdf")
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"[+] Downloaded: {file_path}")
        else:
            print(f"[x] Failed to download: {url}")
    except Exception as e:
        print(f"[!] Error downloading {url}: {e}")


def print_papers(papers, download=False, download_dir="downloads"):
    """Print paper info and optionally download PDFs."""
    is_oa = 0
    is_not_oa = 0
    for paper in papers:
        locations = paper.get('locations') or []
        if not locations:
            is_not_oa += 1
            continue
        location = locations[0]
        is_open = location.get('is_oa', False)
        landing_url = location.get('landing_page_url', 'No URL')

        if is_open:
            is_oa += 1
            print("[OA]", landing_url)
            if download:
                pdf_url = location.get('pdf_url')
                if pdf_url:
                    download_pdf(pdf_url, paper['title'], download_dir)
                else:
                    print(f"[!] No direct PDF URL for: {paper['title']}")
        else:
            is_not_oa += 1
            print("[Closed]", landing_url)

    return is_oa, is_not_oa


def main():
    parser = argparse.ArgumentParser(description="Fetch citing papers from OpenAlex")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--id', help='OpenAlex ID of the paper (e.g., W2088676066)')
    group.add_argument('--doi', help='DOI of the paper')
    group.add_argument('--title', help='Title of the paper')

    parser.add_argument('--mode', choices=['citing_papers', 'summary', 'references', 'related', 'authors', 'concepts'],
                        action='append', default=[],
                        help='What to fetch: citing_papers (default), summary only, references (bibliography), related works, authors, or concepts/topics')
    parser.add_argument('--download', action='store_true', help='Download available OA PDFs (only for citing_papers mode)')
    parser.add_argument('--max-citations', type=str, default="50", dest='max_citations',
                        help="Max papers to fetch or 'all' (default: 50)")
    parser.add_argument('--output-dir', default='.', help='Directory to save output files')
    parser.add_argument('--email', default=None, help='Email for OpenAlex polite pool (faster API access)')
    parser.add_argument('--open-access', choices=['any', 'true', 'false'], default='any',
                        help='Filter works by Open Access status')
    parser.add_argument('--publication-year-from', type=int, default=None,
                        help='Only include works published in or after this year')
    parser.add_argument('--publication-year-to', type=int, default=None,
                        help='Only include works published in or before this year')
    parser.add_argument('--work-type', default='any',
                        choices=['any', 'article', 'book', 'book-chapter', 'dataset', 'dissertation', 'editorial', 'letter', 'paratext', 'peer-review', 'preprint', 'reference-entry', 'report', 'review'],
                        help='Filter works by OpenAlex work type')
    parser.add_argument('--sort', default='none',
                        choices=['none', 'cited_by_count:desc', 'publication_year:desc', 'publication_year:asc'],
                        help='Sort order for work lists')
    parser.add_argument('--include-xpac', action='store_true',
                        help='Include works from Crossref posted-content archives in OpenAlex list queries')

    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    download_dir = os.path.join(output_dir, "downloads")
    modes = args.mode or ['citing_papers']

    if args.max_citations.lower() == "all":
        max_citations = None
    else:
        max_citations = int(args.max_citations)

    try:
        # Resolve identifier type
        if args.title:
            id_type = "title"
            identifier = args.title
        elif args.doi:
            id_type = "doi"
            identifier = args.doi
        else:
            id_type = "id"
            identifier = args.id

        openalex_id = resolve_openalex_id(identifier, id_type, email=args.email)
        work_data = None
        citing_papers = None

        for mode in modes:
            if mode == 'citing_papers' or mode == 'summary':
                if citing_papers is None:
                    citing_papers = fetch_citing_papers(
                        openalex_id,
                        max_citations=max_citations,
                        email=args.email,
                        args=args,
                    )

                if mode == 'citing_papers':
                    is_oa, is_not_oa = print_papers(citing_papers, download=args.download, download_dir=download_dir)
                    print(f"\nSummary:")
                    print(f"Total citing papers: {len(citing_papers)}")
                    print(f"Open Access papers: {is_oa}")
                    print(f"Closed Access papers: {is_not_oa}")

                    if work_data is None:
                        work_data = fetch_work_details(openalex_id, email=args.email)
                    write_summary(work_data, citing_papers, os.path.join(output_dir, "summary.txt"))
                    write_citing_papers_tsv(citing_papers, os.path.join(output_dir, "citing_papers.tsv"))
                else:
                    # summary only mode
                    if work_data is None:
                        work_data = fetch_work_details(openalex_id, email=args.email)
                    cited_by_count = work_data.get('cited_by_count', 0)
                    print(f"Title: {work_data.get('title', 'N/A')}")
                    print(f"DOI: {work_data.get('doi', 'N/A')}")
                    print(f"Year: {work_data.get('publication_year', 'N/A')}")
                    print(f"Cited by: {cited_by_count}")
                    print(f"Type: {work_data.get('type', 'N/A')}")
                    is_oa = work_data.get('open_access', {}).get('is_oa', False)
                    print(f"Open Access: {is_oa}")
                    write_summary(work_data, citing_papers, os.path.join(output_dir, "summary.txt"))
                    write_citing_papers_tsv(citing_papers, os.path.join(output_dir, "citing_papers.tsv"))

            elif mode == 'references':
                referenced = fetch_referenced_works(openalex_id, max_works=max_citations, email=args.email, args=args)
                print(f"\nFound {len(referenced)} referenced works")
                write_referenced_works_tsv(referenced, os.path.join(output_dir, "referenced_works.tsv"))
                print(f"Saved to referenced_works.tsv")

            elif mode == 'related':
                related = fetch_related_works(openalex_id, max_works=max_citations, email=args.email, args=args)
                print(f"\nFound {len(related)} related works")
                write_related_works_tsv(related, os.path.join(output_dir, "related_works.tsv"))
                print(f"Saved to related_works.tsv")

            elif mode == 'authors':
                authors = fetch_author_details(openalex_id, email=args.email)
                print(f"\nFound {len(authors)} authors")
                write_authors_tsv(authors, os.path.join(output_dir, "authors.tsv"))
                print(f"Saved to authors.tsv")

            elif mode == 'concepts':
                if work_data is None:
                    work_data = fetch_work_details(openalex_id, email=args.email)
                concepts = work_data.get('concepts', [])
                print(f"\nFound {len(concepts)} concepts/topics")
                write_concepts_tsv(work_data, os.path.join(output_dir, "concepts.tsv"))
                print(f"Saved to concepts.tsv")

    except Exception as e:
        print(f"[!] Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
