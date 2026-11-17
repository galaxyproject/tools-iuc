import argparse
import os
import re
import sys
import time
import unicodedata
from urllib.parse import unquote, urlparse

import requests


OPENALEX_API_BASE = "https://api.openalex.org"
MAX_RETRIES = 3
RETRY_DELAY = 2
WORK_SELECT = "id,title,doi,publication_year,publication_date,cited_by_count,type,open_access,locations,authorships"


def make_request(url, params=None, email=None):
    """Make a request to the OpenAlex API with retry logic and optional polite pool."""
    request_params = dict(params or {})
    if email:
        request_params["mailto"] = email

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, params=request_params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError:
            if response.status_code == 429:
                wait_time = RETRY_DELAY * (attempt + 1) * 2
                print(f"[!] Rate limited, waiting {wait_time}s before retry...")
                time.sleep(wait_time)
            elif response.status_code == 404:
                raise ValueError(f"OpenAlex record not found: {url}")
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


def tsv_value(value):
    if value is None:
        return "N/A"
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def short_openalex_id(value):
    if not value:
        return "N/A"
    return value.rstrip("/").split("/")[-1]


def extract_openalex_id(identifier, prefix):
    match = re.search(rf"\b{prefix}\d+\b", identifier, flags=re.IGNORECASE)
    if match:
        return match.group(0).upper()
    return None


def normalize_doi(identifier):
    value = identifier.strip()
    parsed = urlparse(value)
    if parsed.netloc.lower() in {"doi.org", "dx.doi.org"}:
        value = unquote(parsed.path.lstrip("/"))
    if value.lower().startswith("doi:"):
        value = value[4:]
    match = re.search(r"10\.\d{4,9}/\S+", value, flags=re.IGNORECASE)
    if not match:
        return None
    return match.group(0).rstrip(".,;")


def normalize_orcid(identifier):
    value = identifier.strip()
    parsed = urlparse(value)
    if parsed.netloc.lower() == "orcid.org":
        value = parsed.path.lstrip("/")
    if value.lower().startswith("orcid:"):
        value = value[6:]
    match = re.fullmatch(r"\d{4}-\d{4}-\d{4}-\d{3}[\dXx]", value)
    if match:
        return match.group(0).upper()
    return None


def normalize_ror(identifier):
    value = identifier.strip()
    parsed = urlparse(value)
    if parsed.netloc.lower() == "ror.org":
        return f"https://ror.org/{parsed.path.lstrip('/')}"
    if re.fullmatch(r"0[a-hj-km-np-tv-z|0-9]{6}\d{2}", value, flags=re.IGNORECASE):
        return f"https://ror.org/{value}"
    return None


def search_variants(term):
    variants = [term]
    transliterated = (
        term.replace("ä", "ae")
        .replace("ö", "oe")
        .replace("ü", "ue")
        .replace("Ä", "Ae")
        .replace("Ö", "Oe")
        .replace("Ü", "Ue")
        .replace("ß", "ss")
    )
    ascii_folded = unicodedata.normalize("NFKD", term).encode("ascii", "ignore").decode("ascii")
    for variant in (transliterated, ascii_folded):
        if variant and variant not in variants:
            variants.append(variant)
    return variants


def read_identifiers(args):
    if args.identifiers:
        values = args.identifiers.replace("__cn__", "\n").splitlines()
    else:
        with open(args.identifier_file, encoding="utf-8") as handle:
            values = handle.readlines()

    identifiers = []
    for line in values:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            line = line.split("\t", 1)[0].strip()
        if line:
            identifiers.append(line)
    if not identifiers:
        raise ValueError("No identifiers were provided.")
    return identifiers


def with_xpac(args, params=None):
    params = dict(params or {})
    if args.include_xpac:
        params["include_xpac"] = "true"
    return params


def fetch_work_details(openalex_id, email=None):
    url = f"{OPENALEX_API_BASE}/works/{openalex_id}"
    return make_request(url, email=email)


def search_first(endpoint, search_term, args, entity_name):
    for variant in search_variants(search_term):
        params = with_xpac(args, {"search": variant, "per_page": 1})
        data = make_request(f"{OPENALEX_API_BASE}/{endpoint}", params=params, email=args.email)
        results = data.get("results", [])
        if results:
            return results[0]
    raise ValueError(f"No OpenAlex {entity_name} found for: {search_term}")


def resolve_work(identifier, args):
    work_id = extract_openalex_id(identifier, "W")
    if work_id:
        data = fetch_work_details(work_id, email=args.email)
    else:
        doi = normalize_doi(identifier)
        if doi:
            data = make_request(f"{OPENALEX_API_BASE}/works/https://doi.org/{doi}", email=args.email)
        else:
            data = search_first("works", identifier, args, "work")
    return {
        "input": identifier,
        "entity_id": short_openalex_id(data["id"]),
        "entity_name": data.get("title") or "N/A",
        "data": data,
    }


def resolve_author(identifier, args):
    author_id = extract_openalex_id(identifier, "A")
    orcid = normalize_orcid(identifier)
    if author_id:
        data = make_request(f"{OPENALEX_API_BASE}/authors/{author_id}", email=args.email)
        filter_clause = f"author.id:{author_id}"
    elif orcid:
        data = make_request(f"{OPENALEX_API_BASE}/authors/orcid:{orcid}", email=args.email)
        filter_clause = f"author.orcid:{orcid}"
    else:
        data = search_first("authors", identifier, args, "author")
        author_id = short_openalex_id(data["id"])
        filter_clause = f"author.id:{author_id}"
    return {
        "input": identifier,
        "entity_id": short_openalex_id(data["id"]),
        "entity_name": data.get("display_name") or "N/A",
        "filter": filter_clause,
        "data": data,
    }


def resolve_institution(identifier, args):
    institution_id = extract_openalex_id(identifier, "I")
    ror = normalize_ror(identifier)
    if institution_id:
        data = make_request(f"{OPENALEX_API_BASE}/institutions/{institution_id}", email=args.email)
        filter_clause = f"institutions.id:{institution_id}"
    elif ror:
        data = make_request(f"{OPENALEX_API_BASE}/institutions/ror:{ror}", email=args.email)
        filter_clause = f"institutions.ror:{ror}"
    else:
        data = search_first("institutions", identifier, args, "institution")
        institution_id = short_openalex_id(data["id"])
        filter_clause = f"institutions.id:{institution_id}"
    return {
        "input": identifier,
        "entity_id": short_openalex_id(data["id"]),
        "entity_name": data.get("display_name") or "N/A",
        "filter": filter_clause,
        "data": data,
    }


def build_work_filters(base_filter, args):
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
    if sort_by == "none":
        return works
    reverse = sort_by.endswith(":desc")
    field = sort_by.split(":", 1)[0]
    return sorted(works, key=lambda work: work.get(field) or 0, reverse=reverse)


def fetch_work_list(base_filter, max_works, args):
    works = []
    per_page = 100
    params = with_xpac(
        args,
        {
            "filter": build_work_filters(base_filter, args),
            "per_page": per_page,
            "select": WORK_SELECT,
            "cursor": "*",
        },
    )
    if args.sort != "none":
        params["sort"] = args.sort

    while True:
        data = make_request(f"{OPENALEX_API_BASE}/works", params=params, email=args.email)
        results = data.get("results", [])
        if not results:
            break
        works.extend(results)
        if max_works and len(works) >= max_works:
            return works[:max_works]
        next_cursor = data.get("meta", {}).get("next_cursor")
        if not next_cursor or len(results) < per_page:
            break
        params["cursor"] = next_cursor
    return works


def fetch_referenced_works(openalex_id, max_works, args):
    work_data = fetch_work_details(openalex_id, email=args.email)
    referenced_works = work_data.get("referenced_works", [])
    details = []
    for ref_id in referenced_works:
        try:
            ref_data = fetch_work_details(short_openalex_id(ref_id), email=args.email)
        except ValueError:
            continue
        if filter_work(ref_data, args):
            details.append(ref_data)
            if max_works and len(details) >= max_works:
                break
    return sort_works(details, args.sort)


def get_best_pdf_url(work):
    locations = work.get("locations") or []
    for location in locations:
        if location.get("is_oa") and location.get("pdf_url"):
            return location["pdf_url"]
    primary = work.get("primary_location") or {}
    if primary.get("is_oa") and primary.get("pdf_url"):
        return primary["pdf_url"]
    return None


def safe_pdf_name(index, resolved, work=None):
    work = work or resolved["data"]
    title = work.get("title") or resolved["entity_name"] or resolved["entity_id"]
    safe_title = "".join(char for char in title if char.isalnum() or char in " _-").strip()
    safe_title = safe_title[:160] or short_openalex_id(work.get("id")) or resolved["entity_id"]
    return f"{index:03d}_{short_openalex_id(work.get('id'))}_{safe_title}.pdf"


def download_pdf(url, file_path):
    try:
        response = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
        response.raise_for_status()
        with open(file_path, "wb") as handle:
            handle.write(response.content)
        print(f"[+] Downloaded: {file_path}")
    except Exception as error:
        print(f"[!] Error downloading {url}: {error}")


def work_row(resolved, work):
    return [
        resolved["input"],
        resolved["entity_id"],
        resolved["entity_name"],
        short_openalex_id(work.get("id", "N/A")),
        work.get("title"),
        work.get("doi"),
        work.get("publication_year"),
        work.get("type"),
        work.get("cited_by_count", 0),
        work.get("open_access", {}).get("is_oa", False),
    ]


def write_rows(path, header, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(tsv_value(value) for value in row) + "\n")


def write_summary(resolved_works, output_path):
    rows = []
    for resolved in resolved_works:
        work = resolved["data"]
        rows.append(
            [
                resolved["input"],
                resolved["entity_id"],
                resolved["entity_name"],
                work.get("doi"),
                work.get("publication_year"),
                work.get("type"),
                work.get("cited_by_count", 0),
                work.get("open_access", {}).get("is_oa", False),
            ]
        )
    write_rows(
        output_path,
        ["Input", "Resolved_OpenAlex_ID", "Resolved_Title", "DOI", "Year", "Type", "Cited_By_Count", "Is_OA"],
        rows,
    )


def write_work_list(path, rows):
    write_rows(
        path,
        [
            "Input",
            "Resolved_Entity_ID",
            "Resolved_Entity_Name",
            "Work_OpenAlex_ID",
            "Title",
            "DOI",
            "Year",
            "Type",
            "Cited_By_Count",
            "Is_OA",
        ],
        rows,
    )


def write_authors(resolved_works, output_path):
    rows = []
    for resolved in resolved_works:
        for authorship in resolved["data"].get("authorships", []):
            author = authorship.get("author", {})
            institutions = authorship.get("institutions") or [{}]
            institution = institutions[0]
            rows.append(
                [
                    resolved["input"],
                    resolved["entity_id"],
                    resolved["entity_name"],
                    author.get("display_name"),
                    author.get("orcid"),
                    short_openalex_id(author.get("id", "N/A")),
                    institution.get("display_name", "N/A"),
                    institution.get("country_code", "N/A"),
                ]
            )
    write_rows(
        output_path,
        [
            "Input",
            "Resolved_OpenAlex_ID",
            "Resolved_Title",
            "Name",
            "ORCID",
            "Author_OpenAlex_ID",
            "Affiliation",
            "Country",
        ],
        rows,
    )


def write_concepts(resolved_works, output_path):
    rows = []
    for resolved in resolved_works:
        for concept in resolved["data"].get("concepts", []):
            rows.append(
                [
                    resolved["input"],
                    resolved["entity_id"],
                    resolved["entity_name"],
                    concept.get("display_name"),
                    f"{concept.get('score', 0):.4f}",
                    concept.get("wikidata"),
                ]
            )
    write_rows(
        output_path,
        ["Input", "Resolved_OpenAlex_ID", "Resolved_Title", "Concept", "Score", "Wikidata_ID"],
        rows,
    )


def validate_modes(input_entity, modes):
    work_modes = {"summary", "citing_papers", "references", "related", "authors", "concepts", "download_pdfs"}
    author_modes = {"author_works", "download_pdfs"}
    institution_modes = {"institution_works", "download_pdfs"}
    allowed = {"work": work_modes, "author": author_modes, "institution": institution_modes}[input_entity]
    invalid = sorted(set(modes) - allowed)
    if invalid:
        raise ValueError(f"Mode(s) {', '.join(invalid)} are not valid for {input_entity} inputs.")


def parse_args():
    parser = argparse.ArgumentParser(description="Fetch OpenAlex work, author, and institution data")
    parser.add_argument("--input-entity", choices=["work", "author", "institution"], required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--identifiers", help="One identifier/title/entity per line")
    group.add_argument("--identifier-file", help="Text or tabular file with identifiers in the first non-empty column")
    parser.add_argument(
        "--mode",
        choices=[
            "summary",
            "citing_papers",
            "references",
            "related",
            "authors",
            "concepts",
            "download_pdfs",
            "author_works",
            "institution_works",
        ],
        action="append",
        default=[],
    )
    parser.add_argument("--max-citations", type=str, default="50", dest="max_citations")
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--email", default=None)
    parser.add_argument("--open-access", choices=["any", "true", "false"], default="any")
    parser.add_argument("--publication-year-from", type=int, default=None)
    parser.add_argument("--publication-year-to", type=int, default=None)
    parser.add_argument(
        "--work-type",
        default="any",
        choices=[
            "any",
            "article",
            "book",
            "book-chapter",
            "dataset",
            "dissertation",
            "editorial",
            "letter",
            "paratext",
            "peer-review",
            "preprint",
            "reference-entry",
            "report",
            "review",
        ],
    )
    parser.add_argument(
        "--sort",
        default="none",
        choices=["none", "cited_by_count:desc", "publication_year:desc", "publication_year:asc"],
    )
    parser.add_argument("--include-xpac", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    download_dir = os.path.join(args.output_dir, "downloads")
    modes = args.mode or ["citing_papers"]
    validate_modes(args.input_entity, modes)

    if args.max_citations.lower() == "all":
        max_works = None
    else:
        max_works = int(args.max_citations)

    try:
        if args.publication_year_from and args.publication_year_to and args.publication_year_from > args.publication_year_to:
            raise ValueError("Publication year from must be less than or equal to publication year to.")
        identifiers = read_identifiers(args)
        if args.input_entity == "work":
            resolved_works = [resolve_work(identifier, args) for identifier in identifiers]

            if "summary" in modes:
                write_summary(resolved_works, os.path.join(args.output_dir, "summary.tsv"))

            if "citing_papers" in modes:
                rows = []
                for resolved in resolved_works:
                    works = fetch_work_list(f"cites:{resolved['entity_id']}", max_works, args)
                    rows.extend(work_row(resolved, work) for work in works)
                write_work_list(os.path.join(args.output_dir, "citing_papers.tsv"), rows)

            if "references" in modes:
                rows = []
                for resolved in resolved_works:
                    works = fetch_referenced_works(resolved["entity_id"], max_works, args)
                    rows.extend(work_row(resolved, work) for work in works)
                write_work_list(os.path.join(args.output_dir, "referenced_works.tsv"), rows)

            if "related" in modes:
                rows = []
                for resolved in resolved_works:
                    works = fetch_work_list(f"related_to:{resolved['entity_id']}", max_works, args)
                    rows.extend(work_row(resolved, work) for work in works)
                write_work_list(os.path.join(args.output_dir, "related_works.tsv"), rows)

            if "authors" in modes:
                write_authors(resolved_works, os.path.join(args.output_dir, "authors.tsv"))

            if "concepts" in modes:
                write_concepts(resolved_works, os.path.join(args.output_dir, "concepts.tsv"))

            if "download_pdfs" in modes:
                os.makedirs(download_dir, exist_ok=True)
                for index, resolved in enumerate(resolved_works, start=1):
                    pdf_url = get_best_pdf_url(resolved["data"])
                    if pdf_url:
                        download_pdf(pdf_url, os.path.join(download_dir, safe_pdf_name(index, resolved)))
                    else:
                        print(f"[!] No direct OA PDF URL for: {resolved['entity_name']}")

        elif args.input_entity == "author":
            resolved_authors = [resolve_author(identifier, args) for identifier in identifiers]
            rows = []
            for resolved in resolved_authors:
                works = fetch_work_list(resolved["filter"], max_works, args)
                rows.extend(work_row(resolved, work) for work in works)
                if "download_pdfs" in modes:
                    os.makedirs(download_dir, exist_ok=True)
                    for index, work in enumerate(works, start=1):
                        pdf_url = get_best_pdf_url(work)
                        if pdf_url:
                            download_pdf(pdf_url, os.path.join(download_dir, safe_pdf_name(index, resolved, work)))
            write_work_list(os.path.join(args.output_dir, "author_works.tsv"), rows)

        elif args.input_entity == "institution":
            resolved_institutions = [resolve_institution(identifier, args) for identifier in identifiers]
            rows = []
            for resolved in resolved_institutions:
                works = fetch_work_list(resolved["filter"], max_works, args)
                rows.extend(work_row(resolved, work) for work in works)
                if "download_pdfs" in modes:
                    os.makedirs(download_dir, exist_ok=True)
                    for index, work in enumerate(works, start=1):
                        pdf_url = get_best_pdf_url(work)
                        if pdf_url:
                            download_pdf(pdf_url, os.path.join(download_dir, safe_pdf_name(index, resolved, work)))
            write_work_list(os.path.join(args.output_dir, "institution_works.tsv"), rows)

    except Exception as error:
        print(f"[!] Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
