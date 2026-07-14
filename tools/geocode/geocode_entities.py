#!/usr/bin/env python
# Copyright 2006 The Galaxy Project. All rights reserved.
"""
Geocode Named Entities for Galaxy

Extracts location entities from NLP-annotated JSON (spaCy or Stanza)
and geocodes them using Nominatim (OpenStreetMap) via urllib.
No external dependencies — uses only the Python standard library.

Author: Keith Suderman
License: MIT
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter


NOMINATIM_PUBLIC = "nominatim.openstreetmap.org"


# Force CI refresh - GitHub Actions environment issue
def detect_format(data):
    """Detect whether JSON is spaCy or Stanza format."""
    if "sentences" in data and data["sentences"]:
        first_sent = data["sentences"][0]
        if "tokens" in first_sent:
            return "stanza"
    if "tokens" in data:
        return "spacy"
    return "unknown"


def extract_entities_spacy(data, entity_types):
    """Extract named entities from spaCy JSON."""
    entities = data.get("entities", [])
    location_entities = []
    for ent in entities:
        label = ent.get("label", "")
        if label in entity_types:
            location_entities.append({
                "text": ent.get("text", ""),
                "type": label,
            })
    return location_entities


def extract_entities_stanza(data, entity_types):
    """Extract named entities from Stanza JSON."""
    location_entities = []
    for sent in data.get("sentences", []):
        for ent in sent.get("entities", []):
            ent_type = ent.get("type", ent.get("label", ""))
            if ent_type in entity_types:
                location_entities.append({
                    "text": ent.get("text", ""),
                    "type": ent_type,
                })
    return location_entities


def normalize_entity(text):
    """Normalize entity text for geocoding.

    Strips leading articles and extra whitespace that can confuse Nominatim.
    """
    text = text.strip()
    text = re.sub(r'^(?:the|a|an)\s+', '', text, flags=re.IGNORECASE)
    return text.strip()


def nominatim_geocode(query, domain, user_agent, timeout=10):
    """Geocode a query string using the Nominatim API.

    Returns (latitude, longitude, display_name) or None if not found.
    """
    params = urllib.parse.urlencode({
        "q": query,
        "format": "jsonv2",
        "limit": 1,
        "accept-language": "en",
    })
    scheme = "https" if domain == NOMINATIM_PUBLIC else "http"
    url = f"{scheme}://{domain}/search?{params}"

    req = urllib.request.Request(url, headers={"User-Agent": user_agent})

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            results = json.loads(resp.read())
    except urllib.error.URLError as e:
        print(f"  Network error: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  Error: {e}", file=sys.stderr)
        return None

    if not results:
        return None

    hit = results[0]
    try:
        lat = float(hit["lat"])
        lon = float(hit["lon"])
        display_name = hit.get("display_name", "")
        return lat, lon, display_name
    except (KeyError, ValueError):
        return None


def geocode_entities(entities, domain, user_agent, is_public):
    """Geocode unique entity names and return results with counts."""
    mention_counts = Counter()
    entity_type_map = {}
    original_text_map = {}
    for ent in entities:
        text = ent["text"].strip()
        if text:
            normalized = normalize_entity(text)
            mention_counts[normalized] += 1
            entity_type_map[normalized] = ent["type"]
            original_text_map[normalized] = text

    results = []
    cache = {}
    unique_names = list(mention_counts.keys())
    total = len(unique_names)

    print(f"Geocoding {total} unique entity name(s)...")

    for i, name in enumerate(unique_names):
        original = original_text_map[name]
        if name in cache:
            location = cache[name]
        else:
            location = nominatim_geocode(name, domain, user_agent)

            # Retry once on failure
            if location is None:
                print(f"  Retrying '{name}'...")
                time.sleep(1)
                location = nominatim_geocode(name, domain, user_agent, timeout=15)

            cache[name] = location

            # Rate limit for public Nominatim (1 req/sec)
            if is_public and i < total - 1:
                time.sleep(1.1)

        if location:
            lat, lon, display_name = location
            results.append({
                "text": original,
                "type": entity_type_map[name],
                "latitude": lat,
                "longitude": lon,
                "count": mention_counts[name],
                "display_name": display_name,
            })
            print(f"  [{i + 1}/{total}] {original} -> ({lat:.4f}, {lon:.4f})")
        else:
            print(f"  [{i + 1}/{total}] {original} -> not found")

    return results


def write_geojson(results, output_path):
    """Write results as GeoJSON FeatureCollection."""
    features = []
    for r in results:
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [r["longitude"], r["latitude"]],
            },
            "properties": {
                "name": r["text"],
                "entity_type": r["type"],
                "count": r["count"],
                "display_name": r["display_name"],
            },
        }
        features.append(feature)

    geojson = {
        "type": "FeatureCollection",
        "features": features,
    }

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(geojson, f, indent=2, ensure_ascii=False)


def write_tabular(results, output_path):
    """Write results as tab-separated table."""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["entity", "type", "latitude", "longitude", "count", "display_name"])
        for r in sorted(results, key=lambda x: -x["count"]):
            writer.writerow([
                r["text"],
                r["type"],
                f"{r['latitude']:.6f}",
                f"{r['longitude']:.6f}",
                r["count"],
                r["display_name"],
            ])


def main():
    parser = argparse.ArgumentParser(description="Geocode named entities from NLP JSON")
    parser.add_argument("--input", required=True, help="Input NLP-annotated JSON file")
    parser.add_argument("--geojson", required=True, help="Output GeoJSON file")
    parser.add_argument("--tabular", required=True, help="Output tabular file")
    parser.add_argument("--entity-types", default="GPE,LOC,FAC",
                        help="Comma-separated entity types to geocode")
    parser.add_argument("--nominatim-url", default=None,
                        help="Custom Nominatim domain (e.g. nominatim.example.edu)")
    parser.add_argument("--user-agent", default="galaxy-nlp-geocoder",
                        help="User-Agent string for Nominatim requests")

    args = parser.parse_args()

    entity_types = set(args.entity_types.split(","))

    # Read input JSON
    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading input JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Detect format and extract entities
    fmt = detect_format(data)
    print(f"Detected input format: {fmt}")

    if fmt == "stanza":
        entities = extract_entities_stanza(data, entity_types)
    elif fmt == "spacy":
        entities = extract_entities_spacy(data, entity_types)
    else:
        print("Error: Could not detect JSON format.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(entities)} location entity mention(s)")

    if not entities:
        write_geojson([], args.geojson)
        write_tabular([], args.tabular)
        print("No location entities found. Wrote empty outputs.")
        sys.exit(0)

    # Set up geocoding
    domain = args.nominatim_url if args.nominatim_url else NOMINATIM_PUBLIC
    is_public = args.nominatim_url is None

    # Geocode
    results = geocode_entities(entities, domain, args.user_agent, is_public)

    print(f"\nSuccessfully geocoded {len(results)} of "
          f"{len(set(e['text'] for e in entities))} unique entities")

    # Write outputs
    write_geojson(results, args.geojson)
    write_tabular(results, args.tabular)

    print(f"Wrote GeoJSON to {args.geojson}")
    print(f"Wrote tabular summary to {args.tabular}")


if __name__ == "__main__":
    main()
