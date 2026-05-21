# Galaxy Wrapper for Named Entity Geocoding

This Galaxy tool extracts location-type named entities from NLP-annotated JSON and geocodes them to geographic coordinates using Nominatim (OpenStreetMap), enabling interactive map visualization and spatial analysis of text corpora.

## Features

- **Entity type selection**: Geocode GPE, LOC, FAC, and ORG entities from NER output
- **Dual output formats**: GeoJSON for map visualization and TSV for analysis
- **Interactive maps**: GeoJSON viewable directly in Galaxy via OpenLayers plugin
- **Flexible geocoding**: Public OpenStreetMap API or custom Nominatim instances
- **Deduplication**: Automatically handles duplicate entity mentions
- **No dependencies**: Pure Python implementation using only `urllib.request`

## Requirements

- **Input**: JSON output from Galaxy NLP tools (spaCy, Stanza, or CoreNLP) with NER annotations
- **No dependencies**: No external Python packages or model downloads required

## Entity Types Supported

| Type | Description | Geocoding Success Rate |
|---|---|---|
| **GPE** | Geo-Political Entities (countries, cities, states) | High - most reliable |
| **LOC** | Non-GPE locations (mountains, rivers, regions) | Good - geographic features |
| **FAC** | Facilities (buildings, airports, highways) | Variable - depends on specificity |
| **ORG** | Organizations (may resolve if place-named) | Low - use with caution |

## Input Format

The tool expects JSON input with NER annotations:
```json
{
  "sentences": [
    {
      "entities": [
        {
          "text": "Paris",
          "label": "GPE",
          "start_char": 10,
          "end_char": 15
        }
      ]
    }
  ]
}
```

## Output Formats

### GeoJSON Output
Standard GeoJSON FeatureCollection with Point geometries. Each feature includes:
- `name`: Entity text as found in the document
- `entity_type`: NER label (GPE, LOC, FAC, ORG)
- `count`: Number of mentions in the text
- `display_name`: Full address resolved by Nominatim
- `coordinates`: [longitude, latitude] for map plotting

**Interactive Visualization**: Viewable directly in Galaxy using the **OpenLayers** plugin for interactive maps.

### Tabular Summary (TSV)
Tab-separated file with columns:
- `entity`: Entity name
- `type`: Entity type (GPE, LOC, FAC, ORG)
- `latitude`: Decimal degrees
- `longitude`: Decimal degrees  
- `count`: Number of mentions
- `display_name`: Full geocoded address

Results are sorted by mention count (most frequent entities first).

## Nominatim Configuration

### Public OpenStreetMap (Default)
- Uses the public Nominatim API at `nominatim.openstreetmap.org`
- **Rate limited**: 1 request per second
- Suitable for small to medium datasets
- No setup required

### Custom Nominatim Instance
- Configure your own Nominatim server for high-volume geocoding
- No rate limiting
- Better performance for large corpora
- Self-hosting options:
  - Docker: `mediagis/nominatim`
  - Regional extracts from Geofabrik to reduce storage
  - Local OpenStreetMap data processing

## Example Use Cases

- **Literary geography**: Map fictional and real places mentioned in novels
- **Historical research**: Visualize geographic references in archival documents
- **Travel writing analysis**: Track routes and destinations in travel narratives
- **News analysis**: Map geographic focus of news articles over time
- **Digital humanities**: Spatial analysis of large text corpora

## Example Workflow

```
Upload corpus → spaCy NLP (NER, JSON output) → 
Geocode Entities (GPE + LOC) → GeoJSON → 
OpenLayers visualization (interactive map)
```

## Installation

Install this tool from the Galaxy Toolshed: `geocode_entities`

No additional setup required - ready to use with the public Nominatim API.

## Geocoding Notes

- **Disambiguation**: Nominatim returns the most likely match for ambiguous place names
- **Accuracy**: Results depend on entity recognition quality and Nominatim's coverage
- **Rate limits**: Public API is limited to 1 req/sec; consider self-hosting for large datasets
- **Coverage**: Best for well-known places; may fail on very local or historical place names

## Citation

This tool uses OpenStreetMap data via the Nominatim geocoding service. Please cite:

```
OpenStreetMap contributors. OpenStreetMap. 2024. 
https://www.openstreetmap.org
```

