import os
import json
import argparse

def main():
    parser = argparse.ArgumentParser(description="Font Data Manager")
    parser.add_argument("--output", required=True, help="Output JSON file")
    args = parser.parse_args()

    font_data = [
        {
            "value": "SAKURATA",
            "name": "SAKURATA Font",
            "version": "1.0",
            "path": "fonts/SAKURATA.ttf"
        },
        {
            "value": "GermaniaOne-Regular",
            "name": "Germania One Font",
            "version": "1.0",
            "path": "fonts/GermaniaOne-Regular.ttf"
        },
        {
            "value": "Hunting",
            "name": "Hunting Font",
            "version": "1.0",
            "path": "fonts/Hunting.ttf"
        },
        {
            "value": "IranNastaliq",
            "name": "IranNastaliq Font",
            "version": "1.0",
            "path": "fonts/IranNastaliq.ttf"
        },
        {
            "value": "NotoKufiArabic",
            "name": "Noto Kufi Arabic Font",
            "version": "1.0",
            "path": "fonts/NotoKufiArabic-VariableFont_wght.ttf"
        }
    ]

    output_data = {
        "data_tables": {
            "font_data": font_data
        }
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=4)

if __name__ == "__main__":
    main()
