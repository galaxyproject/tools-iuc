# make_EAR_glxy.py
# CAUTION: This is for the Galaxy version!
# by Diego De Panis
# ERGA Sequencing and Assembly Committee
EAR_version = "v24.05.20_glxy_beta"

import sys
import argparse
import logging
import re
import yaml
import os
import pytz
import requests
import json
import glob
import math
from math import ceil
from datetime import datetime
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import cm


def make_report(yaml_file):
    logging.basicConfig(filename='EAR.log', level=logging.INFO)
    # Read the content from EAR.yaml file
    with open(yaml_file, "r") as file:
        yaml_data = yaml.safe_load(file)

    # FUNCTIONS ###################################################################################

    def format_number(value):
        try:
            value_float = float(value)
            if value_float.is_integer():
                # format as an integer if no decimal part
                return f'{int(value_float):,}'
            else:
                # format as a float
                return f'{value_float:,}'
        except ValueError:
            # return the original value if it can't be converted to a float
            return value


    # extract gfastats values
    def extract_gfastats_values(content, keys):
        return [re.findall(f"{key}: (.+)", content)[0] for key in keys]

    keys = [
        "Total scaffold length",
        "GC content %",
        "# gaps in scaffolds",
        "Total gap length in scaffolds",
        "# scaffolds",
        "Scaffold N50",
        "Scaffold L50",
        "Scaffold L90",
        "# contigs",
        "Contig N50",
        "Contig L50",
        "Contig L90",
    ]

    display_names = keys.copy()
    display_names[display_names.index("Total scaffold length")] = "Total bp"
    total_length_index = keys.index("Total scaffold length")
    display_names[display_names.index("GC content %")] = "GC %"
    display_names[display_names.index("Total gap length in scaffolds")] = "Total gap bp"
    display_names[display_names.index("# scaffolds")] = "Scaffolds"
    display_names[display_names.index("# contigs")] = "Contigs"

    gaps_index = keys.index("# gaps in scaffolds")
    exclusion_list = ["# gaps in scaffolds"]


    # extract Total bp from gfastats report
    def extract_total_bp_from_gfastats(gfastats_path):
        with open(gfastats_path, "r") as f:
            content = f.read()
        total_bp = re.search(r"Total scaffold length: (.+)", content).group(1)
        total_bp = int(total_bp.replace(',', ''))
        return "{:,}".format(total_bp)


    # compute EBP quality metric
    def compute_ebp_metric(haplotype, gfastats_path, qv_value):
        keys_needed = ["Contig N50", "Scaffold N50"]
        content = ''
        with open(gfastats_path, "r") as f:
            content = f.read()

        values = extract_gfastats_values(content, keys_needed)
        contig_n50_log = math.floor(math.log10(int(values[0].replace(',', ''))))
        scaffold_n50_log = math.floor(math.log10(int(values[1].replace(',', ''))))
    
        return f"Obtained EBP quality metric for {haplotype}: {contig_n50_log}.{scaffold_n50_log}.Q{math.floor(float(qv_value))}"


    # extract qv values
    def get_qv_value(file_path, order, tool, haplotype):
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                if len(lines) > order and (len(lines) == 1 or lines[2].split('\t')[0].strip() == "Both"):
                    target_line = lines[order]
                    fourth_column_value = target_line.split('\t')[3]
                    return fourth_column_value
        except Exception as e:
            logging.error(f"Error reading {file_path} for tool {tool} and haplotype {haplotype}: {str(e)}")
        return ''


    # extract Kmer completeness values       
    def get_completeness_value(file_path, order, tool, haplotype):
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                if len(lines) > order:
                    target_line = lines[order]
                    fifth_column_value = target_line.split('\t')[4].strip()
                    return fifth_column_value
        except Exception as e:
            logging.warning(f"Error reading {dir_path}: {str(e)}")
            return ''


    # Getting kmer plots for curated asm
    def get_png_files(dir_path):
        png_files = glob.glob(f"{dir_path}/*.ln.png")
        if len(png_files) < 4:
            logging.warning(f"Warning: Less than 4 png files found in {dir_path}. If this is diploid, some images may be missing.")
            # fill missing with None
            while len(png_files) < 4:
                png_files.append(None)
        return png_files[:4]


    # get unique part in file names
    def find_unique_parts(file1, file2):
        # Split filenames into parts
        parts1 = file1.split('.')
        parts2 = file2.split('.')
        # Find unique parts
        unique_parts1 = [part for part in parts1 if part not in parts2]
        unique_parts2 = [part for part in parts2 if part not in parts1]
    
        return ' '.join(unique_parts1), ' '.join(unique_parts2)

    
    # extract BUSCO values
    def extract_busco_values(file_path):
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                results_line = re.findall(r"C:.*n:\d+", content)[0]
                s_value = re.findall(r"S:(\d+\.\d+%)", results_line)[0]
                d_value = re.findall(r"D:(\d+\.\d+%)", results_line)[0]
                f_value = re.findall(r"F:(\d+\.\d+%)", results_line)[0]
                m_value = re.findall(r"M:(\d+\.\d+%)", results_line)[0]
                return s_value, d_value, f_value, m_value
        except Exception as e:
            logging.warning(f"Error reading {file_path}: {str(e)}")
            return '', '', '', ''


    # extract BUSCO info
    def extract_busco_info(file_path):
        busco_version = None
        lineage_info = None
    
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                version_match = re.search(r"# BUSCO version is: ([\d.]+)", content)
                if version_match:
                    busco_version = version_match.group(1)
                lineage_match = re.search(r"The lineage dataset is: (.*?) \(Creation date:.*?, number of genomes: (\d+), number of BUSCOs: (\d+)\)", content)
                if lineage_match:
                    lineage_info = lineage_match.groups()
                if not lineage_info:
                    lineage_match = re.search(r"The lineage dataset is: (.*?) \(Creation date:.*?, number of species: (\d+), number of BUSCOs: (\d+)\)", content)
                    if lineage_match:
                        lineage_info = lineage_match.groups()
    
        except Exception as e:
            logging.warning(f"Error reading {file_path}: {str(e)}")
    
        return busco_version, lineage_info


    # Function to check and generate warning messages
    def generate_warning_paragraphs(expected, observed, trait):
        paragraphs = []
        try:
            if trait == "Haploid size (bp)":
                expected_val = int(expected.replace(',', ''))
                observed_val = int(observed.replace(',', ''))
                if abs(expected_val - observed_val) / expected_val > 0.20:
                    message = f". Observed {trait} has >20% difference with Expected"
                    paragraphs.append(Paragraph(message, styles["midiStyle"]))
            elif trait in ["Haploid Number", "Ploidy"]:
                # Ensure both values are integers for comparison
                expected_val = int(expected)
                observed_val = int(observed)
                if expected_val != observed_val:
                    message = f". Observed {trait} is different from Expected"
                    paragraphs.append(Paragraph(message, styles["midiStyle"]))
            elif trait == "Sample Sex":
                # Compare case-insensitive and trimmed strings
                if expected.strip().lower() != observed.strip().lower():
                    message = f". Observed sex is different from Sample sex"
                    paragraphs.append(Paragraph(message, styles["midiStyle"]))
        except Exception as e:
            logging.warning(f"Error in generating warning for {trait}: {str(e)}")
    
        return paragraphs


    # Generate warnings for curated haplotypes (qv, kcomp, busco)
    def generate_curated_warnings(haplotype, qv_value, completeness_value, busco_scores):
        paragraphs = []
        try:
            # Ensure values are correctly interpreted as floats
            qv_val = float(qv_value)
            completeness_val = float(completeness_value)
            s_value = float(busco_scores[0].rstrip('%'))
            d_value = float(busco_scores[1].rstrip('%'))

            # Check QV value
            if qv_val < 40:
                message = f". QV value is less than 40 for {haplotype}"
                paragraphs.append(Paragraph(message, styles["midiStyle"]))

            # Check Kmer completeness value
            if completeness_val < 90:
                message = f". Kmer completeness value is less than 90 for {haplotype}"
                paragraphs.append(Paragraph(message, styles["midiStyle"]))

            # Check BUSCO s_value
            if s_value < 90:
                message = f". BUSCO single copy value is less than 90% for {haplotype}"
                paragraphs.append(Paragraph(message, styles["midiStyle"]))

            # Check BUSCO d_value
            if d_value > 5:
                message = f". BUSCO duplicated value is more than 5% for {haplotype}"
                paragraphs.append(Paragraph(message, styles["midiStyle"]))

        except Exception as e:
            logging.warning(f"Error in generating warnings for {haplotype}: {str(e)}")

        return paragraphs


    # Generate warnings for curated haplotypes (loss, gaps, 90inChrom) 
    def generate_assembly_warnings(asm_data, gaps_per_gbp_data, obs_haploid_num):
        warnings = []

        # Iterate over haplotypes and generate warnings based on the criteria
        for haplotype in asm_stages:
            pre_curation_bp = extract_total_bp_from_gfastats(asm_data['Pre-curation'][haplotype]['gfastats--nstar-report_txt'])
            curated_bp = extract_total_bp_from_gfastats(asm_data['Curated'][haplotype]['gfastats--nstar-report_txt'])
            scaffold_l90 = float(gfastats_data[('Curated', haplotype)][display_names.index('Scaffold L90')].replace(',', ''))

            # Check for assembly length loss > 3%
            if pre_curation_bp and curated_bp:
                loss_percentage = (float(pre_curation_bp.replace(',', '')) - float(curated_bp.replace(',', ''))) / float(pre_curation_bp.replace(',', '')) * 100
                if loss_percentage > 3:
                    warnings.append(Paragraph(f". Assembly length loss > 3% for {haplotype}", styles["midiStyle"]))

            # Check for more than 1000 gaps/Gbp
            gaps_gbp = gaps_per_gbp_data.get(('Curated', haplotype), 0)
            if gaps_gbp > 1000:
                warnings.append(Paragraph(f". More than 1000 gaps/Gbp for {haplotype}", styles["midiStyle"]))

            # Check if Scaffold L90 value is more than Observed Haploid number
            if scaffold_l90 > float(obs_haploid_num):
                warnings.append(Paragraph(f". Not 90% of assembly in chromosomes for {haplotype}", styles["midiStyle"]))

        return warnings

    # Parse pipeline and generate "tree"
    def generate_pipeline_tree(pipeline_data):
        tree_lines = []
        indent = "&nbsp;" * 2  # Adjust indent spacing as needed

        for tool_version_param in pipeline_data:
            parts = tool_version_param.split('|')
            tool_version = parts[0]
            tool, version = tool_version.split('_v') if '_v' in tool_version else (tool_version, "NA")
            
            # Handle parameters: join all but the first (which is tool_version) with ', '
            param_text = ', '.join(parts[1:]) if len(parts) > 1 else "NA"

            # Tool line
            tool_line = f"- <b>{tool}</b>"
            tree_lines.append(tool_line)

            # Version line
            version_line = f"{indent*2}|_ <i>ver:</i> {version}"
            tree_lines.append(version_line)

            # Param line(s)
            if param_text != "NA":
                for param in param_text.split(','):
                    param = param.strip()
                    param_line = f"{indent*2}|_ <i>key param:</i> {param if param else 'NA'}"
                    tree_lines.append(param_line)
            else:
                param_line = f"{indent*2}|_ <i>key param:</i> NA"
                tree_lines.append(param_line)

        # Join lines with HTML break for paragraph
        tree_diagram = "<br/>".join(tree_lines)
        return tree_diagram

    

    # Reading SAMPLE INFORMATION section from yaml ################################################

    # Check for required fields
    required_fields = ["ToLID", "Species", "Sex", "Submitter", "Affiliation", "Tags"]
    missing_fields = [field for field in required_fields if field not in yaml_data or not yaml_data[field]]

    if missing_fields:
        logging.error(f"# GENERAL INFORMATION section in the yaml file is missing or empty for the following information: {', '.join(missing_fields)}")
        sys.exit(1)

    # Check that "Species" field is a string
    if not isinstance(yaml_data["Species"], str):
        logging.error(f"# GENERAL INFORMATION section in the yaml file contains incorrect data type for 'Species'. Expected 'str' but got '{type(yaml_data['Species']).__name__}'.")
        sys.exit(1)


    ## Get data for Header, ToLID table and submitter
    tol_id = yaml_data["ToLID"]
    species = yaml_data["Species"]
    sex = yaml_data["Sex"]
    submitter = yaml_data["Submitter"]
    affiliation = yaml_data["Affiliation"]
    tags = yaml_data["Tags"]

    # Check if tag is valid
    valid_tags = ["ERGA-BGE", "ERGA-Pilot", "ERGA-Satellite"]
    if tags not in valid_tags:
        tags += "[INVALID TAG]"
        logging.warning(f"# SAMPLE INFORMATION section in the yaml file contains an invalid tag. Valid tags are ERGA-BGE, ERGA-Pilot and ERGA-Satellite")


    # Get data from GoaT based on species name
    # urllib.parse.quote to handle special characters and spaces in the species name
    species_name = requests.utils.quote(species)

    # Get stuff from GoaT
    goat_response = requests.get(f'https://goat.genomehubs.org/api/v2/search?query=tax_name%28{species_name}%29&result=taxon')
    goat_data = goat_response.json() # convert json to dict
    
    taxon_number = goat_data['results'][0]['result']['taxon_id']
    
    goat_results = goat_data['results']
    
    class_name = 'NA'
    order_name = 'NA'
    haploid_number = 'NA'
    haploid_source = 'NA'
    ploidy = 'NA'
    ploidy_source = 'NA'

    for result in goat_results:
        lineage = result['result']['lineage']
        for node in lineage:
            if node['taxon_rank'] == 'class':
                class_name = node['scientific_name']
            if node['taxon_rank'] == 'order':
                order_name = node['scientific_name']


    goat_second_response = requests.get(f'https://goat.genomehubs.org/api/v2/record?recordId={taxon_number}&result=taxon&taxonomy=ncbi')
    goat_second_data = goat_second_response.json()

    ploidy_info = goat_second_data['records'][0]['record']['attributes']['ploidy']

    ploidy = ploidy_info['value']
    ploidy_source = ploidy_info['aggregation_source']

    haploid_info = goat_second_data['records'][0]['record']['attributes']['haploid_number']

    haploid_number = haploid_info['value']
    haploid_source = haploid_info['aggregation_source']

    sp_data = [
        ["TxID", "ToLID", "Species", "Class", "Order"],
        [taxon_number, tol_id, species, class_name, order_name]
    ]

    # Transpose the data
    transposed_sp_data = list(map(list, zip(*sp_data)))


    # Reading SEQUENCING DATA section from yaml ###################################################

    # get DATA section from yaml
    data_list = yaml_data.get('DATA', [])

    # Prepare headers
    headers = ['Data']
    data_values = ['Coverage']

    # Extract data from YAML and format it for the table
    for item in data_list:
        for technology, coverage in item.items():
            headers.append(technology)
            data_values.append('NA' if not coverage else coverage)

    # Create a list of lists for the table
    table_data = [headers, data_values]


    # Extract pipeline data from 'Pre-curation' category
    asm_pipeline_data = yaml_data.get('ASSEMBLIES', {}).get('Pre-curation', {}).get('pipeline', [])
    asm_pipeline_tree = generate_pipeline_tree(asm_pipeline_data)


    # Extract pipeline data from 'Curated' category
    curation_pipeline_data = yaml_data.get('ASSEMBLIES', {}).get('Curated', {}).get('pipeline', [])
    curation_pipeline_tree = generate_pipeline_tree(curation_pipeline_data)


    # Reading GENOME PROFILING DATA section from yaml #############################################

    profiling_data = yaml_data.get('PROFILING')

    # Check if profiling_data is available
    if not profiling_data:
        logging.error('Error: No profiling data found in the YAML file.')
        sys.exit(1)

    # Handle GenomeScope specific processing
    genomescope_data = profiling_data.get('GenomeScope')
    if genomescope_data:
        summary_file = genomescope_data.get('genomescope_summary_txt')
        if summary_file and os.path.exists(summary_file):
            with open(summary_file, "r") as f:
                summary_txt = f.read()
            genome_haploid_length = re.search(r"Genome Haploid Length\s+([\d,]+) bp", summary_txt).group(1)
            proposed_ploidy_match = re.search(r"p = (\d+)", summary_txt)
            proposed_ploidy = proposed_ploidy_match.group(1) if proposed_ploidy_match else 'NA'
        else:
            logging.error(f"File {summary_file} not found for GenomeScope.")
            sys.exit(1)
    else:
        logging.error("GenomeScope data is missing in the PROFILING section.")
        sys.exit(1)

    # Handle Smudgeplot specific processing
    smudgeplot_data = profiling_data.get('Smudgeplot')
    if smudgeplot_data:
        verbose_summary_file = smudgeplot_data.get('smudgeplot_verbose_summary_txt')
        if verbose_summary_file and os.path.exists(verbose_summary_file):
            with open(verbose_summary_file, "r") as f:
                smud_summary_txt = f.readlines()
            for line in smud_summary_txt:
                if line.startswith("* Proposed ploidy"):
                    proposed_ploidy = line.split(":")[1].strip()
                    break
        else:
            logging.warning(f"Verbose summary file {verbose_summary_file} not found for Smudgeplot; skipping detailed Smudgeplot analysis.")
    else:
        logging.warning("Smudgeplot data is missing in the PROFILING section; skipping Smudgeplot analysis.")


    # Reading ASSEMBLY DATA section from yaml #####################################################
    
    asm_data = yaml_data.get('ASSEMBLIES', {})
    
    # make a list from the assemblies available in asm_data
    asm_stages = []
    for asm_stage, stage_properties in asm_data.items():
        for haplotypes in stage_properties.keys():
            if haplotypes != 'pipeline' and haplotypes not in asm_stages:
                asm_stages.append(haplotypes)
    
    # get gfastats-based data
    gfastats_data = {}
    for asm_stage, stage_properties in asm_data.items():
        for haplotypes, haplotype_properties in stage_properties.items():
            if isinstance(haplotype_properties, dict):
                if 'gfastats--nstar-report_txt' in haplotype_properties:
                    file_path = haplotype_properties['gfastats--nstar-report_txt']
                    with open(file_path, 'r') as file:
                        content = file.read()
                    gfastats_data[(asm_stage, haplotypes)] = extract_gfastats_values(content, keys)
    
    gaps_per_gbp_data = {}
    for (asm_stage, haplotypes), values in gfastats_data.items():
        try:
            gaps = float(values[gaps_index])
            total_length = float(values[total_length_index])
            gaps_per_gbp = round((gaps / total_length * 1_000_000_000), 2)
            gaps_per_gbp_data[(asm_stage, haplotypes)] = gaps_per_gbp
        except (ValueError, ZeroDivisionError):
            gaps_per_gbp_data[(asm_stage, haplotypes)] = ''
    
    # Define the contigging table (column names) DON'T MOVE THIS AGAIN!!!!!!!
    asm_table_data = [["Metrics"] + [f'{asm_stage} \n {haplotypes}' for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]]]
    
    # Fill the table with the gfastats data
    for i in range(len(display_names)):
        metric = display_names[i]
        if metric not in exclusion_list:
            asm_table_data.append([metric] + [format_number(gfastats_data.get((asm_stage, haplotypes), [''])[i]) if (asm_stage, haplotypes) in gfastats_data else '' for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]])
    
    # Add the gaps/gbp in between
    gc_index =  display_names.index("GC %")
    asm_table_data.insert(gaps_index + 1, ['Gaps/Gbp'] + [format_number(gaps_per_gbp_data.get((asm_stage, haplotypes), '')) for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]])
    
    # get QV, Kmer completeness and BUSCO data
    qv_data = {}
    completeness_data = {}
    busco_data = {metric: {} for metric in ['BUSCO sing.', 'BUSCO dupl.', 'BUSCO frag.', 'BUSCO miss.']}
    for asm_stage, stage_properties in asm_data.items():
        asm_stage_elements = [element for element in stage_properties.keys() if element != 'pipeline']
        for i, haplotypes in enumerate(asm_stage_elements):
            haplotype_properties = stage_properties[haplotypes]
            if isinstance(haplotype_properties, dict):
                if 'merqury_qv' in haplotype_properties:
                    qv_data[(asm_stage, haplotypes)] = get_qv_value(haplotype_properties['merqury_qv'], i, asm_stage, haplotypes)
                if 'merqury_completeness_stats' in haplotype_properties:
                    completeness_data[(asm_stage, haplotypes)] = get_completeness_value(haplotype_properties['merqury_completeness_stats'], i, asm_stage, haplotypes)
                if 'busco_short_summary_txt' in haplotype_properties:
                    s_value, d_value, f_value, m_value = extract_busco_values(haplotype_properties['busco_short_summary_txt'])
                    busco_data['BUSCO sing.'].update({(asm_stage, haplotypes): s_value})
                    busco_data['BUSCO dupl.'].update({(asm_stage, haplotypes): d_value})
                    busco_data['BUSCO frag.'].update({(asm_stage, haplotypes): f_value})
                    busco_data['BUSCO miss.'].update({(asm_stage, haplotypes): m_value})
    
    # Fill the table with the QV data
    asm_table_data.append(['QV'] + [qv_data.get((asm_stage, haplotypes), '') for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]])
    
    # Fill the table with the Kmer completeness data
    asm_table_data.append(['Kmer compl.'] + [completeness_data.get((asm_stage, haplotypes), '') for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]])
    
    # Fill the table with the BUSCO data
    for metric in ['BUSCO sing.', 'BUSCO dupl.', 'BUSCO frag.', 'BUSCO miss.']:
        asm_table_data.append([metric] + [busco_data[metric].get((asm_stage, haplotypes), '') for asm_stage in asm_data for haplotypes in asm_stages if haplotypes in asm_data[asm_stage]])
        

    # Reading CURATION NOTES section from yaml ####################################################

    obs_haploid_num = yaml_data.get("NOTES", {}).get("Obs_Haploid_num", "NA")
    obs_sex = yaml_data.get("NOTES", {}).get("Obs_Sex", "NA")    
    interventions_per_gb = yaml_data.get("NOTES", {}).get("Interventions_per_Gb", "NA")
    contamination_notes = yaml_data.get("NOTES", {}).get("Contamination_notes", "NA")
    other_notes = yaml_data.get("NOTES", {}).get("Other_notes", "NA")

    # Extract Total bp for each haplotype and find the maximum
    curated_assemblies = yaml_data.get('ASSEMBLIES', {}).get('Curated', {})
    total_bp_values = []
    for haplotype, properties in curated_assemblies.items():
        if 'gfastats--nstar-report_txt' in properties:
            total_bp = extract_total_bp_from_gfastats(properties['gfastats--nstar-report_txt'])
            total_bp_values.append(total_bp)
    
    max_total_bp = max(total_bp_values, default='NA')
    
    # Create table data
    genome_traits_table_data = [
        ["Genome Traits", "Expected", "Observed"],
        ["Haploid size (bp)", genome_haploid_length, f"{max_total_bp}"],
        ["Haploid Number", f"{haploid_number} (source: {haploid_source})", obs_haploid_num],
        ["Ploidy", f"{ploidy} (source: {ploidy_source})", proposed_ploidy],
        ["Sample Sex", sex, obs_sex]
    ]

    # Get curator notes
    curator_notes_text = (
        f". Interventions/Gb: {interventions_per_gb}<br/>"
        f". Contamination notes: &quot;{contamination_notes}&quot;<br/>"
        f". Other observations: &quot;{other_notes}&quot;"
    )
    
    

    # PDF CONSTRUCTION ############################################################################

    # Set up the PDF file
    pdf_filename = "EAR.pdf"
    margin = 0.5 * 72  # 0.5 inch in points (normal margin is 1 inch)
    pdf = SimpleDocTemplate(pdf_filename,
                            pagesize=A4,
                            leftMargin=margin,
                            rightMargin=margin,
                            topMargin=margin,
                            bottomMargin=margin)
    elements = []

    # Set all the styles
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='TitleStyle', fontName='Courier', fontSize=20))
    styles.add(ParagraphStyle(name='subTitleStyle', fontName='Courier', fontSize=16))
    styles.add(ParagraphStyle(name='normalStyle', fontName='Courier', fontSize=12))
    styles.add(ParagraphStyle(name='midiStyle', fontName='Courier', fontSize=10))
    styles.add(ParagraphStyle(name='LinkStyle', fontName='Courier', fontSize=10, textColor='blue', underline=True))
    styles.add(ParagraphStyle(name='treeStyle', fontName='Courier', fontSize=10, leftIndent=12))
    styles.add(ParagraphStyle(name='miniStyle', fontName='Courier', fontSize=8))
    styles.add(ParagraphStyle(name='FileNameStyle', fontName='Courier', fontSize=6))
    

    # PDF SECTION 1 -------------------------------------------------------------------------------
    
    # Add the title
    title = Paragraph("ERGA Assembly Report", styles['TitleStyle'])
    elements.append(title)

    # Spacer
    elements.append(Spacer(1, 12))

    # Add version 
    ver_paragraph = Paragraph(EAR_version, styles['normalStyle'])
    elements.append(ver_paragraph)

    # Spacer
    elements.append(Spacer(1, 12))

    # Add tags
    tags_paragraph = Paragraph(f"Tags: {tags}", styles['normalStyle'])
    elements.append(tags_paragraph)

    # Spacer
    elements.append(Spacer(1, 24))

    # Create the SPECIES DATA table with the transposed data
    sp_data_table = Table(transposed_sp_data)

    # Style the table
    sp_data_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, -1), '#e7e7e7'),  # Grey background for column 1
        ("BACKGROUND", (1, 0), (1, -1), colors.white),  # White background for column 2
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ('FONTNAME', (0, 0), (0, 0), 'Courier'),  # Regular font for row1, col1
        ('FONTNAME', (1, 0), (1, 0), 'Courier'),
        ('FONTNAME', (0, 1), (-1, -1), 'Courier'),  # Regular font for the rest of the table
        ('FONTNAME', (1, 1), (1, 1), 'Courier-Bold'),  # Bold font for row1, col2
        ("FONTSIZE", (0, 0), (-1, -1), 14),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black)
    ]))

    # Add SPECIES DATA table
    elements.append(sp_data_table)

    # Spacer
    elements.append(Spacer(1, 32))

    # Create the GENOME TRAITS table
    genome_traits_table = Table(genome_traits_table_data)
    
    # Style the table
    genome_traits_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (0, -1), '#e7e7e7'),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, -1), 'Courier'),
        ('FONTSIZE', (0, 0), (-1, -1), 12),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black)
    ]))

    # Add GENOME TRAITS table
    elements.append(genome_traits_table)

    # Spacer
    elements.append(Spacer(1, 28))
    
    # Add EBP METRICS SECTION subtitle
    subtitle = Paragraph("EBP metrics summary and curation notes", styles['subTitleStyle']) 
    elements.append(subtitle)

    # Spacer
    elements.append(Spacer(1, 24))

    # Iterate over haplotypes in the Curated category to get data for EBP metrics
    curated_assemblies = yaml_data.get('ASSEMBLIES', {}).get('Curated', {})
    haplotype_names = [key for key in curated_assemblies.keys() if key != 'pipeline']

    for haplotype in haplotype_names:
        properties = curated_assemblies[haplotype]
        if 'gfastats--nstar-report_txt' in properties and 'merqury_qv' in properties:
            gfastats_path = properties['gfastats--nstar-report_txt']
            order = haplotype_names.index(haplotype) # Determine the order based on the position of the haplotype in the list
            qv_value = get_qv_value(properties['merqury_qv'], order, 'Curated', haplotype)
        
            ebp_quality_metric = compute_ebp_metric(haplotype, gfastats_path, qv_value)
            EBP_metric_paragraph = Paragraph(ebp_quality_metric, styles["midiStyle"])
            
            # Add the EBP quality metric paragraph to elements
            elements.append(EBP_metric_paragraph)

    # Spacer        
    elements.append(Spacer(1, 8))

    # Add sentence
    Textline = Paragraph("The following metrics were automatically flagged as below EBP recommended standards or different from expected:", styles['midiStyle'])
    elements.append(Textline)

    # Spacer        
    elements.append(Spacer(1, 4))

    # Apply checks and add warning paragraphs to elements
    elements += generate_warning_paragraphs(genome_haploid_length, max_total_bp, "Haploid size (bp)")
    elements += generate_warning_paragraphs(haploid_number, obs_haploid_num, "Haploid Number")
    elements += generate_warning_paragraphs(proposed_ploidy, ploidy, "Ploidy")
    elements += generate_warning_paragraphs(sex, obs_sex, "Sample Sex")

    # Spacer        
    elements.append(Spacer(1, 4))

    # Iterate over haplotypes in the Curated category and apply checks
    for haplotype in haplotype_names:
        properties = curated_assemblies[haplotype]
        if isinstance(properties, dict) and 'merqury_qv' in properties and 'merqury_completeness_stats' in properties and 'busco_short_summary_txt' in properties:
            order = haplotype_names.index(haplotype)
            qv_value = get_qv_value(properties['merqury_qv'], order, "Curated", haplotype)
            completeness_value = get_completeness_value(properties['merqury_completeness_stats'], order, "Curated", haplotype)
            busco_scores = extract_busco_values(properties['busco_short_summary_txt'])

            warnings = generate_curated_warnings(haplotype, qv_value, completeness_value, busco_scores)
            elements += warnings

    assembly_warnings = generate_assembly_warnings(asm_data, gaps_per_gbp_data, obs_haploid_num)
    elements.extend(assembly_warnings)

    # Spacer
    elements.append(Spacer(1, 24))

    # Add small subtitle for Curator notes    
    subtitle = Paragraph("Curator notes", styles['normalStyle']) 
    elements.append(subtitle)

    # Spacer
    elements.append(Spacer(1, 8))

    # Curator notes
    curator_notes_paragraph = Paragraph(curator_notes_text, styles["midiStyle"])
    elements.append(curator_notes_paragraph)
        
    # Page break
    elements.append(PageBreak())


    # PDF SECTION 2 -------------------------------------------------------------------------------

    # Add quality metrics section subtitle
    subtitle = Paragraph("Quality metrics table", styles['TitleStyle'])
    elements.append(subtitle)

    # Spacer 
    elements.append(Spacer(1, 48))

    # create QUALITY METRICS table
    asm_table = Table(asm_table_data)

    # Style the table
    asm_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), '#e7e7e7'),  # grey background for the header
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),       # center alignment
        ('FONTNAME', (0, 0), (-1, -1), 'Courier'),  # bold font for the header
        ('FONTSIZE', (0, 0), (-1, -1), 11),           # font size
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black)
    ]))

    # Add QUALITY METRICS table
    elements.append(asm_table)

    # Spacer
    elements.append(Spacer(1, 5))

    # Store BUSCO version and lineage information from each file in list
    busco_info_list = []
    for asm_stages, stage_properties in asm_data.items():
        for haplotype_keys, haplotype_properties in stage_properties.items():
            if isinstance(haplotype_properties, dict):
                if 'busco_short_summary_txt' in haplotype_properties:
                    busco_version, lineage_info = extract_busco_info(haplotype_properties['busco_short_summary_txt'])
                    if busco_version and lineage_info:
                        busco_info_list.append((busco_version, lineage_info))
    
    # Checking if all elements in the list are identical
    if all(info == busco_info_list[0] for info in busco_info_list):
        busco_version, (lineage_name, num_genomes, num_buscos) = busco_info_list[0]
        elements.append(Paragraph(f"BUSCO {busco_version} Lineage: {lineage_name} (genomes:{num_genomes}, BUSCOs:{num_buscos})", styles['miniStyle']))
    else:
        elements.append(Paragraph("Warning: BUSCO versions or lineage datasets are not the same across results", styles['miniStyle']))
        logging.warning(f"WARNING!!! BUSCO versions or lineage datasets are not the same across results")
   
    # Page break
    elements.append(PageBreak())


    # PDF SECTION 3 -------------------------------------------------------------------------------

    # Add hic maps section subtitle
    subtitle = Paragraph("HiC contact map of curated assembly", styles['TitleStyle'])
    elements.append(subtitle)

    # Spacer 
    elements.append(Spacer(1, 36))
    
    # Initialize counter
    tool_count = 0

    # Add title and images for each step
    for idx, (asm_stages, stage_properties) in enumerate(asm_data.items(), 1):
        if asm_stages == 'Curated':
            tool_elements = [element for element in stage_properties.keys() if element != 'pipeline']
        
            images_with_names = []

            for haplotype in tool_elements:
                haplotype_properties = stage_properties[haplotype]
                
                # Check if there is an image and/or a link
                png_file = haplotype_properties.get('hic_FullMap_png', '')
                link = haplotype_properties.get('hic_FullMap_link', '')

                # Prepare paragraphs for the image and link
                if png_file:
                    # Create image object
                    img = Image(png_file, width=11 * cm, height=11 * cm)
                    images_with_names.append([img])
                else:
                    # Add paragraph for missing image
                    missing_png_paragraph = Paragraph(f"<b>{haplotype}</b> HiC PNG is missing!", styles["midiStyle"])
                    images_with_names.append([missing_png_paragraph])
                
                # Add paragraph for the link
                if link:
                    link_html = f'<b>{haplotype}</b> <link href="{link}" color="blue">[LINK]</link>'
                else:
                    link_html = f'<b>{haplotype}</b> File link is missing!'

                link_paragraph = Paragraph(link_html, styles["midiStyle"])
                images_with_names.append([link_paragraph])
                
                # Append a spacer only if the next element is an image
                if len(tool_elements) > 1 and tool_elements.index(haplotype) < len(tool_elements) - 1:
                    images_with_names.append([Spacer(1, 12)]) 

            # Add images and names to the elements in pairs
            for i in range(0, len(images_with_names), 4):  # Process two images (and their names) at a time
                elements_to_add = images_with_names[i:i+4]
        
                # Create table for the images and names
                table = Table(elements_to_add)
                table.hAlign = 'CENTER'
                elements.append(table)
        
                # Add a page break conditionally
                next_elements_start = i + 4
                if next_elements_start < len(images_with_names):
                    if len(images_with_names[next_elements_start]) > 0 and isinstance(images_with_names[next_elements_start][0], Image):
                        elements.append(PageBreak())
        
            tool_count += 1

    elements.append(PageBreak())


    # PDF SECTION 4 -------------------------------------------------------------------------------

    # Add kmer spectra section subtitle
    subtitle = Paragraph("K-mer spectra of curated assembly", styles['TitleStyle'])
    elements.append(subtitle)

    # Spacer 
    elements.append(Spacer(1, 48))

    # Initialize counter
    counter = 0

    # Iterate over haplotypes in the Curated category to get K-mer spectra images
    curated_assemblies = yaml_data.get('ASSEMBLIES', {}).get('Curated', {})
    haplotype_names = [key for key in curated_assemblies.keys() if key != 'pipeline']

    # Get paths for spectra files
    spectra_files = {
        'hap1': {
            'spectra_cn_png': curated_assemblies.get('hap1', {}).get('merqury_hap_spectra_cn_png', None),
        },
        'hap2': {
            'spectra_cn_png': curated_assemblies.get('hap2', {}).get('merqury_hap_spectra_cn_png', None),
        },
        'common': {
            'spectra_cn_png': curated_assemblies.get('hap1', {}).get('merqury_spectra_cn_png', None),
            'spectra_asm_png': curated_assemblies.get('hap1', {}).get('merqury_spectra_asm_png', None),
        }
    }

    # Filter out None values and empty strings
    spectra_files = {k: {sk: v for sk, v in sv.items() if v} for k, sv in spectra_files.items()}

    # Determine the number of spectra-cn files and assign unique names if needed
    spectra_cn_files = [
        spectra_files['common'].get('spectra_cn_png', None),
        spectra_files['hap1'].get('spectra_cn_png', None),
        spectra_files['hap2'].get('spectra_cn_png', None)
    ]
    spectra_cn_files = [f for f in spectra_cn_files if f]  # Filter out None values

    if len(spectra_cn_files) == 3:
        # For 3 spectra-cn files
        shortest_spectra_cn_file = min(spectra_cn_files, key=lambda f: len(os.path.basename(f)), default=None)
        similar_files = [f for f in spectra_cn_files if f != shortest_spectra_cn_file]
        if similar_files:
            unique_name1, unique_name2 = find_unique_parts(os.path.basename(similar_files[0]), os.path.basename(similar_files[1]))
    else:
        shortest_spectra_cn_file = spectra_cn_files[0] if spectra_cn_files else None
        unique_name1 = unique_name2 = None

    # Create image objects and add filename below each image
    images = []

    for label, file_dict in spectra_files.items():
        for key, png_file in file_dict.items():
            if png_file:
                image = Image(png_file, width=8.4 * cm, height=7 * cm)
                filename = os.path.basename(png_file)

                if filename.endswith("spectra-asm.ln.png"):
                    text = "Distribution of k-mer counts coloured by their presence in reads/assemblies"
                elif filename.endswith("spectra-cn.ln.png"):
                    if len(spectra_cn_files) == 3:
                        # For 3 spectra-cn files use particular text
                        if png_file == shortest_spectra_cn_file:
                            text = "Distribution of k-mer counts per copy numbers found in asm (dipl.)"
                        else:
                            if png_file == spectra_files['hap1'].get('spectra_cn_png', None):
                                text = f"Distribution of k-mer counts per copy numbers found in <b>{unique_name1}</b> (hapl.)"
                            elif png_file == spectra_files['hap2'].get('spectra_cn_png', None):
                                text = f"Distribution of k-mer counts per copy numbers found in <b>{unique_name2}</b> (hapl.)"
                            else:
                                text = "Distribution of k-mer counts per copy numbers found in asm"
                    else:
                        # For 2 spectra-cn files use same text
                        text = "Distribution of k-mer counts per copy numbers found in asm"
                else:
                    text = filename
                
                images.append([image, Paragraph(text, styles["midiStyle"])])

    # Filter None values
    images = [img for img in images if img[0] is not None]

    # Get number of rows and columns for the table
    num_rows = (len(images) + 1) // 2  # +1 to handle odd numbers of images
    num_columns = 2

    # Create the table with dynamic size
    image_table_data = [[images[i * num_columns + j] if i * num_columns + j < len(images) else [] for j in range(num_columns)] for i in range(num_rows)]
    image_table = Table(image_table_data)

    # Style the "table"
    table_style = TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 20),  #20 here is a spacer between rows
    ])

    # Set the style
    image_table.setStyle(table_style)

    # Add image table to elements
    elements.append(image_table)

    # Increase counter by the number of PNGs added
    counter += len(images)

    # If counter is a multiple of 4, insert a page break and reset counter
    if counter % 4 == 0:
        elements.append(PageBreak())

    # Add spacer
    elements.append(Spacer(1, 12))

    # If we have processed all haps and the last page does not contain exactly 4 images, insert a page break
    if counter % 4 != 0:
        elements.append(PageBreak())
 
    
    # PDF SECTION 5 -------------------------------------------------------------------------------

    # Add contamination section subtitle
    subtitle = Paragraph("Post-curation contamination screening", styles['TitleStyle'])
    elements.append(subtitle)

    # Spacer 
    elements.append(Spacer(1, 36))
    
    # Initialize counter
    tool_count = 0

    # Add title and images for each step
    for idx, (asm_stages, stage_properties) in enumerate(asm_data.items(), 1):
        if asm_stages == 'Curated':  # Check if the current stage is 'Curated'
            tool_elements = [element for element in stage_properties.keys() if element != 'pipeline']

            for haplotype in tool_elements:
                haplotype_properties = stage_properties[haplotype]
                if isinstance(haplotype_properties, dict) and 'blobplot_cont_png' in haplotype_properties:
                    # Get image path
                    png_file = haplotype_properties['blobplot_cont_png']

                    # If png_file is not empty, display it
                    if png_file:
                        # Create image object
                        img = Image(png_file, width=20 * cm, height=20 * cm)
                        elements.append(img)

                        # Create paragraph for filename with haplotype name
                        blob_text = f"<b>{haplotype}.</b> Bubble plot circles are scaled by sequence length, positioned by coverage and GC proportion, and coloured by taxonomy. Histograms show total assembly length distribution on each axis."
                        blob_paragraph = Paragraph(blob_text, styles["midiStyle"])
                        elements.append(blob_paragraph)
                    else:
                        # Add paragraph for missing image
                        missing_png_paragraph = Paragraph(f"<b>{haplotype}</b> PNG is missing!", styles["midiStyle"])
                        elements.append(missing_png_paragraph)

                    # Add a page break after each image and its description
                    elements.append(PageBreak())

            tool_count += 1


    # SECTION 6 -----------------------------------------------------------------------------------

    # Add data profile section subtitle
    subtitle = Paragraph("Data profile", styles['TitleStyle']) 
    elements.append(subtitle)

    # Spacer
    elements.append(Spacer(1, 24))

    # Create the DATA PROFILE table
    data_table = Table(table_data)

    # Style the table
    data_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (0, -1), '#e7e7e7'),  # grey background for the first column
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),         # center alignment
        ('FONTNAME', (0, 0), (-1, -1), 'Courier'),      # remove bold font
        ('FONTSIZE', (0, 0), (-1, -1), 12),             # font size for the header
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.black)
    ]))

    # Add DATA PROFILE table
    elements.append(data_table)

    # Spacer
    elements.append(Spacer(1, 32))

    # Add assembly pipeline section subtitle
    subtitle = Paragraph("Assembly pipeline", styles['TitleStyle']) 
    elements.append(subtitle)
 
    # Spacer
    elements.append(Spacer(1, 24))

    # Add ASM PIPELINE tree
    elements.append(Paragraph(asm_pipeline_tree, styles['treeStyle']))

    # Spacer
    elements.append(Spacer(1, 32))

    # Add curation pipeline section subtitle
    subtitle = Paragraph("Curation pipeline", styles['TitleStyle'])  
    elements.append(subtitle)
 
    # Spacer
    elements.append(Spacer(1, 24))

    # Add CURATION PIPELINE tree
    elements.append(Paragraph(curation_pipeline_tree, styles['treeStyle']))
    
    # Spacer
    elements.append(Spacer(1, 48))

    # Add submitter, affiliation
    submitter_paragraph_style = ParagraphStyle(name='SubmitterStyle', fontName='Courier', fontSize=10)
    elements.append(Paragraph(f"Submitter: {submitter}", submitter_paragraph_style))
    elements.append(Paragraph(f"Affiliation: {affiliation}", submitter_paragraph_style))

    # Spacer
    elements.append(Spacer(1, 8))

    # Add the date and time (CET) of the document creation
    cet = pytz.timezone("CET")
    current_datetime = datetime.now(cet)
    formatted_datetime = current_datetime.strftime("%Y-%m-%d %H:%M:%S %Z")
    elements.append(Paragraph(f"Date and time: {formatted_datetime}", submitter_paragraph_style))



    # Build the PDF ###############################################################################
    pdf.build(elements)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create an ERGA Assembly Report (EAR) from a YAML file. Visit https://github.com/ERGA-consortium/EARs for more information')
    parser.add_argument('yaml_file', type=str, help='Path to the YAML file')
    args = parser.parse_args()
    
    make_report(args.yaml_file)