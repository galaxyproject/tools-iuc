import shutil
import sys
import tempfile


def parse_gff_attributes(attr_str):
    """
    Parses a GFF/GTF attribute string and returns a dictionary of name-value
    pairs. The general format for a GFF3 attributes string is

        name1=value1;name2=value2

    The general format for a GTF attribute string is

        name1 "value1" ; name2 "value2"

    The general format for a GFF attribute string is a single string that
    denotes the interval's group; in this case, method returns a dictionary
    with a single key-value pair, and key name is 'group'
    """
    attributes_list = attr_str.split(";")
    attributes = {}
    for name_value_pair in attributes_list:
        # Try splitting by '=' (GFF3) first because spaces are allowed in GFF3
        # attribute; next, try double quotes for GTF.
        pair = name_value_pair.strip().split("=")
        if len(pair) == 1:
            pair = name_value_pair.strip().split("\"")
        if len(pair) == 1:
            # Could not split for some reason -- raise exception?
            continue
        if pair == '':
            continue
        name = pair[0].strip()
        if name == '':
            continue
        # Need to strip double quote from values
        value = pair[1].strip(" \"")
        attributes[name] = value

    if len(attributes) == 0:
        # Could not split attributes string, so entire string must be
        # 'group' attribute. This is the case for strictly GFF files.
        attributes['group'] = attr_str
    return attributes


def gff_attributes_to_str(attrs, gff_format):
    """
    Convert GFF attributes to string. Supported formats are GFF3, GTF.
    """
    if gff_format == 'GTF':
        format_string = '%s "%s"'
        # Convert group (GFF) and ID, parent (GFF3) attributes to transcript_id, gene_id
        id_attr = None
        if 'group' in attrs:
            id_attr = 'group'
        elif 'ID' in attrs:
            id_attr = 'ID'
        elif 'Parent' in attrs:
            id_attr = 'Parent'
        if id_attr:
            attrs['transcript_id'] = attrs['gene_id'] = attrs[id_attr]
    elif gff_format == 'GFF3':
        format_string = '%s=%s'
    attrs_strs = []
    for name, value in attrs.items():
        attrs_strs.append(format_string % (name, value))
    return " ; ".join(attrs_strs)


stderr = sys.argv[1]
global_model_file_name = sys.argv[2]
transcripts = sys.argv[3]

# Read standard error to get total map/upper quartile mass.
total_map_mass = -1
with open(stderr, 'r') as tmp_stderr2:
    for line in tmp_stderr2:
        if line.lower().find("map mass") >= 0 or line.lower().find("upper quartile") >= 0:
            total_map_mass = float(line.split(":")[1].strip())
            break

if global_model_file_name != "None":
    # Global model is simply total map mass from original run.
    with open(global_model_file_name, 'r') as global_model_file:
        global_model_total_map_mass = float(global_model_file.readline())

    # Ratio of global model's total map mass to original run's map mass is
    # factor used to adjust FPKM.
    fpkm_map_mass_ratio = total_map_mass / global_model_total_map_mass

    # Update FPKM values in transcripts.gtf file.
    with open(transcripts, 'r') as transcripts_file:
        with tempfile.NamedTemporaryFile(dir=".", delete=False) as new_transcripts_file:
            for line in transcripts_file:
                fields = line.split('\t')
                attrs = parse_gff_attributes(fields[8])
                attrs["FPKM"] = str(float(attrs["FPKM"]) * fpkm_map_mass_ratio)
                attrs["conf_lo"] = str(float(attrs["conf_lo"]) * fpkm_map_mass_ratio)
                attrs["conf_hi"] = str(float(attrs["conf_hi"]) * fpkm_map_mass_ratio)
                fields[8] = gff_attributes_to_str(attrs, "GTF")
                new_transcripts_file.write("%s\n" % '\t'.join(fields))
    shutil.move(new_transcripts_file.name, transcripts)

if total_map_mass > -1:
    with open("global_model.txt", 'w') as f:
        f.write("%f\n" % total_map_mass)
