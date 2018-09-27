#!/usr/bin/env python

import inspect
import io
from jinja2 import Template
import os
from planemo import tool_builder
import re
import scanpy.api as sc
from lxml import etree


ADATA_INPUT_TEMPLATE = """
#if $input.format == 'csv'
    adata = sc.read_csv($input.adata, )
# else if $input.format == 'loom'
    adata = sc.read_hdf(
        $input.adata,
        sparse=$input.sparse,
        cleanup=$input.cleanup,
        X_name=$input.x_name,
        obs_names=$input.obs_names,
        var_names=$input.var_names)
# else if $input.format == 'h5'
    adata = sc.read_loom(
        $input.adata,
        key=$input.key)
"""

CMD_TEMPLATE = Template("""
sc.{{ name }}({{ params }})
""")

ADATA_OUTPUT_TEMPLATE = """
adata.write_loom($csv_output)
adata.write_csv($loom_output)
"""

SCRIPT_TEMPLATE = Template("""
import scanpy.api as sc
{{ settings }}
{{ inputs }}
{{ script }}
{{ outputs }}
""")


MOD_MACRO_FILENAME = Template("{{ mod_name }}_macros.xml")

parser = etree.XMLParser(remove_blank_text=True)

def init_param(name, el_type, desc=''):
    '''Create a etree Element describing a parameter'''
    # create etree Element
    param = etree.Element(el_type, name=name)
    # extract format
    type_regex = re.search('^`(?P<type>.+?)`', desc)
    if type_regex:
        param.set('type', type_regex.group("type"))
    else:
        param.set('type', 'data')
        param.set('format', desc.split('(')[0].split(',')[0])
    # extract default value
    if el_type == 'param':
        value_regex = re.search('default: `(?P<value>.+?)`', desc)
        value = ''
        if value_regex:
            value = value_regex.group("value")
        param.set('value', value)
    # extract if optional
    if 'optional' in desc:
        param.set('optional', 'true')
    # set label
    if el_type == 'param':
        param.set('label', name)
    else:
        param.set('label', '${tool.name} on ${on_string}: %s' % name)
    # set help
    if el_type == 'param':
        param.set('help', '')
    return param


def get_formatted_params(string, el_type):
    '''Extract formatted params from a string'''
    format_params = []
    for l in string.split('\n'):
        if l == '':
            continue
        regex = "^(?P<name>.+?) : (?P<format>.+)"
        param_regex = re.search(regex, l)
        if param_regex:
            format_params.append(init_param(
                param_regex.group("name"),
                desc=param_regex.group("format"),
                el_type=el_type))
        elif l.startswith(' '):
            prev_param = format_params[-1]
            if 'help' in prev_param.attrib:
                prev_param.set('help', prev_param.get('help') + l.rstrip())
        else:
            format_params.append(init_param(l.rstrip(), el_type=el_type))
    return format_params


def extract_doc(met):
    '''Extract from the documentation of the method'''
    doc = inspect.getdoc(met)
    help_section = ''
    short_desc = ''
    params = []
    returns = []
    if doc:
        # keep full documentation for the help
        help_section = doc
        # keep first sentence for the description
        short_desc = doc.split('\n')[0]
        # extract the parameter and output description
        doc_regex = re.search(
            r'(?P<desc>[\s\S]*)\nParameters\n-+\n(?P<followup>[\s\S]*)',
            doc)
        if doc_regex:
            param_followup = doc_regex.group("followup")
            param_str = param_followup
        else:
            param_followup = doc
            param_str = None
        
        return_regex = re.search(
            r'(?P<params>[\s\S]*)\nReturns\n-+\n(?P<followup>[\s\S]*)',
            param_followup)
        if return_regex:
            if param_str:
                params = get_formatted_params(return_regex.group("params"), el_type='param')
            return_followup = return_regex.group("followup")
            notes_regex = re.search(r'(?P<returns>[\s\S]*)\nNotes', return_followup)
            examples_regex = re.search(r'(?P<returns>[\s\S]*)\nExamples', return_followup)
            if notes_regex:
                returns = get_formatted_params(notes_regex.group('returns'), el_type='data')
            elif examples_regex:
                returns = get_formatted_params(examples_regex.group('returns'), el_type='data')
            else:   
                returns = get_formatted_params(return_followup, el_type='data')
        elif param_str:
            params = get_formatted_params(param_str, el_type='param')                    
    return short_desc, help_section, params, returns   
    

def generate_xml(met, from_mod):
    '''Generate the XML for the method'''
    # Get description and parameter for the method
    met_name = met.__name__
    if from_mod:
        met_name = "%s.%s" % (from_mod, met_name)
    short_desc, help_section, params, returns = extract_doc(met)
    xml_path = "%s.xml" % met_name
    # Generate the XML skeleton
    kwds = {
        'id': met_name,
        'name': met_name,
        'tool': xml_path,
        'description': short_desc,
        'command': 'python $script_file',
        'help_text': help_section,
        'help_from_command': None,
        'doi': ['10.1186/s13059-017-1382-0'],
        'cite_url': [],
        'test_case': None,
        'version': '1.3.1+galaxy1',
        'requirement': ['scanpy@1.3.1'],
        'container': None,
        'macros': True,
        'force': True
    }
    ctx = {}
    tool_description = tool_builder.build(**kwds)
    tool_builder.write_tool_description(ctx, tool_description, **kwds)
    # Get XML structure
    tool_xml = etree.parse(xml_path, parser=parser)
    root = tool_xml.getroot()
    for element in root.iter():
        element.tail = None
    # Create an empty test
    tests = etree.Element("tests")
    test_xml = etree.SubElement(tests, "test")
    # Format the inputs and fill the test
    param_py_str = ''
    adata_inputs = ''
    input_xml = root.find('inputs')
    for p in params:
        if p.get('type') == 'data' and 'AnnData' in p.get('format'):
            adata_inputs += '@CMD_read_inputs'
            n_p = etree.SubElement(input_xml, 'expand', macro = 'inputs_anndata')
            n_p.tail = "\n    "
            param_py_str += "\n   %s = '$adata'," %(p.get('name'))
        else:
            n_p = etree.SubElement(input_xml, p.tag, p.attrib)
            n_p.tail = "\n    "
            param_py_str += "\n   %s = '$%s'," %(n_p.get('name'), n_p.get('name'))    
            test = etree.SubElement(test_xml, "param")
            test.set('name', n_p.get('name'))
            test.set('value', '')
    if met_name.startswith('pl'):
        pl_settings = etree.SubElement(input_xml, 'section')
        pl_settings.set('name', 'pl_settings')
        pl_settings.set('title', 'Plot settings')
        pl_settings.set('expanded', 'false')
        etree.SubElement(pl_settings, 'expand', macro = 'inputs_set_figure_params')
        test = etree.SubElement(test_xml, 'expand', macro = 'tests_set_figure_params')
    # Format the output and fill the test
    output_xml = root.find('outputs')
    adata_outputs = ''
    for o in returns:
        n_o = etree.SubElement(output_xml, o.tag, o.attrib)
        n_o.tail = "\n    "
        if 'adata' in n_o.get('name'):
            n_o.set('name', "csv_output")
            n_o.set('type', "data")
            n_o.set('format', 'csv')
            n_o.set('label', "${tool.name} on ${on_string}: Annotated matrix (csv)")
            loom_el = etree.Element('data')
            loom_el.set('name', "loom_output")
            loom_el.set('type', "data")
            loom_el.set('format', 'loom')
            loom_el.set('label', "${tool.name} on ${on_string}: Annotated matrix (loom)")
            loom_el.tail = "\n    "
            returns.append(loom_el)
            test_loom = etree.SubElement(test_xml, "output")
            test_loom.set('name', loom_el.get('name'))
            test_loom.set('file', '')
            adata_outputs += ADATA_OUTPUT_TEMPLATE
        test = etree.SubElement(test_xml, "output")
        test.set('name', n_o.get('name'))
        test.set('file', '')
    # Add the test
    root.insert(6, tests)
    # Get the script for the config section
    cmd = CMD_TEMPLATE.render(**{'name': met_name, 'params': param_py_str})
    if met_name.startswith('pl'):
        settings = '@CMD_set_figure_params'
    else:
        settings = ''
    script = SCRIPT_TEMPLATE.render(**{
        'settings': settings,
        'inputs': adata_inputs,
        'script': cmd,
        'outputs': adata_outputs})
    # Add configfile
    configfiles = etree.Element("configfiles")
    configfile = etree.SubElement(configfiles, "configfile", name="script_file")
    configfile.text = etree.CDATA(script)
    root.insert(4, configfiles)
    # transform help into CDATA
    help_xml = root.find('help')
    help_xml.text = etree.CDATA(help_xml.text)
    # transform command into CDATA
    cmd_xml = root.find('command')
    cmd_xml.text = etree.CDATA(cmd_xml.text)
    # overwrite file
    with open(xml_path, 'w') as f:
        f.write(etree.tostring(tool_xml, pretty_print=True, encoding='unicode'))


def expand_macros(met, macro_root):
    '''Add cmd, inputs and outputs of met to macros file'''
    # Get description and parameter for the method
    met_name = met.__name__
    short_desc, help_section, params, returns = extract_doc(met)
    # add tests
    test_xml = etree.SubElement(macro_root, "xml", name='tests_%s' % met_name)
    # add inputs
    param_py_str = ''
    input_xml = etree.SubElement(macro_root, "xml", name='inputs_%s' % met_name)
    for p in params:
        n_p = etree.SubElement(input_xml, p.tag, p.attrib)
        n_p.tail = "\n    "
        test = etree.SubElement(test_xml, "param")
        test.set('name', n_p.get('name'))
        test.set('value', '')
        param_py_str += "\n   %s = '$pl_settings.%s'," %(n_p.get('name'), n_p.get('name'))
    # add outputs
    output_xml = etree.SubElement(macro_root, "xml", name='outputs_%s' % met_name)
    for o in returns:
        n_o = etree.SubElement(output_xml, o.tag, o.attrib)
        n_o.tail = "\n    "
    # add cmd
    cmd_token = etree.SubElement(macro_root, "token", name='@CMD_%s@' % met_name)
    cmd_token.text = etree.CDATA(CMD_TEMPLATE.render(
        **{'name': met_name, 'params': param_py_str}))


def create_option(xml_parent, value):
    '''Create an XML option subelement'''
    option = etree.SubElement(xml_parent, "option", value=value)
    option.text = value


def create_when(xml_parent, value='csv'):
    '''Create an XML when subelement with input parameter'''
    when = etree.SubElement(xml_parent, "when", value=value)
    when_param = etree.SubElement(when, "param", name='adata')
    when_param.set('type', 'data')
    when_param.set('format', value)
    when_param.set('label', 'Annotated data matrix')
    return when


def write_element_xml(mod, macro_root, from_root = None):
    '''Write XML for functions in a module and its submodules'''
    # modules
    if not from_root:
        for smod in inspect.getmembers(mod, inspect.ismodule):
            if 'logging' in smod[0] or 'settings' in smod[0]:
                continue
            print(smod)
            write_element_xml(smod[1], macro_root, from_root = smod[0])
    # functions
    for func in inspect.getmembers(mod, inspect.isfunction):
        if func[0].startswith("_") or func[0].startswith("read") or func[0].startswith("write"):
            continue
        print(func)
        el = func[1]
        if func[0] == 'set_figure_params':
            expand_macros(el, macro_root)
        else:
            generate_xml(el, from_mod=from_root)
    # methods
    for met in inspect.getmembers(mod, inspect.ismethod):
        if met[0].startswith("read") or met[0].startswith("write"):
            continue
        print(met)
        el = met[1]
        generate_xml(el, from_mod=from_root)


if __name__ == '__main__':
    # create macro file
    macro_root = etree.Element("macros")
    # add requirement
    requirements_macros = etree.SubElement(macro_root, "xml", name="requirements")
    requirements = etree.SubElement(requirements_macros, "requirements")
    req = etree.SubElement(requirements, "requirement")
    req.set('type',"package")
    req.set('version',"1.3.1")
    req.text = 'scanpy'
    # add citation
    citations_macros = etree.SubElement(macro_root, "xml", name="citations")
    citations = etree.SubElement(citations_macros, "citations")
    req = etree.SubElement(citations, "citation", type = "doi")
    req.text = '10.1186/s13059-017-1382-0'    
    # add macros for input and its format
    input_xml = etree.SubElement(macro_root, "xml", name='inputs_anndata')
    input_cond = etree.SubElement(input_xml, "conditional", name='input')
    format_select = etree.SubElement(input_cond, "param", name='format')
    format_select.set('type', 'select')
    format_select.set('label', 'Format for the annotated data matrix')
    create_option(format_select, value='csv')
    create_option(format_select, value='loom')
    create_option(format_select, value='h5')
    # add when for csv
    create_when(input_cond, value='csv')
    # add when for loom
    when = create_when(input_cond, value='loom')
    subparam = etree.SubElement(when, "param", name='sparse')
    subparam.set('type', 'boolean')
    subparam.set('truevalue', 'True')
    subparam.set('falsevalue', 'False')
    subparam.set('checked', 'true')
    subparam.set('label', 'Is the data matrix to read sparse?')
    subparam = etree.SubElement(when, "param", name='cleanup')
    subparam.set('type', 'boolean')
    subparam.set('truevalue', 'True')
    subparam.set('falsevalue', 'False')
    subparam.set('checked', 'false')
    subparam.set('label', 'Cleanup?')
    subparam = etree.SubElement(when, "param", name='x_name')
    subparam.set('type', 'text')
    subparam.set('value', 'spliced')
    subparam.set('label', 'X_name')
    subparam = etree.SubElement(when, "param", name='obs_names')
    subparam.set('type', 'text')
    subparam.set('value', 'CellID')
    subparam.set('label', 'obs_names')
    subparam = etree.SubElement(when, "param", name='var_names')
    subparam.set('type', 'text')
    subparam.set('value', 'Gene')
    subparam.set('label', 'var_names')
    # add when for h5
    when = create_when(input_cond, value='h5')
    subparam = etree.SubElement(when, "param", name='key')
    subparam.set('type', 'text')
    subparam.set('value', '')
    subparam.set('label', 'Name of dataset in the file')
    # add token with command
    cmd_token = etree.SubElement(macro_root, "token", name='@CMD_read_inputs')
    cmd_token.text = etree.CDATA(ADATA_INPUT_TEMPLATE)
    # parse all elements
    write_element_xml(sc, macro_root)
    # write macros file
    macros_file = 'macros.xml'
    with open(macros_file, 'w') as f:
        f.write(etree.tostring(macro_root, pretty_print=True, encoding='unicode'))
