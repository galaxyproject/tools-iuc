#!/usr/bin/env python
# Dan Blankenberg

import json
import logging
import optparse
import os
import shutil
import tempfile
from xml.etree.ElementTree import tostring

try:
    # For Python 3.0 and later
    from shutil import unpack_archive
except ImportError:
    # Fall back to Python 2 import
    from setuptools.archive_util import unpack_archive

try:
    # For Python 3.0 and later
    from urllib.request import urlretrieve
    from urllib.parse import urlsplit
except ImportError:
    # Fall back to Python 2 imports
    from urllib import urlretrieve
    from urlparse import urlsplit

_log_name = __name__
if _log_name == '__builtin__':
    _log_name = 'toolshed.installed.manual.data.manager'
log = logging.getLogger(_log_name)


# --- These methods are called by/within the Galaxy Application
def exec_before_job(app, inp_data, out_data, param_dict, tool=None, **kwd):
    # Look for any data tables that haven't been defined for this data manager before and dynamically add them to Galaxy
    param_dict = dict(**param_dict)
    data_tables_param = param_dict.get('data_tables', [])
    if not isinstance(data_tables_param, list):
        data_tables_param = [data_tables_param]
    if tool:
        tool_shed_repository = tool.tool_shed_repository
    else:
        tool_shed_repository = None
    tdtm = None
    data_manager = app.data_managers.get_manager(tool.data_manager_id, None)
    for data_table_param in data_tables_param:
        data_table_name = data_table_param.get('data_table_name')
        if data_table_name:
            # the 'data_table_name' value in data_table_param is a SelectToolParameter,
            # to get the selected value we need to cast data_table_name to string
            data_table_name = str(data_table_name)
            # get data table managed by this data Manager
            data_table = app.tool_data_tables.get_tables().get(data_table_name)
            if data_table:
                data_table_filename = data_table.get_filename_for_source(data_manager, None)
                if not data_table_filename:
                    if tdtm is None:
                        from tool_shed.tools import data_table_manager
                        tdtm = data_table_manager.ToolDataTableManager(app)
                        target_dir, tool_path, relative_target_dir = tdtm.get_target_install_dir(tool_shed_repository)
                    # Dynamically add this data table
                    log.debug("Attempting to dynamically create a missing Tool Data Table named %s." % data_table_name)
                    repo_info = tdtm.generate_repository_info_elem_from_repository(tool_shed_repository, parent_elem=None)
                    if repo_info is not None:
                        repo_info = tostring(repo_info)
                    tmp_file = tempfile.NamedTemporaryFile(mode="w")
                    tmp_file.write(__get_new_xml_definition(app, data_table, data_manager, repo_info, target_dir))
                    tmp_file.flush()
                    app.tool_data_tables.add_new_entries_from_config_file(tmp_file.name, None, app.config.shed_tool_data_table_config, persist=True)
                    tmp_file.close()


def __get_new_xml_definition(app, data_table, data_manager, repo_info=None, location_file_dir=None):
    sub_dict = {'table_name': data_table.name, 'comment_char': '', 'columns': '', 'file_path': ''}
    sub_dict.update(data_manager.get_tool_shed_repository_info_dict())
    if data_table.comment_char:
        sub_dict['comment_char'] = 'comment_char="%s"' % (data_table.comment_char)
    for i, name in enumerate(data_table.get_column_name_list()):
        if name is not None:
            sub_dict['columns'] = "%s\n%s" % (sub_dict['columns'], '<column name="%s" index="%s" />' % (name, i))
    location_file_dir = location_file_dir or app.config.galaxy_data_manager_data_path
    for filename in data_table.filenames.keys():
        sub_dict['file_path'] = os.path.basename(filename)
        sub_dict['file_path'] = os.path.join(location_file_dir, sub_dict['file_path'])  # os.path.abspath?
        if not os.path.exists(sub_dict['file_path']):
            # Create empty file
            log.debug("Attempting to create a missing location file %s." % sub_dict['file_path'])
            open(sub_dict['file_path'], 'wb+').close()
        break
    sub_dict['repo_info'] = repo_info or ''
    return """
            <tables><table name="%(table_name)s" %(comment_char)s>
                %(columns)s
                <file path="%(file_path)s" />
                %(repo_info)s
            </table></tables>
           """ % sub_dict


def galaxy_code_get_available_data_tables(trans):
    # list of data tables
    return [(x, x, False) for x in trans.app.tool_data_tables.get_tables().keys()]


def galaxy_code_get_available_data_table_columns(trans, data_table_name):
    return [(x, x, True) for x in trans.app.tool_data_tables.get(data_table_name).get_column_name_list()]
# --- End Galaxy called Methods ---


def get_data_table_entries(params, galaxy_data_manager_data_path):
    rval = {}
    data_tables = params.get('data_tables', [])
    for data_table in data_tables:
        entry_dict = {}
        for column in data_table.get('columns', []):
            value = column.get('data_table_column_value', '')
            if column.get('is_path', {}).get('is_path_selector') == 'yes' and column.get('is_path', {}).get('abspath') == 'abspath':
                value = os.path.abspath(os.path.join(galaxy_data_manager_data_path, value))
            entry_dict[column.get('data_table_column_name', '')] = value
        data_table_name = data_table['data_table_name']
        rval[data_table_name] = rval.get(data_table_name, [])
        rval[data_table_name].append(entry_dict)
    return rval


def get_file_content(params, target_directory):
    directory_content = params.get('directory_content', [])
    for content in directory_content:
        target_path = os.path.join(target_directory, content.get('subdir', ''))
        try:
            os.makedirs(target_path)
        except OSError:
            pass
        if content.get('file_source', {}).get('file_source_selector') == 'URL':
            (filename, headers) = urlretrieve(content.get('file_source', {}).get('file_URL'))
            try:
                bname = headers['Content-Disposition']
            except KeyError:
                bname = os.path.basename(urlsplit(content.get('file_source', {}).get('file_URL')).path)
        else:
            filename = content.get('file_source', {}).get('file_history')
            bname = os.path.basename(filename)
        file_action = content.get('file_action', {}).get('file_action_selector')
        if file_action == 'unpack':
            unpack_archive(filename, target_path)
        else:
            filename_override = content.get('file_action', {}).get('filename_override')
            if filename_override:
                target_path = os.path.join(target_path, filename_override)
            else:
                target_path = os.path.join(target_path, bname)
            shutil.copyfile(filename, target_path)
    return len(directory_content)


def main():
    parser = optparse.OptionParser()
    parser.add_option('', '--galaxy_data_manager_data_path', dest='galaxy_data_manager_data_path', default='', help='Root path for galaxy_data_manager_data_path')
    (options, args) = parser.parse_args()

    filename = args[0]

    with open(filename) as fh:
        params = json.loads(fh.read())
    target_directory = params['output_data'][0]['extra_files_path']

    data_table_entries = get_data_table_entries(params['param_dict'], options.galaxy_data_manager_data_path)

    # save info to json file
    with open(filename, 'w') as fh:
        fh.write(json.dumps({"data_tables": data_table_entries}, sort_keys=True))

    get_file_content(params['param_dict'], target_directory)


if __name__ == "__main__":
    main()
