#!/usr/bin/env python
# Dan Blankenberg
from __future__ import print_function

import datetime
import logging
import optparse
import os
import shutil
import subprocess
import tempfile
from json import (
    dumps,
    loads
)
from os.path import basename
from xml.etree.ElementTree import tostring
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import urlopen

_log_name = __name__
if _log_name == '__builtin__':
    _log_name = 'toolshed.installed.g2.rsync.data.manager'
log = logging.getLogger( _log_name )

# Get the Data from the Galaxy Project rsync server
RSYNC_CMD = 'rsync'
RSYNC_SERVER = "rsync://datacache.g2.bx.psu.edu/"
LOCATION_DIR = "location"
INDEX_DIR = "indexes"

# Pull the Tool Data Table files from github
# FIXME: These files should be accessible from the rsync server directly.
TOOL_DATA_TABLE_CONF_XML_URLS = {'main': "https://raw.githubusercontent.com/galaxyproject/usegalaxy-playbook/master/env/main/files/galaxy/config/tool_data_table_conf.xml",
                                 'test': "https://raw.githubusercontent.com/galaxyproject/usegalaxy-playbook/master/env/test/files/galaxy/config/tool_data_table_conf.xml" }

# Replace data table source entries with local temporary location
GALAXY_DATA_CANONICAL_PATH = "/galaxy/data/"
TOOL_DATA_TABLE_CONF_XML_REPLACE_SOURCE = '<file path="%slocation/' % ( GALAXY_DATA_CANONICAL_PATH )
TOOL_DATA_TABLE_CONF_XML_REPLACE_TARGET = '<file path="%s/'

# Some basic Caching, so we don't have to reload and download everything every time
CACHE_TIME = datetime.timedelta( minutes=10 )
TOOL_DATA_TABLES_LOADED_BY_URL = {}

# Entries will not be selected by default
DEFAULT_SELECTED = False

# Exclude data managers without 'path' column or that are in the manual exclude list
PATH_COLUMN_NAMES = ['path']
EXCLUDE_DATA_TABLES = []
# TODO: Make additional handler actions available for tables that can't fit into the the basic
# "take the value of path" as a dir and copy contents.
# e.g. mafs. Although this maf table is goofy and doesn't have path defined in <table> def,
# but it does exit in the .loc.


# --- These methods are called by/within the Galaxy Application
def exec_before_job( app, inp_data, out_data, param_dict, tool=None, **kwd ):
    # Look for any data tables that haven't been defined for this data manager before and dynamically add them to Galaxy
    param_dict = dict( **param_dict )
    param_dict['data_table_entries'] = param_dict.get( 'data_table_entries', [] )
    if not isinstance( param_dict['data_table_entries'], list ):
        param_dict['data_table_entries'] = [param_dict['data_table_entries']]
    param_dict['data_table_entries'] = ",".join( param_dict['data_table_entries'] )
    if tool:
        tool_shed_repository = tool.tool_shed_repository
    else:
        tool_shed_repository = None
    tdtm = None
    data_manager = app.data_managers.get_manager( tool.data_manager_id, None )
    data_table_entries = get_data_table_entries( param_dict )
    data_tables = load_data_tables_from_url( data_table_class=app.tool_data_tables.__class__ ).get( 'data_tables' )
    for data_table_name, entries in data_table_entries.items():
        # get data table managed by this data Manager
        has_data_table = app.tool_data_tables.get_tables().get( data_table_name )
        if has_data_table:
            has_data_table = bool( has_data_table.get_filename_for_source( data_manager, None ) )
        if not has_data_table:
            if tdtm is None:
                from tool_shed.tools import data_table_manager
                tdtm = data_table_manager.ToolDataTableManager( app )
                target_dir, tool_path, relative_target_dir = tdtm.get_target_install_dir( tool_shed_repository )
            # Dynamically add this data table
            log.debug( "Attempting to dynamically create a missing Tool Data Table named %s." % data_table_name )
            data_table = data_tables[data_table_name]
            repo_info = tdtm.generate_repository_info_elem_from_repository( tool_shed_repository, parent_elem=None )
            if repo_info is not None:
                repo_info = tostring( repo_info )
            tmp_file = tempfile.NamedTemporaryFile()
            tmp_file.write( get_new_xml_definition( app, data_table, data_manager, repo_info, target_dir ) )
            tmp_file.flush()
            app.tool_data_tables.add_new_entries_from_config_file( tmp_file.name, None, app.config.shed_tool_data_table_config, persist=True )
            tmp_file.close()


def galaxy_code_get_available_data_tables( trans ):
    # list of data tables
    found_tables = get_available_tables( trans )
    rval = [ ( x, x, DEFAULT_SELECTED ) for x in found_tables]
    return rval


def galaxy_code_get_available_data_tables_entries( trans, dbkey, data_table_names ):
    # available entries, optionally filtered by dbkey and table names
    if dbkey in [ None, '', '?' ]:
        dbkey = None
    if data_table_names in [ None, '', '?' ]:
        data_table_names = None
    found_tables = get_available_tables_for_dbkey( trans, dbkey, data_table_names )
    dbkey_text = '(%s) ' % ( dbkey ) if dbkey else ''
    rval = [( "%s%s" % ( dbkey_text, x[0] ), dumps( dict( name=x[0].split( ': ' )[0], entry=x[1] ) ).encode( 'base64' ).rstrip(), DEFAULT_SELECTED ) for x in found_tables.items()]
    return rval

# --- End Galaxy called Methods ---


def rsync_urljoin( base, url ):
    # urlparse.urljoin doesn't work correctly for our use-case
    # probably because it doesn't recognize the rsync scheme
    base = base.rstrip( '/' )
    url = url.lstrip( '/' )
    return "%s/%s" % ( base, url )


def rsync_list_dir( server, dir=None, skip_names=[] ):
    # drwxr-xr-x          50 2014/05/16 20:58:11 .
    if dir:
        dir = rsync_urljoin( server, dir )
    else:
        dir = server
    rsync_response = tempfile.NamedTemporaryFile()
    rsync_stderr = tempfile.NamedTemporaryFile()
    rsync_cmd = [ RSYNC_CMD, '--list-only', dir ]
    return_code = subprocess.call( rsync_cmd, stdout=rsync_response, stderr=rsync_stderr )
    rsync_response.flush()
    rsync_response.seek(0)
    rsync_stderr.flush()
    rsync_stderr.seek(0)
    if return_code:
        msg = "stdout:\n%s\nstderr:\n%s" % ( rsync_response.read(), rsync_stderr.read() )
        rsync_response.close()
        rsync_stderr.close()
        raise Exception( 'Failed to execute rsync command (%s), returncode=%s. Rsync_output:\n%s' % (  rsync_cmd, return_code, msg ) )
    rsync_stderr.close()
    rval = {}
    for line in rsync_response:
        perms, line = line.split( None, 1 )
        line = line.strip()
        size, line = line.split( None, 1 )
        line = line.strip()
        date, line = line.split( None, 1 )
        line = line.strip()
        time, line = line.split( None, 1 )
        name = line.strip()
        if name in skip_names:
            continue
        size = line.strip()
        rval[ name ] = dict( name=name, permissions=perms, bytes=size, date=date, time=time )
    rsync_response.close()
    return rval


def rsync_sync_to_dir( source, target ):
    rsync_response = tempfile.NamedTemporaryFile()
    rsync_stderr = tempfile.NamedTemporaryFile()
    rsync_cmd = [ RSYNC_CMD, '-avzP', source, target ]
    return_code = subprocess.call( rsync_cmd, stdout=rsync_response, stderr=rsync_stderr )
    rsync_response.flush()
    rsync_response.seek(0)
    rsync_stderr.flush()
    rsync_stderr.seek(0)
    if return_code:
        msg = "stdout:\n%s\nstderr:\n%s" % ( rsync_response.read(), rsync_stderr.read() )
        rsync_response.close()
        rsync_stderr.close()
        raise Exception( 'Failed to execute rsync command (%s), returncode=%s. Rsync_output:\n%s' % (  rsync_cmd, return_code, msg ) )
    rsync_response.close()
    rsync_stderr.close()
    return return_code


def data_table_needs_refresh( cached_data_table, url ):
    if cached_data_table is None:
        return True, {}
    if datetime.datetime.now() - cached_data_table.get( 'time_loaded' ) > CACHE_TIME:
        data_table_text = urlopen( url ).read()
        if cached_data_table.get( 'data_table_text', None ) != data_table_text:
            return True, {'data_table_text': data_table_text}
        loc_file_attrs = rsync_list_dir( RSYNC_SERVER, LOCATION_DIR )
        if cached_data_table.get( 'loc_file_attrs', None ) != loc_file_attrs:
            return True, {'loc_file_attrs': loc_file_attrs}
    return False, {}


def load_data_tables_from_url( url=None, site='main', data_table_class=None  ):
    if not url:
        url = TOOL_DATA_TABLE_CONF_XML_URLS.get( site, None )
    assert url, ValueError( 'You must provide either a URL or a site=name.' )

    cached_data_table = TOOL_DATA_TABLES_LOADED_BY_URL.get( url, None )
    refresh, attribs = data_table_needs_refresh( cached_data_table, url )
    if refresh:
        data_table_text = attribs.get( 'data_table_text' )or urlopen( url ).read()
        loc_file_attrs = attribs.get( 'loc_file_attrs' ) or rsync_list_dir( RSYNC_SERVER, LOCATION_DIR )

        tmp_dir = tempfile.mkdtemp( prefix='rsync_g2_' )
        tmp_loc_dir = os.path.join( tmp_dir, 'location' )
        os.mkdir( tmp_loc_dir, mode=0o755 )
        rsync_sync_to_dir( rsync_urljoin( RSYNC_SERVER, LOCATION_DIR ), os.path.abspath( tmp_loc_dir ) )

        new_data_table_text = data_table_text.replace( TOOL_DATA_TABLE_CONF_XML_REPLACE_SOURCE, TOOL_DATA_TABLE_CONF_XML_REPLACE_TARGET % ( tmp_loc_dir ) )
        data_table_fh = tempfile.NamedTemporaryFile( dir=tmp_dir, prefix='rysnc_data_manager_data_table_conf_' )
        data_table_fh.write( new_data_table_text )
        data_table_fh.flush()
        tmp_data_dir = os.path.join( tmp_dir, 'tool-data' )
        os.mkdir( tmp_data_dir, mode=0o755 )
        data_tables = data_table_class( tmp_data_dir, config_filename=data_table_fh.name )
        for name, data_table in list(data_tables.data_tables.items()):
            if name in EXCLUDE_DATA_TABLES or not data_table_has_path_column( data_table ):
                log.debug( 'Removing data table "%s" because it is excluded by name or does not have a defined "path" column.', name )
                del data_tables.data_tables[name]
        cached_data_table = { 'data_tables': data_tables, 'tmp_dir': tmp_dir, 'data_table_text': data_table_text, 'tmp_loc_dir': tmp_loc_dir, 'loc_file_attrs': loc_file_attrs, 'time_loaded': datetime.datetime.now() }
        TOOL_DATA_TABLES_LOADED_BY_URL[ url ] = cached_data_table
        # delete the files
        data_table_fh.close()
        cleanup_before_exit( tmp_dir )
    return cached_data_table


def data_table_has_path_column( data_table ):
    col_names = data_table.get_column_name_list()
    for name in PATH_COLUMN_NAMES:
        if name in col_names:
            return True
    return False


def get_available_tables( trans ):
    # list of data tables
    data_tables = load_data_tables_from_url( data_table_class=trans.app.tool_data_tables.__class__ )
    return list(data_tables.get( 'data_tables' ).get_tables().keys())


def get_new_xml_definition( app, data_table, data_manager, repo_info=None, location_file_dir=None ):
    sub_dict = { 'table_name': data_table.name, 'comment_char': '', 'columns': '', 'file_path': '' }
    sub_dict.update( data_manager.get_tool_shed_repository_info_dict() )
    if data_table.comment_char:
        sub_dict['comment_char'] = 'comment_char="%s"' % ( data_table.comment_char )
    for i, name in enumerate( data_table.get_column_name_list() ):
        if name is not None:
            sub_dict['columns'] = "%s\n%s" % ( sub_dict['columns'], '<column name="%s" index="%s" />' % ( name, i ) )
    location_file_dir = location_file_dir or app.config.galaxy_data_manager_data_path
    for filename in data_table.filenames.keys():
        sub_dict['file_path'] = basename( filename )
        sub_dict['file_path'] = os.path.join( location_file_dir, sub_dict['file_path'] )  # os.path.abspath?
        if not os.path.exists( sub_dict['file_path'] ):
            # Create empty file
            open( sub_dict['file_path'], 'wb+' ).close()
        break
    sub_dict[ 'repo_info' ] = repo_info or ''
    return """
            <tables><table name="%(table_name)s" %(comment_char)s>
                %(columns)s
                <file path="%(file_path)s" />
                %(repo_info)s
            </table></tables>
           """ % sub_dict


def get_available_tables_for_dbkey( trans, dbkey, data_table_names ):
    data_tables = load_data_tables_from_url( data_table_class=trans.app.tool_data_tables.__class__ )
    rval = {}
    for name, data_table in data_tables.get( 'data_tables' ).get_tables().items():
        if ( not data_table_names or name in data_table_names ):
            # TODO: check that columns are similiar
            if not dbkey:
                entry_getter = data_table.get_named_fields_list()
            else:
                entry_getter = data_table.get_entries( 'dbkey', dbkey, None, default=[] )
            for entry in entry_getter:
                name = "%s: %s" % ( data_table.name, dumps( entry ) )
                rval[name] = entry
    return rval


def split_path_all( path ):
    rval = []
    path = path.rstrip( '/' )
    while True:
        head, tail = os.path.split( path )
        if tail:
            rval.append( tail )
            path = head
        elif head:
            rval.append( head )
            break
        else:
            break
    rval.reverse()
    return rval


def get_data_for_path( path, data_root_dir ):
    # We list dir with a /, but copy data without
    # listing with / gives a . entry when its a dir
    # cloning without the / will copy that whole directory into the target,
    # instead of just that target's contents
    if path.startswith( GALAXY_DATA_CANONICAL_PATH ):
        path = path[ len( GALAXY_DATA_CANONICAL_PATH ):]
    make_path = path
    rsync_source = rsync_urljoin( rsync_urljoin( RSYNC_SERVER, INDEX_DIR ), path )
    if rsync_source.endswith( '/' ):
        rsync_source = rsync_source[:-1]
    try:
        dir_list = rsync_list_dir( rsync_source + "/" )
    except Exception:
        dir_list = None
    while not dir_list or '.' not in dir_list:
        head, tail = os.path.split( make_path )
        if not head:
            head = tail
        make_path = head
        rsync_source = rsync_urljoin( rsync_urljoin( RSYNC_SERVER, INDEX_DIR ), head )  # if we error here, likely due to a connection issue
        if rsync_source.endswith( '/' ):
            rsync_source = rsync_source[:-1]
        dir_list = rsync_list_dir( rsync_source + "/" )
    split_path = split_path_all( make_path )
    target_path = data_root_dir
    for p in split_path[:-1]:
        target_path = os.path.join( target_path, p )
        if not os.path.exists( target_path ):
            os.mkdir( target_path, mode=0o755 )
    rsync_sync_to_dir( rsync_source, target_path )
    return path


def get_data_and_munge_path( data_table_name, data_table_entry, data_root_dir ):
    path_cols = []
    for key, value in data_table_entry.items():
        if key in PATH_COLUMN_NAMES:
            path_cols.append( ( key, value ) )
    if path_cols:
        for col_name, value in path_cols:
            if value.startswith( GALAXY_DATA_CANONICAL_PATH ):
                data_table_entry[col_name] = get_data_for_path( value, data_root_dir )
            else:
                print('unable to determine location of rsync data for', data_table_name, data_table_entry)
    return data_table_entry


def fulfill_data_table_entries( data_table_entries, data_manager_dict, data_root_dir ):
    for data_table_name, entries in data_table_entries.items():
        for entry in entries:
            entry = get_data_and_munge_path( data_table_name, entry, data_root_dir )
            _add_data_table_entry( data_manager_dict, data_table_name, entry )
    return data_manager_dict


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][data_table_name].append( data_table_entry )
    return data_manager_dict


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )


def get_data_table_entries( params ):
    rval = {}
    data_table_entries = params.get( 'data_table_entries', None )
    if data_table_entries:
        for entry_text in data_table_entries.split( ',' ):
            entry_text = entry_text.strip().decode( 'base64' )
            entry_dict = loads( entry_text )
            data_table_name = entry_dict['name']
            data_table_entry = entry_dict['entry']
            rval[ data_table_name ] = rval.get( data_table_name, [] )
            rval[ data_table_name ].append( data_table_entry )
    return rval


def main():
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    filename = args[0]

    params = loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory, mode=0o755 )
    data_manager_dict = {}

    data_table_entries = get_data_table_entries( params['param_dict'] )

    # Populate the data Tables
    data_manager_dict = fulfill_data_table_entries( data_table_entries, data_manager_dict, target_directory )

    # save info to json file
    open( filename, 'wb' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
