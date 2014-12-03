from galaxy.datatypes.metadata import MetadataElement, MetadataParameter

from galaxy.datatypes.binary import SQlite, Binary, sqlite, data, dataproviders 

#@dataproviders.decorators.has_dataproviders
class GeminiSQlite ( SQlite ):
    """Class describing a Gemini Sqlite database """
    MetadataElement( name="gemini_version", default='0.10.0' , param=MetadataParameter, desc="Gemini Version", 
                     readonly=True, visible=True, no_value='0.10.0' )
    file_ext = "gemini.sqlite"

    def set_meta( self, dataset, overwrite = True, **kwd ):
        super( GeminiSQlite, self ).set_meta( dataset, overwrite = overwrite, **kwd )
        try:
            conn = sqlite.connect( dataset.file_name )
            c = conn.cursor()
            tables_query = "SELECT version FROM version"
            result = c.execute( tables_query ).fetchall()
            for version, in result:
                dataset.metadata.gemini_version = version
            # TODO: Can/should we detect even more attributes, such as use of PED file, what was input annotation type, etc.
        except Exception, exc:
            log.warn( '%s, set_meta Exception: %s', self, exc )

    def sniff( self, filename ):
        if super( GeminiSQLite, self ).sniff( filename):
            gemini_table_names = [ "gene_detailed", "gene_summary", "resources", "sample_genotype_counts", "sample_genotypes", "samples", 
                                  "variant_impacts", "variants", "version" ]
            try:
                conn = sqlite.connect( filename )
                c = conn.cursor()
                tables_query = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
                result = c.execute( tables_query ).fetchall()
                result = map( lambda x: x[0] for x in result )
                for table_name in gemini_table_names:
                    if table_name not in result:
                        return False
                return True
            except Exception, exc:
                log.warn( '%s, set_meta Exception: %s', self, exc )
        return False

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek  = "Gemini SQLite Database, version %s" % ( dataset.metadata.gemini_version or 'unknown' )
            dataset.blurb = data.nice_size( dataset.get_size() )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def display_peek( self, dataset ):
        try:
            return dataset.peek
        except:
            return "Gemini SQLite Database, version %s" % ( dataset.metadata.gemini_version or 'unknown' )

Binary.register_sniffable_binary_format( "gemini.sqlite", "gemini.sqlite", GeminiSQlite )
