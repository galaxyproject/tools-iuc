from galaxy.datatypes.data import Data

import logging

log = logging.getLogger(__name__)

class Hmmer3( Data ):
    """Class for hmmpress database files."""
    file_ext = 'hmmpress'
    allow_datatype_change = False
    composite_type = 'basic'

    def set_peek( self, dataset, is_multi_byte=False ):
        """Set the peek and blurb text."""
        if not dataset.dataset.purged:
            dataset.peek  = "HMMER Binary database"
            dataset.blurb = "HMMER Binary database"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def display_peek( self, dataset ):
        """Create HTML content, used for displaying peek."""
        try:
            return dataset.peek
        except:
            return "HMMER3 database (multiple files)"


    def __init__(self, **kwd):
        Data.__init__(self, **kwd)
        self.add_composite_file('model.hmm.h3m', is_binary=True) # Binary model
        self.add_composite_file('model.hmm.h3i', is_binary=True) # SSI index for binary model
        self.add_composite_file('model.hmm.h3f', is_binary=True) # Profiles (MSV part)
        self.add_composite_file('model.hmm.h3p', is_binary=True) # Profiles (remained)
