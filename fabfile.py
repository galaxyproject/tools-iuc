from fabric.api import local
from fabric.operations import prompt
import os


def package():
    '''
    Pacakge for upload to toolshed
      packaging 'test' version (default) uses local directory
      otherwise, specify a mercurial tag to package
    '''
    package_dir = 'package'
    base_filename = os.path.join(package_dir, 'htseq-count')
    version = prompt("Enter version number for package [test]:")
    revision_option = ''
    if version != '':
        revision_option = '-r "%s"' % version
    else:
        version = 'test'
    version_filename = '%s_%s.tar.gz' % (base_filename, version)
    local('mkdir -p %s' % package_dir)
    local('rm -f %s' % version_filename)
    if version == 'test':
        local('tar czvf %s --exclude "fabfile.*" --exclude "%s" --exclude ".hg*" --exclude ".DS_Store" --exclude "*.pyc" --exclude "*.swp" *' % (version_filename, package_dir))
    else:
        local('hg archive -t tgz %s -X "fabfile.*" -X "package" -X ".hg*" -p . "%s"' % (revision_option, version_filename))
