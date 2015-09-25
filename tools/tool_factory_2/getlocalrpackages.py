import os
import subprocess
import sys


def find_packages(prefix="package_r_"):
    """
    """
    #locate env.sh | grep -i package_r_
    #/data/extended/galaxyJune14_2014/tool_dependency/readline/6.2/devteam/package_r_2_15_0/8ab0d08a3da1/env.sh
    #/data/home/rlazarus/galaxy/tool_dependency_dir/R_3_1_1/3.1.1/fubar/package_r_3_1_1/5f1b8d22140a/env.sh
    #/data/home/rlazarus/galaxy/tool_dependency_dir/R_3_1_1/3.1.1/fubar/package_r_3_1_1/d9964efbfbe3/env.sh
    #/data/home/rlazarus/galtest/tool_dependency_dir/R_3_1_1/3.1.1/fubar/package_r_3_1_1/63cdb9b2234c/env.sh
    eprefix = prefix
    if prefix.find('/') != -1:
        eprefix = prefix.replace('/','\/') # for grep
    path = '.'
    # fails on nitesh's recent mac - locate not working 
    # cl = ['locate env.sh | grep -i %s' % eprefix,]
    cl = ['find %s -iname "env.sh" | grep -i %s' % (path,eprefix),]
    p = subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    out, err = p.communicate()
    fpaths = out.split('\n')
    fpaths = [x for x in fpaths if len(x) > 1]
    fver = [x.split(os.path.sep)[-4:-1] for x in fpaths]
    # >>> foo.split(os.path.sep)[-4:-1]
    # ['fubar', 'package_r_3_1_1', '63cdb9b2234c']
    if len(fpaths) > 0:
        res = [['%s rev %s owner %s' % (x[1],x[2],x[0]),fpaths[i],False] for i,x in enumerate(fver)]
        res[0][2] = True # selected if more than one
        res.insert(0,['Use default (system) interpreter','system',False])
    else:
        res = [['Use default (system) interpreter','system',True],
        ['**WARNING** NO package env.sh files found - is the "find" system command working? Are any interpreters installed?','system',False]]
    # return a triplet - user_sees,value,selected - all unselected if False
    return res

def testapi():
    host_url = 'http://localhost:8080'
    new_path = [ os.path.join( os.getcwd(), "lib" ) ]
    new_path.extend( sys.path[1:] ) # remove scripts/ from the path
    sys.path = new_path
    from galaxy import config
    aconfig = config.Configuration( )
    M_A_K = aconfig.master_api_key
    tooldeps = aconfig.tool_dependency_dir
    gi = GalaxyInstance(url=host_url, key=M_A_K)


if __name__ == "__main__":
    print(find_packages())

