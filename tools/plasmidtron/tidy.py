#!/usr/bin/env python

import sys
import os

def make_list(tmp_file, lst_file):
    dividers = ["_1", "_F", "_R1", "_2", "_R", "_R2"]
    with open(tmp_file, "r") as fh:
        lines = fh.readlines()
        _lines = []
        processed_list = []
        for i in range(0,len(lines)):
            if (os.path.splitext(lines[i].strip())[-1] == ".fasta") or (os.path.splitext(lines[i].strip())[-1] == ".fa"):
                _lines.append(lines[i].strip())
            else:
                file_name = os.path.basename(lines[i])
                if file_name not in processed_list:
                    detected_devider = [div for div in dividers if div in file_name]
                    if len(detected_devider) > 0: 
                        detected_devider = detected_devider[0]
                        new_file_name = file_name.split(detected_devider)[0]
                        for j in range(i+1, len(lines)):
                            if new_file_name in lines[j]:
                                paired = "{},{}".format(lines[i].strip(), lines[j].strip())
                                _lines.append(paired)
                                processed_list.append(file_name)
    
    with open(lst_file,"w") as fh:
        for _line in _lines:
            fh.writelines("{}\n".format(_line))
                    
if __name__ == "__main__":
    tmp_file = sys.argv[1]
    lst_file = sys.argv[2]
    make_list(tmp_file, lst_file)
