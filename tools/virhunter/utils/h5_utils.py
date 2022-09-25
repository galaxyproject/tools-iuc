#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from pathlib import Path

import h5py
import numpy as np


def add_large_kernel(file_path: str):
    filename_wo_ext = os.path.basename(file_path).split('.h5')[0]
    dirname = os.path.dirname(file_path)
    if filename_wo_ext in ['model_10_500', 'model_10_1000']:
        test_array = np.zeros((1024, 512))
        with h5py.File(dirname + '/' + filename_wo_ext + "_kernel_0_1.h5", "r") as split1, \
                h5py.File(dirname + '/' + filename_wo_ext + "_kernel_0_2.h5", "r") as split2, \
                h5py.File(dirname + '/' + filename_wo_ext + "_kernel_0_3.h5", "r") as split3, \
                h5py.File(dirname + '/' + filename_wo_ext + "_kernel_0_4.h5", "r") as split4, \
                h5py.File(file_path, "a") as final_dest:
            stem = 'dense_' + str(int(final_dest.attrs['layer_names'][2].decode(encoding='utf-8', ).split("_")[1]) * 2)
            test_array[0:512, 0:256] = np.array(split1.get(stem + '/' + stem + '/kernel:0:1'))
            test_array[512:1024, 0:256] = np.array(split2.get(stem + '/' + stem + '/kernel:0:2'))
            test_array[0:512, 256:512] = np.array(split3.get(stem + '/' + stem + '/kernel:0:3'))
            test_array[512:1024, 256:512] = np.array(split4.get(stem + '/' + stem + '/kernel:0:4'))
            current_group = final_dest[stem + '/' + stem]
            current_group.create_dataset('kernel:0', data=test_array, dtype='float32')
            final_dest.close()


def remove_large_kernel(file_path: str):
    # Delete dataset - no problems cause it is the last added dataset
    filename_wo_ext = os.path.basename(file_path).split('.h5')[0]
    dirname = os.path.dirname(file_path)
    stem = dirname + '/' + filename_wo_ext
    if filename_wo_ext in ['model_10_500', 'model_10_1000']:
        dest_final = h5py.File(Path(stem + "_light.h5"), "w")
        with h5py.File(file_path, "r") as src:
            for k in src.attrs.keys():
                dest_final.attrs[k] = src.attrs[k]
            for group in src:
                cg = dest_final.create_group(group)
                for group_k in src[group].attrs.keys():
                    cg.attrs[group_k] = src[group].attrs[group_k]
                for subgroup in src[group]:
                    subcg = cg.create_group(subgroup)
                    for subsubgroup in src[group][subgroup]:
                        path = group + '/' + subgroup + '/' + subsubgroup
                        dataset = src.get(path)
                        dataset_arr = np.array(dataset)
                        if dataset_arr.ndim == 2:
                            if dataset_arr.shape[0] != 1024 and dataset_arr.shape[1] != 512:
                                subcg.create_dataset(subsubgroup, data=dataset)
                        else:
                            subcg.create_dataset(subsubgroup, data=dataset)
            dest_final.close()

        os.remove(file_path)
        os.rename(Path(stem + "_light.h5"), file_path)


fn = "/Users/benjamin/virhunter_iuc/tools-iuc/tools/virhunter/tool-data/weights/peach/model_10_500.h5"
#add_large_kernel(fn)
remove_large_kernel(fn)