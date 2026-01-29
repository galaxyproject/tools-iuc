# INPUT PARAMETERS:
# - path to the zarr dataset
# - downscale in x,y,z, default value 1,1,1
# - channel and other axes/dimensions that would need to be pinned/fixed
# - time points span (to be able to test on a short range)
#
# OUTPUT PARAMETERS:
# - name of the .csv file into which the tracking would be saved
#   (the .csv together with the original input zarr can be opened Mastodon tracking software)
# - OPTIONAL! path to zarr into which segmentation will be saved
#
# Note: The intermediate segmentation results will not be saved for now...

import numpy as np

default_tracking_options = {
    'downscale_factor_x' : 1.0,
    'downscale_factor_y' : 1.0,
    'downscale_factor_z' : 1.0,
    'start_from_tp'      : 0,
    'end_at_tp'          : -1,
    'segmentation_model' : 'cyto3',
    'tracking_model'     : 'ctc'
}



def flag_error_and_quit(error_msg):
    import sys
    print(f"ERROR: {error_msg}", file=sys.stderr)
    sys.exit(1)


def obtain_lazy_view_from_the_zarr_path(input_path, scale_level, list_of_coords_for_non_tzyx_dims):
    """
    'scale_level' = 0 means the finest/highest (spatial) resolution, the "bottom of a pyramid"
    """
    import ngff_zarr as nz
    zarr_handle = nz.from_ngff_zarr(input_path)

    if scale_level < 0 or scale_level >= len(zarr_handle.images):
        flag_error_and_quit("scale index negative or larger than number(-1) of available resolutions that the zarr dataset offers")

    zarr_image = zarr_handle.images[scale_level]
    #zarr_image.data.shape
    #zarr_image.dims

    axes_known = []
    axes_unknown = []
    curr_axis_idx = 0 #aka dim number
    for d in zarr_image.dims:
        if not d in "tzyx":
            # dimension to be "moved" to the front
            axes_unknown.append(curr_axis_idx)
        else:
            axes_known.append(curr_axis_idx)
        curr_axis_idx += 1

    if len(axes_unknown) != len(list_of_coords_for_non_tzyx_dims):
        flag_error_and_quit(f"found {len(axes_unknown)} non_tzyx dimensions but different number ({len(list_of_coords_for_non_tzyx_dims)}) of values for them")

    axes_permutation = [*axes_unknown, *axes_known]
    # NB: TODO, would be great to check the order in the 'axes_known' and possibly adjust it...
    view = zarr_image.data.transpose(axes_permutation)[*list_of_coords_for_non_tzyx_dims]

    if 'z' not in zarr_image.dims:
        # assuming then tyx, thus injecting 'z':
        view = np.reshape(view, (view.shape[0],1,view.shape[1],view.shape[2]))

    if len(view.shape) != 4:
        #ds = [ zarr_image.dims[n] for n in axes_permutation[-4:] ]
        flag_error_and_quit(f"after fixing non_tzyx dimensions, tzyx (4) dimensions were supposed to be left; instead {len(view.shape)} dimensions are available")

    return view


def segmentation(view_into_raw_data, tracking_options = default_tracking_options):
    """
    Input ('view_into_raw_data') must be t,z,y,x even for 2D+t images.
    Output is a (possibly very large!) numpy with segmentation masks,
    and the corresponding (view into) into the 'view_into_raw_data'.
    Both outputs are (down-)scaled (given 'tracking_options') already!

    Check the 'default_tracking_options' dictionary to see what all keys are supported.
    """
    from cellpose import models as cp3_models
    from skimage.transform import resize
    from math import ceil

    m1 = tracking_options.get('segmentation_model','cyto3')
    seg_model = cp3_models.CellposeModel(model_type=m1)

    # just FYI
    do_3D = view_into_raw_data.shape[1] > 1
    print(f"seg model initiated, going to do 3D: {do_3D}")

    # figure out the (possibly) downscaled spatial size (zyx axes)
    down_scale_factors = [ \
        tracking_options.get('downscale_factor_z',1), \
        tracking_options.get('downscale_factor_y',1), \
        tracking_options.get('downscale_factor_x',1) ]
    new_spatial_size = [ ceil(size/scale) for size,scale in zip(view_into_raw_data[0].shape, down_scale_factors) ]
    #
    do_scaling = min(down_scale_factors) != max(down_scale_factors) != 1
    print(f"seg, going to scale images: {do_scaling}")

    # trim (along the time axis) the input data
    t_from = tracking_options.get('start_from_tp', 0)
    t_to = tracking_options.get('end_at_tp', -1)
    if t_to == -1: t_to = view_into_raw_data.shape[0]-1
    view_into_raw_data = view_into_raw_data[t_from:t_to+1]

    # 'all_masks' will be in the new downscaled size, and the trimmed length!
    print("memory allocation for segmentation results started...")
    all_masks = np.empty((view_into_raw_data.shape[0],*new_spatial_size), dtype='uint16')
    #
    print("memory allocation for raw images started...")
    all_raws = np.empty((view_into_raw_data.shape[0],*new_spatial_size), dtype=view_into_raw_data.dtype)

    print("segmenting started...")
    for t in range(view_into_raw_data.shape[0]):
        img = np.array( resize(view_into_raw_data[t], new_spatial_size, preserve_range=True) ) if do_scaling \
              else np.array(view_into_raw_data[t], dtype=view_into_raw_data.dtype)
        masks,_,_ = seg_model.eval([img], channels=[0,0], z_axis=0, do_3D=do_3D, normalize=True)
        print(f"done segmenting frame {t}, input image size was {img.shape}")

        # btw, it is possible to re-use the memory into which the original zarr data landed
        #img[:] = masks[0,:]
        all_masks[t] = masks[0]
        all_raws[t] = img
    print("segmenting done")

    return all_raws, all_masks


def resize(view_into_raw_data, view_into_seg_data, tracking_options = default_tracking_options):
    from skimage.transform import resize
    from math import ceil

    # figure out the (possibly) downscaled spatial size (zyx axes)
    down_scale_factors = [ \
        tracking_options.get('downscale_factor_z',1), \
        tracking_options.get('downscale_factor_y',1), \
        tracking_options.get('downscale_factor_x',1) ]
    new_spatial_size = [ ceil(size/scale) for size,scale in zip(view_into_raw_data[0].shape, down_scale_factors) ]
    #
    do_scaling = min(down_scale_factors) != max(down_scale_factors) != 1
    print(f"resizing, going to scale images: {do_scaling}")

    # trim (along the time axis) the input data
    t_from = tracking_options.get('start_from_tp', 0)
    t_to = tracking_options.get('end_at_tp', -1)
    if t_to == -1: t_to = view_into_raw_data.shape[0]-1
    view_into_raw_data = view_into_raw_data[t_from:t_to+1]
    view_into_seg_data = view_into_seg_data[t_from:t_to+1]

    # 'all_masks' will be in the new downscaled size, and the trimmed length!
    print("memory allocation for segmentation results started...")
    all_masks = np.empty((view_into_raw_data.shape[0],*new_spatial_size), dtype=view_into_seg_data.dtype)
    #
    print("memory allocation for raw images started...")
    all_raws = np.empty((view_into_raw_data.shape[0],*new_spatial_size), dtype=view_into_raw_data.dtype)

    print("resizing started...")
    for t in range(view_into_raw_data.shape[0]):
        all_raws[t] = np.array( resize(view_into_raw_data[t], new_spatial_size, preserve_range=True) ) if do_scaling \
              else np.array(view_into_raw_data[t], dtype=view_into_raw_data.dtype)
        all_masks[t] = np.array( resize(view_into_seg_data[t], new_spatial_size, preserve_range=True, order=0) ) if do_scaling \
              else np.array(view_into_seg_data[t], dtype=view_into_seg_data.dtype)
        print(f"done resizing frame {t}, target image size was {all_raws[t].shape} (source size was {view_into_raw_data[t].shape})")
    print("resizing done")

    return all_raws, all_masks


def tracking(view_into_raw_data, seg_data, tracking_options = default_tracking_options):
    """
    Both inputs ('view_into_raw_data' and 'seg_data') must be t,z,y,x even for 2D+t images,
    and of the same shapes.
    Output is that of Trackastra, and napari tracks; both with possibly downscaled spatial
    coordinates (depending on the 'tracking_options').

    Check the 'default_tracking_options' dictionary to see what all keys are supported.
    """
    from trackastra.model import Trackastra
    from trackastra.tracking import graph_to_napari_tracks

    m2 = tracking_options.get('tracking_model','ctc')
    tra_model = Trackastra.from_pretrained(m2)

    print("tracking started...")
    track_graph = tra_model.track(view_into_raw_data, seg_data, mode="greedy")  # or mode="ilp", or "greedy_nodiv"
    print("tracking done")

    # TODO: upscale the coordinates in zyx axes
    # consuider also tracking_options.start_from_tp to offset the 0-based time coordinate of the 'view_into_data'
    return track_graph, graph_to_napari_tracks(track_graph)


def upscale_trackastra_graph(track_graph, tracking_options = default_tracking_options):
    down_scale_factors = [ \
        tracking_options.get('downscale_factor_z',1), \
        tracking_options.get('downscale_factor_y',1), \
        tracking_options.get('downscale_factor_x',1) ]

    nodes = track_graph.nodes()
    d = nodes.data()
    for idx in nodes.keys():
        node = d[int(idx)]
        orig_coords = node['coords']
        new_coords = ( \
            orig_coords[0] * down_scale_factors[0], \
            orig_coords[1] * down_scale_factors[1], \
            orig_coords[2] * down_scale_factors[2] )
        node['coords'] = new_coords

    return track_graph


def upscale_napari_tracks(ntracks, tracking_options = default_tracking_options):
    down_scale_factors = [ \
        tracking_options.get('downscale_factor_z',1), \
        tracking_options.get('downscale_factor_y',1), \
        tracking_options.get('downscale_factor_x',1) ]

    for i in range(len(ntracks[0])):
        tracking_marker = ntracks[0][i]
        tracking_marker[2] *= down_scale_factors[0]
        tracking_marker[3] *= down_scale_factors[1]
        tracking_marker[4] *= down_scale_factors[2]

    return ntracks


def postprocess_and_save_tracking(t, tracking_options = default_tracking_options):
    """
    't' should be the touple that's returned from the tracking() function.
    """
    # just unpack
    track_graph, ntracks = t

    track_graph = upscale_trackastra_graph(track_graph, tracking_options)
    ntracks = upscale_napari_tracks(ntracks, tracking_options)

    import pickle
    f = open("trackastra_trackgraph.pickle.dat","wb")
    pickle.dump(track_graph, f)
    f.close()

    f = open("napari_tracks.pickle.dat","wb")
    pickle.dump(ntracks, f)
    f.close()


def segment_and_track_entry(zarr_path: str, scale_level: int,
                            list_of_coords_for_non_tzyx_dims_to_reach_raw_channel: list[int],
                            tracking_options = default_tracking_options):
    """
    Check the 'default_tracking_options' dictionary to see what all keys are supported.
    It is worthwhile to downscale in x,y,z if the input images are 500+ pixels per dimension.
    """
    raw_data_view = obtain_lazy_view_from_the_zarr_path(zarr_path, scale_level, list_of_coords_for_non_tzyx_dims_to_reach_raw_channel)
    # NB: now the data_view is guaranteed to be order as: tzyx
    #     and it is truly an unmodified view, not scaled, not trimmed
    #
    raw,seg = segmentation(raw_data_view, tracking_options)
    t = tracking(raw,seg, tracking_options)
    #
    postprocess_and_save_tracking(t, tracking_options)


def track_entry(zarr_path: str, scale_level: int,
                list_of_coords_for_non_tzyx_dims_to_reach_raw_channel: list[int],
                list_of_coords_for_non_tzyx_dims_to_reach_seg_channel: list[int],
                tracking_options = default_tracking_options):
    """
    Check the 'default_tracking_options' dictionary to see what all keys are supported.
    It is worthwhile to downscale in x,y,z if the input images are 500+ pixels per dimension
    (yes, even for the tracking itself!).
    """
    raw_data_view = obtain_lazy_view_from_the_zarr_path(zarr_path, scale_level, list_of_coords_for_non_tzyx_dims_to_reach_raw_channel)
    seg_data_view = obtain_lazy_view_from_the_zarr_path(zarr_path, scale_level, list_of_coords_for_non_tzyx_dims_to_reach_seg_channel)
    # NB: now the data_view is guaranteed to be order as: tzyx
    #     and it is truly an unmodified view, not scaled, not trimmed
    #
    raw,seg = resize(raw_data_view, seg_data_view, tracking_options)
    t = tracking(raw,seg, tracking_options)
    #
    postprocess_and_save_tracking(t, tracking_options)



def example():
    testing_zarr_path = 'https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.5/idr0051/180712_H2B_22ss_Courtney1_20180712-163837_p00_c00_preview.zarr/0'
    #
    tracking_options = default_tracking_options.copy()
    tracking_options['downscale_factor_x'] = 2.0
    tracking_options['downscale_factor_y'] = 2.0
    tracking_options['downscale_factor_z'] = 1.5
    tracking_options['start_from_tp'] = 0
    tracking_options['end_at_tp'] = 5
    #
    # axes of the 'testing_zarr_path' are: 't', 'c', 'z', 'y', 'x'; and just one channel ('c')
    list_of_coords_for_non_tzyx_dims_to_reach_raw_channel = [0] # to choose the channel
    segment_and_track_entry(testing_zarr_path,0, list_of_coords_for_non_tzyx_dims_to_reach_raw_channel, tracking_options)


if __name__ == '__main__':
    import argparse
    import json
    import sys
    
    parser = argparse.ArgumentParser(
        description='Trackastra: track cell instances in time-lapse movies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Segment and track from a zarr dataset
  python trackastra_wrapper.py segment_and_track \\
    --zarr_path /path/to/data.zarr
    
  # Track with pre-existing segmentation
  python trackastra_wrapper.py track \\
    --zarr_path /path/to/data.zarr
        '''
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Segment and track command
    seg_track_parser = subparsers.add_parser('segment_and_track',
                                              help='Segment cells and perform tracking')
    seg_track_parser.add_argument('--zarr_path', required=True,
                                  help='Path to the zarr dataset (URL, s3://, or local path)')
    seg_track_parser.add_argument('--scale_level', type=int, default=0,
                                  help='Pyramid scale level to use (default: 0 = highest resolution)')
    seg_track_parser.add_argument('--channel_coords', type=str, default='0',
                                  help='Coordinates to select non-tzyx dimensions (comma/space-separated, default: 0)')
    seg_track_parser.add_argument('--downscale_x', type=float, default=1.0,
                                  help='Downscale factor for X axis (default: 1.0)')
    seg_track_parser.add_argument('--downscale_y', type=float, default=1.0,
                                  help='Downscale factor for Y axis (default: 1.0)')
    seg_track_parser.add_argument('--downscale_z', type=float, default=1.0,
                                  help='Downscale factor for Z axis (default: 1.0)')
    seg_track_parser.add_argument('--start_tp', type=int, default=0,
                                  help='Starting time point (default: 0)')
    seg_track_parser.add_argument('--end_tp', type=int, default=-1,
                                  help='Ending time point (default: -1 = all)')
    seg_track_parser.add_argument('--segmentation_model', type=str, default='cyto3',
                                  choices=['cyto3', 'cyto2', 'nuclei'],
                                  help='Cellpose segmentation model to use (default: cyto3)')
    seg_track_parser.add_argument('--tracking_model', type=str, default='ctc',
                                  help='Trackastra tracking model to use (default: ctc)')
    
    # Track only command
    track_parser = subparsers.add_parser('track',
                                         help='Perform tracking on pre-segmented data')
    track_parser.add_argument('--zarr_path', required=True,
                              help='Path to the zarr dataset (URL, s3://, or local path)')
    track_parser.add_argument('--scale_level', type=int, default=0,
                              help='Pyramid scale level to use (default: 0 = highest resolution)')
    track_parser.add_argument('--raw_channel_coords', type=str, default='0',
                              help='Coordinates to select raw data non-tzyx dimensions (comma/space-separated, default: 0)')
    track_parser.add_argument('--seg_channel_coords', type=str, default='0',
                              help='Coordinates to select segmentation non-tzyx dimensions (comma/space-separated, default: 0)')
    track_parser.add_argument('--downscale_x', type=float, default=1.0,
                              help='Downscale factor for X axis (default: 1.0)')
    track_parser.add_argument('--downscale_y', type=float, default=1.0,
                              help='Downscale factor for Y axis (default: 1.0)')
    track_parser.add_argument('--downscale_z', type=float, default=1.0,
                              help='Downscale factor for Z axis (default: 1.0)')
    track_parser.add_argument('--start_tp', type=int, default=0,
                              help='Starting time point (default: 0)')
    track_parser.add_argument('--end_tp', type=int, default=-1,
                              help='Ending time point (default: -1 = all)')
    track_parser.add_argument('--tracking_model', type=str, default='ctc',
                              help='Trackastra tracking model to use (default: ctc)')
    
    args = parser.parse_args()
    
    def parse_coords(coord_string):
        """Parse comma or space-separated coordinate string to list of ints"""
        if isinstance(coord_string, str):
            # Replace commas with spaces and split
            coords = coord_string.replace(',', ' ').split()
            return [int(c) for c in coords if c.strip()]
        return [int(coord_string)] if coord_string else [0]
    
    if args.command == 'segment_and_track':
        tracking_options = default_tracking_options.copy()
        tracking_options['downscale_factor_x'] = args.downscale_x
        tracking_options['downscale_factor_y'] = args.downscale_y
        tracking_options['downscale_factor_z'] = args.downscale_z
        tracking_options['start_from_tp'] = args.start_tp
        tracking_options['end_at_tp'] = args.end_tp
        tracking_options['segmentation_model'] = args.segmentation_model
        tracking_options['tracking_model'] = args.tracking_model
        
        channel_coords = parse_coords(args.channel_coords)
        segment_and_track_entry(args.zarr_path, args.scale_level,
                               channel_coords, tracking_options)
        print("Tracking completed successfully")
    
    elif args.command == 'track':
        tracking_options = default_tracking_options.copy()
        tracking_options['downscale_factor_x'] = args.downscale_x
        tracking_options['downscale_factor_y'] = args.downscale_y
        tracking_options['downscale_factor_z'] = args.downscale_z
        tracking_options['start_from_tp'] = args.start_tp
        tracking_options['end_at_tp'] = args.end_tp
        tracking_options['tracking_model'] = args.tracking_model
        
        raw_coords = parse_coords(args.raw_channel_coords)
        seg_coords = parse_coords(args.seg_channel_coords)
        track_entry(args.zarr_path, args.scale_level,
                   raw_coords, seg_coords, tracking_options)
        print("Tracking completed successfully")
    
    else:
        parser.print_help()
        sys.exit(1)

