import scanpy as sc
import anndata as ad
import squidpy as sq
import matplotlib.pyplot as plt
# import matplotlib.image
import argparse
import pandas as pd

def cluster_features(features: pd.DataFrame, like=None):
    """
    Calculate leiden clustering of features.

    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs["leiden"]

parser = argparse.ArgumentParser(description="Command parameters for squidpy", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--annData", help="annData object")
parser.add_argument("-i", "--library_id", help="name of library id")
args = parser.parse_args()

datapath = args.annData
library_id=args.library_id
# Load the pre-processed dataset
adata=ad.read_h5ad(datapath)

# /////////////////////////////////////////////////////////////////////////////

scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
img = sq.im.ImageContainer(adata.uns['spatial'][library_id]['images']['hires'],
                           scale=scale, library_id=library_id)
img.show(channelwise=True) #ADD IMG.SHOW TO PDF FILE

# array = img.compute('image')

# matplotlib.image.imsave('name.png', array)
# # Image Processing for Segmentation
sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0, chunks=1000)

# Plot the resulting segmentation
fig, ax = plt.subplots(1, 2)
img_crop = img.crop_corner(0,0, size=500)
img_crop.show(layer="image", channel=0, ax=ax[0])
img_crop.show(
    layer="segmented_watershed",
    channel=0,
    ax=ax[1],
)

# define image layer to use for segmentation
features_kwargs = {"segmentation": {"label_layer": "segmented_watershed"}}
# calculate segmentation features
library_id = 'None'
sq.im.calculate_image_features(
    adata,
    img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    n_jobs=1,
    features_kwargs=features_kwargs,
)

#Create list of features to be looped through
features = adata.obsm["features_segmentation"]
features = [i for i in features if "intensity_mean" in i]
features.append("segmentation_label")
features.append("clusters")
features

# plot results and compare with gene-space clustering using above features list
sq.pl.spatial_scatter(
    sq.pl.extract(adata, "features_segmentation"),
    color=features,
    frameon=False,
    ncols=2,
)

# define different feature calculation combinations
scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_lowres_scalef']
params = {
    # all features, corresponding only to tissue underneath spot
    # "features_orig": {
    #     "features": ["summary", "texture", "histogram"],
    #     "scale": 1.0,
    #     "mask_circle": True,
    # },
    # summary and histogram features with a bit more context, original resolution
    "features_context": {"features": ["summary", "histogram"], "scale": 1.0},
    # summary and histogram features with more context and at lower resolution
    "features_lowres": {"features": ["summary", "histogram"], "scale" : 0.25},
        #Scale is not working with anything but 1.0 due to xarray issues therefore 
        # we will not use it in this iteration but leave it here for future updates
}

for feature_name, cur_params in params.items():
    # features will be saved in `adata.obsm[feature_name]`
    sq.im.calculate_image_features(adata, img, layer="image", key_added=feature_name, n_jobs=1, **cur_params)

# # combine features in one dataframe
# adata.obsm["features"] = pd.concat([adata.obsm[f] for f in params.keys()], axis="columns")
# # make sure that we have no duplicated feature names in the combined table
# adata.obsm["features"].columns = ad.utils.make_index_unique(adata.obsm["features"].columns)

# #Calculte feature clusters useing different features and compare them to gene clusters
# adata.obs["features_summary_cluster"] = cluster_features(adata.obsm["features"], like="summary")
# adata.obs["features_histogram_cluster"] = cluster_features(adata.obsm["features"], like="histogram")
# adata.obs["features_texture_cluster"] = cluster_features(adata.obsm["features"], like="texture")

# # Plot cluster features
# sc.set_figure_params(facecolor="white", figsize=(8, 8))
# sq.pl.spatial_scatter(
#     adata,
#     color=[
#         "features_summary_cluster",
#         "features_histogram_cluster",
#         "features_texture_cluster",
#         "clusters",
#     ],
#     ncols=3,
# )

# # View areas of similarity using heatmap of neighborhood enrichment
# sq.gr.spatial_neighbors(adata)
# sq.gr.nhood_enrichment(adata, cluster_key="clusters")
# sq.pl.nhood_enrichment(adata, cluster_key="clusters")


# # Show all clusters and their co-occurance acrouss spatial dimensions
# clusters_list = adata.obs['clusters'].cat.categories
# for i in clusters_list:
#     sq.gr.co_occurrence(adata, cluster_key="clusters")
#     sq.pl.co_occurrence(
#         adata,
#         cluster_key="clusters",
#         clusters=i,
#         figsize=(8, 4),
#     )
