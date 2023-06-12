import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Command parameters for squidpy", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--visium", help="Visium data")
parser.add_argument("-g", "--genome", help="Optional Metadata from visium")
parser.add_argument("-c", "--count", help="Count file found in command line")
parser.add_argument("-l", "--library", help="Optional visium dtata library ID")
parser.add_argument("-s", "--source-image", help="Path to high-resolution tissue images (*/outs/spatial/tissue_hires_image.png)")
parser.add_argument("-x", "--prefix", help="Prefix label by which data will be separated (ex. 'MT-')")
parser.add_argument("-m", "--min-counts", help="Minimum cells to filter for")
parser.add_argument("-M", "--max-counts", help="Maximum cells to filter for")
parser.add_argument("-p", "--pct-counts", help="Percent of cells meeting minimum counts")
parser.add_argument("-e", "--min-cells", help="Minimum cells requred to pass")
parser.add_argument("-n", "--n-top-genes", help="Keep top N genes")
args = parser.parse_args()

datapath = args.visium
genome_meta = args.genome
countfile = args.count
libid = args.library
imgs = args.source_image
min_counts = int(args.min_counts)
max_counts = int(args.max_counts)
pct_counts = float(args.pct_counts)
min_cells = int(args.min_cells)
n_top_genes = int(args.n_top_genes)
prefix = args.prefix

visium_readin = sc.read_visium(path=datapath, 
                               genome=genome_meta, 
                               count_file=countfile, 
                               library_id=libid, 
                               source_image_path=imgs)
visium_readin.var_names_make_unique()
visium_readin.var[prefix] = visium_readin.var_names.str.startswith(prefix)
sc.pp.calculate_qc_metrics(visium_readin, qc_vars=[prefix], inplace=True)
sc.pp.filter_cells(visium_readin, min_counts=min_counts)
sc.pp.filter_cells(visium_readin, max_counts=max_counts)
counts_var = "pct_counts_" + str(prefix)
visium_readin = visium_readin[visium_readin.obs[counts_var] < pct_counts]
sc.pp.filter_genes(visium_readin, min_cells=min_cells)
sc.pp.normalize_total(visium_readin, inplace=True)
sc.pp.log1p(visium_readin)
sc.pp.highly_variable_genes(visium_readin, flavor="seurat", n_top_genes=n_top_genes)
sc.pp.pca(visium_readin)
sc.pp.neighbors(visium_readin)
sc.tl.umap(visium_readin)
sc.tl.leiden(visium_readin, key_added="clusters")
visium_readin.write_h5ad(filename="squidpy_input.h5ad")