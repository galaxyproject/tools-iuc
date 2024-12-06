
            # rctd script
            # This file is used to specify the parameters for the rctd from spacexr package

            # Load the spacexr library
            library('spacexr')
            library('Matrix')

            counts <- read.table(file = 'inputs/sc_count.tabular', row.names = 1, sep = '\t', header = T)
            metadata <- read.table(file = 'inputs/metadata.tabular', sep = '\t', header = T)

            # create cell_types named list
            cell_types <- metadata[,"annotation"]; names(cell_types) <- metadata[,"barcode"]
            
            # convert to factor data type
            cell_types <- as.factor(cell_types)

                # create nUMI named list
                nUMI <- metadata[, "nUMI"]; names(nUMI) <- metadata[,"barcode"]

            # Create reference object
            reference <- Reference(
                                counts = counts,
                                cell_types = cell_types,
                                    nUMI = nUMI,
                                n_max_cells = 10000,
                                min_UMI = 100
                                )

            counts <- read.table(file = 'inputs/st_count.tabular', row.names = 1, sep = '\t', header = T)
            coords <- read.table(file = 'inputs/coords.tabular', row.names = 1, sep = '\t', header = T)

            nUMI <- colSums(counts) # In tutorials it is always the sum of counts

            # Create SpatialRNA object
            puck <- SpatialRNA(
                            coords = coords,
                            counts = counts,
                            nUMI= nUMI,
                            )

            # provide a basic plot of the nUMI of each pixel on the plot:
            pdf('figures/nUMI_plot.pdf')
            plot_puck_continuous(
                puck = puck,
                barcodes = colnames(puck@counts),
                plot_val = puck@nUMI,
                ylimit = c(0,round(quantile(puck@nUMI,0.9))),
                title ='plot of nUMI')
            dev.off()

            myRCTD <- create.RCTD(
                            spatialRNA = puck,
                            reference = reference,
                            gene_cutoff = 0.000125,
                            fc_cutoff = 0.5,
                            gene_cutoff_reg = 0.0002,
                            fc_cutoff_reg = 0.75,
                            UMI_min = 100,
                            UMI_max = 20000000,
                            counts_MIN = 10,
                            UMI_min_sigma = 300,
                            class_df = NULL, # set as default
                            CELL_MIN_INSTANCE = 25,
                            MAX_MULTI_TYPES = 4,
                            keep_reference = F, # set as default
                            cell_type_profiles = NULL, # set as default
                            CONFIDENCE_THRESHOLD = 5,
                            DOUBLET_THRESHOLD = 20,)

            myRCTD <- run.RCTD(
                            myRCTD,
                            doublet_mode = "doublet")


            # save results
                results <- myRCTD@results

                # save the data frame
                result_df <- results["results_df"]
                write.table(result_df, file = 'results/doublet_results_df.tabular', sep = '\t', quote = F, row.names = T)

                # RCTD plots
                # normalize the cell type proportions to sum to 1.
                norm_weights <- normalize_weights(results[["weights"]])
                cell_type_names <- myRCTD@cell_type_info[["info"]][[2]] #list of cell type names
                spatialRNA <- myRCTD@spatialRNA

                resultsdir <- 'figures'
                
                # make the plots
                # Plots the confident weights for each cell type as in full_mode (saved as 'figures/cell_type_weights.pdf')
                plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
                
                # Plots all weights for each cell type as in full_mode. (saved as 'figures/cell_type_weights_unthreshold.pdf')
                plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
                
                # Plots the weights for each cell type as in doublet_mode. (saved as 'figures/cell_type_weights_doublets.pdf')
                plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results[["weights_doublet"]],results[["results_df"]])

                # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'figures/cell_type_occur.pdf')
                plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

                # makes a map of all cell types, (saved as 'figures/all_cell_types.pdf')
                plot_all_cell_types(results[["results_df"]], spatialRNA@coords, cell_type_names, resultsdir) 

                # doublets
                #obtain a dataframe of only doublets
                doublets <- results[["results_df"]][results[["results_df"]][["spot_class"]] == "doublet_certain",] 
                
                # Plots all doublets in space (saved as 'figures/all_doublets.pdf')
                plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names)
            
                # Plots all doublets in space for each cell type (saved as 'figures/all_doublets_type.pdf')
                plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
                
                # a table of frequency of doublet pairs 
                doub_occur <- table(doublets[["second_type"]], doublets[["first_type"]]) 
                # Plots a stacked bar plot of doublet ocurrences (saved as 'figures/doublet_stacked_bar.pdf')
                plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
            


            # save rds file
                saveRDS(myRCTD, file = 'results/rctd_results.rds')
        