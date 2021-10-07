suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to estimate the effectiveness of proportion methods

## Load Conf
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

## Estimate cell type proportions
est_prop <- music_prop(
    bulk.eset = bulk_eset, sc.eset = scrna_eset,
    clusters = celltypes_label,
    samples = samples_label, select.ct = celltypes, verbose = T)

estimated_music_props <- est_prop$Est.prop.weighted
estimated_nnls_props <- est_prop$Est.prop.allgene

## Show different in estimation methods
## Jitter plot of estimated cell type proportions
jitter.fig <- Jitter_Est(
    list(data.matrix(estimated_music_props),
         data.matrix(estimated_nnls_props)),
    method.name = methods, title = "Jitter plot of Est Proportions")


## Make a Plot
## A more sophisticated jitter plot is provided as below. We separated
## the T2D subjects and normal subjects by their Disease_Factor levels.
estimated_music_props_flat = melt(estimated_music_props)
estimated_nnls_props_flat = melt(estimated_nnls_props)

m_prop <- rbind(estimated_music_props_flat,
                estimated_nnls_props_flat)
colnames(m_prop) <- c("Sub", "CellType", "Prop")

if (is.null(celltypes)){
    celltypes = levels(m_prop$CellType)
    print("No celltypes declared, using:")
    print(celltypes)
}

if (phenotype_target_threshold == -99){
    phenotype_target_threshold = -Inf
    print("phenotype target threshold set to -Inf")
}

if (is.null(phenotype_factors)){
    phenotype_factors = colnames(pData(bulk_eset))
}
## filter out unwanted factors like "sampleID" and "subjectName"
phenotype_factors = phenotype_factors[!(phenotype_factors %in% phenotype_factors_always_exclude)]
print("Phenotype Factors to use:")
print(phenotype_factors)


m_prop$CellType <- factor(m_prop$CellType, levels = celltypes) # nolint
m_prop$Method <- factor(rep(methods, each = nrow(estimated_music_props_flat)), # nolint
                        levels = methods)
m_prop$Disease_factor <- rep(bulk_eset[[phenotype_target]], 2 * length(celltypes)) # nolint
m_prop <- m_prop[!is.na(m_prop$Disease_factor), ]
## Generate a TRUE/FALSE table of Normal == 1 and Disease == 2
sample_groups = c("Normal", sample_disease_group)
m_prop$Disease <- factor(sample_groups[(m_prop$Disease_factor > phenotype_target_threshold) + 1], # nolint
                         levels = sample_groups)

## Binary to scale: e.g. TRUE / 5 = 0.2
m_prop$D <- (m_prop$Disease ==   # nolint
             sample_disease_group) / sample_disease_group_scale
## NA's are not included in the comparison below
m_prop <- rbind(subset(m_prop, Disease != sample_disease_group),
               subset(m_prop, Disease == sample_disease_group))

jitter.new <- ggplot(m_prop, aes(Method, Prop)) +
    geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease),
               size = 2, alpha = 0.7,
               position = position_jitter(width = 0.25, height = 0)) +
    facet_wrap(~ CellType, scales = "free") +
    scale_colour_manual(values = c("white", "gray20")) +
    scale_shape_manual(values = c(21, 24)) + theme_minimal()

## Plot to compare method effectiveness
## Create dataframe for beta cell proportions and Disease_factor levels
m_prop_ana <- data.frame(pData(bulk_eset)[rep(1:nrow(estimated_music_props), 2),
                                          phenotype_factors],
                        ct.prop = c(estimated_music_props[, 2],
                                    estimated_nnls_props[, 2]),
                        Method = factor(rep(methods, each = nrow(estimated_music_props)),
                                        levels = methods))
colnames(m_prop_ana)[1:length(phenotype_factors)] <- phenotype_factors
m_prop_ana <- subset(m_prop_ana, !is.na(m_prop_ana[phenotype_target]))
m_prop_ana$Disease <- factor(sample_groups[(  # nolint
    m_prop_ana[phenotype_target] > phenotype_target_threshold) + 1], sample_groups)
m_prop_ana$D <- (m_prop_ana$Disease ==        # nolint
                 sample_disease_group) / sample_disease_group_scale

jitt_compare <- ggplot(m_prop_ana, aes_string(phenotype_target, "ct.prop")) +
    geom_smooth(method = "lm",  se = FALSE, col = "black", lwd = 0.25) +
    geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease),
               size = 2, alpha = 0.7) +  facet_wrap(~ Method) +
    ggtitle(compare_title) + theme_minimal() +
    scale_colour_manual(values = c("white", "gray20")) +
    scale_shape_manual(values = c(21, 24))


pdf(file = outfile_pdf, width = 8, height = 8)
plot_grid(jitter.fig, jitter.new, labels = "auto", ncol = 1, nrow = 2)
jitt_compare
dev.off()

## Summary table
for (meth in methods) {
    ##lm_beta_meth = lm(ct.prop ~ age + bmi + hba1c + gender, data =
    sub_data = subset(m_prop_ana, Method == meth)
    ## We can only do regression where there are more than 1 factors
    ## so we must find and exclude the ones which are not
    gt1_facts = sapply(phenotype_factors, function(facname){
        return(length(unique(sort(sub_data[[facname]])))==1)
    })
    form_factors = phenotype_factors
    exclude_facts = names(gt1_facts)[gt1_facts]
    if (length(exclude_facts) > 0){
        print("Factors with only one level will be excluded:")
        print(exclude_facts)
        form_factors = phenotype_factors[!(phenotype_factors %in% exclude_facts)]
    }
    lm_beta_meth <- lm(as.formula(
        paste("ct.prop", paste(form_factors, collapse = " + "),
              sep = " ~ ")), data = sub_data)
    print(paste0("Summary: ", meth))
    capture.output(summary(lm_beta_meth),
                   file = paste0("report_data/summ_", meth, ".txt"))
}
