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


## Show different in estimation methods
## Jitter plot of estimated cell type proportions
jitter.fig <- Jitter_Est(
    list(data.matrix(est_prop$Est.prop.weighted),
         data.matrix(est_prop$Est.prop.allgene)),
    method.name = methods, title = "Jitter plot of Est Proportions")


## Make a Plot
## A more sophisticated jitter plot is provided as below. We separated
## the T2D subjects and normal subjects by their HbA1c levels.
m_prop <- rbind(melt(est_prop$Est.prop.weighted),
               melt(est_prop$Est.prop.allgene))

colnames(m_prop) <- c("Sub", "CellType", "Prop")

m_prop$CellType <- factor(m_prop$CellType, levels = celltypes) # nolint
m_prop$Method <- factor(rep(methods, each = 89 * 6), levels = methods) # nolint
m_prop$HbA1c <- rep(bulk_eset$hba1c, 2 * 6) # nolint
m_prop <- m_prop[!is.na(m_prop$HbA1c), ]
m_prop$Disease <- factor(sample_groups[(m_prop$HbA1c > 6.5) + 1], # nolint
                         levels = sample_groups)

m_prop$D <- (m_prop$Disease ==   # nolint
             sample_disease_group) / sample_disease_group_scale
m_prop <- rbind(subset(m_prop, Disease == healthy_phenotype),
               subset(m_prop, Disease != healthy_phenotype))

jitter.new <- ggplot(m_prop, aes(Method, Prop)) +
    geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease),
               size = 2, alpha = 0.7,
               position = position_jitter(width = 0.25, height = 0)) +
    facet_wrap(~ CellType, scales = "free") +
    scale_colour_manual(values = c("white", "gray20")) +
    scale_shape_manual(values = c(21, 24)) + theme_minimal()

## Plot to compare method effectiveness
## Create dataframe for beta cell proportions and HbA1c levels
m_prop_ana <- data.frame(pData(bulk_eset)[rep(1:89, 2), phenotype_factors],
                        ct.prop = c(est_prop$Est.prop.weighted[, 2],
                                    est_prop$Est.prop.allgene[, 2]),
                        Method = factor(rep(methods, each = 89),
                                        levels = methods))
colnames(m_prop_ana)[1:4] <- phenotype_factors
m_prop_ana <- subset(m_prop_ana, !is.na(m_prop_ana[phenotype_gene]))
m_prop_ana$Disease <- factor(sample_groups[(  # nolint
    m_prop_ana[phenotype_gene] > 6.5) + 1], sample_groups)
m_prop_ana$D <- (m_prop_ana$Disease ==        # nolint
                 sample_disease_group) / sample_disease_group_scale

jitt_compare <- ggplot(m_prop_ana, aes_string(phenotype_gene, "ct.prop")) +
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
    ##subset(m_prop_ana, Method == meth))
    lm_beta_meth <- lm(as.formula(
        paste("ct.prop", paste(phenotype_factors, collapse = " + "),
              sep = " ~ ")),
        data = subset(m_prop_ana, Method == meth))
    print(paste0("Summary: ", meth))
    capture.output(summary(lm_beta_meth),
                   file = paste0("report_data/summ_", meth, ".txt"))
}
