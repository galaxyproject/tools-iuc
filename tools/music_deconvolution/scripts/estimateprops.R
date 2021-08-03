suppressPackageStartupMessages(library(xbioc))
suppressPackageStartupMessages(library(MuSiC))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
## We use this script to estimate the effectiveness of proportion methods

## Load Conf
args = commandArgs(trailingOnly=TRUE)
source(args[1])

##bulk_eset = readRDS('test-data/GSE50244bulkeset.rds')
##scrna_eset = readRDS('test-data/EMTABesethealthy.rds')

## clusters_label = 'cellType'
## samples_label = 'sampleID'
## phenotype_factors = c('age', 'bmi', 'hba1c', 'gender')
## celltypes = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal')
## methods = c('MuSiC', 'NNLS')
## phenotype_gene = 'hba1c' ## *must* be in list of phenotype_factors
## sample_groups = c('Normal', 'T2D')
## sample_disease_group = "T2D"
## sample_disease_group_scale = 5
## healthy_phenotype = "Normal"
## compare_title=""

print(bulk_eset)
print(scrna_eset)

## Estimate cell type proportions
est.prop = music_prop(bulk.eset = bulk_eset, sc.eset = scrna_eset,
                      clusters = clusters_label,
                      samples = samples_label, select.ct = celltypes, verbose = T)


## Show different in estimation methods
## Jitter plot of estimated cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(est.prop$Est.prop.weighted),
                             data.matrix(est.prop$Est.prop.allgene)),
                        method.name = methods, title = 'Jitter plot of Est Proportions')


## Make a Plot
## A more sophisticated jitter plot is provided as below. We seperated the T2D subjects and normal
#subjects by their HbA1c levels.
m.prop = rbind(melt(est.prop$Est.prop.weighted),
                        melt(est.prop$Est.prop.allgene))

colnames(m.prop) = c('Sub', 'CellType', 'Prop')

m.prop$CellType = factor(m.prop$CellType, levels = celltypes)

m.prop$Method = factor(rep(methods, each = 89*6), levels = methods)
m.prop$HbA1c = rep(bulk_eset$hba1c, 2*6)
m.prop = m.prop[!is.na(m.prop$HbA1c), ]
m.prop$Disease = factor(sample_groups[(m.prop$HbA1c > 6.5)+1], levels = sample_groups)

m.prop$D = (m.prop$Disease == sample_disease_group)/sample_disease_group_scale
m.prop = rbind(subset(m.prop, Disease == healthy_phenotype),
               subset(m.prop, Disease != healthy_phenotype))

jitter.new = ggplot(m.prop, aes(Method, Prop)) +
  geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease),
             size = 2, alpha = 0.7, position = position_jitter(width=0.25, height=0)) +
  facet_wrap(~ CellType, scales = 'free') + scale_colour_manual( values = c('white', "gray20")) +
  scale_shape_manual(values = c(21, 24))+ theme_minimal()

## Plot to compare method effectiveness
## Create dataframe for beta cell proportions and HbA1c levels
m.prop.ana = data.frame(pData(bulk_eset)[rep(1:89, 2), phenotype_factors],
                        ct.prop = c(est.prop$Est.prop.weighted[, 2],
                                    est.prop$Est.prop.allgene[, 2]),
                        Method = factor(rep(methods, each = 89),
                                        levels = methods))
colnames(m.prop.ana)[1:4] = phenotype_factors
m.prop.ana = subset(m.prop.ana, !is.na(m.prop.ana[phenotype_gene]))
m.prop.ana$Disease = factor(sample_groups[(m.prop.ana[phenotype_gene] > 6.5) + 1], sample_groups)
m.prop.ana$D = (m.prop.ana$Disease == sample_disease_group)/sample_disease_group_scale

jitt.compare <- ggplot(m.prop.ana, aes_string(phenotype_gene, "ct.prop")) +
    geom_smooth(method = 'lm',  se = FALSE, col = 'black', lwd = 0.25) +
    geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease), size = 2, alpha = 0.7) +  facet_wrap(~ Method) +
    ggtitle(compare_title) + theme_minimal() +
    scale_colour_manual( values = c('white', "gray20")) +
    scale_shape_manual(values = c(21, 24))


pdf(file=outfile_pdf, width=8, height=8)
plot_grid(jitter.fig, jitter.new, labels = 'auto', ncol=1, nrow=2)
jitt.compare
dev.off()

## Summary table
for (meth in methods){
    ##lm.beta.meth = lm(ct.prop ~ age + bmi + hba1c + gender, data = subset(m.prop.ana, Method == meth))
    lm.beta.meth = lm(as.formula(paste("ct.prop", paste(phenotype_factors, collapse=" + "), sep=" ~ ")),
                      data = subset(m.prop.ana, Method == meth))
    print(paste0("Summary: ", meth))
    capture.output(summary(lm.beta.meth), file = paste0(meth, ".txt"))
}