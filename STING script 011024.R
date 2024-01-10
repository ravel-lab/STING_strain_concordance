### STING urogenital strain concordance analysis
### Kayla A. Carter
### 1 / 10 / 24






#############
### setup ###
#############

library(DataCombine)
library(data.table)
library(dbplyr)
library(dtplyr)
library(epiDisplay)
library(gmodels)
library(plyr)
library(stringi)
library(reshape2)
library(tidyselect)
library(tidyverse)
library(ggnetwork)
library(ggnewscale)
library(network)
library(sna)
library(ggraph)
library(plotly)
library(tidygraph)
library(colorRamps)
library(viridis)
library(grDevices)
library(ggpattern)
library(patchwork)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(GUniFrac)

# set working directory or add directory path to file names
data <- readRDS("STING_inStrain_compare_genomeWide_compare.rds")
# data should have 2730918 obs of 10 vars
mag <- readRDS('STING MAG inventory.rds')
# mag should have 507 obs of 7 vars
part <- readRDS("STING merge partnership level.rds")
# part should have 42 obs of 325 vars
person <- readRDS("STING merge participant level.rds")
# person should have 138 obs of 217 vars
sites <- readRDS("STING comparable sites.rds")
# sites should have 3613240 obs of 5 vars

# create function that is opposite of %in%
'%!in%' <- function(x, y)!('%in%'(x, y))






#########################################################################
### make trimmed, cleaned inStrain dataset of just concordant strains ###
#########################################################################

# add taxonomy variable to data by merging with mag
data <- left_join(data, mag[, colnames(mag) %in% c("bin_id", "taxonomy")],
                  by = c("genome" = "bin_id"))
data <- rename(data, bug = taxonomy)
# data should have 2730918 obs of 11 vars

# make indicator of percent of genomes compared >=50%
data$per_comp_50 <- ifelse(data$percent_compared >= 0.5, 1, 0)
# data should have 2730918 obs of 12 vars

# make indicator of population ANI >=0.9999
data$popANI_9999 <- ifelse(data$popANI >= 0.9999, 1, 0)
# data should have 2730918 obs of 13 vars

# make indicator of concordant strains with percent of genomes compared >=50%
# AND population ANI >=0.9999, restrict to these comparisons
data$share <- ifelse(is.na(data$per_comp_50) == T, NA,
                     ifelse(is.na(data$popANI_9999) ==  T, NA,
                            ifelse((data$per_comp_50 == 1 & data$popANI_9999 == 1), 1, 0)))
# data should have 2730918 obs of 14 vars

# make variable indicating whether comparison includes a control
data$control <- 0
data$control[grep("Negative", data$name1)] <- 1
data$control[grep("Negative", data$name2)] <- 1
data$control[grep("Zymo", data$name1)] <- 1
# data should have 2730918 obs of 15 vars

# save data including control samples
all_comp <- data
# all_comp should have 2730918 obs of 15 vars

# save version of data that is only controls
controls <- data[data$control == 1,]
# controls should have 78855 obs of 15 vars

# save version of data excluding comparisons with controls
data <- data[data$control == 0,]
# data should have 2652063 obs of 15 vars

# restrict data to concordance events
data_trim <- data[is.na(data$share) == F & data$share == 1,]
# data_trim should have 347 obs of 15 vars

# make pop SNP/Mbp variable
data_trim$popSNP_mbp <- (data_trim$population_SNPs / data_trim$compared_bases_count) * 10^6
# data_trim should have 347 obs of 16 vars






#####################################
### merge metadata into data_trim ###
#####################################

# make concatenated barcode variable in the partnership dataset to merge with 
# data_trim, concatenate in both orders
part$part_barcode_12 <- paste0(part$barcode_1, part$barcode_2)
part$part_barcode_21 <- paste0(part$barcode_2, part$barcode_1)
# part should have 42 obs of 327 vars

# make concatenated barcode variable in data_trim to merge with
data_trim$part_barcode <- paste0(data_trim$name1, data_trim$name2)
# data_trim should have 347 obs of 17 vars

# add partnership variable to data_trim, make factor version
data_trim$partner <- 0
data_trim$partner[data_trim$part_barcode %in% c(part$part_barcode_12, part$part_barcode_21)] <- 1
data_trim$partner_f <- as.factor(data_trim$partner)
# data_trim should have 347 obs of 19 vars

# merge in sex data, drop observations without sex data for both people
data_trim <- left_join(data_trim, person[, colnames(person) %in% c('barcode', 'sex')],
                       by = c('name1' = 'barcode'))
data_trim <- rename(data_trim, sex_1 = sex)
data_trim <- left_join(data_trim, person[, colnames(person) %in% c('barcode', 'sex')],
                       by = c('name2' = 'barcode'))
data_trim <- rename(data_trim, sex_2 = sex)
data_trim <- data_trim[!is.na(data_trim$sex_1),]
data_trim <- data_trim[!is.na(data_trim$sex_2),]
# data_trim should have 334 obs of 21 vars

# merge in other metadata of interest, need to do separately for each contact, 
# will need to rename new cols to indicate which contact
data_trim <- left_join(data_trim, 
                       person[, colnames(person) %in% c('barcode', 'study_day', 'site', 'age', 'n_part', 'cst', 'sub_cst')], 
                       by = c('name1' = 'barcode'))
data_trim <- rename(data_trim, c(study_day_1 = study_day,
                                 site_1 = site,
                                 age_1 = age,
                                 n_part_1 = n_part,
                                 cst_1 = cst,
                                 sub_cst_1 = sub_cst))
data_trim <- left_join(data_trim, 
                       person[, colnames(person) %in% c('barcode', 'study_day', 'site', 'age', 'n_part', 'cst', 'sub_cst')], 
                       by = c('name2' = 'barcode'))
data_trim <- rename(data_trim, c(study_day_2 = study_day,
                                 site_2 = site,
                                 age_2 = age,
                                 n_part_2 = n_part,
                                 cst_2 = cst,
                                 sub_cst_2 = sub_cst))
# data_trim should have 334 obs of 33 vars

# make variables for time difference between enrollment dates, which are coded as
# study_day, which equals 1 for the date of the first enrollment and increases
# in increments of 1 over consecutive days until the final enrollment. so difference
# in days is just simple subtraction, don't need to use date functions
data_trim$day_diff <- abs(data_trim$study_day_1 - data_trim$study_day_2)
data_trim$week_diff <- data_trim$day_diff / 7
# data_trim should have 334 obs of 35 vars

# add indicator of same CST
data_trim$same_cst <- ifelse(data_trim$cst_1 == data_trim$cst_2, 1, 0)
# data_trim should have 334 obs of 36 vars

# restrict to unique observations based on samples and concordant taxon because 
# multiple concordance events of same taxon between the same people likely means 
# the concordant strain was called concordant according to multiple MAGs for the 
# taxon. will assume these redundant events represent a single concordance event
# and retain the event with greatest number of bases compared.
# make var that is concatenated name1, name2, bug to group by
data_trim$part_bug <- paste0(data_trim$part_barcode, data_trim$bug)
# data_trim should have 334 obs of 37 vars
# order data_trim by part_bug and compared_bases_count (descending)
data_trim <- data_trim[order(data_trim$part_bug, -data_trim$compared_bases_count),]
# restrict to unique concordance events within pairs
share_u <- do.call(rbind, by(data_trim, list(data_trim$part_bug), 
                             FUN=function(x) head(x, 1)))
# share_u should have 119 obs of 37 vars

# assuming all concordance between contacts with <5 SNP/Mbp is sexually transmitted
# make indicator of transmission
share_u$trans <- ifelse(share_u$partner == 1 & share_u$popSNP_mbp < 5, 1, 0)
# share_u should have 119 obs of 38 vars

# restrict to whatever was dropped from share when making share_u
share_drop <- anti_join(data_trim, share_u)
# share_drop should have 215 obs of 37 vars






#########################################
### identify and remove contamination ###
#########################################

# restrict controls dataframe to concordance events using same criteria as for 
# non-control samples
control_trim <- controls[is.na(controls$share) == F & controls$share == 1,]
# control_trim should have 19 obs of 15 vars

# confirm all in control_trim are with negative control, not positive
tab1(control_trim$name2)

# make variable in control_trim that is concatenated genome and name1
control_trim$genome_name <- paste0(control_trim$genome, control_trim$name1)
# control_trim should have 19 obs of 16 vars

# make variables in share_u that are concatenated genome and name1, genome and name2
share_u$genome_name1 <- paste0(share_u$genome, share_u$name1)
share_u$genome_name2 <- paste0(share_u$genome, share_u$name2)
# share_u should have 119 obs of 40 vars

# make indicator of contamination in share_u. contamination defined as, for a 
# given concordance event between mon-control samples, both non-control samples 
# are also concordant for the same strain with a control sample
share_u$name1_control <- ifelse(share_u$genome_name1 %in% control_trim$genome_name,
                                1, 0)
share_u$name2_control <- ifelse(share_u$genome_name2 %in% control_trim$genome_name,
                                1, 0)
share_u$contam <- ifelse(share_u$name1_control == 1 & share_u$name2_control == 1,
                         1, 0)
# share_u should have 119 obs of 43 vars

# save dataset of just contamination events, restrict share_u to non-contamination
contam <- share_u[share_u$contam == 1,]
# contam should have 4 obs of 43 vars
share_u <- share_u[share_u$contam == 0,]
# share_u should have 115 obs of 43 vars






#########################################################
### Figure 3 - differential abundance analysis by sex ###
#########################################################

### prepare relative abundance matrix, metadata dataframe for ZicoSeq

# generate relative abundance matrix/array with samples as columns, rename BVAB1
ra <- t(person[, c(32:165)])
colnames(ra) <- person$barcode
rownames(ra)[rownames(ra) == 'BVAB1'] <- 'Ca_Lachnocurva_vaginae'
# ra should be num [1:134, 1:138]

# drop chlamydia trachomatis from relative abundances because artificially
# inflated in male samples due to lower total bacterial load
ra <- ra[row.names(ra) != 'Chlamydia_trachomatis',]
# ra should be num [1:133, 1:138]

# make metadata dataframe, needs to have samples as rows and row names same as 
# col names in relative abundance matrix
meta <- person[, colnames(person) %in% c('barcode', 'sex')]
rownames(meta) <- meta$barcode
# meta should have 138 obs of 2 vars

# zicoseq differential abundance analysis
zico <- ZicoSeq(meta, ra, grp.name = 'sex', feature.dat.type = 'proportion',
                prev.filter = 0.02, # filtering things detected in <3 samples
                mean.abund.filter = 0.00001, 
                max.abund.filter = 0.0001, 
                is.winsor = T, outlier.pct = 0.03, winsor.end = 'top',
                link.func = list(function (x) x^0.5), stats.combine.func = max,
                perm.no = 999, strata = NULL,
                ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                is.fwer = T, verbose = T, return.feature.dat = T)

# ZicoSeq.plot does not pass text.size argument to geom_text_repel labels and 
# the labels are way too big, so going to edit the function to make smaller
# labels and repel parameters
# trace(ZicoSeq.plot, edit=TRUE)
# edit: ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf, color = 'black')
# to: ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf, color = 'black', size = 2, force = 2.5, force_pull = 0.5)

# ZicoSeq.plot labels all taxa with p< specified cutoff, but want to label all
# taxa that are presented in stacked bar plots, network diagram, so going to edit
# function to add these
# trace(ZicoSeq.plot, edit=TRUE)
# edit:
# if (length(unique(meta.dat[, grp.name])) > 2 & !is.numeric(meta.dat[,grp.name])) {
#   ...
# }
# else {
#   signs <- t(sign(coefs))[t(ZicoSeq.obj$R2 == R2)]
#   signs[signs == 0] <- 1
#   plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
#                           prevalence = prevalence, 
#                           abundance = abundance, 
#                           R2 = R2 * signs, 
#                           taxa = rownames(ZicoSeq.obj$R2))
#   plot.data[plot.data$pvals > cutoff, 'taxa'] <- ''
# }
# to: 
# tax_interest <- c('Ca_Lachnocurva_vaginae', 'g_Corynebacterium', 'Atopobium_vaginae', 
#                   'Gardnerella_vaginalis', 'Lactobacillus_crispatus', 
#                   'Lactobacillus_iners', 'Lactobacillus_jensenii', 
#                   'g_Megasphaera', 'g_Neisseria', 'Prevotella_amnii', 
#                   'Prevotella_bivia', 'g_Propionimicrobium', 'Sneathia_amnii', 
#                   'g_Staphylococcus', 'g_Streptococcus', 'g_Ureaplasma', 
#                   'g_Mobiluncus', 'Prevotella_timonensis', 'Sneathia_sanguinegens')
# '%!in%' <- function(x, y)!('%in%'(x, y))
# if (length(unique(meta.dat[, grp.name])) > 2 & !is.numeric(meta.dat[,grp.name])) {
#   ...
# }
# else {
#   signs <- t(sign(coefs))[t(ZicoSeq.obj$R2 == R2)]
#   signs[signs == 0] <- 1
#   plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
#                           prevalence = prevalence, 
#                           abundance = abundance, 
#                           R2 = R2 * signs, 
#                           taxa = rownames(ZicoSeq.obj$R2))
#   plot.data[(plot.data$pvals > cutoff & plot.data$taxa %!in% tax_interest), 'taxa'] <- ''
# }

plot_zico <- ZicoSeq.plot(zico, meta, pvalue.type = 'p.adj.fwer', cutoff = 0.01) +
  ggtitle('') +
  ylab(expression(paste(log[10], '(FWER-adjusted p)'))) +
  annotate(geom = 'text', x = -0.043, y = 3.75, 
           label = 'Enriched in\nvaginal', size = 3, fontface = 'bold') +
  annotate(geom = 'text', x = 0.1, y = 3.75, 
           label = 'Enriched in penile', size = 3, fontface = 'bold') +
  theme_bw(base_size = 10)
# plot_zico






###################################################################################
### Figure S1 - Differential abundance by CT/NG status among penile metagenomes ###
###################################################################################

### prepare relative abundance matrix, metadata dataframe for ZicoSeq

# generate relative abundance matrix/array with samples as columns, rename BVAB1
ra_m <- t(person[person$sex == 'm', c(32:165)])
colnames(ra_m) <- person$barcode[person$sex == 'm']
rownames(ra_m)[rownames(ra_m) == 'BVAB1'] <- 'Ca_Lachnocurva_vaginae'
# ra_m should be num [1:134, 1:64]

# drop chlamydia trachomatis, neisseria genus because they're basically the exposure
ra_m <- ra_m[row.names(ra_m) %!in% c('Chlamydia_trachomatis', 'g_Neisseria'),]
# ra_m should be num [1:132, 1:64]

# make metadata dataframe, needs to have samples as rows and row names same as 
# col names in relative abundance matrix
CrossTable(person$ct[person$sex == 'm'], person$gc[person$sex == 'm'])
# all men with gonorrhea also have chlamydia so can use the ct variable as exposure
meta_m <- person[person$sex == 'm', colnames(person) %in% c('barcode', 'ct')]
rownames(meta_m) <- meta_m$barcode
# meta_m should have 64 obs of 2 vars

# restrict meta_m and ra_m to participants with non-missing ct data
meta_m <- meta_m[!is.na(meta_m$ct),]
# meta_m should have 62 obs of 2 vars
ra_m <- ra_m[, colnames(ra_m) %in% meta_m$barcode]
# ra_m should be num [1:132, 1:62] 

# drop taxa that are all 0's from ra_m
ra_m <- ra_m[!rowSums(ra_m) == 0,]
# ra_m should be num [1:122, 1:62] 

# zicoseq differential abundance analysis
zico_m <- ZicoSeq(meta_m, ra_m, grp.name = 'ct', feature.dat.type = 'proportion',
                  prev.filter = 0.05, # filtering things detected in <3 samples
                  mean.abund.filter = 0.00001, 
                  max.abund.filter = 0.0001, 
                  is.winsor = T, outlier.pct = 0.03, winsor.end = 'top',
                  link.func = list(function (x) x^0.5), stats.combine.func = max,
                  perm.no = 999, strata = NULL,
                  ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                  is.fwer = T, verbose = T, return.feature.dat = T)

# ZicoSeq.plot does not pass text.size argument to geom_text_repel labels and 
# the labels are too big, so going to edit the function to make smaller labels
# trace(ZicoSeq.plot, edit=TRUE)
# edit: ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf, color = 'black')
# to: ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf, color = 'black', size = 3.25, , fontface = 'italic')

# plotting FWER-adjusted p values with cutoff at fwer-adj p < 0.1
plot_zico_m <- ZicoSeq.plot(zico_m, pvalue.type = 'p.adj.fwer', cutoff = 0.1) +
  ggtitle('') +
  ylab(expression(paste(log[10], '(FWER-adjusted p)'))) +
  annotate(geom = 'text', x = -0.1, y = 2.75, 
           label = 'Enriched in\nSTI-', size = 3.5, fontface = 'bold') +
  annotate(geom = 'text', x = 0.045, y = 2.75, 
           label = 'Enriched in\nSTI+', size = 3.5, fontface = 'bold') +
  theme_bw(base_size = 10)
# plot_zico_m






#####################################################
### Figure 4 - strain concordance network diagram ###
#####################################################

# rename KA00274_sp902373515, UBA629_sp005465875, Bifidobacterium
share_u$bug[share_u$bug == 'KA00274_sp902373515'] <- 'BVAB2'
share_u$bug[share_u$bug == 'UBA629_sp005465875'] <- 'Ca_Lachnocurva_vaginae'
share_u$bug[share_u$bug == 'Bifidobacterium'] <- 'Gardnerella_sp'

# make network graph object to plot
mat <- as.matrix(share_u[, colnames(share_u) %in% c('name1', 'name2')])
names(mat) <- c('to', 'from')
net <- network(mat, directed = F, hyper = F, loops = T, multiple = T, matrix.type = 'edgelist')

# add taxon, partnership data as edge attributes
net <- as_tbl_graph(net) %>%
  activate(edges) %>%
  mutate(bug = share_u$bug) %>%
  activate(edges) %>%
  mutate(partner_f = share_u$partner_f) 

# add sex data as node attribute
# make dataframe of barcodes in same order as in net, add sex, then add sex to net
sex <- as.data.frame(names(as.list(net[[1:54]])))
colnames(sex) <- 'barcode'
sex <- left_join(sex, unique(person[, colnames(person) %in% c('barcode', 'sex')]))
# sex should have 54 obs of 2 vars
net <- net %>%
  activate(nodes) %>%
  mutate(sex = sex$sex)

# define taxon colors for network diagram
net_colors <- c('#ff0000', 
                '#73FCD6', 
                '#ff8c00', 
                '#333333', 
                '#c7aa8f', 
                '#0000cd', 
                '#3388ff', 
                '#ff738a', 
                '#b0b0b0', 
                '#ffc0cb', 
                '#800080', 
                '#f8fca9', 
                '#c1ffc1', 
                '#BFCE6D', 
                '#ac8fc7', 
                '#6963ff', 
                '#225151', 
                '#20b2aa', 
                '#981752', 
                '#B5F46E', 
                '#97C0DD', 
                '#AFCCB7', 
                '#1C792B', 
                '#7A4300', 
                '#FDAA8B') 

names(net_colors) <- c('Lactobacillus_crispatus', 
                       'Gardnerella_vaginalis', 
                       'Lactobacillus_iners', 
                       'Lactobacillus_jensenii',  
                       'Sneathia_sanguinegens', 
                       'Fannyhessea_vaginae', 
                       'Fannyhessea', 
                       'Streptococcus_agalactiae', 
                       'Prevotella_sp000758925', 
                       'Streptococcus_anginosus', 
                       'Staphylococcus_epidermidis', 
                       'Prevotella_timonensis', 
                       'Gardnerella_sp', 
                       'Prevotella',  
                       'Prevotella_amnii', 
                       'Sneathia_vaginalis',  
                       'Gardnerella_leopoldii',
                       'Gardnerella_swidsinkii', 
                       'Lactobacillus_mulieris', 
                       'Limosilactobacillus_portuensis', 
                       'Streptococcus_mitis', 
                       'Gardnerella_piotii', 
                       'Gardnerella_vaginalis_A', 
                       'Mobiluncus_curtisii', 
                       'Megasphaera_lornae') 

net_colors <- net_colors[order(factor(names(net_colors)))]

# plot
plot_net <- ggraph(net, layout = 'nicely') + 
  geom_edge_fan(aes(edge_color = bug, edge_width = partner_f, alpha = partner_f),
                strength = 4.5)  +
  scale_edge_width_manual(values = c(1, 2.9),
                          labels = c('Non-contacts', 'Sexual contacts'),
                          'Concordance between:') +
  scale_edge_alpha_manual(values = c(0.75, 1),
                          labels = c('Non-contacts', 'Sexual contacts'),
                          'Concordance between:') +
  scale_edge_color_manual(values = net_colors,
                          labels = c(expression(paste(italic('Fannyhessea'), ' spp.')),
                                     expression(paste(italic('Fannyhessea vaginae'), '')),
                                     expression(paste(italic('Gardnerella leopoldii'), '')),
                                     expression(paste(italic('Gardnerella piotii'), '')),
                                     expression(paste(italic('Gardnerella'), ' spp.')),
                                     expression(paste(italic('Gardnerella swidsinskii'), '')),
                                     expression(paste(italic('Gardnerella vaginalis'), '')),
                                     expression(paste(italic('Gardnerella vaginalis'), ' A')),
                                     expression(paste(italic('Lactobacillus crispatus'), '')),
                                     expression(paste(italic('Lactobacillus iners'), '')),
                                     expression(paste(italic('Lactobacillus jensenii'), '')),
                                     expression(paste(italic('Lactobacillus mulieris'), '')),
                                     expression(paste(italic('Limosilactobacillus portuensis'), '       ')),
                                     expression(paste(italic('Megasphaera lornae'), '')),
                                     expression(paste(italic('Mobiluncus curtisii'), '')),
                                     expression(paste(italic('Prevotella'), ' spp.')),
                                     expression(paste(italic('Prevotella amnii'), '')),
                                     expression(paste(italic('Prevotella'), ' sp000758925')),
                                     expression(paste(italic('Prevotella timonensis'), '')),
                                     expression(paste(italic('Sneathia sanguinegens'), '')),
                                     expression(paste(italic('Sneathia vaginalis'), '')),
                                     expression(paste(italic('Staphylococcus epidermidis   '), '')),
                                     expression(paste(italic('Streptococcus agalactiae'), '')),
                                     expression(paste(italic('Streptococcus anginosus'), '')),
                                     expression(paste(italic('Streptococcus mitis'), ''))),
                          'Taxon') +
  geom_node_point(aes(shape = sex, size = sex), fill = 'black') +
  scale_shape_manual(values = c(21, 23),
                     labels = c('F', 'M'),
                     'Sex') +
  scale_size_manual(values = c(4, 4),
                    labels = c('F', 'M'),
                    'Sex') +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(size = c(3, 3))),
         size = guide_legend(order = 1,
                             override.aes = list(size = c(3, 3))),
         edge_alpha = guide_legend(order = 2,
                                   reverse = T),
         edge_width = guide_legend(order = 2,
                                   reverse = T),
         edge_color = guide_legend(ncol = 1,
                                   order = 3,
                                   override.aes = list(edge_width = 2))) +
  theme_blank(base_size = 10) +
  theme(legend.text.align = 0,
        legend.key.size = unit(1, 'char'), 
        legend.key.height = unit(1, 'char'),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.margin = margin(0, 0, 0, 0))
# plot_net



# summarize network characteristics
net <- network(mat, directed = F, hyper = F, loops = T, multiple = T, 
               matrix.type = 'edgelist')

network.edgecount(net)
# 115 edges

network.size(net)
# 54 nodes






##########################################
### Figure S5 - sexual network diagram ###
##########################################

# make network graph object to plot
mat <- as.matrix(part[, colnames(part) %in% c('sting_id_1', 'sting_id_2')])
no_contacts <- matrix(nrow =length(person$sting_id[person$sting_id %!in% unique(c(mat[,1], mat[,2]))]),
                      ncol = 2)
no_contacts[, 1] <- person$sting_id[person$sting_id %!in% unique(c(mat[,1], mat[,2]))]
no_contacts[, 2] <- 'drop'
mat <- rbind(mat, no_contacts)
names(mat) <- c('to', 'from')
sex_net <- network(mat, directed = F, hyper = F, loops = T, multiple = F, 
                   matrix.type = 'edgelist') 
sex_net <- delete.vertices(sex_net, vid = c(1))

# add sex data as node attribute
# make dataframe of barcodes in same order as in net, add sex, then add sex to sex_net
sex_net <- as_tbl_graph(sex_net)
attrib <- as.data.frame(names(as.list(sex_net[[1:138]])))
colnames(attrib) <- 'sting_id'
attrib <- left_join(attrib, unique(person[, colnames(person) %in% c('sting_id', 'sex', 'n_part')]))
attrib$n_part <- factor(ifelse(attrib$n_part == 0, 0, 1), levels = c(1, 0))
attrib$sex_part <- factor(ifelse(attrib$sex == 'f' & attrib$n_part == 1, 1,
                                 ifelse(attrib$sex == 'm' & attrib$n_part == 1, 2,
                                        ifelse(attrib$sex == 'f', 3, 4))))
# attrib should have 138 obs of 4 vars
sex_net <- sex_net %>%
  activate(nodes) %>%
  mutate(sex_part = attrib$sex_part)

# network diagram
plot_sex_net <- ggraph(sex_net, layout = 'nicely') + 
  geom_edge_fan(width = 1)  +
  geom_node_point(aes(shape = sex_part, fill = sex_part), stroke = 1, size = 4) +
  scale_shape_manual(values = c(21, 23, 21, 23),
                     labels = c('F in contact dyad', 
                                'M in contact dyad',
                                'F no contact',
                                'M no contact'),
                     'Sex and dyad status     ') +
  scale_fill_manual(values = c('black', 'black', 'white', 'white'),
                    labels = c('F in contact dyad', 
                               'M in contact dyad',
                               'F no contact',
                               'M no contact'),
                    'Sex and dyad status     ') +
  guides(shape = guide_legend(override.aes = list(size = c(3, 3))),
         fill = guide_legend(override.aes = list(size = c(3, 3)))) +
  theme_blank(base_size = 10) +
  theme(legend.text.align = 0,
        legend.key.size = unit(1, 'char'), 
        legend.key.height = unit(1, 'char'),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.margin = margin(0, 0, 0, 0))
# plot_sex_net






################################################
### general descriptives of study population ###
################################################

tab1(person$race)
tab1(person$race[person$sex == 'f'])
tab1(person$race[person$sex == 'm'])

summary(person$age)
summary(person$age[person$sex == 'f'])
summary(person$age[person$sex == 'm'])

summary(person$part_m[!is.na(person$part_m)])
sum(is.na(person$part_m))
summary(person$part_m[!is.na(person$part_m) & person$sex == 'f'])
sum(is.na(person$part_m[person$sex == 'f']))
summary(person$part_m[!is.na(person$part_m) & person$sex == 'm'])
sum(is.na(person$part_m[person$sex == 'm']))

summary(person$part_f[!is.na(person$part_f)])
sum(is.na(person$part_f))
summary(person$part_f[!is.na(person$part_f) & person$sex == 'f'])
sum(is.na(person$part_f[person$sex == 'f']))
summary(person$part_f[!is.na(person$part_f) & person$sex == 'm'])
sum(is.na(person$part_f[person$sex == 'm']))

tab1(person$condom[person$sex == 'f'])
tab1(person$hormonal[person$sex == 'f'])
tab1(person$iud_copper[person$sex == 'f'])

person <- person %>% mutate(other = case_when(foam == 1 ~ 1,
                                              withdrawal == 1 ~ 1,
                                              abstinence == 1 ~ 1,
                                              emergency == 1 ~ 1))
person$other[is.na(person$other) & person$sex == 'f'] <- 0
# person should have 138 obs of 218 vars
tab1(person$other[person$sex == 'f'])

sum(person$wave1 == 1 & person$n_part == 0)
sum(person$wave1 == 1 & person$n_part != 0)
sum(person$wave2 == 1)
sum(person$wave3 == 1)
sum(person$wave4 == 1)

sum(person$wave1[person$sex == 'f'] == 1 & person$n_part[person$sex == 'f'] == 0)
sum(person$wave1[person$sex == 'f'] == 1 & person$n_part[person$sex == 'f'] != 0)
sum(person$wave2[person$sex == 'f'] == 1)
sum(person$wave3[person$sex == 'f'] == 1)
sum(person$wave4[person$sex == 'f'] == 1)

sum(person$wave1[person$sex == 'm'] == 1 & person$n_part[person$sex == 'm'] == 0)
sum(person$wave1[person$sex == 'm'] == 1 & person$n_part[person$sex == 'm'] != 0)
sum(person$wave2[person$sex == 'm'] == 1)
sum(person$wave3[person$sex == 'm'] == 1)
sum(person$wave4[person$sex == 'm'] == 1)

tab1(person$site)
tab1(person$site[person$sex == 'f'])
tab1(person$site[person$sex == 'm'])

tab1(person$cst)
tab1(person$cst[person$sex == 'f'])
tab1(person$cst[person$sex == 'm'])

tab1(person$ct)
tab1(person$ct[person$sex == 'f'])
tab1(person$ct[person$sex == 'm'])

tab1(person$gc)
tab1(person$gc[person$sex == 'f'])
tab1(person$gc[person$sex == 'm'])

tab1(person$pid)

tab1(person$dis_eval_f[person$sex == 'f'])
tab1(person$cervix_dis[person$sex == 'f'])
tab1(person$bv_any[person$sex == 'f'])
tab1(person$mpc[person$sex == 'f'])
tab1(person$yeast[person$sex == 'f'])

tab1(person$dis_eval_m[person$sex == 'm'])
tab1(person$urethritis[person$sex == 'm'])
tab1(person$circumcised[person$sex == 'm'])

tab1(person$symp)
tab1(person$symp[person$sex == 'f'])
tab1(person$symp[person$sex == 'm'])






###########################################
### data analysis of concordant strains ###
###########################################

### descriptives

# summarize n events per taxon
tab1(share_u$bug)
tab1(share_u$bug[share_u$partner == 1])
tab1(share_u$bug[share_u$partner == 0])



### summarize popSNPs per Mbp compared by bug, partnership status

SNP_mbp_bug <- data.frame(t(sapply(unique(share_u$bug),
                                   function(x) summary(share_u$popSNP_mbp[share_u$bug == x]))))
kruskal.test(share_u$popSNP_mbp ~ share_u$bug)

SNP_mbp_bug_part <- data.frame(t(sapply(unique(share_u$bug),
                                        function(x) summary(share_u$popSNP_mbp[share_u$bug == x & share_u$partner == 1]))))
SNP_mbp_bug_no_part <- data.frame(t(sapply(unique(share_u$bug),
                                           function(x) summary(share_u$popSNP_mbp[share_u$bug == x & share_u$partner == 0]))))

SNP_mbp_part <- as.data.frame(rbind(summary(share_u$popSNP_mbp[share_u$partner_f == 1]),
                                    summary(share_u$popSNP_mbp[share_u$partner_f == 0])))
SNP_mbp_part$contact <- c('yes', 'no')
kruskal.test(share_u$popSNP_mbp ~ share_u$partner_f)



### summarize concordance of metadata for concordance events

tab1(share_u$same_cst, graph = F)
CrossTable(share_u$bug, share_u$same_cst)

summary(share_u$week_diff)
week_bug <- data.frame(t(sapply(unique(share_u$bug),
                                function(x) summary(share_u$week_diff[share_u$bug == x]))))

# code variable that indicates sex of both participants with concordance
share_u$sex_sex <- paste0(share_u$sex_1, share_u$sex_2)
share_u$sex_sex[share_u$sex_sex == 'mf'] <- 'fm'
# share_u should have 115 obs of 44 vars
CrossTable(share_u$bug[share_u$partner == 0], share_u$sex_sex[share_u$partner == 0])



### examine relative abundances of concordant strains

# merge both participants' relevant relative abundance data into share_u
share_u <- left_join(share_u, person[, colnames(person) %in% c('barcode', 'Atopobium_vaginae', 'Atopobium_parvulum', 'Gardnerella_vaginalis', 'Lactobacillus_crispatus', 'Lactobacillus_iners', 'Lactobacillus_jensenii', 'Lactobacillus_vaginalis', 'g_Megasphaera', 'g_Mobiluncus', 'Prevotella_amnii', 'Prevotella_timonensis', 'Sneathia_sanguinegens', 'Sneathia_amnii', 'g_Staphylococcus', 'g_Streptococcus')],
                     by = c('name1' = 'barcode'))
colnames(share_u)[45:ncol(share_u)] <- paste0(colnames(share_u)[45:ncol(share_u)], '_1')
share_u$g_Atopobium_1 <- share_u$Atopobium_parvulum_1 + share_u$Atopobium_vaginae_1
share_u <- left_join(share_u, person[, colnames(person) %in% c('barcode', 'Atopobium_vaginae', 'Atopobium_parvulum', 'Gardnerella_vaginalis', 'Lactobacillus_crispatus', 'Lactobacillus_iners', 'Lactobacillus_jensenii', 'Lactobacillus_vaginalis', 'g_Megasphaera', 'g_Mobiluncus', 'Prevotella_amnii', 'Prevotella_timonensis', 'Sneathia_sanguinegens', 'Sneathia_amnii', 'g_Staphylococcus', 'g_Streptococcus')],
                     by = c('name2' = 'barcode'))
colnames(share_u)[61:ncol(share_u)] <- paste0(colnames(share_u)[61:ncol(share_u)], '_2')
share_u$g_Atopobium_2 <- share_u$Atopobium_parvulum_2 + share_u$Atopobium_vaginae_2
# share_u should have 115 obs of 76 vars

# code variables for RA ratio between participants with concordance
share_u <- share_u %>% group_by(part_bug) %>% mutate(fvag_ratio = max((Atopobium_vaginae_1 / Atopobium_vaginae_2),
                                                                      (Atopobium_vaginae_2 / Atopobium_vaginae_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(fannyhessea_ratio = max((g_Atopobium_1 / g_Atopobium_2),
                                                                             (g_Atopobium_2 / g_Atopobium_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(gard_ratio = max((Gardnerella_vaginalis_1 / Gardnerella_vaginalis_2),
                                                                      (Gardnerella_vaginalis_2 / Gardnerella_vaginalis_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(lcrisp_ratio = max((Lactobacillus_crispatus_1 / Lactobacillus_crispatus_2),
                                                                        (Lactobacillus_crispatus_2 / Lactobacillus_crispatus_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(liners_ratio = max((Lactobacillus_iners_1 / Lactobacillus_iners_2),
                                                                        (Lactobacillus_iners_2 / Lactobacillus_iners_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(ljens_ratio = max((Lactobacillus_jensenii_1 / Lactobacillus_jensenii_2),
                                                                       (Lactobacillus_jensenii_2 / Lactobacillus_jensenii_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(lvag_ratio = max((Lactobacillus_vaginalis_1 / Lactobacillus_vaginalis_2),
                                                                      (Lactobacillus_vaginalis_2 / Lactobacillus_vaginalis_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(mega_ratio = max((g_Megasphaera_1 / g_Megasphaera_2),
                                                                      (g_Megasphaera_2 / g_Megasphaera_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(mobil_ratio = max((g_Mobiluncus_1 / g_Mobiluncus_2),
                                                                       (g_Mobiluncus_2 / g_Mobiluncus_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(pamnii_ratio = max((Prevotella_amnii_1 / Prevotella_amnii_2),
                                                                        (Prevotella_amnii_2 / Prevotella_amnii_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(ptim_ratio = max((Prevotella_timonensis_1 / Prevotella_timonensis_2),
                                                                      (Prevotella_timonensis_2 / Prevotella_timonensis_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(ssang_ratio = max((Sneathia_sanguinegens_1 / Sneathia_sanguinegens_2),
                                                                       (Sneathia_sanguinegens_2 / Sneathia_sanguinegens_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(samnii_ratio = max((Sneathia_amnii_1 / Sneathia_amnii_2),
                                                                        (Sneathia_amnii_2 / Sneathia_amnii_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(staph_ratio = max((g_Staphylococcus_1 / g_Staphylococcus_2),
                                                                       (g_Staphylococcus_2 / g_Staphylococcus_1)))

share_u <- share_u %>% group_by(part_bug) %>% mutate(strep_ratio = max((g_Streptococcus_1 / g_Streptococcus_2),
                                                                       (g_Streptococcus_2 / g_Streptococcus_1)))
# share_u should have 115 obs of 91 vars

# summarize relative abundances for concordant taxa within contact dyads
summary(c(share_u$Atopobium_vaginae_1[share_u$partner == 1 & share_u$bug == 'Fannyhessea_vaginae'],
           share_u$Atopobium_vaginae_2[share_u$partner == 1 & share_u$bug == 'Fannyhessea_vaginae']))

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 1 & share_u$bug == 'Gardnerella_vaginalis_A'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 1 & share_u$bug == 'Gardnerella_vaginalis_A']))

summary(c(share_u$Prevotella_amnii_1[share_u$partner == 1 & share_u$bug == 'Prevotella_amnii'],
          share_u$Prevotella_amnii_2[share_u$partner == 1 & share_u$bug == 'Prevotella_amnii']))

summary(c(share_u$Sneathia_sanguinegens_1[share_u$partner == 1 & share_u$bug == 'Sneathia_sanguinegens'],
          share_u$Sneathia_sanguinegens_2[share_u$partner == 1 & share_u$bug == 'Sneathia_sanguinegens']))

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 1 & share_u$bug == 'Gardnerella_leopoldii'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 1 & share_u$bug == 'Gardnerella_leopoldii']))

summary(c(share_u$Sneathia_amnii_1[share_u$partner == 1 & share_u$bug == 'Sneathia_vaginalis'],
          share_u$Sneathia_amnii_2[share_u$partner == 1 & share_u$bug == 'Sneathia_vaginalis']))

summary(c(share_u$Lactobacillus_iners_1[share_u$partner == 1 & share_u$bug == 'Lactobacillus_iners'],
          share_u$Lactobacillus_iners_2[share_u$partner == 1 & share_u$bug == 'Lactobacillus_iners']))

# summarize relative abundances for concordant taxa within non-contact pairs
summary(c(share_u$g_Atopobium_1[share_u$partner == 0 & share_u$bug == 'Fannyhessea'],
          share_u$g_Atopobium_2[share_u$partner == 0 & share_u$bug == 'Fannyhessea']))
summary(share_u$fannyhessea_ratio[share_u$partner == 0 & share_u$bug == 'Fannyhessea'])

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Gardnerella_leopoldii'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Gardnerella_leopoldii']))
summary(share_u$gard_ratio[share_u$partner == 0 & share_u$bug == 'Gardnerella_leopoldii'])

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Gardnerella_piotii'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Gardnerella_piotii']))
summary(share_u$gard_ratio[share_u$partner == 0 & share_u$bug == 'Gardnerella_piotii'])

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Gardnerella_sp'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Gardnerella_sp']))
summary(share_u$gard_ratio[share_u$partner == 0 & share_u$bug == 'Gardnerella_sp'])

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Gardnerella_swidsinkii'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Gardnerella_swidsinkii']))
summary(share_u$gard_ratio[share_u$partner == 0 & share_u$bug == 'Gardnerella_swidsinkii'])

summary(c(share_u$Gardnerella_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Gardnerella_vaginalis'],
          share_u$Gardnerella_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Gardnerella_vaginalis']))
summary(share_u$gard_ratio[share_u$partner == 0 & share_u$bug == 'Gardnerella_vaginalis'])

summary(c(share_u$Lactobacillus_crispatus_1[share_u$partner == 0 & share_u$bug == 'Lactobacillus_crispatus'],
          share_u$Lactobacillus_crispatus_2[share_u$partner == 0 & share_u$bug == 'Lactobacillus_crispatus']))
summary(share_u$lcrisp_ratio[share_u$partner == 0 & share_u$bug == 'Lactobacillus_crispatus'])

summary(c(share_u$Lactobacillus_iners_1[share_u$partner == 0 & share_u$bug == 'Lactobacillus_iners'],
          share_u$Lactobacillus_iners_2[share_u$partner == 0 & share_u$bug == 'Lactobacillus_iners']))
summary(share_u$liners_ratio[share_u$partner == 0 & share_u$bug == 'Lactobacillus_iners'])

summary(c(share_u$Lactobacillus_jensenii_1[share_u$partner == 0 & share_u$bug == 'Lactobacillus_jensenii'],
          share_u$Lactobacillus_jensenii_2[share_u$partner == 0 & share_u$bug == 'Lactobacillus_jensenii']))
summary(share_u$ljens_ratio[share_u$partner == 0 & share_u$bug == 'Lactobacillus_jensenii'])

summary(c(share_u$Lactobacillus_jensenii_1[share_u$partner == 0 & share_u$bug == 'Lactobacillus_mulieris'],
          share_u$Lactobacillus_jensenii_2[share_u$partner == 0 & share_u$bug == 'Lactobacillus_mulieris']))
summary(share_u$ljens_ratio[share_u$partner == 0 & share_u$bug == 'Lactobacillus_mulieris']) 

summary(c(share_u$Lactobacillus_vaginalis_1[share_u$partner == 0 & share_u$bug == 'Limosilactobacillus_portuensis'],
          share_u$Lactobacillus_vaginalis_2[share_u$partner == 0 & share_u$bug == 'Limosilactobacillus_portuensis']))
summary(share_u$lvag_ratio[share_u$partner == 0 & share_u$bug == 'Limosilactobacillus_portuensis'])

summary(c(share_u$g_Megasphaera_1[share_u$partner == 0 & share_u$bug == 'Megasphaera_lornae'],
          share_u$g_Megasphaera_2[share_u$partner == 0 & share_u$bug == 'Megasphaera_lornae']))
summary(share_u$mega_ratio[share_u$partner == 0 & share_u$bug == 'Megasphaera_lornae'])

summary(c(share_u$g_Mobiluncus_1[share_u$partner == 0 & share_u$bug == 'Mobiluncus_curtisii'],
          share_u$g_Mobiluncus_2[share_u$partner == 0 & share_u$bug == 'Mobiluncus_curtisii']))
summary(share_u$mobil_ratio[share_u$partner == 0 & share_u$bug == 'Mobiluncus_curtisii'])

summary(c(share_u$Prevotella_amnii_1[share_u$partner == 0 & share_u$bug == 'Prevotella_amnii'],
          share_u$Prevotella_amnii_2[share_u$partner == 0 & share_u$bug == 'Prevotella_amnii']))
summary(share_u$pamnii_ratio[share_u$partner == 0 & share_u$bug == 'Prevotella_amnii'])

summary(c(share_u$Prevotella_timonensis_1[share_u$partner == 0 & share_u$bug == 'Prevotella_timonensis'],
          share_u$Prevotella_timonensis_2[share_u$partner == 0 & share_u$bug == 'Prevotella_timonensis']))
summary(share_u$ptim_ratio[share_u$partner == 0 & share_u$bug == 'Prevotella_timonensis'])

summary(c(share_u$g_Staphylococcus_1[share_u$partner == 0 & share_u$bug == 'Staphylococcus_epidermidis'],
          share_u$g_Staphylococcus_2[share_u$partner == 0 & share_u$bug == 'Staphylococcus_epidermidis']))
summary(share_u$staph_ratio[share_u$partner == 0 & share_u$bug == 'Staphylococcus_epidermidis'])

summary(c(share_u$g_Streptococcus_1[share_u$partner == 0 & share_u$bug == 'Streptococcus_agalactiae'],
          share_u$g_Streptococcus_2[share_u$partner == 0 & share_u$bug == 'Streptococcus_agalactiae']))
summary(share_u$strep_ratio[share_u$partner == 0 & share_u$bug == 'Streptococcus_agalactiae'])

summary(c(share_u$g_Streptococcus_1[share_u$partner == 0 & share_u$bug == 'Streptococcus_anginosus'],
          share_u$g_Streptococcus_2[share_u$partner == 0 & share_u$bug == 'Streptococcus_anginosus']))
summary(share_u$strep_ratio[share_u$partner == 0 & share_u$bug == 'Streptococcus_anginosus'])

summary(c(share_u$g_Streptococcus_1[share_u$partner == 0 & share_u$bug == 'Streptococcus_mitis'],
          share_u$g_Streptococcus_2[share_u$partner == 0 & share_u$bug == 'Streptococcus_mitis']))
summary(share_u$strep_ratio[share_u$partner == 0 & share_u$bug == 'Streptococcus_mitis'])






#######################################################
### Figure 5 - popANI histogram, popSNP/Mbp boxplot ###
#######################################################

### histogram of popANI colored by taxon

# make categorical factor of popANI >0.9999, >0.99999, >0.999999, >0.9999999
plot <- share_u[, colnames(share_u) %in% c('bug', 'popANI')]
plot <- plot %>% mutate(popANI_f = as.factor(case_when(popANI == 1 ~ 4,
                                                       popANI > 0.999999 ~ 3,
                                                       popANI > 0.99999 ~ 2,
                                                       popANI > 0.9999 ~ 1)))
# plot should have 115 obs of 3 vars

# plot with legend
ani_histo <- ggplot(plot, aes(x = popANI_f, color = bug, fill = bug)) +
  geom_histogram(stat = 'count') +
  scale_color_manual(values = net_colors, 
                     labels = c(expression(paste(italic('Fannyhessea'), ' spp.')),
                                expression(paste(italic('Fannyhessea vaginae'), '')),
                                expression(paste(italic('Gardnerella leopoldii'), '')),
                                expression(paste(italic('Gardnerella piotii'), '')),
                                expression(paste(italic('Gardnerella'), ' spp.')),
                                expression(paste(italic('Gardnerella swidsinskii'), '')),
                                expression(paste(italic('Gardnerella vaginalis'), '')),
                                expression(paste(italic('Gardnerella vaginalis'), ' A')),
                                expression(paste(italic('Lactobacillus crispatus'), '')),
                                expression(paste(italic('Lactobacillus iners'), '')),
                                expression(paste(italic('Lactobacillus jensenii'), '')),
                                expression(paste(italic('Lactobacillus mulieris'), '')),
                                expression(paste(italic('Limosilactobacillus portuensis'), '       ')),
                                expression(paste(italic('Megasphaera lornae'), '')),
                                expression(paste(italic('Mobiluncus curtisii'), '')),
                                expression(paste(italic('Prevotella'), ' spp.')),
                                expression(paste(italic('Prevotella amnii'), '')),
                                expression(paste(italic('Prevotella'), ' sp000758925')),
                                expression(paste(italic('Prevotella timonensis'), '')),
                                expression(paste(italic('Sneathia sanguinegens'), '')),
                                expression(paste(italic('Sneathia vaginalis'), '')),
                                expression(paste(italic('Staphylococcus epidermidis   '), '')),
                                expression(paste(italic('Streptococcus agalactiae'), '')),
                                expression(paste(italic('Streptococcus anginosus'), '')),
                                expression(paste(italic('Streptococcus mitis'), ''))),
                     'Taxon',
                     guide = guide_legend(ncol = 1, byrow = TRUE)) +
  scale_fill_manual(values = net_colors,
                    labels = c(expression(paste(italic('Fannyhessea'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Gardnerella leopoldii'), '')),
                               expression(paste(italic('Gardnerella piotii'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Gardnerella swidsinskii'), '')),
                               expression(paste(italic('Gardnerella vaginalis'), '')),
                               expression(paste(italic('Gardnerella vaginalis'), ' A')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus mulieris'), '')),
                               expression(paste(italic('Limosilactobacillus portuensis'), '       ')),
                               expression(paste(italic('Megasphaera lornae'), '')),
                               expression(paste(italic('Mobiluncus curtisii'), '')),
                               expression(paste(italic('Prevotella'), ' spp.')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Prevotella'), ' sp000758925')),
                               expression(paste(italic('Prevotella timonensis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Staphylococcus epidermidis   '), '')),
                               expression(paste(italic('Streptococcus agalactiae'), '')),
                               expression(paste(italic('Streptococcus anginosus'), '')),
                               expression(paste(italic('Streptococcus mitis'), ''))),
                    'Taxon',
                    guide = guide_legend(ncol = 1, byrow = TRUE))  +
  scale_x_discrete(labels = c('> 99.99%', '> 99.999%', '> 99.9999%', '100%')) +
  xlab('popANI') +
  ylab('Count') +
  theme_light(base_size = 10) +
  theme(legend.text.align = 0, 
        legend.key.size = unit(0.75, 'char'), 
        legend.key.height = unit(0.75, 'char'),
        legend.spacing.y = unit(0.03, 'in'),
        legend.title = element_text(margin = margin(t = 15, b = 6)),
        axis.title.x = element_text(vjust = -0.75),
        axis.title.y = element_text(vjust = 2))
# ani_histo

# plot with no legend for patchwork
ani_histo_no <- ggplot(plot, aes(x = popANI_f, color = bug, fill = bug)) +
  geom_histogram(stat = 'count',
                 show.legend = F) +
  scale_color_manual(values = net_colors) +
  scale_fill_manual(values = net_colors)  +
  scale_x_discrete(labels = c('> 99.99%', '> 99.999%', '> 99.9999%', '100%')) +
  xlab('popANI') +
  ylab('Count') +
  ggtitle('A') +
  theme_light(base_size = 10) +
  theme(axis.title.x = element_text(vjust = -0.75),
        axis.title.y = element_text(vjust = 2),
        plot.title = element_text(size = 15))
# ani_histo_no



### box and whisker, dot plot of population SNPs/Mbp in concordant strains

# plot with legend
plot_SNP_mbp <- ggplot(share_u, aes(x = bug, y = popSNP_mbp)) +
  geom_boxplot(aes(fill = bug), outlier.shape = NA,
               linewidth = c(rep(1.05, 2),
                             0.7,
                             rep(1.05, 2),
                             0.7,
                             rep(1.05, 2),
                             rep(0.7, 4),
                             rep(1.05, 3),
                             rep(0.7, 2),
                             rep(1.05, 8)),
               color = c('#3388ff', # Fannyhessea spp.
                         '#0000cd', # F. vaginae
                         'black',
                         '#AFCCB7', # G. piotii
                         '#c1ffc1', # Gardnerella spp.
                         'black',
                         '#73FCD6', # G. vaginalis
                         '#1C792B', # G. vaginalis A
                         rep('black', 4),
                         '#B5F46E', # L. portuensis
                         '#FDAA8B', # M. lornae
                         '#7A4300', # M. curtisii
                         rep('black', 2),
                         '#b0b0b0', # Prevotella sp000758925
                         '#f8fca9', # P timonensis
                         '#c7aa8f', # S. sanguinegens
                         '#6963ff', # S. vaginalis
                         '#800080', # S. epidermidis
                         '#ff738a', # S. agalactiae
                         '#ffc0cb', # S. anginosus
                         '#97C0DD')) + # S. mitis 
  scale_fill_manual(values = net_colors, 
                    labels = c(expression(paste(italic('Fannyhessea'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Gardnerella leopoldii'), '')),
                               expression(paste(italic('Gardnerella piotii'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Gardnerella swidsinskii'), '')),
                               expression(paste(italic('Gardnerella vaginalis'), '')),
                               expression(paste(italic('Gardnerella vaginalis'), ' A')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus mulieris'), '')),
                               expression(paste(italic('Limosilactobacillus portuensis'), '       ')),
                               expression(paste(italic('Megasphaera lornae'), '')),
                               expression(paste(italic('Mobiluncus curtisii'), '')),
                               expression(paste(italic('Prevotella'), ' spp.')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Prevotella'), ' sp000758925')),
                               expression(paste(italic('Prevotella timonensis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Staphylococcus epidermidis   '), '')),
                               expression(paste(italic('Streptococcus agalactiae'), '')),
                               expression(paste(italic('Streptococcus anginosus'), '')),
                               expression(paste(italic('Streptococcus mitis'), ''))),
                    'Taxon',
                    guide = guide_legend(ncol = 1, order = 1,
                                         override.aes = list(linetype = 0,
                                                             shape = 0))) +
  geom_jitter(aes(shape = partner_f, size = partner_f), 
              color = 'black',
              width = 0.12, height = 0.1) +
  scale_shape_manual(values = c(4, 19),
                     labels = c('Non-contacts', 'Sexual contacts'),
                     'Concordance between:',
                     guide = guide_legend(reverse = T,
                                          override.aes = list(size = c(3, 2.5)),
                                          order = 2)) +
  scale_size_manual(values = c(0.85, 1.25),
                    labels = c('Non-contacts', 'Sexual contacts'),
                    'Concordance between:',
                    guide = guide_legend(reverse = T,
                                         override.aes = list(size = c(4, 3)),
                                         order = 2)) +
  xlab('Taxon') +
  ylab('popSNP/Mbp') +
  theme_light(base_size = 10) +
  theme(legend.text.align = 0, legend.spacing = unit(0, 'in'),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        legend.key.size = unit(1, 'char'), 
        legend.key.height = unit(1, 'char'))
# plot_SNP_mbp

# pull legend for patchwork
legend_ani_snp <- as_ggplot(get_legend(plot_SNP_mbp))

# no legend, for patchwork
plot_SNP_mbp_no <- ggplot(share_u, aes(x = bug, y = popSNP_mbp)) +
  geom_boxplot(aes(fill = bug), outlier.shape = NA,
               linewidth = c(rep(1.05, 2),
                             0.7,
                             rep(1.05, 2),
                             0.7,
                             rep(1.05, 2),
                             rep(0.7, 4),
                             rep(1.05, 3),
                             rep(0.7, 2),
                             rep(1.05, 8)),
               color = c('#3388ff', # Fannyhessea spp.
                         '#0000cd', # F. vaginae
                         'black',
                         '#AFCCB7', # G. piotii
                         '#c1ffc1', # Gardnerella spp.
                         'black',
                         '#73FCD6', # G. vaginalis
                         '#1C792B', # G. vaginalis A
                         rep('black', 4),
                         '#B5F46E', # L. portuensis
                         '#FDAA8B', # M. lornae
                         '#7A4300', # M. curtisii
                         rep('black', 2),
                         '#b0b0b0', # Prevotella sp000758925
                         '#f8fca9', # P timonensis
                         '#c7aa8f', # S. sanguinegens
                         '#6963ff', # S. vaginalis
                         '#800080', # S. epidermidis
                         '#ff738a', # S. agalactiae
                         '#ffc0cb', # S. anginosus
                         '#97C0DD'), # S. mitis 
               show.legend = F) + 
  scale_fill_manual(values = net_colors) +
  geom_jitter(aes(shape = partner_f, size = partner_f), 
              color = 'black',
              width = 0.12, height = 0.1,
              show.legend = F) +
  scale_shape_manual(values = c(4, 19)) +
  scale_size_manual(values = c(0.85, 1.25)) +
  xlab('Taxon') +
  ylab('popSNP/Mbp') +
  ggtitle('B') +
  theme_light(base_size = 10) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        plot.title = element_text(size = 15))
# plot_SNP_mbp_no



### patchwork of ani_histo_no, plot_SNP_mbp_no, legend_ani_snp

# generate plot spacer to use in final patchwork plot, need spacer because end of 
# 'Limosilactobacillus portuensis' in legend is cut off in saved PDF
plot_spacer <- plot_spacer()

# patchwork together ani_histo_no, plot_SNP_mbp_no, legend_ani_snp, plot_spacer
# area(Top, Left, Bottom, Rright)
# (ani_histo_no | plot_SNP_mbp_no | legend_ani_snp | plot_spacer) +
#   plot_layout(design = c(area(1, 1, 98, 177), # ANI histogram
#                          area(102, 1, 200, 177), # SNP box plot
#                          area(1, 188, 200, 260), # legend
#                          area(1, 260, 200, 262))) # spacer






###############################################
### Figure 2 - sample-wise stacked bar plot ###
###############################################

# excluding chlamydia trachomatis from stacked bar plots because it was more likely 
# to be easily visualized on the plot in male CT+ samples then female CT+ samples, 
# likely due penile urethral microbiomes having lower overall bacterial abundance, 
# which would increase the likelihood of recovering chlamydia sequences in penile 
# samples compared to vaginal since STI pathogens are typically present at much 
# lower abundance than members of the microbiota - including it may suggest spurious 
# relationships between sex and chlamydia abundance, so we chose to exclude it

# drop ct
person_noCT <- person[, colnames(person) %!in% c('Chlamydia_trachomatis')]
# person_noCT should have 138 obs 217 vars

# identify taxa with top 10 mean relative abundances among women
mean_f <- colSums(person[person$sex == 'f', 32:164], na.rm = T) / sum(person$sex == 'f')
names(mean_f) <- colnames(person)[32:164]
mean_f <- sort(mean_f, decreasing = T)
# mean_f[1:10]
top10_f <- names(mean_f)[1:10]
# Gardnerella_vaginalis, Lactobacillus_iners, Lactobacillus_crispatus, BVAB1, Atopobium_vaginae
# Lactobacillus_jensenii, g_Streptococcus, Prevotella_bivia, Prevotella_amnii, g_Megasphaera

# identify taxa with top 10 mean relative abundances among men
mean_m <- colSums(person_noCT[person_noCT$sex == 'm', 32:164], na.rm = T) / sum(person_noCT$sex == 'm')
names(mean_m) <- colnames(person_noCT)[32:164]
mean_m <- sort(mean_m, decreasing = T)
# mean_m[1:10]
top10_m <- names(mean_m)[1:10]
# Gardnerella_vaginalis, g_Propionimicrobium, g_Corynebacterium, g_Ureaplasma, g_Streptococcus,
# Lactobacillus_iners, g_Staphylococcus, Sneathia_amnii, g_Neisseria, Atopobium_vaginae

# identify unique set of taxa in female top 10, male top 10
tax <- unique(c(top10_f, top10_m))
# tax should have 16 elements

# add any taxa that have putative sexual transmission but are not in top 10 mean 
# for female or male participants
# unique(share_u$bug[share_u$trans == 1])
# tax
# need to add Sneathia_sanguinegens 
tax <- append(tax, 'Sneathia_sanguinegens')
# tax should have 17 elements

# generate wide dataframe limited to variables of interest, add variable for 'other' rel abundance
plot <- person[, colnames(person) %in% c('barcode', 'sex', 'ct', 'cst', 'wave', tax)]
plot$other <- 1 - rowSums(plot[, 6:ncol(plot)], na.rm = T)
# plot should have 138 obs of 23 vars

# rename some bugs to match names in colors, remove leading 'g_'
plot <- plot %>% rename(Ca_Lachnocurva_vaginae = BVAB1,
                        Fannyhessea_vaginae = Atopobium_vaginae,
                        Sneathia_vaginalis = Sneathia_amnii,
                        g_Gardnerella = Gardnerella_vaginalis)

# convert to long
plot <- reshape2::melt(plot,
                       id.vars = c('barcode', 'sex', 'ct', 'cst', 'wave'),
                       measure.vars = colnames(plot)[colnames(plot) %!in% c('barcode', 'sex', 'ct', 'cst', 'wave')])
plot <- plot %>% rename(bug = variable,
                        ra = value)
# plot should have 2484 obs of 7 vars

# factor bug for desired plotting order for taxa (alphabetical) to color plot by
plot$bug <- factor(plot$bug, levels = c('Ca_Lachnocurva_vaginae',
                                        'g_Corynebacterium',
                                        'Fannyhessea_vaginae', 
                                        'g_Gardnerella',
                                        'Lactobacillus_crispatus',
                                        'Lactobacillus_iners',
                                        'Lactobacillus_jensenii',
                                        'g_Megasphaera',
                                        'g_Neisseria',
                                        'Prevotella_amnii',
                                        'Prevotella_bivia',
                                        'g_Propionimicrobium',
                                        'Sneathia_sanguinegens',
                                        'Sneathia_vaginalis', 
                                        'g_Staphylococcus',
                                        'g_Streptococcus',
                                        'g_Ureaplasma',
                                        'other'))

# make barcode in plot into factor that is ordered according sex, within sex 
# ordered by CST, within CST by ordered by decreasing RA of characteristic taxon 
# so that plot is grouped by sex, CST, RA 
female <- person[person$sex == 'f',]
female <- female[order(female$cst),]
male <- person[person$sex == 'm',]
male <- male[order(male$cst),]

female_I <- female[is.na(female$cst) == F & female$cst == 'I',]
female_I <- female_I %>% arrange(desc(Lactobacillus_crispatus))

female_III <- female[is.na(female$cst) == F & female$cst == 'III',]
female_III <- female_III %>% arrange(desc(Lactobacillus_iners))

female_IVA <- female[is.na(female$cst) == F & female$cst == 'IV-A',]
female_IVA <- female_IVA %>% arrange(desc(BVAB1))

female_IVB <- female[is.na(female$cst) == F & female$cst == 'IV-B',]
female_IVB <- female_IVB %>% arrange(desc(Gardnerella_vaginalis))

female_IVC <- female[is.na(female$cst) == F & female$cst == 'IV-C',]
female_IVC <- female_IVC %>% arrange(desc(g_Corynebacterium))

female_V <- female[is.na(female$cst) == F & female$cst == 'V',]
female_V <- female_V %>% arrange(desc(Lactobacillus_jensenii))

female <- rbind(female_I, female_III, female_IVA, female_IVB, female_IVC, female_V)
# female should have 74 obs of 218 vars

# no male CST I

male_III <- male[is.na(male$cst) == F & male$cst == 'III',]
male_III <- male_III %>% arrange(desc(Lactobacillus_iners))

male_IVA <- male[is.na(male$cst) == F & male$cst == 'IV-A',]
male_IVA <- male_IVA %>% arrange(desc(BVAB1))

male_IVB <- male[is.na(male$cst) == F & male$cst == 'IV-B',]
male_IVB <- male_IVB %>% arrange(desc(Gardnerella_vaginalis))

male_IVC <- male[is.na(male$cst) == F & male$cst == 'IV-C',]
male_IVC <- male_IVC %>% arrange(desc(g_Corynebacterium))

male_V <- male[is.na(male$cst) == F & male$cst == 'V',]
male_V <- male_V %>% arrange(desc(Lactobacillus_jensenii))

male <- rbind(male_III, male_IVA, male_IVB, male_IVC, male_V)
# male should have 64 obs of 218 vars

# going to factor barcode in plot dataset according to order in female and male,
# but adding some extra barcodes between female and male to give some empty space
# between female and male in plot for ease of viewing
bars <- c(female$bar, 'empty_1', 'empty_2', male$barcode)

# relevel plot$barcode as ordered in person
plot$barcode <- factor(plot$barcode, levels = bars)

# some samples missing CT status and plotting at default gray, scale_fill_manual
# won't let me assign a color for NAs so instead changing NA so can color accordingly
plot$ct <- as.numeric(plot$ct)
# tab1(plot$ct)
# values are 1 for no CT, 2 for CT, change NA to 0
plot$ct[is.na(plot$ct) == T] <- 0
# factor
plot$ct <- factor(plot$ct, levels = c(2, 1, 0))

# add empty spacer barcodes
empty_1 <- c('empty_1', NA, NA, NA, NA, NA, 0.00)
empty_2 <- c('empty_2', NA, NA, NA, NA, NA, 0.00)
plot <- rbind(plot, empty_1, empty_2)
plot$ra <- as.numeric(plot$ra)

# define colors for stacked bar plots
stack_colors <- c('#ff0000', 
                  '#ff8c00', 
                  '#333333', 
                  '#0000cd', 
                  '#20b2aa', 
                  '#008B45', 
                  '#bfbfbf', 
                  '#808080', 
                  '#b31900', 
                  '#ffc0cb', 
                  '#800080', 
                  '#ac8fc7', 
                  '#6963ff', 
                  '#D1279C',  
                  '#6AB8F6', 
                  '#B8A404', 
                  '#ffff00',
                  '#c7aa8f') 

names(stack_colors) <- c('Lactobacillus_crispatus', 
                         'Lactobacillus_iners', 
                         'Lactobacillus_jensenii', 
                         'Fannyhessea_vaginae', 
                         'g_Gardnerella', 
                         'g_Megasphaera', 
                         'Prevotella_bivia', 
                         'other', 
                         'Ca_Lachnocurva_vaginae', 
                         'g_Streptococcus', 
                         'g_Staphylococcus', 
                         'Prevotella_amnii', 
                         'Sneathia_vaginalis',  
                         'g_Neisseria', 
                         'g_Propionimicrobium', 
                         'g_Ureaplasma', 
                         'g_Corynebacterium',
                         'Sneathia_sanguinegens') 

stack_colors <- stack_colors[order(factor(names(stack_colors), levels = rev(levels(plot$bug))))]

# plot with all legends
stack_bar_all <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug))) +
  scale_fill_manual(values = stack_colors,  
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.02)) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.02)) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.02)) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stack_bar_all

# plot with taxon legend only
stack_bar_tax <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug))) +
  scale_fill_manual(values = stack_colors,
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  scale_x_discrete(labels = factor(plot$barcode)) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(angle=90),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stack_bar_tax

# pull taxon legend for patchworking
legend_tax <- as_ggplot(get_legend(stack_bar_tax))

# plot with cst legend only
stack_bar_cst <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug)),
           show.legend = F) +
  scale_fill_manual(values = stack_colors, 
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.02)) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stack_bar_cst

# pull cst legend for patchworking
legend_cst <- as_ggplot(get_legend(stack_bar_cst))

# plot with sex legend only
stack_bar_sex <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug)),
           show.legend = F) +
  scale_fill_manual(values = stack_colors, 
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.02)) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stack_bar_sex

# pull taxon legend for patchworking
legend_sex <- as_ggplot(get_legend(stack_bar_sex))

# plot with ct legend only
stack_bar_ct <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug)),
           show.legend = F) +
  scale_fill_manual(values = stack_colors,  
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.02),
            show.legend = F) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.02)) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# stack_bar_ct

# pull ct legend for patchworking
legend_ct <- as_ggplot(get_legend(stack_bar_ct))

# plot with NO legends
stack_bar_no <- ggplot(plot, aes(x = factor(barcode), y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(bug)),
           show.legend = F) +
  scale_fill_manual(values = stack_colors,  
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 1,
                             override.aes = list(size = 2),
                             order = 4)) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.07, fill = cst,
                                height = 0.03),
            show.legend = F) +
  scale_fill_manual(values = c('#FF0000', # I
                               '#FF8C00', # III
                               '#1F8C72', # IV-A
                               '#0044FF', # IV-B
                               '#A2805D', # IV-C
                               '#FFF200'), # V
                    'CST',
                    na.translate = F) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.11, fill = sex,
                                height = 0.03),
            show.legend = F) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  guides(fill = guide_legend(order = 2, override.aes = list(size = 2))) +
  new_scale('fill') +
  geom_tile(plot, mapping = aes(x = factor(barcode), y = 1.03, fill = ct,
                                height = 0.03),
            show.legend = F) +
  scale_fill_manual(values = c('#000000', 
                               'grey38',
                               'grey76'),
                    labels = c('Yes', 
                               'No',
                               'Missing'),
                    'Genital CT',
                    na.translate = F) +
  guides(fill = guide_legend(order = 3, override.aes = list(size = 2))) +
  annotate(geom = 'text', x = 150.1, y = 1.07,
           label = 'CST', hjust = 1,
           color = 'black', size = 2.7) +
  annotate(geom = 'text', x = 149.1, y = 1.11,
           label = 'Sex', hjust = 1,
           color = 'black', size = 2.7) +
  annotate(geom = 'text', x = 147.7, y = 1.03,
           label = 'CT', hjust = 1,
           color = 'black', size = 2.7) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(0.85, 'char'), legend.key.height = unit(0.85, 'char'), legend.text.align = 0,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(0, 10, -20, 10))
# stack_bar_no

# patchwork together stack_bar_no, legend_cst, legend_ct, legend_sex, legend_tax
# area(Top, Left, Bottom, Right)
# (stack_bar_patch <- (stack_bar_no | 
#                        legend_sex | 
#                        legend_cst | 
#                        legend_ct | 
#                        legend_tax |
#                        plot_spacer)) +
#   plot_layout(design = c(area(1, 1, 110, 220), # stacked bar and tiles
#                          area(19, 215, 34, 240), # sex legend
#                          area(48, 208, 70, 240), # CST legend
#                          area(87, 219.5, 102, 240), # CT legend
#                          area(15, 251, 100, 318), # taxon legend
#                          area(1, 318, 110, 319))) # spacer






##############################################################################################
### Figure 6 - stacked bar with panel for each sexual contact dyad with concordant strains ###
##############################################################################################

# generate wide dataframe limited to sexual contacts and variables of interest, 
# add variable for 'other' rel abund
plot <- share_u[share_u$partner == 1, colnames(share_u) %in% c('name1', 'name2', 'bug', 'part_barcode', 'sex_1', 'sex_2')]
# length(unique(c(plot$part_barcode)))
# need 4 panels, one panel for each unique contact dyad with concordant strains
# order plot by part_barcode
plot <- plot %>%
  arrange(part_barcode)
# restrict to first observation for each pair
plot <- do.call(rbind, by(plot, list(plot$part_barcode), 
                          FUN=function(x) head(x, 1)))
# plot should have 4 obs of 6 vars

# split out second member of dyad as own observations and restructure plot dataframe
names <- c(plot$name1, plot$name2)
bugs <- rep(plot$bug, 2)
part_barcodes <- rep(plot$part_barcode, 2)
sexes <- c(plot$sex_1, plot$sex_2)
plot <- as.data.frame(cbind(names, bugs, part_barcodes, sexes))
colnames(plot) <- c('barcode', 'bug', 'part_barcode', 'sex')
# plot should have 8 obs of 4 vars

# merge in taxon relative abundances, add other RA variable
plot <- left_join(plot, person[, colnames(person) %in% c('barcode', tax)])
plot$other <- 1 - rowSums(plot[ , 5:ncol(plot)])
# plot should have 8 obs of 22 vars

# convert to long
plot <- reshape2::melt(plot,
                       id.vars = colnames(plot)[colnames(plot) %!in% c(tax, 'other')],
                       measure.vars = c(tax, 'other'))
plot <- plot %>% rename(tax = variable,
                        ra = value)
# plot should have 144 obs of 6 vars

# for each heterosexual panel, want female on left, male on right. for WSW dyad
# will have 722441 on right (not for any particular reason) but still want x 
# axis labels to indicate sex so going to do three level factor with labels
plot$samp <- ifelse(plot$sex == '1', 1, 2)
plot$samp[plot$barcode == '722441'] <- 3
plot$samp <- factor(plot$samp, labels = c('F', 'M', ' F '))
# plot should have 144 obs of 7 vars

# recode some taxon names to match color names
plot$tax <- as.character(plot$tax)
plot$tax[plot$tax == 'Atopobium_vaginae'] <- 'Fannyhessea_vaginae'
plot$tax[plot$tax == 'BVAB1'] <- 'Ca_Lachnocurva_vaginae'
plot$tax[plot$tax == 'Gardnerella_vaginalis'] <- 'g_Gardnerella'
plot$tax[plot$tax == 'Sneathia_amnii'] <- 'Sneathia_vaginalis'

# factor tax_name for desired plotting order for taxa (alphabetical) to color plot by
plot$tax <- factor(plot$tax, levels = c('Ca_Lachnocurva_vaginae',
                                        'g_Corynebacterium',
                                        'Fannyhessea_vaginae', 
                                        'g_Gardnerella',
                                        'Lactobacillus_crispatus',
                                        'Lactobacillus_iners',
                                        'Lactobacillus_jensenii',
                                        'g_Megasphaera',
                                        'g_Neisseria',
                                        'Prevotella_amnii',
                                        'Prevotella_bivia',
                                        'g_Propionimicrobium',
                                        'Sneathia_sanguinegens',
                                        'Sneathia_vaginalis', 
                                        'g_Staphylococcus',
                                        'g_Streptococcus',
                                        'g_Ureaplasma',
                                        'other'))

# factor bug for desired plotting order of panels 
plot$bug <- factor(plot$bug, levels = c('Fannyhessea_vaginae',
                                        'Gardnerella_leopoldii',
                                        'Sneathia_vaginalis',
                                        'Lactobacillus_iners'))

# factor part_barcode for desired panel plotting order
plot$part_barcode <- factor(plot$part_barcode, levels = c('711740127049',
                                                          '779174126428',
                                                          '830051752480',
                                                          '565548722441'))

# manually generate facet labels that say what was concordant
labs <- c('565548722441' = 'L. iners',
          '711740127049' = 'F. vaginae, P. amnii,\nS. sanguinegens',
          '779174126428' = 'G. leopoldii',
          '830051752480' = 'S. vaginalis')

# plot
stack_bar_concord <- ggplot(plot, aes(x = samp, y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(tax))) +
  scale_fill_manual(values = stack_colors, 
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '   ')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 3, title.position = 'top')) +
  facet_wrap(~part_barcode, ncol = 2, scales = 'free',
             labeller = labeller(part_barcode = labs, multi_line = T)) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_light(base_size = 9) +
  theme(legend.key.size = unit(0.9, 'char'), legend.key.height = unit(0.9, 'char'), 
        legend.text.align = 0, legend.position = 'bottom',
        legend.margin = margin(-10, 0, 0, 0), legend.title.align = 0.5,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'italic', margin = margin(0.05, 0, 0.05, 0, 'cm'), color = 'black'), #size = 7,
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)))
# stack_bar_concord






####################################################################################################
### Figure S2 - stacked bar with panel for each dyad with L. crispatus,  L. jensenii concordance ###
####################################################################################################

# generate wide dataframe limited to these participants and variables of interest, 
# add variable for 'other' rel abund
plot <- share_u[share_u$bug %in% c('Lactobacillus_crispatus', 'Lactobacillus_jensenii'), 
                colnames(share_u) %in% c('name1', 'name2', 'bug', 'part_barcode', 'sex_1', 'sex_2')]
# length(unique(c(plot$part_barcode)))
# need 20 panels because one panel for each unique contact dyad with concordant strains
# order plot by part_barcode
plot <- plot %>%
  arrange(part_barcode)
# restrict to first observation for each pair
plot <- do.call(rbind, by(plot, list(plot$part_barcode), 
                          FUN=function(x) head(x, 1)))
# plot should have 20 obs of 6 vars

# split out second member of dyad as own observations and restructure plot dataframe
names <- c(plot$name1, plot$name2)
bugs <- rep(plot$bug, 2)
part_barcodes <- rep(plot$part_barcode, 2)
sexes <- c(plot$sex_1, plot$sex_2)
plot <- as.data.frame(cbind(names, bugs, part_barcodes, sexes))
colnames(plot) <- c('barcode', 'bug', 'part_barcode', 'sex')
# plot should have 40 obs of 4 vars

# merge in taxon relative abundances, add other RA variable
plot <- left_join(plot, person[, colnames(person) %in% c('barcode', tax)])
plot$other <- 1 - rowSums(plot[ , 5:ncol(plot)])
# plot should have 40 obs of 22 vars

# convert to long
plot <- reshape2::melt(plot,
                       id.vars = colnames(plot)[colnames(plot) %!in% c(tax, 'other')],
                       measure.vars = c(tax, 'other'))
plot <- plot %>% rename(tax = variable,
                        ra = value)
# plot should have 720 obs of 6 vars

# recode some taxon names to match color names
plot$tax <- as.character(plot$tax)
plot$tax[plot$tax == 'Atopobium_vaginae'] <- 'Fannyhessea_vaginae'
plot$tax[plot$tax == 'BVAB1'] <- 'Ca_Lachnocurva_vaginae'
plot$tax[plot$tax == 'Gardnerella_vaginalis'] <- 'g_Gardnerella'
plot$tax[plot$tax == 'Sneathia_amnii'] <- 'Sneathia_vaginalis'

# factor tax_name for desired plotting order for taxa (alphabetical) to color plot by
plot$tax <- factor(plot$tax, levels = c('Ca_Lachnocurva_vaginae',
                                        'g_Corynebacterium',
                                        'Fannyhessea_vaginae', 
                                        'g_Gardnerella',
                                        'Lactobacillus_crispatus',
                                        'Lactobacillus_iners',
                                        'Lactobacillus_jensenii',
                                        'g_Megasphaera',
                                        'g_Neisseria',
                                        'Prevotella_amnii',
                                        'Prevotella_bivia',
                                        'g_Propionimicrobium',
                                        'Sneathia_sanguinegens',
                                        'Sneathia_vaginalis', 
                                        'g_Staphylococcus',
                                        'g_Streptococcus',
                                        'g_Ureaplasma',
                                        'other'))

# factor part_barcode for desired panel plotting order
plot$part_barcode <- factor(plot$part_barcode, levels = c('151881255671',
                                                          '151881284782',
                                                          '151881438567',
                                                          '151881734660',
                                                          '151881876520',
                                                          '151881913319',
                                                          '284782438567',
                                                          '438567255671',
                                                          '438567734660',
                                                          '867207255671',
                                                          '876520255671',
                                                          '876520438567',
                                                          '876520734660',
                                                          '913319255671',
                                                          '913319438567',
                                                          '913319734660',
                                                          '913319876520',
                                                          '255671734660',
                                                          '284782255671',
                                                          '284782734660'))

# manually generate facet labels that say what was concordant
labs <- c('151881255671' = 'L. crispatus',
          '151881284782' = 'L. crispatus',
          '151881438567' = 'L. crispatus',
          '151881734660' = 'L. crispatus',
          '151881876520' = 'L. crispatus',
          '151881913319' = 'L. crispatus',
          '284782438567' = 'L. crispatus',
          '438567255671' = 'L. crispatus',
          '438567734660' = 'L. crispatus',
          '867207255671' = 'L. crispatus',
          '876520255671' = 'L. crispatus',
          '876520438567' = 'L. crispatus',
          '876520734660' = 'L. crispatus',
          '913319255671' = 'L. crispatus',
          '913319438567' = 'L. crispatus',
          '913319734660' = 'L. crispatus',
          '913319876520' = 'L. crispatus',
          '255671734660' = 'L. crispatus,\nL. jensenii',
          '284782255671' = 'L. jensenii',
          '284782734660' = 'L. jensenii')

# plot
stack_bar_lacto <- ggplot(plot, aes(x = barcode, y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(tax))) +
  scale_fill_manual(values = stack_colors, 
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 3, title.position = 'top')) +
  facet_wrap(~part_barcode, ncol = 4, scales = 'free',
             labeller = labeller(part_barcode = labs, multi_line = T)) +
  scale_x_discrete(labels = c('F', 'F')) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_light(base_size = 9) +
  theme(legend.key.size = unit(0.9, 'char'), legend.key.height = unit(0.9, 'char'), 
        legend.text.align = 0, legend.position = 'bottom',
        legend.margin = margin(-10, 0, 0, 0), legend.title.align = 0.5,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'italic', margin = margin(0.05, 0, 0.05, 0, 'cm'), color = 'black'), #size = 7,
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)))
# stack_bar_lacto






########################################################################################
### Figure S3 - stacked bar with panel for each dyad with G. swidsinskii concordance ###
########################################################################################

# generate wide dataframe limited to these participants and variables of interest, 
# add variable for 'other' rel abund
plot <- share_u[share_u$bug %in% c('Gardnerella_swidsinkii'), 
                colnames(share_u) %in% c('name1', 'name2', 'bug', 'part_barcode', 'sex_1', 'sex_2')]
# length(unique(c(plot$part_barcode)))
# need 52 panels because one panel for each unique contact dyad with concordant strains
# order plot by part_barcode
plot <- plot %>%
  arrange(part_barcode)
# restrict to first observation for each pair
plot <- do.call(rbind, by(plot, list(plot$part_barcode), 
                          FUN=function(x) head(x, 1)))
# plot should have 52 obs of 6 vars

# for each F-M panel, want female on left, male on right. want x axis labels for 
# all so going to make three level variable
plot$sex_1 <- as.character(plot$sex_1)
plot$sex_2 <- as.character(plot$sex_2)
plot$sex_2[plot$sex_2 == 'f'] <- ' f '

# split out second member of dyads as own observations and restructure plot dataframe
names <- c(plot$name1, plot$name2)
bugs <- rep(plot$bug, 2)
part_barcodes <- rep(plot$part_barcode, 2)
sexes <- c(plot$sex_1, plot$sex_2)
plot <- as.data.frame(cbind(names, bugs, part_barcodes, sexes))
colnames(plot) <- c('barcode', 'bug', 'part_barcode', 'sex')
plot$sex <- toupper(plot$sex)
# plot should have 104 obs of 4 vars

# want F-M panels at end, factor part_barcode accordingly
fm <- c('084720608644', '189168608644', '289622608644', '374018608644', 
        '565788608644', '608644324554', '608644847464', '722441608644', 
        '999308608644', '999308779003')
ff <- unique(plot$part_barcode[plot$part_barcode %!in% fm])
plot$part_barcode <- factor(plot$part_barcode, levels = c(ff, fm))

# merge in taxon relative abundances, add other RA variable
plot <- left_join(plot, person[, colnames(person) %in% c('barcode', tax)])
plot$other <- 1 - rowSums(plot[ , 5:ncol(plot)])
# plot should have 104 obs of 22 vars

# convert to long
plot <- reshape2::melt(plot,
                       id.vars = colnames(plot)[colnames(plot) %!in% c(tax, 'other')],
                       measure.vars = c(tax, 'other'))
plot <- plot %>% rename(tax = variable,
                        ra = value)
# plot should have 1872 obs of 6 vars

# recode some taxon names to match color names
plot$tax <- as.character(plot$tax)
plot$tax[plot$tax == 'Atopobium_vaginae'] <- 'Fannyhessea_vaginae'
plot$tax[plot$tax == 'BVAB1'] <- 'Ca_Lachnocurva_vaginae'
plot$tax[plot$tax == 'Gardnerella_vaginalis'] <- 'g_Gardnerella'
plot$tax[plot$tax == 'Sneathia_amnii'] <- 'Sneathia_vaginalis'

# factor tax_name for desired plotting order for taxa (alphabetical) to color plot by
plot$tax <- factor(plot$tax, levels = c('Ca_Lachnocurva_vaginae',
                                        'g_Corynebacterium',
                                        'Fannyhessea_vaginae', 
                                        'g_Gardnerella',
                                        'Lactobacillus_crispatus',
                                        'Lactobacillus_iners',
                                        'Lactobacillus_jensenii',
                                        'g_Megasphaera',
                                        'g_Neisseria',
                                        'Prevotella_amnii',
                                        'Prevotella_bivia',
                                        'g_Propionimicrobium',
                                        'Sneathia_sanguinegens',
                                        'Sneathia_vaginalis', 
                                        'g_Staphylococcus',
                                        'g_Streptococcus',
                                        'g_Ureaplasma',
                                        'other'))

# plot
stack_bar_swid <- ggplot(plot, aes(x = sex, y = ra)) +
  geom_bar(stat = 'identity',
           aes(fill = fct_rev(tax))) +
  scale_fill_manual(values = stack_colors, 
                    labels = c('Other',
                               expression(paste(italic('Ureaplasma'), ' spp.')),
                               expression(paste(italic('Streptococcus'), ' spp.')),
                               expression(paste(italic('Staphylococcus'), ' spp.')),
                               expression(paste(italic('Sneathia vaginalis'), '')),
                               expression(paste(italic('Sneathia sanguinegens'), '')),
                               expression(paste(italic('Propionimicrobium'), ' spp.')),
                               expression(paste(italic('Prevotella bivia'), '')),
                               expression(paste(italic('Prevotella amnii'), '')),
                               expression(paste(italic('Neisseria'), ' spp.')),
                               expression(paste(italic('Megasphaera'), ' spp.')),
                               expression(paste(italic('Lactobacillus jensenii'), '')),
                               expression(paste(italic('Lactobacillus iners'), '')),
                               expression(paste(italic('Lactobacillus crispatus'), '')),
                               expression(paste(italic('Gardnerella'), ' spp.')),
                               expression(paste(italic('Fannyhessea vaginae'), '')),
                               expression(paste(italic('Corynebacterium'), ' spp.')),
                               expression(paste(italic('Ca.'), 'Lachnocurva vaginae'))),
                    'Taxon') +
  guides(fill = guide_legend(reverse = T, ncol = 4, title.position = 'top')) +
  facet_wrap(~part_barcode, ncol = 6, scales = 'free') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  labs(x='', y = 'Taxon relative abundance') +
  theme_light(base_size = 9) +
  theme(legend.key.size = unit(0.9, 'char'), legend.key.height = unit(0.9, 'char'), 
        legend.text.align = 0, legend.position = 'bottom',
        legend.margin = margin(-10, 0, 0, 0), legend.title.align = 0.5,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'italic', margin = margin(0.05, 0, 0.05, 0, 'cm'), color = 'black'), #size = 7,
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        strip.background = element_blank(), strip.text.x = element_blank(),
        panel.spacing = unit(0.2, 'cm'))
# stack_bar_swid






##################################################################################################
### Figure S4 - L. crispatus relative abundance histogram among WSW contacts and male contacts ###
##################################################################################################

# prepare data for plotting

part$sex_sex <- paste0(part$sex_1, part$sex_2)
# part should have 42 obs of 328 vars
plot_1 <- part[!is.na(part$Lactobacillus_crispatus_1), colnames(part) %in% c ('sting_id_1', 'sex_1', 'Lactobacillus_crispatus_1')]
plot_2 <- part[!is.na(part$Lactobacillus_crispatus_2), colnames(part) %in% c ('sting_id_2', 'sex_2', 'Lactobacillus_crispatus_2')]
# both should have 42 obs of 3 vars
colnames(plot_1) <- c('sex', 'sting_id', 'crisp')
colnames(plot_2) <- c('sex', 'sting_id', 'crisp')
plot <- rbind(plot_1, plot_2)
plot <- unique(plot)
# plot should have 74 obs of 3 vars
plot$crisp <- plot$crisp*100

# drop F contacts in heterosexual dyads 
drop <- c(part$sting_id_1[part$sex_1 == 'f' & part$sex_2 == 'm'], part$sting_id_2[part$sex_2 == 'f' & part$sex_1 == 'm'])
plot <- plot[plot$sting_id %!in% drop,]
# plot should have 42 obs of 3 vars

# plot
crisp_hist <- ggplot(plot, aes(x = crisp, fill = sex)) +
  geom_histogram(binwidth = 0.25) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex',
                    na.translate = F) +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab(expression(paste(italic('L. crispatus'), ' relative abundance'))) +
  ylab('Count') +
  theme_light(base_size = 10)
# crisp_hist






##############################################################
### Figure S6 - boxplot of processed non-host reads by sex ###
##############################################################

plot_reads <- ggplot(person) +
  geom_boxplot(aes(x = factor(sex, labels = c('F', 'M')), y = (10^log10_reads), 
                   fill = sex),
               show.legend = F, color = 'black', outlier.size = 0.75) +
  geom_jitter(aes(x = factor(sex, labels = c('F', 'M')), y = (10^log10_reads)),
              width = 0.2, height = 0.1, size = 0.75, 
              show.legend = F) +
  scale_fill_manual(values = c('#053061',
                               '#92c5de'),
                    labels = c('Female', 
                               'Male'),
                    'Sex') +
  xlab('') +
  ylab('Processed non-host reads') +
  scale_y_continuous(trans = 'log10') +
  theme_light(base_size = 10) 
# plot_reads






#######################################################################
### Analysis of comparable sites among female and male participants ###
#######################################################################

# summarize proportion of metagenome-MAG mappings with 0 comparable sites

sum(sites$sites_gt5x[sites$sex == 'f'] == 0)
sum(sites$sites_gt5x[sites$sex == 'f'] == 0) / sum(sites$sex == 'f')
# 1463779 (53.3%) metagenome-MAG mappings for female participants have 0 comparable sites

sum(sites$sites_gt5x[sites$sex == 'm'] == 0)
sum(sites$sites_gt5x[sites$sex == 'm'] == 0) / sum(sites$sex == 'm')
# 473513 (54.7%) metagenome-MAG mappings for male participants have 0 comparable sites



# summarize per-MAG comparable sites for mappings with >0 comparable sites
summary(sites$sites_gt5x[sites$sites_gt5x > 0 & sites$sex == 'f'])
summary(sites$sites_gt5x[sites$sites_gt5x > 0 & sites$sex == 'm'])

# test if per-MAG comparable sites differ between female and male participants
kruskal.test(sites$sites_gt5x, sites$sex)





