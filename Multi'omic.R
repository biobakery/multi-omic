# Omnibus tests on mgx and mtx data from the HMP2. Baseline only data.

dir.create("./R_Multi-omics") # Create a new directory
dir.create("./R_Multi-omics/Data") # Create a new directory
setwd("./R_Multi-omics") # Change the current working directory 
getwd()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
install.packages("circlize")

# Load the packages needed
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
#####
# MGX species
#####
# You can either download the file from bitbucket page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/taxonomic_profiles_pcl_week0.csv", "./Data/taxonomic_profiles_pcl_week0.csv") # Download the mgx species data and put it into the data directory

# Read the taxonomic data into R environment
tax = read.csv(file = "./Data/taxonomic_profiles_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(tax) # Inspect the tax
dim(tax) # dimensions of the tax
str(tax) # structure of the tax
names(tax) # column names of tax
row.names(tax) # row.names of tax


# Prepare the data
# Extract the metadata
metadata = data.frame(tax[1:5])
metadata[1:5,] # check the output
str(metadata) #structure of metadata
is.na(metadata) # Check for NAs that will mess with the PERMANOVAS: Age has some
count(is.na(metadata$consent_age)) # Check for how many there are: 6/96
# If this was a discrete variable we could just classify the NAs as Unknown and keep them in the model, 
# but since Age is a continuous variable typically we would either remove those from the data or impute the median. 
# In this case let's impute the median in order to keep samples. 
unique(metadata$consent_age)
metadata$consent_age[is.na(metadata$consent_age)] = median(metadata$consent_age, na.rm = T)
unique(metadata$consent_age) # Check the output: good to go

# Extract species data and transpose the df
species = data.frame(t(tax[6:ncol(tax)]))
str(species) # everything is numeric, good to go
row.names(species)
# species[] = as.data.frame(sapply(species, function(x) as.numeric(as.character(x))))
species[1:8,1:4] # check the output
# subset to species only
# which don't have "t__"
tmp.ind = grep("\\|t__", rownames(species), invert = T) # grep the rows that do not include strain stratifications
tmp.ind # check the output
tmp = species[tmp.ind,]  # Create a new dataframe with only those row numbers
tmp.ind = grep("\\|s__", rownames(tmp)) # grep the rows that only include down to species stratifications
tmp.ind # check the output
species = tmp[tmp.ind,] # Create a new dataframe with only those row numbers
rm(tmp,tmp.ind) # remove temp files to clear up space
row.names(species) # Check the output to make sure that we only have species level stratifications

# trim species names
rownames(species) = gsub(".*\\|", "", rownames(species))
row.names(species) # Check the output, looks great
colSums(species) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1

# filter for beta div (we will keep species as is for alpha diversity)
dim(species)
dim(species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]) 
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]
#Let's transpose it for easier use downstream
species_filt = data.frame(t(species_filt), check.names = F)
species = data.frame(t(species), check.names = F)

#####
# MGX pathways
#####
# You can either download the file from bitbuckert page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/dna_pathabundance_relab_pcl_week0.csv", "./Data/dna_pathabundance_relab_pcl_week0.csv")

# Read the dna pathway data into R environment
dna_path = read.csv(file = "./Data/dna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(dna_path) # Inspect the data
dim(dna_path) # dimensions of the data
str(dna_path) # structure of the data
names(dna_path) # column names of data
row.names(dna_path) # row.names of data
# Remove metadata and keep only pathways and transpose the data
dna_path = data.frame(t(dna_path[6:ncol(dna_path)]))
str(dna_path) # everything is numeric, good to go
row.names(dna_path)

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(dna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
dna_path_unstratified = dna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(dna_path_unstratified) # check the output: looks great
colSums(dna_path_unstratified) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1


# filter for beta div
dim(dna_path_unstratified)
dim(dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]) 
dna_path_unstratified_filt = dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(dna_path_unstratified)), ]
#Let's transpose it for easier use downstream
dna_path_unstratified_filt = data.frame(t(dna_path_unstratified_filt), check.names = F)
dna_path_unstratified = data.frame(t(dna_path_unstratified), check.names = F)

#####
# MTX pathways
#####
# You can either download the file from bitbuckert page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/rna_pathabundance_relab_pcl_week0.csv", "./Data/rna_pathabundance_relab_pcl_week0.csv")

# Read the rna pathway data into R environment
rna_path = read.csv(file = "./Data/rna_pathabundance_relab_pcl_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(rna_path) # Inspect the data
dim(rna_path) # dimensions of the data
str(rna_path) # structure of the data
names(rna_path) # column names of data
row.names(rna_path) # row.names of data
# Remove metadata and keep only pathways and transpose the data
rna_path = data.frame(t(rna_path[6:ncol(rna_path)]))
str(rna_path) # everything is numeric, good to go
row.names(rna_path)
# minimize the metadata to just the samples available in these data
dim(metadata)
list = names(rna_path) # make a list of sample ids to subset on
list # check the output
metadata_rna = subset(metadata, row.names(metadata) %in% list)
dim(metadata_rna)
metadata_rna # check the output

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(rna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
rna_path_unstratified = rna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(rna_path_unstratified) # check the output: looks great
colSums(rna_path_unstratified) # Check the sample sums to make sure they are in proportion format (0-1) and are all ~1


# filter for beta div
dim(rna_path_unstratified)
dim(rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]) 
rna_path_unstratified_filt = rna_path_unstratified[apply(rna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(rna_path_unstratified)), ]
#Let's transpose it for easier use downstream
rna_path_unstratified_filt = data.frame(t(rna_path_unstratified_filt), check.names = F)
rna_path_unstratified = data.frame(t(rna_path_unstratified), check.names = F)



#####
# RNA/DNA pathway ratios
#####
# You can either download the file from bitbuckert page or do it with the below code:
download.file("https://bitbucket.org/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/downloads/rna_dna_path_relative_expression_week0.csv", "./Data/rna_dna_path_relative_expression_week0.csv")

# Read the rna_dna pathway data into R environment
rna_dna_path = read.csv(file = "./Data/rna_dna_path_relative_expression_week0.csv", header = T, row.names = 1, check.names = FALSE)

head(rna_dna_path) # Inspect the data
dim(rna_dna_path) # dimensions of the data
str(rna_dna_path) # structure of the data
names(rna_dna_path) # column names of data
row.names(rna_dna_path) # row.names of data
# Transpose the data
rna_dna_path = data.frame(t(rna_dna_path))
str(rna_dna_path) # everything is numeric, good to go
row.names(rna_dna_path)
# minimize the metadata to just the samples available in these data
dim(metadata)
list = names(rna_dna_path) # make a list of sample ids to subset on
list # check the output
metadata_rna_dna = subset(metadata, row.names(metadata) %in% list)
dim(metadata_rna_dna)
metadata_rna_dna # check the output: we only have one race now, so let's get rid of that column
metadata_rna_dna$race = NULL

# Remove species stratifications
tmp.ind = grep("\\|.*", rownames(rna_dna_path), invert = T) # grep the rows that do not include species stratifications 
tmp.ind # check the output
rna_dna_path_unstratified = rna_dna_path[tmp.ind,] # Create a new dataframe with only those unstratified rows
rm(tmp.ind) # Remove tmp.ind to clear space
row.names(rna_dna_path_unstratified) # check the output: looks great


# filter for beta div
#Only keep RNA/DNA pathways that passed filtering for DNA pathways
#Create a list of the pathway names (col names) for subsetting the data
list = names(dna_path_unstratified_filt)
list
#Check the dimensions to make sure it matches with the DNA numbers before subsetting
dim(dna_path_unstratified_filt)
dim(subset(rna_dna_path_unstratified, row.names(rna_dna_path_unstratified) %in% list))
#subset
rna_dna_path_unstratified_filt = rna_dna_path_unstratified[list,]
dim(rna_dna_path_unstratified_filt)
head(rna_dna_path_unstratified_filt)

#Let's transpose the dataframe for easier use downstream
rna_dna_path_unstratified_filt = data.frame(t(rna_dna_path_unstratified_filt), check.names = F)

# log transform the RNA/DNA ratio
rna_dna_path_unstratified_filt_log = log2(rna_dna_path_unstratified_filt + 1)
head(rna_dna_path_unstratified_filt_log)

#####
#If you have Jeremy's code still loaded into your R console you can start here!
#####
#Compare the taxonomy to the DNA Pathways
#Use Spearman correlations to compare the taxonomy to the dna pathway relative abundance
m.corr.tax.path = cor(species_filt, dna_path_unstratified_filt, method="spearman")

heatmap.obj.tax.path = Heatmap(m.corr.tax.path,
                         name = "Comparison All Taxa to All Pathways",
                         col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         show_column_names=TRUE,
                         column_names_side = "bottom",
                         show_row_names=TRUE,
                         row_names_side = "left",
                         column_names_gp = gpar(fontsize = 4),
                         row_names_gp = gpar(fontsize = 4),
                         #row_title="Top 25 Species",
                         #row_title_side = "right",
                         column_title="DNA pathways",
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 4),
                         clustering_distance_rows = "pearson",
                         clustering_method_rows = "average",
                         clustering_distance_columns = "pearson",
                         clustering_method_columns = "average",
                         heatmap_legend_param = list(at = c(-1,0,1),
                                                     color_bar = "continuous",
                                                     legend_direction="horizontal",
                                                     labels_gp = gpar(fontsize = 12),
                                                     legend_width = unit(10, "cm"),
                                                     title_position = "topcenter",
                                                     title = "Spearman correlation"),
                         show_heatmap_legend = TRUE,
                         width = unit(5, "cm")
)
draw(heatmap.obj.tax.path,
     heatmap_legend_side = "bottom",
     row_title="Top 25 species",
     row_title_side="right",
     column_title = "Top 50 DNA pathways",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 14),
     row_title_gp = gpar(fontsize = 14)
)

#Let's take a look at E. coli
ecoli = subset(m.corr.tax.path, row.names(m.corr.tax.path) == "s__Escherichia_coli")
ecoli

#Use spearman correlation only on the most abundant taxa and pathways
top = names(sort(colMeans(species_filt), decreasing = TRUE))[1:50]
species_top = species_filt[,names(species_filt) %in% top]

#Use spearman correlation only on the most abundant pathways
top_path = names(sort(colMeans(dna_path_unstratified), decreasing = TRUE))[1:50]
dna_path_top = dna_path_unstratified[,names(dna_path_unstratified) %in% top_path]

m.corr.tax.path.top = cor(species_top, dna_path_top, method = "spearman")

heatmap.obj.tax.path.top = Heatmap(m.corr.tax.path.top,
                               name = "Comparison Top Taxa to Top Pathways",
                               col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                               show_row_dend = FALSE,
                               show_column_dend = FALSE,
                               show_column_names=TRUE,
                               column_names_side = "bottom",
                               show_row_names=TRUE,
                               row_names_side = "left",
                               column_names_gp = gpar(fontsize = 4),
                               row_names_gp = gpar(fontsize = 8),
                               #row_title="Top 25 Species",
                               #row_title_side = "right",
                               column_title="DNA pathways",
                               column_title_side = "top",
                               column_title_gp = gpar(fontsize = 10),
                               clustering_distance_rows = "pearson",
                               clustering_method_rows = "average",
                               clustering_distance_columns = "pearson",
                               clustering_method_columns = "average",
                               heatmap_legend_param = list(at = c(-1,0,1),
                                                           color_bar = "continuous",
                                                           legend_direction="horizontal",
                                                           labels_gp = gpar(fontsize = 12),
                                                           legend_width = unit(4, "cm"),
                                                           title_position = "topcenter",
                                                           title = "Spearman correlation"),
                               show_heatmap_legend = TRUE,
                               width = unit(10, "cm")
)
draw(heatmap.obj.tax.path.top,
     heatmap_legend_side = "bottom",
     row_title="Top 50 species",
     row_title_side="right",
     column_title = "Top 50 DNA pathways",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 14),
     row_title_gp = gpar(fontsize = 14)
)

#Here the Abundance of s__Faecalibacterium_prausnitzii appears to be very correlated to the pathway BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched amino acid biosynthesis
#Lets pull them out and graph the abundance of both
#Create an object called scatter with those two abundance profiles to make construction of the graph easier
scatter = data.frame(dna_path_unstratified$`BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched amino acid biosynthesis`, species_filt$s__Faecalibacterium_prausnitzii)
row.names(scatter) = row.names(metadata)
names(scatter) = c("Pathway", "Taxa")

p = ggplot(scatter, aes(x= Taxa, y= Pathway)) + geom_point(size=4) + theme_bw()
p

# What do you think of this plot?


#Final data reducing stratagy- increase the filters that we choose to identify the abundance changes in those bugs
#that are common in the dataset
species = data.frame(t(species))
dim(species)
dim(species[apply(species, 1, function(x) sum(x > 0.0001) > 0.25 * ncol(species)), ]) 
species_filt2 = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.25 * ncol(species)), ]
#Let's transpose it for easier use downstream
species_filt2 = data.frame(t(species_filt2), check.names = F)



#Apply a filter level of 80% prevalence (previous reported in the literature to represent "core" pathways)
dna_path_unstratified = t(dna_path_unstratified)
dim(dna_path_unstratified)
dim(dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.8 * ncol(dna_path_unstratified)), ]) 
dna_path_core = dna_path_unstratified[apply(dna_path_unstratified, 1, function(x) sum(x > 0.0001) > 0.8 * ncol(dna_path_unstratified)), ]
#Let's transpose it for easier use downstream
dna_path_core = data.frame(t(dna_path_core), check.names = F)


m.corr.tax.path.core = cor(species_filt2, dna_path_core, method = "spearman")

heatmap.obj.tax.path.top = Heatmap(m.corr.tax.path.top,
                                   name = "Comparison of Prevalent Taxa to Core Pathways",
                                   col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                                   show_row_dend = FALSE,
                                   show_column_dend = FALSE,
                                   show_column_names=TRUE,
                                   column_names_side = "bottom",
                                   show_row_names=TRUE,
                                   row_names_side = "left",
                                   column_names_gp = gpar(fontsize = 4),
                                   row_names_gp = gpar(fontsize = 8),
                                   #row_title="Top 25 Species",
                                   #row_title_side = "right",
                                   column_title="DNA pathways",
                                   column_title_side = "top",
                                   column_title_gp = gpar(fontsize = 10),
                                   clustering_distance_rows = "pearson",
                                   clustering_method_rows = "average",
                                   clustering_distance_columns = "pearson",
                                   clustering_method_columns = "average",
                                   heatmap_legend_param = list(at = c(-1,0,1),
                                                               color_bar = "continuous",
                                                               legend_direction="horizontal",
                                                               labels_gp = gpar(fontsize = 12),
                                                               legend_width = unit(4, "cm"),
                                                               title_position = "topcenter",
                                                               title = "Spearman correlation"),
                                   show_heatmap_legend = TRUE,
                                   width = unit(10, "cm")
)
draw(heatmap.obj.tax.path.top,
     heatmap_legend_side = "bottom",
     row_title="Prevalent species",
     row_title_side="right",
     column_title = "Core DNA pathways",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 14),
     row_title_gp = gpar(fontsize = 14)
)

#Do you see anything else you would like to pull out here and look at the abundance plots?

#One Final look at the correlations: After this you can play with these codes and the additional omic datasets loaded into this
#R-session. Jeremy wrote codes for Species, DNA Pathways, RNA Pathways and RNA/DNA Ratios
diagnosis = metadata$diagnosis
names(diagnosis) = rownames(metadata)

m.corr.nonIBD = cor(species_top[diagnosis=="nonIBD",],dna_path_top[diagnosis=="nonIBD",],method="spearman")
m.corr.cd = cor(species_top[diagnosis=="CD",],dna_path_top[diagnosis=="CD",],method="spearman")
m.corr.uc = cor(species_top[diagnosis=="UC",],dna_path_top[diagnosis=="UC",],method="spearman")

pdf("diagnosis_spearman_heatmap.pdf", height = 14, width = 21)
heatmap.obj.nonIBD = Heatmap(m.corr.nonIBD,
                         name = "nonIBD",
                         col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         show_column_names=TRUE,
                         column_names_side = "bottom",
                         show_row_names=TRUE,
                         row_names_side = "left",
                         column_names_gp = gpar(fontsize = 4),
                         row_names_gp = gpar(fontsize = 8),
                         #row_title="Top 25 Species",
                         #row_title_side = "right",
                         column_title="DNA pathways",
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 12),
                         clustering_distance_rows = "pearson",
                         clustering_method_rows = "average",
                         clustering_distance_columns = "pearson",
                         clustering_method_columns = "average",
                         heatmap_legend_param = list(at = c(-1,0,1),
                                                     color_bar = "continuous",
                                                     legend_direction="horizontal",
                                                     labels_gp = gpar(fontsize = 12),
                                                     legend_width = unit(10, "cm"),
                                                     title_position = "topcenter",
                                                     title = "Spearman correlation"),
                         show_heatmap_legend = TRUE,
                         width = unit(10, "cm")
)
heatmap.obj.cd = Heatmap(m.corr.cd,
                         name = "CD",
                         col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         show_column_names=TRUE,
                         column_names_side = "bottom",
                         show_row_names=FALSE,
                         #row_names_side = "left",
                         column_names_gp = gpar(fontsize = 4),
                         #row_names_gp = gpar(fontsize = 6),
                         row_title="Top 25 Species",
                         row_title_side = "right",
                         column_title="DNA pathways (CD)",
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 12),
                         row_order = row_order(heatmap.obj.nonIBD),
                         column_order = column_order(heatmap.obj.nonIBD),
                         width = unit(10, "cm"),
                         show_heatmap_legend = FALSE
)
heatmap.obj.uc = Heatmap(m.corr.uc,
                         name = "UC",
                         col=colorRamp2(c(-1,0,1),c("blue","white","red")),
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         show_column_names=TRUE,
                         column_names_side = "bottom",
                         show_row_names=FALSE,
                         #row_names_side = "left",
                         column_names_gp = gpar(fontsize = 4),
                         #row_names_gp = gpar(fontsize = 6),
                         row_title="Top 25 Species",
                         row_title_side = "right",
                         column_title="DNA pathways (UC)",
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 12),
                         row_order = row_order(heatmap.obj.nonIBD),
                         column_order = column_order(heatmap.obj.nonIBD),
                         width = unit(10, "cm"),
                         show_heatmap_legend = FALSE
)

draw(heatmap.obj.nonIBD + heatmap.obj.cd + heatmap.obj.uc,
     heatmap_legend_side = "bottom",
     row_title="Top 25 species",
     row_title_side="right",
     column_title = "Top 50 DNA pathways",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 22),
     row_title_gp = gpar(fontsize = 22)
)
dev.off()

#Try this kind of correlational analysis with the other datasets!
