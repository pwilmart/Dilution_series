
# load the R libraries
library(tidyverse)
library(stringr)
library(edgeR) # used in Parts 2 and 3

# Set up to select the proper columns and label them
keep <- c(1, 2, 5, 6, 9, 10)
labels = c("a_25", "b_20", "c_15", "d_10", "e_5", "f_2.5")

# read in the PSM dataset, keep the 6 columns, and set labels
psm_start <- read_csv("psm_tmt.csv")
psm_start <- psm_start[keep]
colnames(psm_start) <- labels

# PAW Grouped peptides
peptide_start <- read_csv("grouped_peptide_summary_TMT_8.csv")
peptide_start <- peptide_start[keep]
colnames(peptide_start) <- labels

# PAW grouped proteins
protein_start <- read_csv("grouped_protein_summary_TMT_8.csv")
protein_start <- protein_start[keep]
colnames(protein_start) <- labels

# filter out rows with any missing data points and see how many rows remain
psms <- psm_start[apply(psm_start, 1, function(x) all(x > 0)), ] 
print("PSMs (before and after filtering zeros):")
format(c(nrow(psm_start), nrow(psms)), digits = 0, big.mark = ',')
                        
peptides <- peptide_start[apply(peptide_start, 1, function(x) all(x > 0)), ] 
print("Peptides (before and after filtering zeros):")
format(c(nrow(peptide_start), nrow(peptides)), digits = 0, big.mark = ',')
                                
proteins <- protein_start[apply(protein_start, 1, function(x) all(x > 0)), ] 
print("Proteins (before and after filtering zeros):")
format(c(nrow(protein_start), nrow(proteins)), digits = 0, big.mark = ',')

# check column sums for PSMs, peptides, and proteins
print("PSMs column sums:")
format(round(colSums(psms), digits = 0), big.mark = ",")

print("Peptide column sums:")
format(round(colSums(peptides), digits = 0), big.mark = ",")

print("Protein column sums:")
format(round(colSums(proteins), digits = 0), big.mark = ",")

show_densities <- function(df, xmin, xmax, title) {
    # makes density distributions of the log2 intensities
        # df - data frame with TMT intensities
        # xmin, xmax - x axis limits
        # title- plot title
    
    # get data in long form
    df_long <- gather(log2(df), key = "dilution", value = "intensity")

    # make the density plots colored by dilution
    ggplot(df_long, aes(x = intensity, fill = dilution)) +
        geom_density(alpha = 0.2) + 
        ggtitle(title) +
        xlim(xmin, xmax)
}

# compare PSM intensity (peak heights) distributions
show_densities(psms, 5, 25, "PSM data")

# compare peptide intensity (peak heights) distributions
show_densities(peptides, 5, 25, "Peptide data")

# compare protein intensity (peak heights) distributions
show_densities(proteins, 0, 30, "Protein data")

# helper function for scatter plots
show_scatter <- function(df, title) {
    # makes linear scale scatter plots (all series against reference)
        # df - tidy data frame with intensities and average reference vector
        # title - plot title
    
    # add the reference column (average) and put into long form for ggplot
    df$ref <- rowMeans(df)
    gg_df <- gather(df, key = dilution, value = intensity, a_25:f_2.5)
    
    # make the full range scatter plot
    full  <- ggplot(data = gg_df, aes(x = ref, y = intensity)) +
        geom_point(aes(color = dilution, shape = dilution), alpha = 0.5) +
        geom_smooth(aes(color = dilution), method = "lm") +
        ggtitle(title)

    # expanded scale
    expanded <- ggplot(data = gg_df, aes(x = ref, y = intensity)) +
        geom_point(aes(color = dilution, shape = dilution), alpha = 0.5) + 
        geom_smooth(aes(color = dilution), method = "lm") +
        coord_cartesian(xlim = c(0, 500000), ylim = c(0, 500000)) +
        ggtitle(str_c(title, " (expanded scale)"))
    
    print(full) 
    print(expanded)
}

# scater plots for PSMs
show_scatter(psms, "raw PSMs")

# scater plots for peptides
show_scatter(peptides, "Peptides")

# scater plots for proteins
show_scatter(proteins, "Proteins")

# MA style plot for PSMs
show_MA_plot <- function(df, title) {
    # transforms the data and makes MA plots (all versus reference)
        # df - data frame with intensities and reference vector
        # title - plot title
    
    # take ratios relative to reference then log transform
    df$ref <- rowMeans(df)
    log_df <- log2(df[1:6] / df$ref) # these are log 2
    log_df$ref <- log10(df$ref) # this stays in log 10
    
    # we need data in tidy long form
    gg_log_df <- gather(log_df, key = dilution, value = log_ratios, a_25:f_2.5)
    
    # compute the expected horizontal lines
    ratios <- log2(colMeans(df)[1:6] / colMeans(df)[7])
    
    # make the plot
    ggplot(data = gg_log_df, aes(x = ref, y = log_ratios)) +
        geom_point(aes(color = dilution, shape = dilution), alpha = 0.3) +
        geom_hline(yintercept = ratios) +
        xlab("log ref intensity") + 
        ggtitle(title)
}
    
# make some MA plots for PSMs
show_MA_plot(psms, "PSMs (MA plot)")

# make some MA plots for peptides
show_MA_plot(peptides, "Peptides (MA plot)")

# make some MA plots for proteins
show_MA_plot(proteins, "Proteins (MA plot)")

# how close were the measured dilution factors?
tot_int <- sum(colSums(proteins))
unit <- tot_int / sum(c(25, 20, 15, 10, 5, 2.5))
round(colSums(proteins) / unit, 2)

# loading the edgeR library was done at the top of the notebook
# load data into DGEList objects (counts are intensities, group are sample labels)
# the original data frames have added columns so we need just the first 6 cols

# PSMs first
y_psms  <- DGEList(counts = psms[1:6], group = factor(labels))
y_psms$samples

# peptides next
y_peptides  <- DGEList(counts = peptides[1:6], group = factor(labels))
y_peptides$samples

# finally, proteins
y_proteins  <- DGEList(counts = proteins[1:6], group = factor(labels))
y_proteins$samples

# these are the original column sums (they should match the library sizes)
round(colSums(proteins[1:6]), 0)

# dividing by library sizes makes each sample sum to 1.0
scaled_proteins <- sweep(proteins[1:6], 2, y_proteins$samples$lib.size, FUN = "/")
round(head(scaled_proteins), 8)
colSums(scaled_proteins)

# make our own cpm scale
cpm_proteins <- scaled_proteins * 1000000
round(head(cpm_proteins), 0)
colSums(cpm_proteins)

# use the cpm function in edgeR
round(head(cpm(y_proteins)), 0)
colSums(cpm(y_proteins))

# compute TMM factors - they get added to $samples
y_proteins_tmm <- calcNormFactors(y_proteins)
# round(y_proteins_tmm$samples$norm.factors, 6)
y_proteins_tmm$samples

# see what column sums we get after a cpm function call now that there are TMM factors
round(colSums(cpm(y_proteins_tmm)), 0)
round(mean(colSums(cpm(y_proteins_tmm))))

# adjust the library sizes by the TMM factors to get the real library size factors
real_factors <- y_proteins_tmm$samples$lib.size * y_proteins_tmm$samples$norm.factors

# apply the factors to get the TMM normalized data table (fractional scale)
proteins_tmm_norm <- sweep(proteins[1:6], 2, real_factors, FUN = "/")

# put TMM normalized data on a one million scale
by_hand_cpm <- proteins_tmm_norm * 1000000

round(head(by_hand_cpm), 0)
round(colSums(by_hand_cpm), 0)
round(mean(colSums(by_hand_cpm)))

# make multi-dimensional scaling plots to see if series are similar
plotMDS(calcNormFactors(y_psms), main = "PSMs")
plotMDS(calcNormFactors(y_peptides), main = "Peptides")
plotMDS(calcNormFactors(y_proteins), main = "Proteins")

# define a function for normalizing data and computing CVs
get_cv <- function(DGE) {
    # normalizes edgeR data, makes cpm scale, computes CV distribution
        # DGE - edgeR DGEList object
        # returns a CV distribution vector
    
    # get a cmp scale table of normalized data
    CPM <- as.data.frame(cpm(calcNormFactors(DGE)))

    # compute CVs for first 4 channels (drop the 5 and 2.5 dilutions)
    cv <- 100 * apply(CPM[1:4], 1, sd) / rowMeans(CPM[1:4])
}

# compute the CVs and look at the distribution summary numbers
psm_cvs <- get_cv(y_psms)
peptide_cvs <- get_cv(y_peptides)
protein_cvs <- get_cv(y_proteins)

print("PSM CVs:")
summary(psm_cvs)

print("Peptide CVs:")
summary(get_cv(y_peptides))

print("Protein CVs:")
summary(get_cv(y_proteins))

BP <- function(cv, ymax, title) {
    # returns a labeled boxplot
        # cv - vector of CV values
        # ymax - upper limit for plot
        # title - main plot title
    
    # make the box plot and label the median
    bp  <- boxplot(cv, ylim = c(0, ymax), notch = TRUE, main = title)
    label <- text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
    list(bp, label)
}

# compare some box plots of the CV distributions
par(mfrow = c(1, 3))
plot.label <- BP(get_cv(y_psms), 100, "PSMs")
plot.label <- BP(get_cv(y_peptides), 100, "Peptides")
plot.label <- BP(get_cv(y_proteins), 100, "Proteins")
par(mfrow = c(1, 1))

# log the R session (we should always end notebooks with this)
sessionInfo()
