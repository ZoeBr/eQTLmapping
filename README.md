# eQTLmapping

Calculating molecular population statistics
```
salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=150G --job-name=copy --time=24:00:00 --partition=general --account=a_senv_ege srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l
module load anaconda3
conda activate pixy
pixy --stats pi dxy --vcf /scratch/user/s4480088/wgs_dataset_sf1.vcf.gz --populations /scratch/user/s4480088/Populations1.txt --window_size 10000 --n_cores 10 --chromosomes scaffold_1 --output_prefix 'wgs_test1' --output_folder '/scratch/user/s4480088'
```
Calculating logTPM from raw count data
```
exp2$length_kb <- (exp2$end - exp2$start) / 1000
> count_data <- exp2[, 5:(ncol(exp2)-1)]
> rpk <- count_data / exp2$length_kb
> tpm <- apply(rpk, 2, function(x) x / sum(x) * 1e6)
> tpm_df <- cbind(exp2[, 1:4], tpm)
> fwrite(tpm_df, "gene_exp_input_TPM.bed.gz", sep="\t", quote=FALSE, row.names=FALSE, compress = "gzip")
```
Tensorqtl code
For cis-eQTLs
```
#Load required modules
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post
#Load file paths
plink_prefix_path = '/scratch/user/s4480088/plink_prefix_path_eQTL2'
expression_bed = '/scratch/user/s4480088/gene_exp_input_TPM.bed.gz'
covariates_file = '/scratch/user/s4480088/Covariates_subset2.tsv'
prefix = 'TensorQTL_test4'

#run tensorqtl
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df.loc[phenotype_pos_df['chr'] == 'scaffold_1'], phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'scaffold_1'], covariates_df=covariates_df, window=25000, seed=123456)
cis_df = cis_df.reset_index()
cis_df.to_csv("/scratch/user/s4480088/cis_eQTL_results_scaff1_AG.tsv", sep='\t', index=False)
for f in files:
   scaffold = os.path.basename(f).split('.')[0]
   df = pd.read_csv(f, sep='\t')  # or sep=',' if CSV
   df['scaffold'] = scaffold
   all_dfs.append(df)

file_pattern = "cis_eQTL_results_scaff*_GR.tsv"
files = sorted(glob.glob(file_pattern))
all_dfs = []
for f in files:
# Extract scaffold name from filename, e.g., 'scaffold1.txt' → 'scaffold1'
scaffold = os.path.basename(f).split('.')[0]


for f in files:
  scaffold = os.path.basename(f).split('.')[0]
		df = pd.read_csv(f, sep='\t')  # or sep=',' if CSV
		df['scaffold'] = scaffold
		all_dfs.append(df)

combined_df = pd.concat(all_dfs, ignore_index=True)
combined_df.to_csv("combined_tensorqtl_output.tsv", sep='\t', index=False)
df = pd.read_csv("combined_tensorqtl_output.tsv", sep='\t')

filtered_df_AG = df_AG[
    (df['af'] >= 0.05) &
    (df['slope'].abs() <= 4) &
    (df['pval_beta'] <= 0.05)
]

filtered_df_GR = df_GR[
    (df_GR['af'] >= 0.05) &
    (df_GR['slope'].abs() <= 4) &
    (df_GR['pval_beta'] <= 0.05)
]

filtered_df = tensorqtl_results_df[
    (tensorqtl_results_df['af'] >= 0.05) &
    (tensorqtl_results_df['slope'].abs() <= 2) &
    (tensorqtl_results_df['pval_beta'] <= 0.05)
]

filtered_df2 = tensorqtl_results_df[
    (tensorqtl_results_df['af'] >= 0.05) &
    (tensorqtl_results_df['pval_beta'] <= 0.05)
]
```
For Trans eQTLs
```
#Run batch script - too much memory.
#!/bin/bash
#SBATCH --job-name=tensorqtl_trans
#SBATCH --output=tensorqtl_trans_%j.out
#SBATCH --error=tensorqtl_trans_%j.err
#SBATCH --time=24:00:00             
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            
#SBATCH --mem=64G                    
#SBATCH --partition=short            

# Load required modules
module load python/3.12              
source activate qtltools             

# Run TensorQTL
python3 -m tensorqtl \
  /scratch/user/s4480088/plink_prefix_path_eQTL \
  /scratch/user/s4480088/gene_exp_input_logTPM.bed.gz \
  TensoreQTL_trans_output \
  --covariates /scratch/user/s4480088/Covariates_subset.tsv \
  --mode trans \
  --batch_size 1000 \
  --maf_threshold 0.05 \
  --output_text
# Identify stringent associations
import pandas as pd
from statsmodels.stats.multitest import multipletests

# Load TensorQTL results
df = pd.read_csv("TensoreQTL_trans_AG.trans_qtl_pairs.txt.gz", sep="\t", compression="gzip")

# Compute FDR (Benjamini-Hochberg)
df['qval'] = multipletests(df['pval'], method='fdr_bh')[1]

stringent_df2 = df[
    (df['qval'] < 0.001) &
    (df['b'].abs() > 2) &
    (df['af'] > 0.2) &
	(df['pval'] < 1e-12)
]

stringent_df2.to_csv("TensoreQTL_trans_AG_significant.txt.gz", sep="\t", index=False, compression="gzip")

very_stringent_df = df[
    (df['qval'] < 1e-5) &
    (df['pval'] < 1e-15) &
    (df['b'].abs() > 3) &
    (df['af'] > 0.3)
]
```
Identification of overlapping cis and trans eQTLs
```
#Find overlapping cis and trans eQTLs

import pandas as pd

cis_df = pd.read_csv("filtered_tensorqtl_results.tsv", sep=",")  # or sep="\t" if it's tab-delimited
trans_df = pd.read_csv("trans_eqtl_stringent3.txv.gz", sep="\t", compression="gzip")

overlap_variants = pd.merge(cis_df, trans_df, on="variant_id", suffixes=("_cis", "_trans"))

overlap_variants.to_csv("overlapping_variants.tsv", sep="\t", index=False)
	
overlap_phenotypes = pd.merge(cis_df, trans_df, on="phenotype_id", suffixes=("_cis", "_trans"))

overlap_phenotypes.to_csv("overlapping_phenotypes.tsv", sep="\t", index=False)	

# Standardize effect sizes
cis_df["cis_effect_strength"] = cis_df["slope"] / cis_df["slope_se"]
trans_df["trans_effect_strength"] = trans_df["b"] / trans_df["b_se"]

# Get sets of phenotype IDs
cis_genes = set(cis_df["phenotype_id"])
trans_genes = set(trans_df["phenotype_id"])
both_genes = cis_genes & trans_genes
cis_only_genes = cis_genes - both_genes
trans_only_genes = trans_genes - both_genes

# 1. Cis-only
cis_only_df = cis_df[cis_df["phenotype_id"].isin(cis_only_genes)].copy()
cis_only_df["group"] = "cis only"
cis_only_df["effect_strength"] = cis_only_df["cis_effect_strength"]

# 2. Cis in overlapping genes
cis_overlap_df = cis_df[cis_df["phenotype_id"].isin(both_genes)].copy()
cis_overlap_df["group"] = "cis (overlapping)"
cis_overlap_df["effect_strength"] = cis_overlap_df["cis_effect_strength"]

# 3. Trans-only
trans_only_df = trans_df[trans_df["phenotype_id"].isin(trans_only_genes)].copy()
trans_only_df["group"] = "trans only"
trans_only_df["effect_strength"] = trans_only_df["trans_effect_strength"]

# 4. Trans in overlapping genes
trans_overlap_df = trans_df[trans_df["phenotype_id"].isin(both_genes)].copy()
trans_overlap_df["group"] = "trans (overlapping)"
trans_overlap_df["effect_strength"] = trans_overlap_df["trans_effect_strength"]

# Combine all
plot_df = pd.concat([
    cis_only_df[["phenotype_id", "effect_strength", "group"]],
    cis_overlap_df[["phenotype_id", "effect_strength", "group"]],
    trans_only_df[["phenotype_id", "effect_strength", "group"]],
    trans_overlap_df[["phenotype_id", "effect_strength", "group"]],
])

import seaborn as sns
import matplotlib.pyplot as plt

# Order groups as you specified
group_order = ["cis only", "cis (overlapping)", "trans only", "trans (overlapping)"]

plt.figure(figsize=(10, 6))
sns.boxplot(data=plot_df, x="group", y="effect_strength", order=group_order, palette="Set2")

plt.ylabel("Scaled Effect Size (Effect / SE)")
plt.xlabel("Gene Category")
plt.title("Comparison of Scaled Effect Sizes Across eQTL Categories")

plt.tight_layout()
plt.savefig("cis_trans_effect_strength_boxplot.png", dpi=300)
plt.show()

cis_only_df = 1081 interactions
cis_overlap_df = 75 interactions
trans_only_df = 13479 interactions
trans_overlap_df = 229 interactions

df = overlap_phenotypes.copy()

# Remove missing values just in case
df = df.dropna(subset=["slope", "slope_se", "b", "b_se"])

# Compute scaled effect sizes (effect / SE)
df["cis_scaled"] = df["slope"] / df["slope_se"]
df["trans_scaled"] = df["b"] / df["b_se"]

# Classify direction agreement
df["direction"] = df.apply(
    lambda row: "same" if row["cis_scaled"] * row["trans_scaled"] > 0 else "opposite",
    axis=1
)

# Plot
plt.figure(figsize=(7, 7))
sns.scatterplot(
    data=df,
    x="cis_scaled",
    y="trans_scaled",
    hue="direction",
    palette={"same": "green", "opposite": "red"},
    alpha=0.7
)

# Add reference lines at 0
plt.axhline(0, color="gray", linestyle="--", linewidth=1)
plt.axvline(0, color="gray", linestyle="--", linewidth=1)

# Labels and title
plt.xlabel("Scaled Cis Effect Size (slope / SE)")
plt.ylabel("Scaled Trans Effect Size (b / SE)")
plt.title("Scaled Cis vs Trans Effect Sizes for Overlapping Genes")
plt.legend(title="Direction Agreement")

# Save and show
plt.tight_layout()
plt.savefig("cis_vs_trans_scaled_scatter.png", dpi=300)

r, p = pearsonr(df["cis_scaled"], df["trans_scaled"])
>>> print(f"Pearson r = {r:.2f}, p = {p:.2g}")
Pearson r = -0.32, p = 7.9e-07

direction_counts = df["direction"].value_counts()
>>> print(direction_counts)
direction
opposite    217
same         12
```
GO term enrichment analysis
```
import goatools
>>> from goatools.obo_parser import GODag
>>> from goatools.go_enrichment import GOEnrichmentStudy
>>> obodag = GODag("go-basic.obo")
go-basic.obo: fmt(1.2) rel(2025-10-10) 42,666 Terms
>>> gene2go = {}

with open("gene2go.txt") as f:
    for line in f:
        gene, terms = line.strip().split('\t')
        gene2go[gene] = set(terms.split(','))
for g in gene2go:
    gene2go[g] = set(gene2go[g])
goea = GOEnrichmentStudy(population_genes, gene2go, obodag, methods=['fdr_bh'])
results = goea.run_study(study_genes)
goea.print_summary(results, pval=0.05)

results = goea.run_study(study_genes2)
```
Plotting pixy output
```
#In R
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
install.packages('zoo')
library(zoo)

inp<-read.table("wgs_test1_pi.txt",sep="\t",header=T)

chroms <- unique(inp$chromosome)
chrOrder <- sort(chroms)
inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)

if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p)
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Pi for pop", p))+
      labs(x="position of window start", y="Pi")+
      scale_colour_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("piplot_",p, ".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

#Just WAT1 region

if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p & 
                        chromosome == "scaffold_1" & 
                        window_pos_1 >= 35848322 & 
                        window_pos_2 <= 35891354)
    
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color = chromosome)) +
      geom_point() +
      facet_grid(. ~ chromosome) +
      labs(title = paste("Pi for pop", p),
           x = "Position of window start",
           y = "Pi") +
      scale_colour_manual(values = rep(c("black", "gray"), ceiling((length(chrOrder) / 2)))) +
      theme_classic() +
      theme(legend.position = "none")
    
    ggsave(paste("piplot_", p, ".png", sep = ""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

#whole chromosome line graph - this is good code

pixy_to_long <- function(pixy_files){
  pixy_df <- list()
  for(i in 1:length(pixy_files)){
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if(stat_file_type == "pi"){
      df <- read_delim(pixy_files[i], delim ="\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome, key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
      
    } else{
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome, key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_folder <- "wgs_test1"
pixy_files <- "wgs_test1_pi.txt"
#pixy_files <- list.files(pixy_folder, pattern = "\\.txt$", full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

pixy_labeller <- as_labeller(c(avg_pi = "pi"),
                               default = label_parsed)
pixy_df %>%
  filter(chromosome == "scaffold_1") %>%
  filter(statistic %in% c("avg_pi")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/100000) %>%
  ggplot(aes(x = chr_position, y = value, colour = statistic))+
  geom_line(size = 0.25)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistix = pixy_labeller,
                                 value = label_value))+
  xlab("Position on Chromosome (Mb)")+
  ylab("statistic value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_brewer(palette = "Set1")

#Just WAt1 region
pixy_df %>%
  filter(chromosome == "scaffold_1") %>%
  filter(statistic %in% c("avg_pi")) %>%
  filter(window_pos_1 >= 35600000, window_pos_2 <= 35961354) %>%  # 
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1e6) %>%  # Mb for x-axis
  ggplot(aes(x = chr_position, y = value, colour = statistic)) +
  geom_line(size = 0.25) +
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller, value = label_value)) +
  xlab("Position on Chromosome (Mb)") +
  ylab("Statistic Value") +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")


#Attempt to smooth
smoothed_df <- pixy_df %>%
  filter(chromosome == "1", statistic == "avg_pi") %>%  # adjust as needed
  arrange(chr_position) %>%  # make sure data is sorted
  mutate(
    rolling_value = rollmean(value, k = 10, fill = NA, align = "center"),  # window of 10
    chr_position = ((window_pos_1 + window_pos_2)/2)/100000
  )

smoothed_df <- pixy_df %>%
  filter(chromosome == "scaffold_1", statistic == "avg_pi") %>%  # adjust as needed
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/100000) %>%
  arrange(chr_position) %>%  # now this works!
  mutate(
    rolling_value = rollmean(value, k = 10, fill = NA, align = "center")
  )

ggplot(smoothed_df, aes(x = chr_position, y = rolling_value, colour = statistic)) +
  geom_line(size = 0.5) +
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x") +
  xlab("Position on Chromosome (Mb)") +
  ylab("Rolling Mean of Statistic") +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")
```
#Minor Allele Frequency calc
```
 plink2 --pfile plink_prefix_path_eQTL --freq --out minor_allele_freq --allow-extra-chr
 awk 'NR > 1 {maf=$5<0.5?$5:1-$5; print $1, $2, maf}' vcf_allele_freq.afreq > vcf_allele_freq_out.maf
```
