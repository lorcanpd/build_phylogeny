#All analysis

#Meta data file
setwd ("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering")

All_meta_data = read.table("/home/jovyan/md39/fetal_mapping_md39/Meta.Data_Update_6.txt", header = TRUE, sep = "\t")
Blood_coverage = read.table("/home/jovyan/Cancer_Pipeline_Reports_SampleDetails (9).txt", header = TRUE, sep = "\t")
Blood_meta <- read.table("/lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/00_Supplementary_data/Blood_Only2.txt", header = TRUE, sep = "\t")
View(Blood_meta)
View(updated_meta_data_5)
View(Blood_coverage)

#Projects:
#2855
#3137
#3186
#3152
#3187 - blood

# ```{bash}
# #In command line
# cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering
#
# #Add in caveman files and the bam files into the analysis folder in your directory.

#Tim sym link code
# cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering
# project_number=3187
# #CREATE BAM SYMLINKS FOR CGPVAF
# for bam in $(find /nfs/cancer_ref01/nst_links/live/$project_number -iname '*.sample.dupmarked.bam'); do
# ln -s $bam
# done
#
# for bai in $(find /nfs/cancer_ref01/nst_links/live/$project_number -iname '*sample.dupmarked.bam.bai'); do
# ln -s $bai
# done

library(Rsamtools)
library(VariantAnnotation)


# Filter unidirectionally supported variants
bamfiles <- list.files(pattern = "\\.bam$", full.names = FALSE)
sample_names <- sapply(files, function(file) {
    strsplit(file, "\\.", fixed = FALSE)[[1]][1]
})

# create an empty data frame to store the details of the excluded mutations
excluded <- data.frame(
    sample = character(),
    Chrom = character(),
    Pos = numeric(),
    Ref = character(),
    Alt = character(),
    stringsAsFactors = FALSE
)

# loop through each sample
for (sample in sample_names) {
    # read in the vcf for the current sample
    # combine sample name with vcf extension
    vcf_file <- paste0(sample, ".caveman_c.flag.vcf.gz")
    vcf <- readVcf(vcf_file, "hg38")

    bam_file <- paste0(sample, ".sample.dupmarked.bam")


    for (i in seq_along(vcf)) {
        variant <- vcf[i]
        region <- GRanges(
            seqnames(rowRanges(variant)), ranges(rowRanges(variant))
        )
        param <- ScanBamParam(
            which = region, flag = scanBamFlag(isUnmappedQuery = FALSE)
        )
        bam <- scanBam(bam_file, param = param)

        if (length(bam[[1]]$flag) > 0) {
            reads <- bam[[1]]
            flags <- as.integer(reads$flag)
            isReverse <- bitwAnd(flags, 0x10) != 0

            ref_allele <- ref(variant)
            alt_allele <- alt(variant)
            pos <- start(region)

            # get number of forward reads supporting the variant
            forward <- sum(!isReverse & reads$pos == pos & reads$seq == alt_allele)
            # get number of reverse reads supporting the variant
            reverse <- sum(isReverse & reads$pos == pos & reads$seq == alt_allele)
            # get number of reads supporting the reference allele
            ref <- sum(reads$pos == pos & reads$seq == ref_allele)

            # if 95% of reads supporting the variant are in the same direction
            # exclude the variant
            if ((forward / (forward + reverse) > 0.95) |
                (reverse / (forward + reverse) > 0.95)) {
                excluded <- rbind(
                    excluded,
                    data.frame(
                        sample = sample,
                        Chrom = seqnames(region),
                        Pos = start(region),
                        Ref = ref_allele,
                        Alt = alt_allele,
                        stringsAsFactors = FALSE
                    )
                )
                vcf[i] <- NULL
            }



            # Set read quality threshold for reads supporting the variant

        # DEEP SNV -
        }

    }


}




##############
#Tim script to convert VCF to Bed file
#Then run CGPVAF
#Tim script to bypass LCM filter to generate all_snvs.bed file

#Change all_muts*.bed with a number so it doesn't overwrite the previous file. 5 files generated in CGPvaf folder - /lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/02_DNA_processing/CGPvaf
touch all_muts.bed
for sample in $(grep -o 'PD[0-9]\{5\}[a-z]\{0,2\}' /lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/00_Supplementary_data/Blood_Only2.txt
);
do
VCF=/nfs/cancer_ref01/nst_links/live/$project_number/$sample/$sample.caveman_c.flag.vcf.gz
zcat $VCF | grep -v "^#" | grep "PASS" | awk -F$'\t' '{split($0,a,"\t"); n=split(a[8],b,";"); for(i=1; i<=n; i++){if(b[i]~/CLPM/){split(b[i],c,"=")}else if(b[i]~/ASMD/){split(b[i],d,"=")}}; if(c[2] == 0 && d[2] >= 140){print $0}}' | cut -f1,2,4,5 >> all_muts.bed
done

cat all_muts*.bed | sort -u > all_snvs_sorted_unique.bed

#This works on the command line
module load vafcorrect
cgpVaf.pl -d /lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/02_DNA_processing/Variant_calls/Subs/01_caveman -o ./ -a snp -mq 30 -bq 25 -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz -bo 1 -b /lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/02_DNA_processing/Variant_calls/Subs/01_caveman/all_snvs_sorted_unique.bed -nn PDv38is_wgs -tn PD53943aa


cgpVaf.pl -d /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering -o ./ -a snp -mq 30 -bq 25 -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz -bo 1 -b /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/all_snvs_sorted_unique.bed -nn PDv38is_wgs -tn PD53943aa



#Run CGPVAF (input is BAM files and BED files). Runs all chromosomes so no need to concatonate.
#File is cgpvaf_md_persample.sh which is in Tim's script folder

cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/

while read sample; do
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R 'select[mem>3000] rusage[mem=3000]' -M3000 -J "$sample" ./CGPVAF_MD_Final_CH38.sh snp "$sample"
done < /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Blood_Samples.txt


#Concatenate into one sample
while read sample; do
sed '/^##/d' PDv38is*${sample}*.tsv > ${sample}.snp.tsv
done < /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Blood_Samples.txt


options(stringsAsFactors = F)
samples=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Blood_Samples.txt")[,1]
View(samples)
NR=c()
NV=c()
for(sample in samples){
  data = read.table(paste0(sample,".snp.tsv"), comment.char="",header=T)
  NR = cbind(NR,data[,grepl("DEP",colnames(data))&colnames(data)!="PDv38is_wgs_DEP"])
  NV = cbind(NV,data[,grepl("MTR",colnames(data))&colnames(data)!="PDv38is_wgs_MTR"])
}
Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
rownames(NV)=rownames(NR)=Muts
colnames(NR)=colnames(NV)=samples

write.table(NV,"NV_input.txt")
write.table(NR,"NR_input.txt")



#Make sure you have NV_input, NR_input and build.phylogeny.R script
#Submit job to farm
cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R 'select[mem>30000] rusage[mem=30000]' -M30000 -J Treebuilding /software/R-4.1.3/bin/Rscript build_phylogeny.R -r NR_input.txt -v NV_input.txt --mpboot_path /lustre/scratch126/casm/team274sb/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/ --plot_spectra F --max_muts_plot 100000 -m T

cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_False
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R 'select[mem>30000] rusage[mem=30000]' -M30000 -J Treebuilding /software/R-4.1.3/bin/Rscript build_phylogeny.R -r NR_input.txt -v NV_input.txt --mpboot_path /lustre/scratch126/casm/team274sb/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/ --plot_spectra F --max_muts_plot 100000 -m F

# **** check memory **** (not submitted jobs below)
#Plot Spectra Changes - True / Switched on
cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Plot_Spectra_T2
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R 'select[mem>30000] rusage[mem=30000]' -M30000 -J Treebuilding /software/R-4.1.3/bin/Rscript build_phylogeny.R -r NR_input.txt -v NV_input.txt --mpboot_path /lustre/scratch126/casm/team274sb/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/ --plot_spectra T --max_muts_plot 100000 -m T

cd /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_False/Plot_Spectra_T1
bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R 'select[mem>30000] rusage[mem=30000]' -M30000 -J Treebuilding /software/R-4.1.3/bin/Rscript build_phylogeny.R -r NR_input.txt -v NV_input.txt --mpboot_path /lustre/scratch126/casm/team274sb/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/ --plot_spectra T --max_muts_plot 100000 -m F



***
  #If you are looking for a software or folder in the farm, you can use this function
find -name <enter name>
  find -name mpboot-sse-1.1.0-Linux
***

#Mixture modelling True
Variant_depth_MMT <- read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NV_filtered_all.txt", header = TRUE)
Total_depth_MMT <- read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NR_filtered_all.txt", header = TRUE)

VAF_MMT <- Variant_depth_MMT / Total_depth_MMT
VAF_MMF <- Variant_depth_MMF / Total_depth_MMF

setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/")
getwd()
write.table(VAF_MMT, file = "VAF_MMT.txt", sep = "\t", quote = FALSE, col.names = NA)

setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_False/")
getwd()
write.table(VAF_MMF, file = "VAF_MMF.txt", sep = "\t", quote = FALSE, col.names = NA)

View(VAF_MMF)
View(VAF_MMT)

#2 or more reads for every variant -
Mutationburdenvariantdepth <- data.frame(colSums(Variant_depth_MMT > 2))
View(Mutationburdenvariantdepth)
colSums(Mutationburdenvariantdepth) #32016 variants

# Adding row names as a new column
Mutationburdenvariantdepth$SampleName <- rownames(Mutationburdenvariantdepth)
colnames(Mutationburdenvariantdepth) <- c("MutationBurden", "SampleName")
Mutationburdenvariantdepth$SampleName <- as.character(Mutationburdenvariantdepth$SampleName)

# Sort the data frame by MutationBurden in descending order and create an ordered factor
Mutationburdenvariantdepth <- Mutationburdenvariantdepth %>%
  arrange(desc(MutationBurden)) %>%
  mutate(SampleName = factor(SampleName, levels = SampleName))

# Create the scatter plot with a dotted red line at y = 41.9
ggplot(Mutationburdenvariantdepth, aes(x = SampleName, y = MutationBurden)) +
  geom_point() +
  geom_hline(yintercept = 41.9, linetype = "dotted", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x labels for readability
  labs(title = "Scatter Plot - Mixture Modelling True",
       x = "Sample Name",
       y = "Mutation Burden")


dataf4 <- Mutationburdenvariantdepth
View(dataf4)
#Rename columns
colnames(dataf4)[which(names(dataf4) == "colSums.Variant_depth_MMT...2.")] <- "Mutational_Burden"
View(dataf4)

#install.packages("ggplot2")
library(ggplot2)
library(tibble)
df <- tibble::rownames_to_column(dataf4, "Added_Column_MedianVAF")
View(df)
colnames(df)[which(names(df) == "Added_Column_MedianVAF")] <- "Sample"


Blood_meta <- read.table("/lustre/scratch126/casm/team274sb/md39/All_fetal_mapping_md39/00_Supplementary_data/Blood_Only2.txt", header = TRUE, sep = "\t")

metadata = Blood_meta
View(metadata)
df$Histo = NA

for(i in 1:nrow(df)){

  df$Histo[i] = metadata[metadata$Sample_ID == df$Sample[i],]$Histo

}
View(df)

organ_vec <- character(nrow(df))
germ_layer_vec <- character(nrow(df))
histo_vec <- character(nrow(df))
bulk_phenotype_vec <- character(nrow(df))
organ_histo_vec <- character(nrow(df))

# Loop through each row in df
for (i in 1:nrow(df)) {
  # Find the corresponding Sample in metadata
  sample_match <- metadata$Sample_ID == df$Sample[i]

  # If a match is found, populate the new columns
  if (any(sample_match)) {
    organ_vec[i] <- metadata$Organ[sample_match]
    germ_layer_vec[i] <- metadata$Germ_layer[sample_match]
    histo_vec[i] <- metadata$Histo[sample_match]
    bulk_phenotype_vec[i] <- metadata$Bulk_phenotype[sample_match]
    organ_histo_vec[i] <- metadata$Organ_Histo[sample_match]
  } else {
    # If no match is found, you can choose to set default values or leave them as NA
    organ_vec[i] <- NA
    germ_layer_vec[i] <- NA
    histo_vec[i] <- NA
    bulk_phenotype_vec[i] <- NA
    organ_histo_vec[i] <- NA
  }
}

# Add the new columns to df
df$Organ <- organ_vec
df$Germ_layer <- germ_layer_vec
df$Histo <- histo_vec
df$Bulk_phenotype <- bulk_phenotype_vec
df$Organ_Histo <- organ_histo_vec

View(df)
subset_df_blood <- df[df$Histo == "Blood", ]
View(subset_df_blood)
subset_df_blood_final <- subset_df_blood[subset_df_blood$Organ != "Umbilical_cord", ]
View(subset_df_blood_final)

ggplot(subset_df_blood_final, aes(x = reorder(Sample, -Mutational_Burden), y = Mutational_Burden)) +
  geom_point() +
  theme(axis.text.x = element_blank()) +  # This hides the sample names on the x-axis
  labs(x = "Samples", y = "Mutational Burden") +
  theme_minimal()


library(ggplot2)
library(ggbeeswarm)

#Boxplot and dot plot overylayed of all samples
ggplot(subset_df_blood_final) +
  geom_boxplot(mapping = aes(x = Sample, y = Mutational_Burden), outlier.shape = NA) +
  geom_quasirandom(mapping = aes(x = Sample, y = Mutational_Burden)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Overlayed boxpot and dotplot reordered from lowest to highest
ggplot(df, aes(x = reorder(Histo, Mutational_Burden, FUN = median), y = Mutational_Burden)) + geom_boxplot() + geom_quasirandom() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Overlayed boxpot and dotplot reordered from highest to lowest
ggplot(df, aes(x = reorder(Histo, -Mutational_Burden, FUN = median), y = Mutational_Burden)) + geom_boxplot() + geom_quasirandom() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
labs(x = "Mutational Burden", y = "Tissue")

ggplot(df, aes(x = Mutational_Burden, y = reorder(Histo, -Mutational_Burden, FUN = median))) +
  geom_boxplot() +
  geom_quasirandom() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Mutational Burden", y = "Tissue")
///


#Boxplot and dot plot overylayed of all samples
ggplot(df) +
  geom_boxplot(mapping = aes(x = Organ_Histo, y = Mutational_Burden), outlier.shape = NA) +
  geom_quasirandom(mapping = aes(x = Organ_Histo, y = Mutational_Burden)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Overlayed boxpot and dotplot reordered from lowest to highest
ggplot(df, aes(x = reorder(Organ_Histo, Mutational_Burden, FUN = median), y = Mutational_Burden)) + geom_boxplot() + geom_quasirandom() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Overlayed boxpot and dotplot reordered from highest to lowest
ggplot(df, aes(x = reorder(Organ_Histo, -Mutational_Burden, FUN = median), y = Mutational_Burden)) + geom_boxplot() + geom_quasirandom() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#Turn variant depth file into a binary using the >2 mutations to be real cut off

#Variant_depth1 is binary (1 or 0)
Variant_depth1 = Variant_depth_MMT
Variant_depth1[Variant_depth1 < 2] <- 0
Variant_depth1[Variant_depth1 >= 2] <- 1
View(Variant_depth1)
tVariant_depth <- as.data.frame(t(Variant_depth1))
View(tVariant_depth) #1281 rows, 10696 columns
tVariant_depth$Sample_ID <- rownames(tVariant_depth) #extra column added
View(tVariant_depth)


install.packages("MASS")
install.packages("reshape2")
install.packages("reshape")

library(MASS)
library(reshape2)
library(reshape)
library(dplyr)


#Nathan code
tVariant_depth.m <- melt(tVariant_depth, id.vars="Sample_ID")
View(tVariant_depth.m)
tVariant_depth.m <- tVariant_depth.m[tVariant_depth.m$value > 0,]
View(tVariant_depth.m)
colnames(tVariant_depth.m) <- c("Sample_ID", "Variant","Count")
View(tVariant_depth.m)
View(metadata)
#Are there matches? This is to make sure the same samples are in tVariant_depth.m and metadata
check <- intersect(tVariant_depth.m$Sample_ID, metadata$Sample_ID)
View(check) #1281 samples
print(paste(sep="", length(check), " Matches between Variant Depth Dataframe and Metadata"))

#Combine datasets once we confirmed they have the same samples
variants_annotated <- inner_join(tVariant_depth.m, metadata, by="Sample_ID")
View(variants_annotated)

#Check for unique mutations
length(unique(variants_annotated$Variant)) #9855 unique variants
length(unique(variants_annotated$Variant)) == nrow(variants_annotated) #if returns TRUE then there are no duplicates

#Create a table of all mutations with frequency
n_occur <- table(variants_annotated$Variant)
View(n_occur)
class(n_occur)
n_occur <- as.data.frame(n_occur)

#Create a table with mutations in atleast 2 samples
Variantsmorethanone <- variants_annotated[variants_annotated$Variant %in% n_occur$Var1[n_occur$Freq > 1],]
View(variants_annotated)
View(Variantsmorethanone)
nrow(variants_annotated) #Total number of mutations - 61530
nrow(Variantsmorethanone) #This tells you how many mutations are in atleast 2 samples - 55755

#Duplicate rows
variants_annotated[duplicated(variants_annotated$Variant),]

#Numberofduplicaterows
dim(variants_annotated[duplicated(variants_annotated$Variant),])[1]


#Frequency of mutations
Freq <- as.data.frame(table(Variantsmorethanone$Variant))
View(Freq) #Table of mutation frequency

variants_annotated$VAF <- 0
VAF = VAF_MMT
# Loop through each sample in variants_annotated
for (sample in unique(variants_annotated$Sample_ID)) {
  # Check if the sample is a column in the VAF dataframe
  if (sample %in% colnames(VAF)) {
    # Get the indices of the mutations for the current sample
    indices <- match(variants_annotated$Variant[variants_annotated$Sample_ID == sample], rownames(VAF))
    # Extract the VAF values for these mutations
    vaf_values <- VAF[indices, sample]
    # Assign the VAF values to the 'VAF' column for the current sample
    variants_annotated$VAF[variants_annotated$Sample_ID == sample] <- vaf_values
  }
}

# If there are any NAs introduced by mismatches, replace them with 0
variants_annotated$VAF[is.na(variants_annotated$VAF)] <- 0

View(variants_annotated)
# Copy dataframe into a new one and remove specified columns
final_variants_annotated_with_vaf <- subset(variants_annotated, select = -c(PD_ID, Count, Organ_Histo))
View(final_variants_annotated_with_vaf)

#Create a ggplot for samples vs variants with VAF informing the colour
ggplot(final_variants_annotated_with_vaf, aes(x = Sample_ID, y = VAF, color = VAF)) +
  geom_point(alpha = 0.6) + # Use alpha to adjust point transparency if there are many overlapping points
  scale_color_gradient(low = "blue", high = "red") + # Color gradient for VAF values
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + # Rotate x-axis text for readability
  labs(x = "Sample ID", y = "Variant Allele Frequency (VAF)", color = "VAF") +
  ggtitle("VAF of Mutations for Each Sample")

View(final_variants_annotated_with_vaf)

#VAF filter
VAF_filtered <- final_variants_annotated_with_vaf %>%
  filter(VAF > 0.25, VAF < 1)
View(VAF_filtered)

mutational_burden <- VAF_filtered %>%
  group_by(Sample_ID) %>%
  summarise(MutationCount = n())
View(mutational_burden)

colnames(mutational_burden) <- c("SampleName", "MutationBurden")

# Convert factor to character (if necessary)
mutational_burden$SampleName <- as.character(mutational_burden$SampleName)

# Sort the data frame by MutationBurden in descending order
mutational_burden <- mutational_burden %>%
  arrange(desc(MutationBurden)) %>%
  mutate(SampleName = factor(SampleName, levels = SampleName))

View(mutational_burden)
# Create the scatter plot with a dotted red line at y = 41.9
ggplot(mutational_burden, aes(x = SampleName, y = MutationBurden)) +
  geom_point() +
  geom_hline(yintercept = 41.9, linetype = "dotted", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x labels for readability
  labs(title = "Scatter Plot of Mutation Burden by Sample",
       x = "Sample Name",
       y = "Mutation Burden")


variant_frequency_filtered <- VAF_filtered %>%
  group_by(Variant) %>%
  summarise(Frequency = n()) %>%
  arrange(desc(Frequency))

View(variant_frequency_filtered)

Freq$Var1 <- as.character(Freq$Var1)

# Create a new table with the difference in mutations
difference_mutations <- Freq %>%
  filter(!Var1 %in% variant_frequency_filtered$Variant, Freq > 0)

View(difference_mutations) #these mutations are removed by the low VAF filter (0.25 or less and if == 1)


#Generate JBrowse images for 1701 variants = sum(Freq$Freq>1)

setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/")
options(stringsAsFactors = F)
NV=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NV_filtered_all.txt")
dim(NV)
NV = NV[rownames(NV) %in% difference_mutations$Var1, ]
dim(NV)
muts_shared=rownames(NV)[rowSums(NV>3)>1]
Muts_coord=matrix(ncol=4,unlist(strsplit(muts_shared,split="_")),byrow = T) #splitting into 4 columns

intv=50 #interval 50, make new dataframe, where the start is position -50 and the end is position +50 can give coordinates to JBriwse
Muts_bed = data.frame(Chr=Muts_coord[,1],
                      Start=as.numeric(Muts_coord[,2])-intv,
                      End=as.numeric(Muts_coord[,2])+intv)
rownames(Muts_bed)=muts_shared
project_number=3187

samples = colnames(NV)
#Make sure there is a # before JBrowse as it needs this to work
for(mut in muts_shared){
  samples_select=samples[NV[mut,]>3]
  if(length(samples_select)>=5) samples_select=samples_select[1:5] #If 5 or more samples have that mutation, we pick 5 samples (e.g. if we have 10 we would only pick 5),
  if(length(samples_select)<5) samples_select=c(samples_select,sample(samples[!samples %in% samples_select], 2, replace = FALSE))
  header=paste0("# Jbrowse https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F",project_number,"&tracks=DNA%2C",paste(paste0(samples_select,"_bwa"),collapse="%2C"))
  write(header,paste0(mut,".bed.tmp"))
  write.table(Muts_bed[mut,],paste0(mut,".bed.tmp"),row.names = F,col.names = F,quote=F,sep="\t",append=T)
}



#Run in command line:

#Generate a password file
Directory - /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/
cat *.bed.tmp > Jbrowse.bed

chmod +x jbrowse_rasterize_v2.sh
bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R 'select[mem>32000] rusage[mem=32000]' -M32000 -n 10 /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/jbrowse_rasterize_v2.sh

//////////////////////////////////////////////////////////////////////////////////////////////////////
  #Remove artefacts identified in JBRowse

  setwd("/lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/Jbrowse2/Artefacts")
#Create folders - Germline, True_mutations, Artefacts, Unsure
#Inspect each JBrowse image and sort

#Copy files names of all mutations into a single text file and remove.png from the end
ls *.png | sed 's/\.png$//' > All_Artefact_JBrowse_names.txt

#All the mutations are in the NV file
Variant_depth2 <- read.table("/lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/03_DNA_Analysis/Patient_snv_NV_filtered_all.txt", header = TRUE)
Variant_depth2 <- Variantsmorethanone

#We need to remove the artefacts from this
# Read the list of variants to be removed from the text file
variants_to_remove <- readLines("All_Artefact_JBrowse_names.txt") #46 Variants

checkmutation <- Variantsmorethanone$Variant
checkmutation <- rownames(NV)
View(checkmutation)
# Read the 'All_Artefact_JBrowse_names' file
artefact_file <- "All_Artefact_JBrowse_names.txt"
artefact_data <- readLines(artefact_file)

# Create an empty data frame to store the extracted values
extracted_data <- data.frame(chromosome = character(),
                             position = numeric(),
                             base1 = character(),
                             base2 = character(),
                             stringsAsFactors = FALSE)

# Loop through each line in the 'checkmutation' data frame
for (mutation in checkmutation) {
  # Split the mutation into chromosome, position, and bases
  parts <- strsplit(mutation, "_|-")
  chromosome <- parts[[1]][1]
  position <- as.numeric(parts[[1]][2])
  base1 <- parts[[1]][3]
  base2 <- parts[[1]][4]

  # Create a new row in the data frame
  new_row <- data.frame(chromosome = chromosome,
                        position = position,
                        base1 = base1,
                        base2 = base2,
                        stringsAsFactors = FALSE)

  # Add the new row to the extracted data frame
  extracted_data <- rbind(extracted_data, new_row)
}

# Print the extracted data frame
View(extracted_data)


# Read the 'All_Artefact_JBrowse_names' file
artefact_file <- "All_Artefact_JBrowse_names.txt"
artefact_data <- readLines(artefact_file)

# Create an empty data frame to store the extracted values
extracted_data2 <- data.frame(chromosome = character(),
                              start_position = character(),
                              end_position = character(),
                              stringsAsFactors = FALSE)

# Loop through each line in the 'All_Artefact_JBrowse_names' file
for (line in artefact_data) {
  # Split the line into chromosome and position range
  parts <- strsplit(line, "_|-")
  chromosome <- gsub("chr", "", parts[[1]][1])
  start_position <- parts[[1]][2]
  end_position <- parts[[1]][3]

  # Create a new row in the data frame
  new_row <- data.frame(chromosome = chromosome,
                        start_position = start_position,
                        end_position = end_position,
                        stringsAsFactors = FALSE)

  # Add the new row to the extracted data frame
  extracted_data2 <- rbind(extracted_data2, new_row)
}

# Print the extracted data frame
View(extracted_data2)


# Convert chromosome format in extracted_data
extracted_data$chromosome <- gsub("chr", "", extracted_data$chromosome)

# Extract the chromosome number from the chromosome column in extracted_data2
extracted_data2$chromosome_num <- gsub("chr", "", extracted_data2$chromosome) #46 entries here

# Convert position column to numeric in extracted_data and extracted_data2
extracted_data$position <- as.numeric(extracted_data$position)
extracted_data2$start_position <- as.numeric(extracted_data2$start_position)
extracted_data2$end_position <- as.numeric(extracted_data2$end_position)

# Create an empty data frame to store the matched positions
extracted_data3 <- data.frame(chromosome = character(),
                              position = numeric(),
                              stringsAsFactors = FALSE)

# Loop through each unique chromosome in extracted_data
for (chrom in unique(extracted_data$chromosome)) {
  # Get the positions for the current chromosome in extracted_data
  positions <- extracted_data$position[extracted_data$chromosome == chrom]

  # Get the start and end positions for the current chromosome in extracted_data2
  start_positions <- extracted_data2$start_position[extracted_data2$chromosome_num == chrom]
  end_positions <- extracted_data2$end_position[extracted_data2$chromosome_num == chrom]

  # Check if the start_positions and end_positions are non-empty
  if (length(start_positions) > 0 && length(end_positions) > 0) {
    # Create a logical vector to track if a position is within any range
    position_within_range <- logical(length(positions))

    # Loop through each range in extracted_data2
    for (i in 1:length(start_positions)) {
      # Check if each position lies within the current range
      position_within_range <- position_within_range | (positions >= start_positions[i] & positions <= end_positions[i])
    }

    # Filter the positions that lie within any range
    positions_within_range <- positions[position_within_range]

    # Create a new data frame with the matched positions
    matched_positions <- data.frame(chromosome = rep(chrom, length(positions_within_range)),
                                    position = positions_within_range,
                                    stringsAsFactors = FALSE)

    # Add the matched positions to the extracted_data3 data frame
    extracted_data3 <- rbind(extracted_data3, matched_positions)
  }
}

#Remove duplicates
extracted_data3 <- unique(extracted_data3)

# Print the extracted_data3 data frame
View(extracted_data3) #There are 47 variants here none duplicated (it should be 46)


extracted_data3$chromosome <- paste0("chr", extracted_data3$chromosome)

# Filter rows based on matching chromosome and position
extracted_data4 <- extracted_data3[extracted_data3$chromosome %in% Muts_coord$chromosome & extracted_data3$position %in% Muts_coord$position, ]

# Reset row names
rownames(extracted_data4) <- NULL

# View the resulting data frame
View(extracted_data4) #Now there are 46 variants


#Remove artefacts identified in extracted_data4 from Variant_depth2 (rownames across all samples)
# Create a new column in extracted_data4 combining chromosome and position
extracted_data4$variant <- paste0(extracted_data4$chromosome, "_", extracted_data4$position)

# Create Variant_depth3 by removing rows found in extracted_data4 from Variant_depth2
Variant_depth3 <- Variant_depth2

for (variant_extracted in extracted_data4$variant) {
  # Extract the chromosome and position from the extracted variant
  variant_parts <- strsplit(variant_extracted, "_")
  variant_chr <- variant_parts[[1]][1]
  variant_pos <- variant_parts[[1]][2]

  # Find matching rows in Variant_depth3 and remove them
  rows_to_remove <- rownames(Variant_depth3)[grepl(paste0(variant_chr, "_", variant_pos), rownames(Variant_depth3))]
  Variant_depth3 <- Variant_depth3[!(rownames(Variant_depth3) %in% rows_to_remove), ]
}

View(Variant_depth3) #This has 46 less variants than Variant_depth2

dataf5 <- data.frame(colSums(Variant_depth3 > 2))


#Rename columns
colnames(dataf5)[which(names(dataf5) == "colSums.Variant_depth2...2.")] <- "Mutational_Burden"
df2 <- tibble::rownames_to_column(dataf5, "Added_Column_MedianVAF")
colnames(df2)[which(names(df2) == "Added_Column_MedianVAF")] <- "Sample"
colnames(df2)[which(names(df2) == "colSums.Variant_depth3...2.")] <- "Mutational_Burden"
View(df2)

df2$Organ = NA

for(i in 1:nrow(df2)){

  df2$Organ[i] = variants_annotated[variants_annotated$Sample_ID == df2$Sample[i],]$base_organ

}

View(df2)

#Boxplot and dot plot overylayed of all samples
ggplot(df2) +
  geom_boxplot(mapping = aes(x = Organ, y = Mutational_Burden), outlier.shape = NA) +
  geom_quasirandom(mapping = aes(x = Organ, y = Mutational_Burden)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Overlayed boxpot and dotplot reordered from lowest to highest
ggplot(df2, aes(x = reorder(Organ, -Mutational_Burden, FUN = median), y = Mutational_Burden)) +
  geom_boxplot() +
  geom_quasirandom() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Organ", y = "Mutational Burden")


#Use the artefact dataframe to remove mutations from NV matrix and then plot mutational signature




#Generate heatmap of mutations
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

AtVariant_depth <- tVariant_depth
AtVariant_depth$Sample_ID <- rownames(tVariant_depth)
AtVariant_depth.m <- melt(AtVariant_depth, id.vars="Sample_ID")
colnames(AtVariant_depth.m) <- c("Sample_ID", "Variant","Count")

#Are there matches? This is to make sure the same samples are in tVariant_depth.m and metadata (should be 160).
check <- intersect(AtVariant_depth.m$Sample_ID, metadata$Sample_ID)
View(check)
print(paste(sep="", length(check), " Matches between Variant Depth Dataframe and Metadata"))

#Combine datasets once we confirmed they have the same samples
Avariants_annotated <- join(AtVariant_depth.m, metadata, by="Sample_ID")
View(Avariants_annotated)
n_distinct(Avariants_annotated$Variant)

#Aggregate
Avariants_annotated$Count = as.numeric(Avariants_annotated$Count)

View(Avariants_annotated[is.na(Avariants_annotated$Count),])
Avariants_annotated <- Avariants_annotated[Avariants_annotated$Variant != "SampleID",]
Avariants_annotated$Variant = as.character(Avariants_annotated$Variant)
sum(is.na(Avariants_annotated))
library(tidyverse)
Avariants_annotated3 = Avariants_annotated %>% group_by(Variant,Histo,Organ) %>% dplyr::summarise(total_count=sum(Count))

Avariants_annotated3 = Avariants_annotated <- Avariants_annotated[c("Variant", "Count", "Organ", "Histo")]

#The rownames need to be unique of the matrix to plot on heatmap
#Avariants_annotated$Variant2 <- paste0(Avariants_annotated$Variant, "_", rownames(Avariants_annotated))

Avariants_annotated3$Histo2 <- paste0(Avariants_annotated3$Organ, "_", Avariants_annotated3$Histo)

Avariants_annotated_wide <- pivot_wider(data = Avariants_annotated3, names_from = "Histo2", values_from = "Count", id_cols = "Variant")

Avariants_annotated_wide2 <- as.matrix(Avariants_annotated_wide[,-1])

rownames(Avariants_annotated_wide2)=Avariants_annotated_wide$Variant

#Load in metadata with germ layer information
meta_data_germ <- read.table("/lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/00_Supplementary_Data/Right and left Meta Data and germ layer 2855 PD53943 Final MD.txt", header = TRUE)
View(meta_data_germ)

meta_data <- read.table("/lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/00_Supplementary_Data/Right and left Meta Data and germ layer 2855 PD53943 Final MD.txt", header = TRUE)


class(meta_data_germ)
# Import mutation data as a matrix or data frame
mutation_matrix <- meta_data_germ

View(mutation_matrix)
# Define row and column annotations
organ_type <- c("Umbilical_cord", "Head", "Leg_right", "Placenta", "Liver", "Distal_bowel" , "Proximal_bowel", "Bladder", "Skin", "Brain", "Spinal_column", "Thyroid_left", "Arm_left", "Arm_right", "Adrenal_right", "Kidney_left", "Lung_left", "Lung_right", "Gonad", "Stomach", "Kidney_right", "Adrenal_left")
germlayer <- c("Ectoderm", "Endoderm", "Mesoderm")

row_annotation <- data.frame(organ_type, germlayer)
col_annotation <- data.frame(organ_type, germlayer)

# Create heatmap
install.packages("ComplexHeatmap")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
col_annotation = HeatmapAnnotation(df = data.frame(
  organ = Avariants_annotated3$Organ[match(colnames(Avariants_annotated_wide2),Avariants_annotated3$Histo2)],
  germline =  ),
  col = list(germlayer = c('Kidney'='red','Brain'='blue')))
Heatmap(Avariants_annotated_wide2,
        show_row_names = FALSE,
        name = "Mutation",
        #col = circlize::colorRamp2(c(0, 1), c("white", "red")),
        row_dend_reorder = TRUE,
        row_title = "Samples",cluster_rows = F,cluster_columns = F,
        #row_split = , column_split = ,
        bottom_annotation = col_annotation,
        #row_annotation = row_annotation,
        #column_annotation = col_annotation,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")


HeatmapKidneysData1$Variant <- factor(HeatmapKidneysData1$Variant,
                                      levels = unique(HeatmapKidneysData1$Variant))







#Plot mutational spectra
setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True")
NV=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NV_filtered_all.txt", header = TRUE)
dim(NV)
NV = NV[rownames(NV) %in% difference_mutations$Var1, ]
dim(NV)
NR=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NR_filtered_all.txt", header = TRUE)
View(NV)

muts_shared=rownames(NV)[rowSums(NV>3)==0] #If shared mutation keep as 1, if total mutation change to 0
Muts_coord=matrix(ncol=4,unlist(strsplit(muts_shared,split="_")),byrow = T)

colnames(Muts_coord) = c("Chrom","Pos","Ref","Alt")

source ("/lustre/scratch126/casm/team274sb/md39/plot_spectrum_hg38.R")
Muts_coord=matrix(ncol=4,unlist(strsplit(muts_shared,split="_")),byrow = T)
plot_spectrum(Muts_coord,save=NULL)

#Plot signature from branch.txt file (output from build.phylogeny.R script)
NV=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_False/Patient_snv_assigned_to_branches.txt", header = TRUE)
View(NV)
NV <- NV[, !(names(NV) %in% c("Branch", "Patient", "SampleID"))]
colnames(NV) = c("Chrom","Pos","Ref","Alt")
plot_spectrum(NV,save=NULL)

#Compare NV MM T and NV MM F mutations - what does the binomial mixture model filter out? What is the signature of this?
NV_MMT=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_assigned_to_branches.txt", header = TRUE)
NV_MMF=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_False/Patient_snv_assigned_to_branches.txt", header = TRUE)
dim(NV_MMT)
dim(NV_MMF)
View(NV_MMT)

# Creating a unique identifier for each mutation in both data frames
NV_MMF$MutationID <- with(NV_MMF, paste(Chr, Pos, Ref, Alt, sep="_"))
NV_MMT$MutationID <- with(NV_MMT, paste(Chr, Pos, Ref, Alt, sep="_"))

# Finding unique MutationIDs in NV_MMF that are not in NV_MMT
unique_MutationIDs <- setdiff(NV_MMT$MutationID, NV_MMF$MutationID)

View(unique_MutationIDs) #This is a list of mutations the beta-binomial filters out.
dim(NV_MMT)
dim(NV_MMF)

#Jbrowse
setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mixture_Modelling")
options(stringsAsFactors = F)
NV=read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_NV_filtered_all.txt")
Muts_coord=matrix(ncol=4,unlist(strsplit(unique_MutationIDs,split="_")),byrow = T)
View(Muts_coord)

combined_info <- apply(Muts_coord, 1, function(x) paste(x, collapse = "_"))

intv=50 #interval 50, make new dataframe, where the start is position -50 and the end is position +50 can give coordinates to JBriwse
muts_shared=combined_info
Muts_bed = data.frame(Chr=Muts_coord[,1],
                      Start=as.numeric(Muts_coord[,2])-intv,
                      End=as.numeric(Muts_coord[,2])+intv)
rownames(Muts_bed)=muts_shared
project_number=3187

samples = colnames(NV)
#Make sure there is a # before JBrowse as it needs this to work
for(mut in muts_shared){
  samples_select=samples[NV[mut,]>3]
  if(length(samples_select)>=5) samples_select=samples_select[1:5] #If 5 or more samples have that mutation, we pick 5 samples (e.g. if we have 10 we would only pick 5),
  if(length(samples_select)<5) samples_select=c(samples_select,sample(samples[!samples %in% samples_select], 2, replace = FALSE))
  header=paste0("# Jbrowse https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F",project_number,"&tracks=DNA%2C",paste(paste0(samples_select,"_bwa"),collapse="%2C"))
  write(header,paste0(mut,".bed.tmp"))
  write.table(Muts_bed[mut,],paste0(mut,".bed.tmp"),row.names = F,col.names = F,quote=F,sep="\t",append=T)
}

ls | grep "\.bed\.tmp$" | wc -l #1681

#Run in command line:

#Generate a password file
Directory - /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mixture_Modelling
  cat *.bed.tmp > Jbrowse.bed

chmod +x jbrowse_rasterize_v2.sh
bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R 'select[mem>32000] rusage[mem=32000]' -M32000 -n 10 /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mixture_Modelling/jbrowse_rasterize_v2.sh


#Signature
#Plot mutational spectra
Muts_coord=matrix(ncol=4,unlist(strsplit(unique_MutationIDs,split="_")),byrow = T)
View(Muts_coord)
colnames(Muts_coord) = c("Chrom","Pos","Ref","Alt")
source ("/lustre/scratch126/casm/team274sb/md39/plot_spectrum_hg38.R")
plot_spectrum(Muts_coord,save=NULL)

#VAF cut off - low VAF filter - inspect JBrowse for these mutations
setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mutations_removed_by_VAF_filter")
View(NV_MMT)
View(VAF_filtered)

# Convert the MutationID column of NV_MMT and Variant column of VAF_filtered into vectors
NV_MMT_Mutations <- NV_MMT$MutationID
VAF_filtered_Mutations <- VAF_filtered$Variant

# Find the difference: mutations in NV_MMT that are not in VAF_filtered
filtered_out_mutations <- setdiff(NV_MMT_Mutations, VAF_filtered_Mutations)

# Create a dataframe from the filtered out mutations
filtered_out_mutations_df <- data.frame(MutationID = filtered_out_mutations)

unique_MutationIDs <- unique(filtered_out_mutations_df)
View(unique_MutationIDs) #1991

#Muts_coord=matrix(ncol=4,unlist(strsplit(unique_MutationIDs,split="_")),byrow = T)
# Extract the MutationID column as a character vector
mutation_ids <- as.character(unique_MutationIDs$MutationID)

# Use strsplit to split the mutation IDs at underscores
split_ids <- strsplit(mutation_ids, split = "_")

# Unlist and create a matrix
Muts_coord <- matrix(unlist(split_ids), ncol = 4, byrow = TRUE)

# View the resulting matrix
View(Muts_coord)

intv=50 #interval 50, make new dataframe, where the start is position -50 and the end is position +50 can give coordinates to JBriwse
combined_info <- apply(Muts_coord, 1, function(x) paste(x, collapse = "_"))
muts_shared=combined_info
Muts_bed = data.frame(Chr=Muts_coord[,1],
                      Start=as.numeric(Muts_coord[,2])-intv,
                      End=as.numeric(Muts_coord[,2])+intv)

rownames(Muts_bed)=muts_shared
project_number=3187
samples = colnames(NV)

#Make sure there is a # before JBrowse as it needs this to work
for(mut in muts_shared){
  samples_select=samples[NV[mut,]>3]
  if(length(samples_select)>=5) samples_select=samples_select[1:5] #If 5 or more samples have that mutation, we pick 5 samples (e.g. if we have 10 we would only pick 5),
  if(length(samples_select)<5) samples_select=c(samples_select,sample(samples[!samples %in% samples_select], 2, replace = FALSE))
  header=paste0("# Jbrowse https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F",project_number,"&tracks=DNA%2C",paste(paste0(samples_select,"_bwa"),collapse="%2C"))
  write(header,paste0(mut,".bed.tmp"))
  write.table(Muts_bed[mut,],paste0(mut,".bed.tmp"),row.names = F,col.names = F,quote=F,sep="\t",append=T)
}


ls | grep "\.bed\.tmp$" | wc -l #1991

#Run in command line:
cat *.bed.tmp > Jbrowse.bed

chmod +x jbrowse_rasterize_v2.sh
bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R 'select[mem>32000] rusage[mem=32000]' -M32000 -n 10 /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mutations_removed_by_VAF_filter/jbrowse_rasterize_v2.sh


#Signature
#Plot mutational spectra
colnames(Muts_coord) = c("Chrom","Pos","Ref","Alt")
source ("/lustre/scratch126/casm/team274sb/md39/plot_spectrum_hg38.R")
plot_spectrum(Muts_coord,save=NULL)



#Mutations after the low VAF filter - what is left - inspect these mutations and plot the spectra
setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Low_VAF_Filter")
View(NV_MMT)
View(VAF_filtered)

# This keeps only the rows from VAF_filtered where the Variant matches with MutationID in NV_MMT
VAF_Filtered_Mutations <- merge(VAF_filtered, NV_MMT, by.x = "Variant", by.y = "MutationID")

# Retain only the columns from VAF_filtered
VAF_Filtered_Mutations <- VAF_Filtered_Mutations[names(VAF_filtered)]

View(VAF_Filtered_Mutations)
unique_MutationIDs <- data.frame(MutationID = VAF_Filtered_Mutations$Variant)
unique_MutationIDs <- unique(unique_MutationIDs)
View(unique_MutationIDs)

#Muts_coord=matrix(ncol=4,unlist(strsplit(unique_MutationIDs,split="_")),byrow = T)
# Extract the MutationID column as a character vector
mutation_ids <- as.character(unique_MutationIDs$MutationID)

# Use strsplit to split the mutation IDs at underscores
split_ids <- strsplit(mutation_ids, split = "_")

# Unlist and create a matrix
Muts_coord <- matrix(unlist(split_ids), ncol = 4, byrow = TRUE)

# View the resulting matrix
View(Muts_coord)

intv=50 #interval 50, make new dataframe, where the start is position -50 and the end is position +50 can give coordinates to JBriwse
combined_info <- apply(Muts_coord, 1, function(x) paste(x, collapse = "_"))
muts_shared=combined_info
Muts_bed = data.frame(Chr=Muts_coord[,1],
Start=as.numeric(Muts_coord[,2])-intv,
End=as.numeric(Muts_coord[,2])+intv)

rownames(Muts_bed)=muts_shared
project_number=3187
samples = colnames(NV)

#Make sure there is a # before JBrowse as it needs this to work
for(mut in muts_shared){
samples_select=samples[NV[mut,]>3]
if(length(samples_select)>=5) samples_select=samples_select[1:5] #If 5 or more samples have that mutation, we pick 5 samples (e.g. if we have 10 we would only pick 5),
if(length(samples_select)<5) samples_select=c(samples_select,sample(samples[!samples %in% samples_select], 2, replace = FALSE))
header=paste0("# Jbrowse https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F",project_number,"&tracks=DNA%2C",paste(paste0(samples_select,"_bwa"),collapse="%2C"))
write(header,paste0(mut,".bed.tmp"))
write.table(Muts_bed[mut,],paste0(mut,".bed.tmp"),row.names = F,col.names = F,quote=F,sep="\t",append=T)
}


ls | grep "\.bed\.tmp$" | wc -l #4311

#Run in command line:
cat *.bed.tmp > Jbrowse.bed

chmod +x jbrowse_rasterize_v2.sh
bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R 'select[mem>32000] rusage[mem=32000]' -M32000 -n 10 /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Low_VAF_Filter/jbrowse_rasterize_v2.sh


#Signature
#Plot mutational spectra
colnames(Muts_coord) = c("Chrom","Pos","Ref","Alt")
source ("/lustre/scratch126/casm/team274sb/md39/plot_spectrum_hg38.R")
plot_spectrum(Muts_coord,save=NULL)


#Mutations in mutations plotted on the tree
setwd("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mutations_Tree_Branches")
#Read in mutations from branches.txt (after beta binomial over dispersion calculation)
Branches <- read.table("/lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Patient_snv_assigned_to_branches.txt", header = TRUE, sep = "\t")
View(Branches)
Branches$Mutations <- paste(Branches$Chr, Branches$Pos, Branches$Ref, Branches$Alt, sep = "_")

unique_MutationIDs <- data.frame(MutationID = Branches$Mutations) #6302 mutations
unique_MutationIDs <- unique(unique_MutationIDs)
View(unique_MutationIDs) #6302 mutations

#Muts_coord=matrix(ncol=4,unlist(strsplit(unique_MutationIDs,split="_")),byrow = T)
# Extract the MutationID column as a character vector
mutation_ids <- as.character(unique_MutationIDs$MutationID)

# Use strsplit to split the mutation IDs at underscores
split_ids <- strsplit(mutation_ids, split = "_")

# Unlist and create a matrix
Muts_coord <- matrix(unlist(split_ids), ncol = 4, byrow = TRUE)

# View the resulting matrix
View(Muts_coord)

intv=50 #interval 50, make new dataframe, where the start is position -50 and the end is position +50 can give coordinates to JBriwse
combined_info <- apply(Muts_coord, 1, function(x) paste(x, collapse = "_"))
muts_shared=combined_info
Muts_bed = data.frame(Chr=Muts_coord[,1],
                      Start=as.numeric(Muts_coord[,2])-intv,
                      End=as.numeric(Muts_coord[,2])+intv)

rownames(Muts_bed)=muts_shared
project_number=3187
samples = colnames(NV)

#Make sure there is a # before JBrowse as it needs this to work
for(mut in muts_shared){
  samples_select=samples[NV[mut,]>3]
  if(length(samples_select)>=5) samples_select=samples_select[1:5] #If 5 or more samples have that mutation, we pick 5 samples (e.g. if we have 10 we would only pick 5),
  if(length(samples_select)<5) samples_select=c(samples_select,sample(samples[!samples %in% samples_select], 2, replace = FALSE))
  header=paste0("# Jbrowse https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F",project_number,"&tracks=DNA%2C",paste(paste0(samples_select,"_bwa"),collapse="%2C"))
  write(header,paste0(mut,".bed.tmp"))
  write.table(Muts_bed[mut,],paste0(mut,".bed.tmp"),row.names = F,col.names = F,quote=F,sep="\t",append=T)
}


ls | grep "\.bed\.tmp$" | wc -l #5970

#Run in command line:
cat *.bed.tmp > Jbrowse.bed

chmod +x jbrowse_rasterize_v2.sh
bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R 'select[mem>32000] rusage[mem=32000]' -M32000 -n 10 /lustre/scratch126/casm/team274sb/md39/Blood_Analysis_No_Filtering/Mixture_Modelling_True/Mutations_Tree_Branches/jbrowse_rasterize_v2.sh


#Signature
#Plot mutational spectra
colnames(Muts_coord) = c("Chrom","Pos","Ref","Alt")
source ("/lustre/scratch126/casm/team274sb/md39/plot_spectrum_hg38.R")
plot_spectrum(Muts_coord,save=NULL)





#Subset NV matrix for PD53943ey
NV_subset <- NV["PD53943ey"]
View(NV_subset)

#Taryn code to split rownames for input into mutational signature deconstruct sigs package. Add a sample column and change chromosome position in column name
#whichsignatures(), supply output in matrix. - muts_coordT = read.table(textConnection(rownames(Variant_depth3)),sep="_")
#Using GRch38 - use this as a reference
#Add 10% cut off to filter - won't be tobacco smoking in fetus etc.
#Anna script

View(muts_coord)