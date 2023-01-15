#Find Probe SNPs
library(dplyr)
library(tibble)
library("caroline")
library(dplyr)
####must have bedtools installed in /usr/bin/bedtools 
###Infer Probe Sequence region from manifest if using a different array
### (this has already been done for EPIC, 450k, and EPIC_450k_overlap; if using one of these
## skip to part 2)

###PART 1####
# ARRAY <- "450k"
# library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
# library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# library(minfi)
# if(ARRAY == "EPIC"){
#     data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19");
#     annoObj = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)}

# if(ARRAY == "450k"){
#     data("IlluminaHumanMethylation450kanno.ilmn12.hg19");
#     annoObj = getAnnotation("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# }
# annot <- as.data.frame(annoObj)
# #subset manifest for the Cytosine position
# cross_ref <- annot[,c(1:4,9)]
###coordinates of the probe 5' end (smaller coordinate than pos for - strand target and greater coordiante for + strand target)
###this step takes about half an hour####
# cross_ref$CpG_probe_5_prime_start <- 0
# i <- 1
# while (i <= nrow(cross_ref)){     
#     if (cross_ref[i, "strand"] == "+"){ 
#         #if + and type 1
#         if(cross_ref[i, "Type"] == "I"){
#          cross_ref[i, "CpG_probe_5_prime_start"] <- cross_ref[i, "pos"] + 49};
#         #if + and type 2
#         if (cross_ref[i, "Type"] == "II"){
#         cross_ref[i, "CpG_probe_5_prime_start"] <- cross_ref[i, "pos"] + 50}
#         };
#     if (cross_ref[i, "strand"] == "-"){
#         #if - and type 1
#         if(cross_ref[i, "Type"] == "I"){
#         cross_ref[i, "CpG_probe_5_prime_start"] <- cross_ref[i, "pos"] - 48};
#         #if - and type 2
#         if (cross_ref[i, "Type"] == "II"){
#         cross_ref[i, "CpG_probe_5_prime_start"] <- cross_ref[i, "pos"] - 49}
#         };
#     i <- i + 1
#     }
# write.csv(cross_ref, paste0(ARRAY,"_cross_ref.csv"))
###create plus and minus targets bed files for intersecting (plus target probes have their 5' start at a higher coordinate 
###than the target C whereas minus strand targets have their 5' start at a lower cooridinate than the target C
# cross_ref$end_of_snp_search <- NA
# cross_ref_plus <- cross_ref[cross_ref$strand == "+",]
# cross_ref_minus <- cross_ref[cross_ref$strand == "-",]
###plus and minus targeted probes have different probe sequence definitions 
# cross_ref_plus <- cross_ref_plus %>% mutate(end_of_snp_search = ifelse(Type == "I", pos - 1, ifelse(Type == "II", pos, end_of_snp_search)))
# cross_ref_minus <- cross_ref_minus %>% mutate(end_of_snp_search = ifelse(Type == "I", pos + 2, ifelse(Type == "II", pos + 1, end_of_snp_search)))
###take out correct columns, put columns in right correct, make the probe region open closed, and write out as a bed file
# cross_ref_plus_bed <- cross_ref_plus[,c("chr", "end_of_snp_search", "CpG_probe_5_prime_start", "Name", "strand", "Type", "pos")]
# cross_ref_minus_bed <- cross_ref_minus[,c("chr", "CpG_probe_5_prime_start", "end_of_snp_search", "Name" , "strand", "Type", "pos")]
# cross_ref_plus_bed <- cross_ref_plus_bed %>% rename("end_of_snp_search" = "start" , "CpG_probe_5_prime_start" = "end")
# cross_ref_minus_bed <- cross_ref_minus_bed %>% rename("CpG_probe_5_prime_start" = "start", "end_of_snp_search" = "end")
# cross_ref_bed <- rbind(cross_ref_plus_bed, cross_ref_minus_bed)
###make it open close coordinates
# cross_ref_bed <-  cross_ref_bed %>% mutate(start = start - 1)
# name = paste0(ARRAY, "_cross_ref.bed")
# write.table(cross_ref_bed, file = name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

##PART 2###Start here if you used 450k, EPIC, or EPIC_450k_overlap (put one of those three prefixes or the prefix of your cross_ref.bed file)

ARRAY <- "EPIC_450k_overlap"
###probe sequence is defined by the the sequence before the target position for - strand targets and after for + strand targets
###see diagram for exact sequence coordinate diagrams as - and + strand target probes differ in their exact probe sequend delimitations
###We also search the SBE position for type 1 and type2 probes, so we search within a 51bp region for all probes

##Input your bed file of SNP coordinates you'd like to compare
## Must be in UCSC BED format (chr, start, stop, id, QUAL) https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/vcf2bed.html
##VCF2bed can be used to convert vcfs to the first 5 column format. File should not have headers.
##or the following commands: 
##1) grep -v '^#' input.vcf | awk -F '\t' '{print $1,($2-1),$2,$3,$6,$4,$5}' > output_fix.ucsc.bed
##2) sed 's/ /\t/g' output_fix.ucsc.bed > correct_format.bed  
##3) sed 's/^/chr/' correct_format.bed > correct_format_with_chr.bed 
##the QUAL column is not used but there must be a 5th column which contains something (can just be .)
##in bedtools regions are delimited with (] notation.  So for SNPs, the stop position in your population_bed file should be the snp coordinate
##column 6 and column 7 must be the ref and alt alleles (these are needed to check if SBE snps are color-channel switching

##SNP_coordinates.bed is the full 1000genomes phase 3 bed file
population_bed <- "SNP_coordinates"
array_cross_ref <- paste0(ARRAY,"_cross_ref.bed")
##only took a few minutes for 78 million snp bed file
system(paste("/usr/bin/bedtools intersect -a " ,population_bed, "-b", array_cross_ref," -wa -wb  > ", file_name))
bedtools_intersect <- read.delim(paste0(ARRAY,"_bed_intersect"), header=FALSE)
##distance to pos (target Cytosine) from snp
bedtools_intersect$distances <- bedtools_intersect$V3 - bedtools_intersect$V14
##rename columns
bedtools_intersect <- bedtools_intersect %>% rename( "chr" = "V1" , "rsids" = "V4", "refs" = "V6" , "alts" = "V7", "probes"="V11" , "pos" = "V14" )
##each row corresponds to a snp in a probe with its associated CpG and distance away from target cytosine
bedtools_intersect$chr_pos <- paste0(bedtools_intersect$chr,":",bedtools_intersect$V3)

###PART 3
###Checks if SBE type 1 snps actually are color-channel switching.  If not, then remove them from the dataframe###
bedtools_intersect <- bedtools_intersect %>% tibble::rownames_to_column('row_num') 

##define the SBE type 1 snps
SBE_1_snps <- rbind(bedtools_intersect %>% filter(probes %in% plus_type1_probes & distances == -1), bedtools_intersect %>% filter(probes %in% minus_type1_probes & distances == 2))
bedtools_intersect[which(bedtools_intersect$probes %in% SBE_1_snps$probes),]
SBE_1_snps
non_CCS_snps <- c()
##check if SBE snps are actually color-channel switching (AKA C<->G or A<->T snps are not CCS)
##if they are not CCS then remove them from the bedtools_intersect
if(nrow(SBE_1_snps)>0){
for (i in 1:nrow(SBE_1_snps)){
    snp <- SBE_1_snps[i,];
    print(snp$row_num);
    #if not a CCS remove from bedtools_intersect
        if((snp$refs == "A" && snp$alts == "T") | (snp$refs =="T" && snp$alts =="A") | (snp$refs == "C" && snp$alts== "G") | (snp$refs =="G" && snp$alts =="C")){
            row <- as.numeric(snp$row_num);
            non_CCS_snps <- c(non_CCS_snps, row)
            };
    i <- i + 1
}
print(non_CCS_snps)
if(length(non_CCS_snps)>0){
    post_CCS_bedtools_intersect <- bedtools_intersect[-c(non_CCS_snps),]}
if(length(non_CCS_snps)==0){
    post_CCS_bedtools_intersect <- bedtools_intersect
}
}
if(nrow(SBE_1_snps)==0){
post_CCS_bedtools_intersect <- bedtools_intersect
}
###now bedtools_intersect contains the probes to be removed and the distances to each snp!!!!! 
write.delim(post_CCS_bedtools_intersect, file=paste0(ARRAY, "snp_probes"))





