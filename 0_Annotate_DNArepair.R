library(tidyverse)
library(splitstackshape)

#all.wgs.txt do mutated annotation first by ANOVAR ( RUN in UPPMAX)
# 
# module load  bioinfo-tools annovar annovar_data
# awk -F '\t' '{print $6, $7, $8, $9, $10, $1, $3}' all.wgs.txt | sed 1d > all.wgs.avinput
# annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 all.wgs.avinput $ANNOVAR_DATA/humandb  # downlaod exonicVariants



## 20230116 Togenther with Pan, we updated HR gene list and add new pathways -- ONCE
#genelist.DNArepair = read_tsv('~/OneDrive - KI.SE/Mac/Project/202205_JEM_reproduce/_Ref/DNArepair.genelist.tsv')
# genelist.DNArepair = readxl::read_excel('~/OneDrive - KI.SE/Mac/Project/202206_Kataegis/_Ref/DNA repair gene list-marked_v2-20220113-pan.230116-Hui.xlsx', sheet = 'DDRgene_20230116')
genelist.DNArepair = readxl::read_excel('_ref/DDR gene and pathways.20230118.xlsx', sheet = 1)
# check whether all genes could be found in our cohort 
annot_all = read_tsv('Results/SBS_Signature/all.wgs.avinput.variant_function', col_names = c('Region', 'Gene', 'Chr','Start', 'End', 'Ref', 'Alt', "Sample"))
#annot_all = read_tsv('../20230120/1_Signature/Combine_202301/Annovar/all.wgs.avinput.variant_function', col_names = c('Region', 'Gene', 'Chr','Start', 'End', 'Ref', 'Alt', "Sample")) # missed sample supply

annot_exonic = read_tsv('Results/SBS_Signature/all.wgs.avinput.exonic_variant_function', col_names = c('Line', 'Func.Alt', 'AAChange.refGene', 'Input')) 
#annot_exonic = read_tsv('../20230120/1_Signature/Combine_202301/Annovar/all.wgs.avinput.exonic_variant_function', col_names = c('Line', 'Func.Alt', 'AAChange.refGene', 'Input')) 


setdiff(genelist.DNArepair$Gene, annot_all$Gene) # ERCC5 (in introgenic region); MRE11A (checked Genecard, should modify to MRE11)



## run 
genelist.DNArepair = readxl::read_excel('_ref/DDR gene and pathways.20230118.xlsx', sheet = 1) %>% 
  cSplit('Pathways', '/', 'long') %>% select(Gene , Pathway = Pathways)
table(genelist.DNArepair$Pathway)

mut.func.splicing = annot_all %>% select(Func.Gene = Gene, everything()) %>% 
  #separate(Input, c("Chr", "Start","End", "Ref","Alt", "Sample"), sep = ' ') %>%
  filter(Region %in% c("splicing", "exonic;splicing" )) %>%
  mutate(Func.Gene = str_remove(Func.Gene, '.*;')) %>% 
  cSplit("Func.Gene",sep =  '(', direction =  'wide') %>%
  mutate(Func.Gene_2 = str_remove(Func.Gene_2, '\\)')) %>% 
  cSplit("Func.Gene_2",  ',', 'long') %>% 
  mutate(Func.Gene_2 = str_remove(Func.Gene_2, '.*:')) %>% 
  filter(!grepl('UTR3|5$', Func.Gene_2), !is.na(Func.Gene_2)) %>% # not include UTR region
  mutate(Func.Gene  = str_c(Func.Gene_1, Func.Gene_2, sep = ':' ), Gene = Func.Gene_1) %>% 
  inner_join(genelist.DNArepair) %>% ## filter gene %in% DDR gene list
  group_by(Sample,  Pathway) %>% 
  summarise(AAChange.refGene = str_c(unique(Func.Gene), collapse = ';')) 
  


df_annot = annot_exonic %>% 
  filter(!Func.Alt %in% c('synonymous SNV', 'unknown')) %>% 
  separate(Input, c("Chr", "Start","End", "Ref","Alt", "Sample"), sep = ' ') %>% 
  cSplit("AAChange.refGene", ',', 'long') %>% 
  mutate(AAChange.refGene = str_replace_all(AAChange.refGene, ':.*:', ':')) %>% 
  mutate(Gene = str_remove(AAChange.refGene, '.p.*')) %>% 
  inner_join(genelist.DNArepair) %>% 
  group_by(Sample,  Pathway) %>% 
  summarise(AAChange.refGene = str_c(unique(AAChange.refGene), collapse = ';')) %>% 
  bind_rows(mut.func.splicing) %>%   
  group_by(Sample,  Pathway) %>% 
  summarise(AAChange.refGene = str_c(unique(AAChange.refGene), collapse = ';')) %>% 
  spread(Pathway , AAChange.refGene)

write_tsv(df_annot , '0.data/Sample.DNArepair.tsv', na = '')
write_tsv(df_annot %>% filter(Sample %in% our.samples ), '0.data/OurSample.DNArepair.xls', na = '')


# write_tsv(df_annot , '../Results_new202301/Sample_statistic/Sample.DNArepair.tsv', na = '')
# write_tsv(df_annot %>% filter(Sample %in% our.samples ), '../Results_new202301/Sample_statistic/OurSample.DNArepair.xls', na = '')

# how many somatic SNVs in IGHV regions
df_ighv = read_tsv('1_Signature/Combine/Matrix/input/all.wgs.txt') %>% 
  #filter(chrom  == 'chr14', pos_start > 106411066, pos_start < 107283226) %>% 
  filter(chrom  == 'chr14', pos_start > 106405840, pos_start < 107288051) %>% # check the IGHV position in Ig_segment.hg19.txt
  filter(mut_type != 'INDEL' ) %>% 
  #cSplict('ref',sep = '', direction = 'long') %>% 
  group_by(Sample) %>% summarise(N_mutation_IGHV = sum(mut_type == 'SNV') + 2*sum(mut_type == 'DINUC') )
write_tsv(df_ighv, 'Results/Sample_statistic/N_mutations_IGHV.tsv')

phe = read_tsv('0.data/Sample_clinicInfo.tsv') %>% 
  left_join(df_annot) %>% 
  left_join(df_ighv) %>%
  mutate(N_mutation_IGHV = if_else(is.na(N_mutation_IGHV) , as.integer(0), N_mutation_IGHV))
 
write_tsv(phe, '0.data/Sample_clinicInfo.addDNArepair.tsv', na = '')
