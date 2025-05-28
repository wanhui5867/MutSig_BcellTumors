
### remove kataegis with high APOBEC pct.
#APOBEC with high C/G mutation at the TCN motif
library(tidyverse)
library(Biostrings) # getSeq
library(BSgenome.Hsapiens.UCSC.hg19)


#order_disease = c('MCL','BL', 'FL', 'DLBCL', 'un-BNHL','CLL')
#setwd("/Users/huwan/OneDrive - KI.SE/Mac/Project/202206_Kataegis//")


df_kat_all = read_tsv('2_Kataegis/kataegis_sep.wgs.txt') %>% mutate(ID = as.character(ID)) %>% 
  bind_rows(read_tsv('~/OneDrive - KI.SE/Mac/Project/Kataegis/0.data/2013_breast_cancer/Matrix/2013_breast_kataegis.wgs.txt') ) # Read kataegis from breast cancer as a bechmark



df_kat = df_kat_all %>% filter(!is.na(ref), !is.na(alt)) %>% rowwise() %>% 
  mutate(Content = toupper(getSeq(Hsapiens,chrom, pos_start-1, pos_start+1)),
         Content5 = toupper(getSeq(Hsapiens,chrom, pos_start-2, pos_start+2)),
         type = if_else(ref %in% c('A','G'), toupper(complement(DNAString(str_c(ref, alt, sep = '-')))) , str_c(ref, alt, sep = '-')),
         subtype = if_else(ref %in% c('A','G'), toupper(reverseComplement(DNAString(Content))), Content),
         Content5_complement = if_else(ref %in% c('A','G'), toupper(reverseComplement(DNAString(Content5))), Content5),
         MutationType = str_c(str_sub(subtype, 1,1), '[', str_replace(type, '-', '>' ), ']', str_sub(subtype, -1,-1)))

#write_tsv(df_kat, 'kataegis_sep.wgs.subtype.txt')

write_tsv(df_kat, '2_Kataegis/kataegis_sep.wgs.subtype.txt')
#str_match(df_kat$MutationType, 'T\\[C.*' )


# calculat TCN percent
df_kat_apo_pct = df_kat %>% group_by(MutationType, Project, kataegis = Sample ) %>%  
  summarize(value = n()) %>% 
  mutate(APOBEC = if_else(str_detect(MutationType, 'T\\[C.*' ), 'TCN', 'other') ) %>%
  ungroup() %>% group_by(kataegis, APOBEC) %>%
  summarise(N_mut = sum(value)) %>%
  mutate(pct_mut = N_mut/sum(N_mut) * 100) %>%
  filter(APOBEC == 'TCN') %>% 
  mutate(Sample = str_remove_all(kataegis, '\\|.*')) %>% 
  left_join(df_kat %>% select(Project, Sample) %>% mutate(Sample = str_remove_all(Sample, '\\|.*')) %>% distinct())


# setting breast cancer' outliers as cutoff
outliers_breast = ceiling(max(boxplot.stats(df_kat_apo_pct %>% filter(Project == '2013_breast_kataegis') %>% pull(pct_mut))$out))

pb_apobec <- ggpubr::ggviolin(df_kat_apo_pct  , x = 'Project', y = 'pct_mut', ylab = 'percentage of TCN mutations', 
                      #add.params = list(point.size = 0.1),
          #order = c(order_disease, "Breast"),
          add = 'boxplot') + 
  geom_hline(yintercept = outliers_breast , color = 'red')

ggsave(pb_apobec, filename = '2_Kataegis/apobec_dectect_violine.pdf', w = 6, h = 4)
#ggsave(pb_apobec, filename = '2_Kataegis/all_sep_kataegis/apobec_dectect_violine.nobox.pdf', w = 6, h = 4)




# save apobec's  kataegis
dir.create('2_Kataegis/apobec_kataegis/')
dir.create('2_Kataegis/apobec_kataegis/Matrix/')
k_apobec_list = df_kat_apo_pct %>% filter(pct_mut >= outliers_breast, Project != '2013_breast_kataegis') %>% pull(kataegis)
write_tsv( df_kat_all %>% filter(Sample %in% k_apobec_list), file = '2_Kataegis/apobec_kataegis/Matrix/kataegis_apobec.wgs.txt')

# list samples who have kataegis
df_apobec = df_kat_all %>% filter(Sample %in% k_apobec_list) %>% 
  mutate(kataegis_id = Sample, Sample = str_remove(Sample, '\\|.*$')) %>% 
  group_by(Sample) %>% summarise(N = length(unique(kataegis_id)))
write_tsv(df_apobec, '2_Kataegis/apobec_kataegis/Summary_samples.apobecKataegis.tsv')


df_mut_apobec_cplt = df_kat %>% filter(Sample %in% k_apobec_list)
write_tsv( df_mut_apobec_cplt, file = '2_Kataegis/apobec_kataegis/Matrix/kataegis_apobec.wgs.subtype.txt')

sig_apobec = df_mut_apobec_cplt %>% group_by(MutationType) %>% 
  summarise(APOBEC = n()) %>% 
  mutate(APOBEC = APOBEC/sum(APOBEC)) 
write_tsv( sig_apobec, file = '2_Kataegis/apobec_kataegis/Matrix/kataegis_apobec_SBS96_Signatures.txt')


# save the kataegis without apobec
dir.create('2_Kataegis/rm_apobec_kataegis/')
dir.create('2_Kataegis/rm_apobec_kataegis/Matrix/')
write_tsv( df_kat_all  %>% filter(!Sample %in% k_apobec_list,  Project != '2013_breast_kataegis'),
           file = '2_Kataegis/rm_apobec_kataegis/Matrix/kataegis_sep.wgs.txt')
write_tsv( df_kat  %>% filter(!Sample %in% k_apobec_list,  Project != '2013_breast_kataegis'),
           file = '2_Kataegis/rm_apobec_kataegis/kataegis_sep.wgs.subtype.txt')



# check number of kataegis 
df_kat_all  %>% filter( Project != '2013_breast_kataegis') %>% pull(Sample) %>% unique() %>% length()
df_kat_all  %>% filter(!Sample %in% k_apobec_list,  Project != '2013_breast_kataegis') %>% pull(Sample) %>% unique() %>% length()
df_kat_all  %>% filter(Sample %in% k_apobec_list) %>% pull(Sample) %>% unique() %>% length()

# check mutations number of kataegis
df_kat_all  %>% filter( Project != '2013_breast_kataegis') %>% nrow()
df_kat_all  %>% filter(!Sample %in% k_apobec_list,  Project != '2013_breast_kataegis') %>% nrow()
df_kat_all  %>% filter(Sample %in% k_apobec_list) %>% nrow()

