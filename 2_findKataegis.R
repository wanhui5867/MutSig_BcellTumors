### This program is to plot inter-mutational distance && dectect kataegis 
### Flow: 1.Calculate inter-mutational distance(IMD) ;
###       2.Use IMD to plot rainfall;
###       3.Calculate abnormal distance line(ADL);
###       4.Shower1:Calculate adjacent mutations and fisher test;
###       5.Shower2:Merge adjacent shower and fisher test & summarise & plot
###       6.zooming IGH kataegis.
### Reference:  2014 Cell: B Cell Sup3.er-Enhancers and Regulatory Clusters
###             Recruit AID Tumorigenic Activity
### hui.wan-------  -----------#
### V1: 2020-05-23

library(optparse)
opt_list  = list(
  make_option(opt_str = c('-i','--inputdir'), action = 'store', 
              help = 'SNV files path'),
  make_option(opt_str = c('-o','--outdir'), action = 'store',
              help = 'output files dirname')
)
parser = OptionParser(usage = "usage:%prog[options]", option_list = opt_list)
args = parse_args(parser, positional_arguments = T)
opt = args$options
files = args$args


library(tidyverse)
library(data.table)
library(gtools)
library(egg) #set_panel_size
library(tibbletime)  ## rollify
#library(GenomicVis)          ###kataegis one-step
#library(VariantAnnotation)   ###signature + SVs + CNV
#library(ClusteredMutations)  ###learn kataegis step-wise


###--------- GOLABLE Define: data load and functions ----------

chr24 = str_c('chr', c(1:22, 'X', 'Y'))
ref_hg19 = read_tsv('~/project/Kataegis/_Ref/hg19.chrom.sizes', col_names = c('Chr', 'size')) %>% filter(Chr %in% chr24) %>% 
  arrange(match(Chr, chr24)) %>% 
  mutate(ChrPos = cumsum(size),
         ChrDif = c(0, cumsum(size)[-24])) 

df_igh = read_tsv('~/project/Kataegis/_Ref/IGH.region.forUCSV.tsv', skip = 3, col_names = F)


col_sig6 = structure(c('skyblue', "black", "firebrick1", "gray", "darkolivegreen2", "pink"),
                     names =  c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))  # define colors for 6 mutation types



num_mut = 10     #shower: adjacent mutations' number
reg_leg = 10000  #shower: region length of adjacent mutations
p_cut = 0.0001   #shower: fisher test cutoff

##FUN: batch read file
read.single <- function(filename,  path){
  samplename = str_remove(filename, '_SNV.txt')
  print(paste0('Reading ', samplename))
  
  read_tsv(str_c(path, filename,sep='/')) %>% 
    select(Chr, Position = Start, Ref , Alt ) %>% 
    mutate(Sample = samplename)
}


##FUN: convert base A/G to T/C in rainplot
reverse <- function(base){
  if(base=='A'){
    rbase = 'T'
  } else if(base == 'T'){
    rbase = 'A'
  } else if(base == 'C'){
    rbase = 'G'
  } else if(base == 'G'){
    rbase = 'C'
  } else {
    rbase = 'N'
  }
  return(rbase)
}

##FUN: roll in shower window size : create list of adjacent mutation position and inter-distance
pos_list_roll <- rollify(~c(.), window = num_mut , unlist = F)
imd_list_roll <- rollify(~diff(.), window = num_mut , unlist = F)

##FUN: fisher test: background_less_ADL      |  interested_less_ADL
##                  backgroud_more_equl_ADL  |  interested_more_equl_ADL
my.fisher.test.p.value <- function(...){
  obj <- try(fisher.test(matrix(c(..1, ..2,..3,..4), nrow = 2)), silent = T)
  if (is(obj,"try-error")) return(NA) else return(obj$p.value)
}


##------- 0.read data: NEED MODIFY ----


# read from forSigprofile file
#myfile = '~/project//1_Signature/Matrix/all.wgs.txt'
myfile = opt$inputdir
df_forSig = read_tsv(myfile, col_names = T)  

# DINUC split to two SNVs
df_dinuc1 = df_forSig %>% 
  filter(mut_type == 'DINUC') %>% 
  mutate(mut_type = 'SNV', ref = str_sub(ref, 1, 1), alt = str_sub(alt, 1,1)) %>% 
  filter(ref != alt)

df_dinuc2 = df_forSig %>% 
  filter(mut_type == 'DINUC') %>% 
  mutate(mut_type = 'SNV', ref = str_sub(ref, 2, 2), alt = str_sub(alt, 2, 2),
         pos_start = pos_start+1) %>% 
  filter(ref != alt)

# prepare each samples'  SNVs  ---
df_all = df_forSig %>% 
  filter(mut_type == 'SNV') %>% 
  bind_rows(df_dinuc1, df_dinuc2) %>% 
 # mutate(Chr = str_c('chr', chrom)) %>%  # use if chrom has no 'chr' prefix
  select(Chr = chrom , Position = pos_start, Ref = ref, Alt = alt, Sample)

  

n_sample = length(unique(df_all$Sample))

##---------1. calculate IMD =  mutation_position - next_mut_posi  ----- 
#dir1 = './2_Kataegis//'
dir1 = opt$outdir
dir.create(dir1)

df_imd = df_all %>%  
  #---sort---#
  arrange(match(Chr, chr24), Position) %>% 
  #---IMD---#
  group_by(Sample, Chr) %>% #all samples
  mutate(distance = c(diff(Position), NA)) %>% 
  ungroup() %>% 
  filter(!is.na(distance)) %>% 
  #---covert---#
  rowwise() %>%  
  mutate(Type = if_else(Ref %in% c('A','G'), 
                        str_c(reverse(Ref),reverse(Alt), sep = '>'), 
                        str_c(Ref, Alt, sep = '>')),
         log10Dis = round(log10(distance),1)) %>%
  ungroup() %>% 
  #---arrange chromsome for plot---#
  left_join(ref_hg19 , by = 'Chr') %>% 
  mutate(GenomicPosition = Position + ChrDif,
         chrmax = ChrPos) %>% 
  select(-size, -ChrPos, - ChrDif)

write_tsv(df_imd, paste0(dir1,'1.all_IMD.tsv'))

##----------2.plot rainfall ---------#
# p <- ggplot(df_imd, aes(x = as.numeric(GenomicPosition), y = as.numeric(log10Dis), color = Type)) +
#   geom_point(size =.1) + 
#   geom_vline(aes(xintercept = as.numeric(chrmax)), size = .3, linetype = "dotted", color = "grey50") + 
#   scale_color_manual(values = col_sig6) +
#   facet_wrap(~ Sample, ncol = 10) +  ## all_samples
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), panel.grid = element_blank()) + 
#   labs(x = 'Genomic Position', y = "log10 Intermutation Distance (bp)")
# pl <- grid.arrange(set_panel_size(p,  width=unit(3,"cm"), height = unit(1, "cm")))
# 
# ggsave(pl, filename = paste0(dir1,'2.rainfall_allsamples.pdf'), w = 35, h =  ceiling(n_sample/10)*1.6 + 3, units = 'cm'  )


###--------------3. calculate ADL = averange(IMD)/10 per sample ------ 

df_adl = df_imd %>% group_by(Sample) %>% 
  summarise(ADL = mean(distance, na.rm = T)/ num_mut,
            All_N_less_ADL = sum(distance < ADL, na.rm =T),
            All_N_more_ADL = sum(distance >= ADL, na.rm =T))

###-------------4. make shower: 10 adjacent mutations to a raw shower ----------

df_shower_raw = df_all %>%
  #---sort by position----#
  arrange(Sample, match(Chr, chr24), Position) %>% 
  #---filter mutation number/per chr/ < adjacent number ---#
  group_by(Sample, Chr) %>% #all samples 
  mutate(mut_n = n()) %>% 
  filter(mut_n > num_mut) %>% 
  #---10 mutations to one shower---#
  mutate(Pos_end = Position,
         Pos_start = lag(Position, num_mut -1),
         list_pos = pos_list_roll(Position),
         distance = Pos_end - Pos_start + 1) %>% 
  ungroup() %>% 
  #---filter shower distance < defined region length---#
  filter(!is.na(distance), distance <= reg_leg) %>% 
  select(- mut_n, -Position)

###-----------5. merge adjacent showers and fisher test ------#
df_shower_merge = df_shower_raw %>% 
  #--- merge shower: next_position_start <= previous_position_end  -----#
  group_by(Sample, Chr) %>% 
  mutate(index = c(0, cumsum( lead(Pos_start) > Pos_end)[-n()])) %>% 
  group_by(Sample, Chr, index) %>% 
  summarize(Pos_start = first(Pos_start),
            Pos_end = last(Pos_end),
            distance = Pos_end - Pos_start + 1,
            list_pos = list(unique(unlist(list_pos))), 
            list_dis = list(diff(sort(unique(unlist(list_pos)))))) %>% 
  ungroup() %>% rowwise() %>% 
  #--- add ADL info ---#
  left_join(df_adl)  %>%
  #--- fisher test ---#
  mutate(mutation_number = length(unlist(list_dis)) + 1,
         N_less_ADL = sum(unlist(list_dis) < ADL),
         N_more_ADL = sum(unlist(list_dis) >= ADL),
         p = my.fisher.test.p.value(All_N_less_ADL, All_N_more_ADL,N_less_ADL, N_more_ADL))  %>% 
  filter(p <= p_cut)

write_tsv(df_shower_merge %>% select(-list_dis, -list_pos), paste0(dir1, '3.kataegis.shower.tsv'))



##----------- 6. output and plot kataegis region --------------
df_kataegis = df_shower_merge %>% 
  select(Sample, Chr, index, list_pos, mutation_number) %>%
  unnest(list_pos) %>% 
  rename(Position = list_pos) %>% 
  right_join(df_imd)  %>% 
  mutate(kataegis = if_else(is.na(index), 0, 1),
         Type = ifelse(kataegis == 1, Type, NA))

df_ktg_Mut = df_kataegis %>% filter(kataegis == 1) %>% 
  select(Sample:Alt)

write_tsv(df_ktg_Mut, paste0(dir1, '4.kataegis.MutType.tsv'))

# p <- ggplot(df_kataegis, aes(x = as.numeric(GenomicPosition), y = as.numeric(log10Dis), color = Type)) +
#   geom_point(size =.1) + 
#   geom_vline(aes(xintercept = as.numeric(chrmax)), size = .3, linetype = "dotted", color = "grey50") + 
#   scale_color_manual(values = col_sig6, na.value = 'grey75') +
#   scale_alpha_manual(values = 1, na.value = 0.3) + 
#   facet_wrap(~ Sample, ncol = 10) +  ## all_samples
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), panel.grid = element_blank()) + 
#   labs(x = 'Genomic Position', y = "log10 Intermutation Distance (bp)") 
# pl <- grid.arrange(set_panel_size(p,  width=unit(3,"cm"), height = unit(1,"cm")))
# 
# ggsave(pl, filename = paste0(dir1,'4.rainfall_kataegis.pdf'), w = 35, h = ceiling(n_sample/10)*1.6 + 3 , units = 'cm')

###-------------- 7. find kataegis in IGH -------##

find_igh <- function(position){
  region = df_igh %>% filter(between(position, X2, X3)) %>% pull(X4)
  ifelse(length(region) == 0, 'Intergenic',region)
}

df_chr14 = df_kataegis %>% 
  filter(Chr == 'chr14', kataegis != 0, 
         between(Position, min(df_igh$X2, df_igh$X3), max(df_igh$X2, df_igh$X3)))  %>% 
  select(Chr, Sample, Position, Ref, Alt ) %>% 
  rowwise() %>% 
  mutate(IGHregion = find_igh(Position))

# statistic ---
df_chr14_sum = df_chr14 %>% group_by(IGHregion) %>% 
  summarise(n = n())
write_tsv(df_chr14_sum, paste0(dir1, '5.IGHkataegis.statistic.tsv'))

## prepare for UCSC ----
df_chr14_o = df_chr14 %>% 
  select(Chr, start = Position, end = Position, IGHregion, start2 = Position, end2 = Position, Color = Sample,  Ref, Alt)  %>% 
  mutate(chr =  Chr,  X1 = 0, X2 = '-') %>% 
  select(chr, start:IGHregion, X1:X2, start2:Color, Ref:Alt) %>% 
  add_row(chr = c('browser position chr14:106053274-107288051',
                  'browser hide all;',
                  'track name="IGH kataegis region" description="IGH regions - CSR" visibility=2 itemRgb="On"'),
          .before = 1  )
write_tsv(df_chr14_o, paste0(dir1,'5.IGH.zooming.tsv'), col_names = F, na = '')

df_IGH_forUSCS = df_igh %>% filter(X4 %in% df_chr14_sum$IGHregion)  %>% 
  add_row(X1 = c('browser position chr14:106053274-107288051',
                 'browser hide all;',
                 'track name="IGH region" description="IGH regions " visibility=2 itemRgb="On"'),
          .before = 1  )
write_tsv(df_IGH_forUSCS, paste0(dir1,'5.IGH.select.forUCSC.tsv'), col_names = F, na='')


## prepare for forSigProfiler
df_kataegis_forSigProfile = df_ktg_Mut %>% 
  left_join(df_forSig %>% select(Project, Sample) %>% distinct()) %>% 
  mutate( Sample = str_c(Sample , Chr, index, sep = '|') ) %>% 
  transmute(Project = Project, Sample = Sample, ID = index, 
            Genome = 'GRCh37', mut_type = 'SNV',  chrom = Chr, 
            pos_start = Position, pos_end = Position, 
            ref = Ref, alt = Alt, Type = 'SOMATIC') 

write_tsv(df_kataegis_forSigProfile, paste0(dir1,'kataegis_sep.wgs.txt'))  

