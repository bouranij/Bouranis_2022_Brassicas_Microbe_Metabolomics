library(tidyverse)
library(magrittr)
library(here)
here::i_am('./ChemRichAnalysis.R')

#Source ChemRich...hopefully one day they actually make into a package...
source("https://raw.githubusercontent.com/barupal/ChemRICH/master/chemrich_minimum_analysis.R ")
library(RCurl)
library(pacman)
pacman::p_load(GGally)
pacman::p_load(DT)
pacman::p_load(RCurl)
pacman::p_load(RJSONIO)
pacman::p_load(ape)
pacman::p_load(devEMF)
pacman::p_load(dynamicTreeCut)
pacman::p_load(extrafont)
pacman::p_load(ggplot2)
pacman::p_load(ggpubr)
pacman::p_load(ggrepel)
pacman::p_load(grid)
pacman::p_load(htmlwidgets)
pacman::p_load(igraph)
pacman::p_load(magrittr)
pacman::p_load(network)
pacman::p_load(officer)
pacman::p_load(openxlsx)
pacman::p_load(phytools)
pacman::p_load(plotly)
pacman::p_load(plotrix)
pacman::p_load(rcdk)
pacman::p_load(readxl)
pacman::p_load(rvg)
pacman::p_load(sna)
pacman::p_load(visNetwork)

#Step 1: Read in the mgf file that was used for Canopus and parse it out extract the necessary information
l <- readLines(here('Data/ChemRich_Data/PosClean.mgf'))
ti <- which(str_detect(l, 'FEATURE_ID'))
ck <- data.frame()
for(i in ti){
  if(str_detect(l[i + 3], 'ION')){
    tmp <- data.frame(canopus_id = str_split(l[i], '=')[[1]][2],
                      ionMass = str_split(l[i + 2], '=')[[1]][2],
                      RTseconds = str_split(l[i + 6], '=')[[1]][2],
                      progenesis_id = str_split(l[i + 5], '=')[[1]][2])
  }else{
    tmp <- data.frame(canopus_id = str_split(l[i], '=')[[1]][2],
                      ionMass = str_split(l[i + 2], '=')[[1]][2],
                      RTseconds = str_split(l[i + 5], '=')[[1]][2],
                      progenesis_id = str_split(l[i + 4], '=')[[1]][2])
    
  }
  ck <- rbind(ck, tmp)
}

canopusout <- read_tsv(here('Data/ChemRich_Data/compound_identifications.tsv')) %>%
  dplyr::select(ionMass, retentionTimeInSeconds, id)

#Now split the canopus_id to give a single value to match what we extracted from the mgf:
csplit <- canopusout %>%
  mutate(splitid = str_split(id, '_')) %>%
  mutate(idmgf = map_chr(splitid, function(x) x[3]))

#All match to the mgf!
sum(csplit$idmgf %in% ck$canopus_id)

#Create the final key linking CanopusIDs back to ProgenesisIDs
finalkey <- inner_join(csplit, ck, by = c('idmgf' = 'canopus_id')) %>%
  dplyr::select(id, idmgf, progenesis_id) %>%
  rename('canopus_id' = id)

#Create the sample type key:
mdata_pos <- read_csv(here('Data/Metabolomics/metadata_pos.csv')) %>%
  mutate(group = ifelse(group == 'A', 'C_Type', 'E_type'))

#Load in the canopus classifications
canopus_classes <- read_tsv(here('Data/ChemRich_Data/canopus_summary.tsv'))

#Load in progenesis data and calculate foldchanges
NCvBrocFull <- read_csv(here('Data/ChemRich_Data/NCvBroc_Pos.csv'), skip = 2) 
NCvBroc <- NCvBrocFull %>%
  #Select columns
  dplyr::select(1,8,9, 16:35) %>%
  #Rename to only those that are necessary
  rename('pvalue' = 'Anova (p)', 'qvalue' = 'q Value') %>%
  #Make tidy
  pivot_longer(cols = starts_with('pos'), names_to = 'sample') %>%
  #Fix names
  mutate(namefix = str_split(sample, '\\.\\.\\.')) %>%
  mutate(sample = map_chr(namefix, function(x) x[1])) %>%
  dplyr::select(-namefix) %>%
  #Join metadata
  left_join(mdata_pos) %>%
  #Group to necessary
  group_by(Compound, pvalue, qvalue, treatment) %>%
  #Grab the median value
  summarise(medianval = median(value) + 0.001) %>%
  ungroup() %>%
  #Make untidy for manipulation
  pivot_wider(names_from = 'treatment', values_from = 'medianval') %>%
  #Calculate the fold change
  mutate(foldchange = broc/no_veg) %>%
  #Join in the link between Progenesis ID and Canopus ID
  inner_join(., finalkey, by = c('Compound' = 'progenesis_id')) %>%
  #Join in Canopus information
  inner_join(., canopus_classes, by = c('canopus_id' = 'name'))

NCvBrusFull <- read_csv(here('Data/ChemRich_Data/NCvBrus_Pos.csv'), skip = 2) 
NCvBrus <- NCvBrusFull %>%
  dplyr::select(1,8,9, 16:35) %>%
  rename('pvalue' = 'Anova (p)', 'qvalue' = 'q Value') %>%
  pivot_longer(cols = starts_with('pos'), names_to = 'sample') %>%
  mutate(namefix = str_split(sample, '\\.\\.\\.')) %>%
  mutate(sample = map_chr(namefix, function(x) x[1])) %>%
  dplyr::select(-namefix) %>%
  left_join(mdata_pos) %>%
  group_by(Compound, pvalue, qvalue, treatment) %>%
  summarise(medianval = median(value) + 0.001) %>%
  ungroup() %>%
  pivot_wider(names_from = 'treatment', values_from = 'medianval') %>%
  mutate(foldchange = brus/no_veg) %>%
  inner_join(., finalkey, by = c('Compound' = 'progenesis_id')) %>%
  inner_join(., canopus_classes, by = c('canopus_id' = 'name'))

NCvComboFull <- read_csv(here('Data/ChemRich_Data/NCvCombo_Pos.csv'), skip = 2) 
NCvCombo <- NCvComboFull %>%
  dplyr::select(1,8,9, 16:35) %>%
  rename('pvalue' = 'Anova (p)', 'qvalue' = 'q Value') %>%
  pivot_longer(cols = starts_with('pos'), names_to = 'sample') %>%
  mutate(namefix = str_split(sample, '\\.\\.\\.')) %>%
  mutate(sample = map_chr(namefix, function(x) x[1])) %>%
  dplyr::select(-namefix) %>%
  left_join(mdata_pos) %>%
  group_by(Compound, pvalue, qvalue, treatment) %>%
  summarise(medianval = median(value) + 0.001) %>%
  ungroup() %>%
  pivot_wider(names_from = 'treatment', values_from = 'medianval') %>%
  mutate(foldchange = combo/no_veg)  %>%
  inner_join(., finalkey, by = c('Compound' = 'progenesis_id')) %>%
  inner_join(., canopus_classes, by = c('canopus_id' = 'name'))

#Pull all the classnames
classnames <- canopus_classes %>%
  dplyr::select('level 5') %>%
  distinct() %>%
  inset('order', value = sample(1:nrow(.)))

#Trim the file to match the needs of ChemRich
#I recommend making and descending into sub-directories before runing this
#ChemRICH overwrite the results each time otherwise
NCvBrocOut <- NCvBroc %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) 
#Write the excel output - ChemRich's desired input
write.xlsx(NCvBrocOut, './NCvBroc_CR.xlsx')
#Run chemrich
run_chemrich_basic(inputfile = './NCvBroc_CR.xlsx')

#Remove the peptides, unannotated compounds, and amino acids for simplicity
NCvBrocOut_nopep <- NCvBroc %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) %>%
  filter(!set %in% c('Peptides', NA, 'Amino acids and derivatives'))

write.xlsx(NCvBrocOut_nopep, './NCvBroc_CR_nopep.xlsx')
run_chemrich_basic(inputfile = './NCvBroc_CR_nopep.xlsx')

#Repeat for each comparison
NCvBrusOut <- NCvBrus %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) 

write.xlsx(NCvBrusOut, './NCvBrus_CR.xlsx')
run_chemrich_basic(inputfile = './NCvBrus_CR.xlsx')

NCvBrusOut_nopep <- NCvBrus %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) %>%
  filter(!set %in% c('Peptides', NA, 'Amino acids and derivatives'))

write.xlsx(NCvBrusOut_nopep, './NCvBrus_CR_nopep.xlsx')
run_chemrich_basic(inputfile = './NCvBrus_CR_nopep.xlsx')

NCvComboOut <- NCvCombo %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) 
write.xlsx(NCvComboOut, './NCvCombo_CR.xlsx')
run_chemrich_basic(inputfile = './NCvCombo_CR.xlsx')

NCvComboOut_nopep <- NCvCombo %>%
  filter(qvalue <= 0.05) %>%
  left_join(., classnames) %>%
  dplyr::select(Compound, order, pvalue, foldchange, `level 5`) %>%
  rename('compound_name' = Compound, 'effect_size' = foldchange, 'set' = `level 5`) %>%
  filter(!set %in% c('Peptides', NA, 'Amino acids and derivatives'))

write.xlsx(NCvComboOut_nopep, './NCvCombo_CR_nopep.xlsx')
run_chemrich_basic(inputfile = './NCvCombo_CR_nopep.xlsx')

