#Figures for PoopSoup2 Manuscript
#Script by John "Yanni" Bouranis

# Environment -------------------------------------------------------------
library(phyloseq)
library(cowplot)
library(tidyverse)
library(magrittr)
library(cowplot)
library(here)
here::i_am('./FigureScript.R')

#Set the palette for the figures
pal <- RColorBrewer::brewer.pal(4, 'Dark2')
names(pal) <- c('broc', 'brus', 'combo', 'NC')


## Helper Functions --------------------------------------------------------
#Log Helper
log_helper <- function(x, min.val){
  log2((x + sqrt(x^2 + min.val^2))/2)
}

#Pareto Scaling:
#Pareto helper
PS_helper <- function(x){
  (x - mean(x))/sqrt(sd(x, na.rm = T))
}	

#Transformation Functions:
#Log Scaling:
log_transform <- function(mtb){
  mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
  min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
  mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
  return(mtb_log_trans)
}

#Pareto Scaling:
pareto_scale <- function(mtb){
  mtb_scaled <- apply(mtb, 2, PS_helper) 
  return(mtb_scaled)
}

#Centered Log Ratio
#Geometric Mean
geoMean <- function(x) exp(mean(log(x)))
#CLR helper function
centerHelper <- function(x) log(x/geoMean(x))
#Actual CLR function
clr <- function(x){
  #Replace zero counts with a count of 1
  nz <- data.frame(apply(x, 2, function(y) replace(y, which(y == 0), 1)))
  #Complete the clr transformation
  data.frame(apply(nz, 2, centerHelper))
}

## Microbiome --------------------------------------------------------------
#Load in teh data
asvtab <- readRDS(here('Data/Microbiome/asv_tab.RDS'))
taxatab <- readRDS(here('Data/Microbiome/tax_tab.RDS'))
metadata <- read.csv(here('Data/Microbiome/microbiome_metadata.csv'))

#Fix the metadata
metadata %<>% mutate(group = ifelse(group == 'A', 'C_Type', 'E_type'))

#Factorize and set the levels on the metadata
metadata$treatment %<>% factor(levels = c('fecal_stock', 'no_veg', 'broc', 'brus', 'combo', 'control_digest'))
metadata$fecal_sample %<>% factor(levels = c('T5631','T5632','T6260','T6291','T4669','T1995','T5627','T5717','T5854','T6382')) 

rownames(metadata) <- metadata$sample
rownames(asvtab) <-metadata$sample
ps_raw <- phyloseq(otu_table(asvtab, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxatab))

#Give arbitrary names to the taxa as opposed to keeping as just DNA-sequences which identify them
taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))

#Fill in missing genus names:
renames <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Genus'])])
taxdf <- tax_table(ps_raw)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(ps_raw)[renames, 'Genus'] <- renamed_genus

#Remove the control digests, these are not relevant to our analysis
ps_raw <- ps_raw %>% subset_samples(treatment != 'control_digest')
#Agglomerate to the genus level
ps_genera <- ps_raw %>% tax_glom(taxrank = "Genus")
#Remove taxa not seen more than 3 times in at least 20% of the samples
ps_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE)
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) )
#Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)

ps_final <- ps %>% subset_samples(treatment != 'fecal_stock')
#Create the count data
ps_final_c <- ps_counts %>% 
  filter_taxa(function(x) mean(x / sum(x)) > 1e-5, TRUE) %>%
  subset_samples(treatment != 'fecal_stock')

microdata <- cbind(sample = data.frame(sample_data(ps_final))$metabolomics_neg_sample, data.frame(otu_table(ps_final))) %>%
  mutate(sample = gsub('neg', 'ms', sample))

microdata_c <- cbind(sample = data.frame(sample_data(ps_final_c))$metabolomics_neg_sample, data.frame(otu_table(ps_final_c))) %>%
  mutate(sample = gsub('neg', 'ms', sample))

#Generate the data to plot
fps <- ps_final %>%
  #Agglomerate taxa t the family level
  tax_glom('Family') %>%
  #Make tidy
  psmelt() %>%
  #Adjust treatment to a character so it can be modified
  modify_at('treatment', as.character) %>%
  #Modify no_veg to be NC 
  mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment)) %>%
  #Combine treatment and sample together for plotting help
  mutate(samptreat = paste0(treatment, '_', fecal_sample)) %>%
  arrange(Family, fecal_sample)

#Relevel to make the plot nice and clean
goodlevels <- c(c(fps$samptreat %>% unique() %>% grep('NC',., value = T)), 
                c(fps$samptreat %>% unique() %>% grep('broc',., value = T)),
                c(fps$samptreat %>% unique() %>% grep('brus',., value = T)),
                c(fps$samptreat %>% unique() %>% grep('combo',., value = T)))
#Set it as a factor
fps$samptreat %<>% factor(., levels = goodlevels)
#Split between the treatment and the fecal sample and pull the fecal sample for plotting
fss <- sapply(goodlevels, function(x) str_split(x, '_')[[1]][2])

#Build the insane palette that will be used to color in this figure
c20 <- c(
  "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "palegreen2","dodgerblue2", "gray70", "khaki2",
  "yellow1", "green1","maroon", "blue1", "orchid1", "steelblue4",
  "darkturquoise", "skyblue2", "yellow3","deeppink1" ,
  "#CAB2D6", # lt purple
  "#FDBF6F" # lt orange
)
#Make a pie chart of the colors to make sure theyre not too wild
pie(rep(1, 20), col = c20)


# Figure 1: Save 3200x1769 ------------------------------------------------
#Plot the relative abundances of each microbe!
h1 <- ggplot(fps,aes(y = Abundance, x = samptreat, fill = Family, group = treatment)) + 
  #Set the borders of the bars to be black
  geom_bar(stat = 'identity', color = 'black') + 
  #Add the label of the fecal samples
  scale_x_discrete(labels = fss) +
  #Cowplot to make pretty
  theme_cowplot() +
  #Adjust axes to scale up properly
  theme(axis.text.x = element_text(angle = 90, size = 28, vjust = 0.5),
        axis.text.y = element_text(size = 28, vjust = 0.8),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        axis.title.x = element_text(size = 35, vjust = 5), 
        axis.title.y = element_text(size = 35 , vjust = 1),
        legend.key.size = unit(3, 'line')) +
  #Make the y-axis percentages
  scale_y_continuous(expand = c(0,0), labels = scales::label_percent()) +
  #Adjust colors
  scale_fill_manual(values = c20) +
  labs(x = '\n Subject ID',
       y = 'Microbial Abundance') 

#Get the legend alone
leg <- plot_grid(NULL, get_legend(h1), ncol = 2, rel_widths = c(01,10))

#Make the top layer to sit ontop of the bar plots
h2 <- ggplot(fps) + 
  geom_bar(aes(x = samptreat, y = 1, fill = treatment), stat = 'identity', width = 1) +
  theme_void() +
  #Remove the legend
  theme(legend.position = 'none') + 
  #Add the text
  annotate('text', x = 5.5, y = 10, label = 'Negative Control', size = 12) +
  annotate('text', x = 15.5, y = 10, label = 'Broccoli', size = 12) +
  annotate('text', x = 25.5, y = 10, label = 'Brussels', size = 12) +
  annotate('text', x = 35.5, y = 10, label = 'Combo', size = 12) +
  #Make it pretty using our theme colors
  #scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, name = 'Dark2')) +
  scale_fill_manual(values = pal)

#Remove the legend rom the original plot
h1 <- h1 + theme(legend.position = 'none')
#Stack-em
plot <- plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.7, 15))
#Combine with the legend to make the final plot
plot_final <- plot_grid(plot, leg, ncol = 2, rel_widths = c(100, 25))
plot_final

# Figure 2: Save 3000x2000 ------------------------------------------------
#Run the PCoA
ordBCall <- ordinate(ps_final, method = 'PCoA', distance = 'bray')

#Alternate color names
palmb <- pal
names(palmb) <- c('broc', 'brus', 'combo', 'no_veg')


#Plot the PCoA
mb <- plot_ordination(ps_final, ordBCall, color = 'treatment') +
  stat_ellipse(aes(group = treatment, color = treatment, fill = treatment), geom = 'polygon', alpha = 0.2) +
  geom_point(size = 10) +
  theme_minimal_grid() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = palmb, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  scale_color_manual(values = palmb, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme(axis.text.x = element_text(angle = 90, size = 28, vjust = 0.5),
        axis.text.y = element_text(size = 28, vjust = 0.8),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        axis.title.x = element_text(size = 35), 
        axis.title.y = element_text(size = 35 , vjust = 1)) +
  coord_fixed()

#Load the GLM output from the analysis
load(here('Data/Plot_Data/GLMOutput.RData'))

glmplot <- glmres %>%
  #Filter to only ASVs which are significant for at least 1 comparison
  filter(OTU %in% unique(glmnbout$OTU)) %>%
  #Label significant hits
  mutate(sp = ifelse(padj <= 0.05, 'Sig',
                     ifelse(sign(l2fc) == 1, 'NSP', 'NSN'))) %>%
  #Label fold changes
  mutate(fc = ifelse(sign(l2fc) == 1, 'Positive', 'Negative')) %>%
  #Make names nice
  mutate(test = ifelse(test == 'NVvBroc', 'NC vs Broc',
                       ifelse(test == 'NVvBrus', 'NC vs Brus',
                              ifelse(test == 'NVvCombo', 'NC vs Combo', 
                                     ifelse(test == 'BrocvCombo', 'Broc vs Combo',
                                            ifelse(test == 'BrusvCombo', 'Brus vs Combo', 'Broc vs Brus'))))))
#Relevel factors plotting
glmplot$sp %<>% factor(levels = c('Sig', 'NSP', 'NSN'))
glmplot$fc %<>% factor(levels = c('Positive', 'Negative'))

gplot <- ggplot(glmplot, aes(y = test, x = OTU, fill = fc)) +
  #Add tiles
  geom_tile(color = 'black', position = position_nudge(x = -0.5, y = -0.5), size = 2) +
  #Move the points so theyre in the center of the block
  geom_point(position = position_nudge(x = -0.5, y = -0.5), size = 10, aes(color = sp)) +
  #Make squares
  coord_fixed() +
  #Make pretty with cowplot
  cowplot::theme_cowplot() +
  #Adjust axes to make it look nice
  theme(axis.text.x = element_text(vjust = -0.3, angle = 90, hjust = 1),
        axis.text.y = element_text(vjust = 1.2),
        legend.position = 'none') +
  #Fill in the squares
  scale_fill_manual(values = c('darkorchid', 'goldenrod3', 'grey74')) +
  #Fill in the dots
  scale_color_manual(values = c('black', 'darkorchid', 'goldenrod3'))  +
  #Add labels
  ylab('Comparison') +
  xlab('ASV') +
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 28),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))

#Make dummy plot to start to form the legend
dplot <- glmplot %>%
  filter(sp == 'Sig') %>%
  #Adjust names so they line up correctly
  mutate(test = ifelse(test == 'NVvBroc', 'NC vs Broc',
                       ifelse(test == 'NVvBrus', 'NC vs Brus',
                              ifelse(test == 'NVvComb', 'NC vs Combo', 
                                     ifelse(test == 'BrocvCombo', 'Broc vs Combo',
                                            ifelse(test == 'BrusvCombo', 'Brus vs Combo', 'Broc vs Brus'))))))

#Make the dummy plot
fcs <- ggplot(dplot, aes(y = test, x = OTU)) +
  geom_point(position = position_nudge(x = -0.5, y = -0.5), size = 4, aes(color = sp)) +
  scale_fill_manual(values = c('darkorchid', 'goldenrod3', 'grey74')) +
  scale_color_manual(values = c('black', 'darkorchid', 'goldenrod'),
                     labels = c('Significant (q ≤ 0.05)'))   +
  labs(color = 'Significance')  +
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 34),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))

#Make a plot showing just LFC
fcp <- ggplot(glmplot, aes(y = test, x = OTU, fill = fc)) +
  geom_tile() +
  scale_fill_manual(values = c('darkorchid', 'goldenrod', 'grey74')) +
  labs(fill = 'Fold Change') +
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 34),
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30))


#Grab legends
fillLeg <- get_legend(fcp)
sigLeg <- get_legend(fcs)

#Combine the legends into one image
leg <- plot_grid(fillLeg, NULL, sigLeg, ncol = 1, align = 'hv', rel_heights = c(1,-0.8,1))
#Combine the legend and the plot to make the final product
ahp <-  plot_grid(gplot, NULL, leg, rel_widths = c(10,-.55, 3), nrow = 1)

#Combine to make the final plot
plot_grid(NULL, mb, ahp, ncol = 1, rel_heights = c(0.1,1,1), labels = c('', 'A', 'B'), label_size = 40)


# Figure 3: Save 3000x2000 ------------------------------------------------

#Load in the data:
negugly <- read_csv(here('Data/Metabolomics/neg_metabolome_full.csv'))
#Fix the metadata
mdata_neg <- read_csv(here('Data/Metabolomics/metadata_neg.csv')) %>%
  mutate(group = ifelse(group == 'A', 'C_Type', 'E_type')) %>%
  mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment)) %>%
  modify_at('treatment', factor, levels = c('NC', 'broc', 'brus', 'combo'))

#Clean it up
neg_full <- negugly %>%
  #Remove compounds with a CV > 50 in the QC
  filter(is.na(cv_g_50)) %>%
  #Keep the columns we want, discard the rest
  dplyr::select(compound, starts_with('neg')) %>%
  #Add neg to the end of all compounds so we know where it came from
  mutate(compound = paste0(compound, '_neg')) %>%
  #Make the data tidy
  pivot_longer(starts_with('neg'), names_to = 'sample', values_to = 'intensity') %>%
  #Make the data messy again
  pivot_wider(names_from = 'compound', values_from = 'intensity') %>%
  #Make it completely numeric
  column_to_rownames('sample')

#Correctly order the metadata
mneg <- left_join(data.frame(sample = rownames(neg_full)), mdata_neg)

#Grabbing data for plotting
sinapicData <- neg_full %>%
  #Grab the feature of interest
  dplyr::select(`20.20_223.0608_neg`) %>%
  #Make the rownames into a column
  rownames_to_column('sample') %>%
  #Join in the metadata
  left_join(mneg) %>%
  #Rename to intensity for easier plotting
  rename('intensity' = `20.20_223.0608_neg`)

sinapicPlot <- ggplot(sinapicData, aes(x = treatment, y = intensity, color = treatment)) +
  #Add errorbars
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  #Add bars for mean
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  #Add dots to show individuals
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  #Fix the labels to be capitals
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  #Make pretty
  theme_cowplot() +
  #Add title
  ggtitle('Sinapic Acid') +
  #Drop the legend and center the title
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  #Remove the gap and increase the y-axis to capture all data comfortably
  scale_y_continuous(limits = c(0, max(sinapicData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  #Fix the label names
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL)  +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))
sinapicPlot

IAAData <- neg_full %>%
  dplyr::select(`18.65_204.0665_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `18.65_204.0665_neg`)

IAAPlot <- ggplot(IAAData, aes(x = treatment, y = intensity, color = treatment)) +
  #Add errorbars
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  #Add bars for mean
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  #Add dots to show individuals
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  #Fix the labels to be capitals
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  #Make pretty
  theme_cowplot() +
  #Add title
  ggtitle('Indole Acetic Acid') +
  #Drop the legend and center the title
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  #Remove the gap and increase the y-axis to capture all data comfortably
  scale_y_continuous(limits = c(0, max(IAAData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  #Fix the label names
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

subericData <- neg_full %>%
  dplyr::select(`18.77_174.0891_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `18.77_174.0891_neg`)

subericPlot <- ggplot(subericData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Suberic Acid') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(subericData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

kaempData <- neg_full %>%
  dplyr::select(`12.42_775.1730_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `12.42_775.1730_neg`)

kaempPlot <- ggplot(kaempData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  ggtitle('Kaempferol Glucoside*') +
  theme_cowplot() +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(kaempData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

aksData <- neg_full %>%
  dplyr::select(`22.78_514.2831_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `22.78_514.2831_neg`) %>%
  filter(sample != 'neg_15')

aksPlot <- ggplot(aksData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Alkanesulfonic Acid \n (m/z: 514.2831, RT: 22.78)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(aksData$intensity)+100), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL)  +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))
aksPlot

dcpData <- neg_full %>%
  dplyr::select(`24.68_606.3078_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `24.68_606.3078_neg`) 

dcpPlot <- ggplot(dcpData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Cyclic Depsipeptide \n (m/z: 606.3078, RT: 24.68)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(dcpData$intensity)+100), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL)  +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))
dcpPlot

fig3ae <- plot_grid(subericPlot, sinapicPlot, kaempPlot, dcpPlot, IAAPlot, nrow = 1, labels = c(LETTERS[1:5]), 
          rel_widths = rep(1, 5), label_size = 20)

#Run the edited ChemRichFunction to return the appropriate data
source(here('./EditedChemRICH.R'))
CRData <- run_chemrich_basic(inputfile = here('Data/Plot_Data/NCvCombo_CR_nopep.xlsx') )
#Remove the categories we don't care about
CRPlot <- CRData %>%
  filter(!name %in% c('Hydroxy bile acids, alcohols and derivatives', 'Glycinated bile acids and derivatives', 
                   'Amino acids and derivatives', 'Peptides'))


#Code Taken from the ChemRICH Sourecode
CHRPlot <- ggplot(CRPlot,aes(x=order,y=-log(pvalues))) +
  geom_point(aes(size=csize, color=upratio)) +
  #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
  scale_color_gradient(low = "blue", high = "red", limits=c(0,1))+
  scale_size(range = c(10, 40)) +
  scale_y_continuous("-log (p-value)",limits = c(0, max(-log(CRPlot$pvalues))+4  )) +
  scale_x_continuous(" Chemical Set Order (As Provided) ") +
  theme_bw() +
  #labs(title = "ChemRICH cluster impact plot") +
  geom_label_repel(aes(label = name), color = "gray20",family="Arial",data=subset(CRPlot, csize>2),
                   force = 12, size = 7, force_pull = 4)+
  theme(text=element_text(family="Arial Black"))+
  theme(
    plot.title = element_text(face="bold", size=30,hjust = 0.5),
    axis.title.x = element_text(face="bold", size=27),
    axis.title.y = element_text(face="bold", size=27, angle=90),
    panel.grid.major = element_blank(), # switch off major gridlines
    panel.grid.minor = element_blank(), # switch off minor gridlines
    legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=12),
    legend.key.size = unit(1.5, "lines"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.spacing = unit(.05, "cm"),
    axis.text.x = element_text(size=20,angle = 0, hjust = 1),
    axis.text.y = element_text(size=25,angle = 0, hjust = 1)
  )

#Final Plot for Figure 3: Save 3000x2000
plot_grid(fig3ae, CHRPlot, ncol = 1, rel_widths = c(1,1), labels = c('','F'), label_size = 23)


# Figure 4 ----------------------------------------------------------------
#Load in the positive mode data:
posugly <- read_csv(here('Data/Metabolomics/pos_metabolome_full.csv'))

#Metadata fix
mdata_pos <- read_csv(here('Data/Metabolomics/metadata_pos.csv')) %>%
  mutate(group = ifelse(group == 'A', 'C_Type', 'E_type')) %>%
  mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment)) %>%
  modify_at('treatment', factor, levels = c('NC', 'broc', 'brus', 'combo'))
  
#Scale and clean-up appropriately
pos_full <- posugly %>%
  filter(is.na(cv_g_50)) %>%
  dplyr::select(compound, starts_with('pos')) %>%
  mutate(compound = paste0(compound, '_pos')) %>%
  pivot_longer(starts_with('pos'), names_to = 'sample', values_to = 'intensity') %>%
  pivot_wider(names_from = 'compound', values_from = 'intensity') %>%
  column_to_rownames('sample')

#Make sure the metadata is good:
mpos <- left_join(data.frame(sample = rownames(pos_full)), mdata_pos)

#Load in the Diablo Data:
diablo <- readRDS(here('Data/Plot_Data/diablo_final.RDS'))

#Generate the consensus data by finding the average between the microbiome and metabolomics components
conDiablo <- data.frame(Reduce('+', diablo$variates)/2) %>%
  #Bind to the metadata
  cbind(diablo$Y) %>%
  #Rename it to be easier to work with  
  rename('treatment' = `diablo$Y`) %>%
  #Make a character for manipulation
  modify_at('treatment', as.character) %>%
  #Update the name of no_veg to NC for plotting
  mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment))
  
#Plot it nice
diPlot <- ggplot(conDiablo, aes(x = comp1, y = comp2, color = treatment)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = treatment, fill = treatment), geom = 'polygon', alpha = 0.4)  +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.title = element_blank()) +
  cowplot::theme_cowplot() +
  xlab('Component 1') +
  ylab('Component 2') +
  #theme(aspect.ratio = 1) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  scale_fill_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme(axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        #aspect.ratio = 1,
        legend.position = c(0.75,0.2)) +
  coord_fixed()

myrData <- pos_full %>%
  dplyr::select(`21.39_286.2365_pos`) %>%
  rownames_to_column('sample') %>%
  left_join(mpos) %>%
  rename('intensity' = `21.39_286.2365_pos`)

myrPlot <- ggplot(myrData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Myristoylglycine') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(myrData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

azaData <- neg_full %>%
  dplyr::select(`20.29_209.0795_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `20.29_209.0795_neg`)

azaPlot <- ggplot(azaData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Azelaic Acid') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(azaData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

triData <- neg_full %>%
  dplyr::select(`23.43_329.2331_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `23.43_329.2331_neg`)

triPlot <- ggplot(triData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('TriHOME') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(triData$intensity)+15), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

hpdData <- neg_full %>%
  dplyr::select(`24.27_311.2225_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `24.27_311.2225_neg`)

hpdPlot <- ggplot(hpdData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('HpODE') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(hpdData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

pinData <- neg_full %>%
  dplyr::select(`23.76_285.2071_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `23.76_285.2071_neg`)

pinPlot <- ggplot(pinData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Pinellic Acid') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(pinData$intensity)+25), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

#Combine it all into one plot
fig4dots <- plot_grid(myrPlot, pinPlot, azaPlot, hpdPlot, triPlot, nrow = 1, labels = c(LETTERS[2:6]), 
          rel_widths = rep(1, 5), label_size = 20)

#Remove the rownames 
rownames(microdata_c) <- NULL
micro_clr <- microdata_c %>%
  #Make numeric
  column_to_rownames('sample') %>%
  #CLR Transform the data
  clr() %>%
  #Add back in sample data
  rownames_to_column('sample')

#Pull the data of interest
hpd_lachData <- hpdData %>%
  #Make sample names generic
  mutate(sample = gsub('neg', 'ms', sample)) %>%
  #Join in the microbiome data
  left_join(micro_clr) %>%
  #Pull what we want
  dplyr::select(treatment, intensity, ASV256)

#Make teh plot
hlPlot <- ggplot(hpd_lachData, aes(x = ASV256, y = intensity, color = treatment)) +
  geom_smooth(method = 'lm', se = F, size = 1.5) +
  #geom_smooth(method = 'lm', inherit.aes = F, aes(x = ASV256, y = intensity), se = F, color = 'black', size = 1.5) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Lachnospiraceae Anaerostipes vs HpODE') +
  theme(legend.position = c(0.75,0.15),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(hpd_lachData$intensity)+25), expand = c(0,0), name = 'Intensity') +
  scale_x_continuous(expand = c(0,0), limits = c(min(hpd_lachData$ASV256)-1, max(hpd_lachData$ASV256+1)),
                     name = 'CLR-Transfromed Abundance') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'),
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 27)) 
hlPlot

myr_lachData <- myrData %>%
  mutate(sample = gsub('pos', 'ms', sample)) %>%
  left_join(micro_clr) %>%
  dplyr::select(treatment, intensity, ASV151)

mlPlot <- ggplot(myr_lachData, aes(x = ASV151, y = intensity, color = treatment)) +
  geom_smooth(method = 'lm', se = F, size = 1.5) +
  #geom_smooth(method = 'lm', inherit.aes = F, aes(x = ASV151, y = intensity), se = F, color = 'black', size = 1.5) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Lachnospiraceae Blautia vs Myristoylglycine') +
  theme(legend.position = c(0.1, 0.2),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(-1, max(myr_lachData$intensity)+2), expand = c(0,0), name = 'Intensity') +
  scale_x_continuous(expand = c(0,0), limits = c(min(myr_lachData$ASV151)-1, max(myr_lachData$ASV151+1)),
                     name = 'CLR-Transfromed Abundance') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'),
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 27)) 
mlPlot

#Comcbine the plots
crPlots <- plot_grid(hlPlot, mlPlot, rel_widths = rep(1,2), labels = c('G', 'H'), label_size = 20)

#Ok now actually combine them in a way we like
ugh1 <- plot_grid(diPlot, myrPlot, pinPlot, azaPlot, nrow = 1, rel_widths = c(2,1,1,1), 
                  labels = c('A', 'B', 'C', 'D'), label_size = 20)
ugh2 <- plot_grid( hpdPlot, triPlot, crPlots, nrow = 1, rel_widths = c(1,1,3),
                   labels = c('E', 'F', ''), label_size = 20)

#Make the final plot
plot_grid(ugh1, ugh2, nrow = 2, rel_heights = c(1,1))


# Figure 5 ----------------------------------------------------------------
#Pull the data we like
bnzData <- neg_full %>%
  dplyr::select(`23.72_605.2251_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `23.72_605.2251_neg`)

#Make the plot
bnzPlot <- ggplot(bnzData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  #geom_jitter(size = 2, width = 0.1, height = 0) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Benzenoid \n (m/z: 605.2251, RT: 23.72)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(bnzData$intensity)+25), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

aliData <- neg_full %>%
  dplyr::select(`23.48_745.2021_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `23.48_745.2021_neg`)

aliPlot <- ggplot(aliData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('3-Alkylindole \n (m/z: 745.2021, RT: 23.48)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(aliData$intensity)+25), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

flavData <- neg_full %>%
  dplyr::select(`22.77_705.3253_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `22.77_705.3253_neg`) %>%
  filter(sample != 'neg_40')

flavPlot <- ggplot(flavData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Flavonoid-3-O-Glycoside \n (m/z: 705.3253, RT: 22.77)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(flavData$intensity)+2), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

terpData <- neg_full %>%
  dplyr::select(`23.83_627.2378_neg`) %>%
  rownames_to_column('sample') %>%
  left_join(mneg) %>%
  rename('intensity' = `23.83_627.2378_neg`) 

terpPlot <- ggplot(terpData, aes(x = treatment, y = intensity, color = treatment)) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', color = 'black', size = 1.25, width = 0.3) +
  stat_summary(fun = 'mean', geom = 'bar', fill = 'white', color = 'black', size = 1.25, width = 0.5) +
  geom_point(size = 4) +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  theme_cowplot() +
  ggtitle('Triterpenoid \n (m/z: 627.2378, RT: 23.83)') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, max(terpData$intensity)+5), expand = c(0,0), name = 'Intensity') +
  scale_x_discrete(labels = c('NC', 'Broc', 'Brus', 'Combo'), name = NULL) +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 25, face = 'plain'))

#Combine them all together
fig5dots <- plot_grid(bnzPlot, aliPlot, flavPlot, terpPlot, nrow = 1, labels = c('A', 'B', 'C', 'D'), label_size = 20)

#Load in the correlation matrix
corData <- read_csv(here('Data/Plot_Data/FullCorrelationTable.csv'))
#Load in teh class data and clean
classData <- read_csv(here('Data/Plot_Data/MicrobialFilterTable.csv')) %>%
  dplyr::select(Compound, `most specific class`, molecularFormula, Polarity, `level 5`, superclass) %>%
  rename('class' = `most specific class`) %>%
  mutate(feature = ifelse(Polarity == 'Negative', paste0(Compound, '_neg'), paste0(Compound, '_pos')))

#Load in the taxonomic information
taxData <- as.data.frame(tax_table(ps_final)) %>%
  rownames_to_column('ASV')

#Pull only the compounds we find interesting
hmCompounds <- classData %>%
  filter(class %in% c('Azoles', 'Thiodioxopiperazines', 'Benzenoids', '1-hydroxy-2-unsubstituted benzenoids', 'Diterpene glycosides',
                      'Medium-chain fatty acids', 'Long-chain fatty acids', 'Alkanesulfonic acids', 'Acyclic monoterpenoids',
                      'Heteroaromatic compounds', 'Eicosanoids', 'Benzene and substituted derivatives', 'Triterpenoids',
                      'Limonoids', '3-alkylindoles', 'Cyclic depsipeptides', 'Indoles', 'Benzoic acids and derivatives',
                      'Flavonoid O-glycosides', '2,4-disubstituted thiazoles', 'Flavonoid-3-O-glycosides', 'Diterpenoids'))
  
#Prep the data for the heatmap
hmData <- corData %>%
  inner_join(., hmCompounds) %>%
  #Format the names proeprly
  mutate(feature = paste0(class, ' (', feature, ')')) %>%
  dplyr::select(starts_with('ASV'), feature) %>%
  pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'cor') %>%
  left_join(taxData) %>%
  dplyr::select(Genus, feature, cor) %>%
  pivot_wider(names_from = 'Genus', values_from = 'cor') %>%
  column_to_rownames('feature')

#Make the heatmap
hmap2 <- pheatmap(hmData, cellwidth = 21, cellheight = 21, fontsize = 20, fontsize_row = 20, fontsize_col = 20)

#Put it all together
plot_grid(fig5dots, hmap2$gtable, ncol = 1, rel_heights = c(1,2), labels = c('', 'E'), label_size = 23)


# Supplemental Figures ----------------------------------------------------


# Supplemental Figure 1 ---------------------------------------------------
#Load Vegan for more analysis
library(vegan)
#Subset our data appropriately
dps <- ps_counts %>%
  subset_samples(treatment != c('fecal_stock') )

#Make the rarefaction curves
d <- dps %>%  
  otu_table()

S <- specnumber(d) # observed number of species
(raremax <- min(rowSums(d)))
Srare <- rarefy(d, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
out <- rarecurve(d, step = 20, sample = raremax, col = "blue", cex = 0.6)

rare <- lapply(out, function(x){b <- as.data.frame(x) 
b <- data.frame(clase = b[,1], raw.read = rownames(b)) 
b$raw.read <- as.numeric(gsub("N", "",  b$raw.read)) 
return(b)})

#convert to data frame:
rare <- map_dfr(rare, function(x){
  z <- data.frame(x) 
  return(z)
}, .id = "Sample") 

#Pull metadata
md <- data.frame(sample_data(dps)) %>%
  dplyr::select(treatment, fecal_sample, sample)

rare %<>% modify_at('Sample', as.integer)

#Pull out the data to plot
uh <- rare %>% 
  dplyr::group_by(Sample) %>% 
  nest() %>%
  cbind(md) %>%
  unnest(cols = 'data') %>%
  modify_at('treatment', as.character) %>%
  mutate(treatment = ifelse(treatment == 'no_veg', 'NC', treatment)) 

#Plot that beautiful data
ggplot(uh, aes(x=raw.read, y=clase, color = treatment, group = Sample))+
  geom_path(size = 1.25) + 
  labs(x = 'Sample Size', y = 'Species') +
  cowplot::theme_cowplot() +
  ggtitle('Rarefaction Curves - Species Level') +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment')


# Supplemental Figure 2 ---------------------------------------------------
#Scale the data
pos_scaled <- pos_full %>% 
  #column_to_rownames('sample') %>%
  log_transform() %>%
  pareto_scale()%>%
  as.data.frame() 

neg_scaled <- neg_full %>% 
  #column_to_rownames('sample') %>%
  log_transform() %>%
  pareto_scale()%>%
  as.data.frame() 

#PCA Analysis:
#Load GGFortify to make plotting easy
library(ggfortify)
#Run the PCA
pca_pos <- prcomp(pos_scaled, center = TRUE)
#Make it pretty and plot it
pcpPlot <- autoplot(pca_pos, data = mpos, colour = 'treatment', size = 3) + ggtitle('Positive Ionization Mode') +
  stat_ellipse(aes(group = treatment, color = treatment, fill = treatment), geom = 'polygon', alpha = 0.2)  +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  cowplot::theme_minimal_grid() +
  scale_fill_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment')  +
  coord_fixed()

pca_neg <- prcomp(neg_scaled, center = TRUE)
pcnPlot <- autoplot(pca_neg, data = mneg, colour = 'treatment', size = 3) + ggtitle('Negative Ionization Mode')  +
  stat_ellipse(aes(group = treatment, color = treatment, fill = treatment), geom = 'polygon', alpha = 0.2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  cowplot::theme_minimal_grid() +
  scale_fill_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment') +
  scale_color_manual(values = pal, labels = c('Broc', 'Brus', 'Combo', 'NC'), name = 'Treatment')  +
  coord_fixed()

#Put them together into one image
plot_grid(pcpPlot, pcnPlot, nrow = 1, labels = c('A', 'B'))


# Supplemental Figure 3 ---------------------------------------------------
library(pROC)

#Run leave one out cross validation:
dtestloo <- perf(diablo, validation = 'loo', auc = TRUE, progressBar = FALSE)

#Extract out the data to make the ROC curves
metab1 <- as.data.frame(dtestloo$predict$nrep1$metab$comp1)
metab2 <- as.data.frame(dtestloo$predict$nrep1$metab$comp2)
micro1 <- as.data.frame(dtestloo$predict$nrep1$micro$comp1)
micro2 <- as.data.frame(dtestloo$predict$nrep1$micro$comp2)

#Make the class matrix
classmat <- data.frame(class = diablo$Y) %>%
  mutate(broc = ifelse(class == 'broc', 1, 0)) %>%
  mutate(brus = ifelse(class == 'brus', 1, 0)) %>%
  mutate(combo = ifelse(class == 'combo', 1, 0)) %>%
  mutate(no_veg = ifelse(class == 'no_veg', 1, 0)) %>%
  dplyr::select(-class)

#Build a list of the outputs 
outputs <- list(metab_comp1 = metab1, metab_comp2 = metab2, micro_comp1 = micro1, micro_comp2 = micro2)

#Custom function to pull the classes and see how they do
pull_roc <- function(prediction_list, classes){
  purrr::map(prediction_list, function(x){
    map2(classes, x, roc)
  })
}

#Pull out the ROCS
rocs <- pull_roc(outputs, classmat)
#Make them into nice plots using ggroc
rocplots <- purrr::map(rocs, ggroc)
#Extract out the AUCs 
aucs <- purrr::map(rocs, function(x) purrr::map(x, auc))

#Make the plots pretty
p1 <- rocplots[[1]] +
  geom_path(size = 2) + 
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = palmb,
                     name = 'Outcome',
                     labels = c(paste0('Broc vs All: ', round(aucs[[1]][[1]], 3)),
                                paste0('Brus vs All: ', round(aucs[[1]][[2]], 3)),
                                paste0('Combo vs All: ', round(aucs[[1]][[3]], 3)),
                                paste0('NC vs All: ', round(aucs[[1]][[4]], 3)))) +
  ggtitle('Metabolome - Component 1')  +
  theme(legend.position = c(0.5,0.25))

p2 <- rocplots[[2]] +
  geom_path(size = 2) + 
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = palmb,
                     name = 'Outcome',
                     labels = c(paste0('Broc vs All: ', round(aucs[[2]][[1]], 3)),
                                paste0('Brus vs All: ', round(aucs[[2]][[2]], 3)),
                                paste0('Combo vs All: ', round(aucs[[2]][[3]], 3)),
                                paste0('NC vs All: ', round(aucs[[2]][[4]], 3)))) +
  ggtitle('Metabolome - Component 2')  +
  theme(legend.position = c(0.5,0.25))

p3 <- rocplots[[3]] +
  geom_path(size = 2) + 
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = palmb,
                     name = 'Outcome',
                     labels = c(paste0('Broc vs All: ', round(aucs[[3]][[1]], 3)),
                                paste0('Brus vs All: ', round(aucs[[3]][[2]], 3)),
                                paste0('Combo vs All: ', round(aucs[[3]][[3]], 3)),
                                paste0('NC vs All: ', round(aucs[[3]][[4]], 3)))) +
  ggtitle('Microbiome - Component 1')  +
  theme(legend.position = c(0.5,0.25))

p4 <- rocplots[[4]] +
  geom_path(size = 2) + 
  geom_abline(intercept = 1) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = palmb,
                     name = 'Outcome',
                     labels = c(paste0('Broc vs All: ', round(aucs[[4]][[1]], 3)),
                                paste0('Brus vs All: ', round(aucs[[4]][[2]], 3)),
                                paste0('Combo vs All: ', round(aucs[[4]][[3]], 3)),
                                paste0('NC vs All: ', round(aucs[[4]][[4]], 3)))) +
  ggtitle('Microbiome - Component 2')  +
  theme(legend.position = c(0.5,0.25))

#Put them all together
allroc <- plot_grid(p1, p2, p3, p4, ncol = 2)

#Make the Correlation Circle:

#Load in the diablo features
dfeatures_neg <- read_csv(here('Data/Plot_Data/diablo_out_neg.csv')) %>%
  dplyr::select(feature, GroupContrib, importance) %>%
  rename('gc_metab' = 'GroupContrib')

dfeatures_pos <- read_csv(here('Data/Plot_Data/diablo_out_pos.csv')) %>%
  dplyr::select(feature, GroupContrib, importance) %>%
  rename('gc_metab' = 'GroupContrib')

dfeatures <- rbind(dfeatures_neg, dfeatures_pos) %>%
  arrange(desc(abs(importance)))
  
dasv_c1 <- read_csv(here('Data/Plot_Data/DiabloLoadings_comp1.csv')) %>%
  dplyr::select(name, GroupContrib, comp1) %>%
  rename('gc_micro' = 'GroupContrib', 'ASV' = name, 'importance' = comp1)

dasv_c2 <- read_csv(here('Data/Plot_Data/DiabloLoadings_comp2.csv')) %>%
  dplyr::select(name, GroupContrib, comp2) %>%
  rename('gc_micro' = 'GroupContrib', 'ASV' = name, 'importance' = comp2)

#Pull otu the importance of each
dasv <- rbind(dasv_c1, dasv_c2) %>%
  unique() %>%
  arrange(desc(abs(importance)))

#Pull the Correlation Circle data and filter for plotting
dcir <- plotVar(diablo, style = 'ggplot2', var.names = T, plot = F) %>%
  filter(names %in% c(dfeatures$feature, dasv$ASV)) 

#Circle function to make plotting easier
circleFun <- function(center = c(-1,1),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#Make the circle
circle <- circleFun(c(0,0), 2, npoints = 100)

#Grab importance data
impdat <- rbind(dasv %>% rename('gc' = gc_micro, 'names' = ASV),
                dfeatures %>% rename('gc' = gc_metab, 'names' = feature))

#Put it all together
dcir2 <- dcir %>%
  left_join(impdat)

#Make a new palette to inlcue 'Tie'
tcol <- RColorBrewer::brewer.pal(n = 6, name = 'Dark2')[6]
names(tcol) <- 'tie'
paltie <- c(palmb, tcol)

#Plot it
corcirc <- ggplot(dcir2, aes(x, y, color = gc, shape = Block)) + 
  geom_point(size = 4) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  xlab('Component 1') +
  ylab('Component 2') +
  ggtitle('Correlation Circle')  +
  #Make the circle independent of other aesthetics
  geom_path(aes(x,y), data = circle, inherit.aes = F)  +
  scale_color_manual(values = paltie, labels = c('Broc', 'Brus', 'Combo', 'NC', 'Tie'), name = 'Group \n Contribution')

#Put it all together into one plot.
plot_grid(allroc, corcirc, labels = c('A', 'B'), ncol = 2)






