options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/nd_trials/')
}
start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...'); flush.console()

suppressMessages(library(tidyverse))
suppressMessages(library(DiagrammeR))
suppressMessages(library(DiagrammeRsvg))
suppressMessages(library(magrittr))
suppressMessages(library(rsvg))
suppressMessages(library(MASS)); select=dplyr::select; summarize=dplyr::summarize;
suppressMessages(library(eulerr))
suppressMessages(library(magick))

##############
# OUTPUT FILE
##############

cat(file=stderr(), 'done!\nCreating output stream...'); flush.console()

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T

###############
#### FUNCTIONS & CONSTANTS
###############

cat(file=stderr(), 'done!\nLoading functions and constants...'); flush.console()

percent = function(proportion,format='fg',digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format=format),"%",sep="") ) )
}

signed_percent = function(proportion) {
  result = gsub(' ','',paste(ifelse(proportion > 0, '+', ''),round(proportion*100,0),"%",sep="") )
  result[is.na(proportion)] = '-'
  return ( result )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

blend_hex_atomic = function(col1, col2, where) {
  rgb1 = strtoi(paste0('0x',c(substr(col1,2,3), substr(col1,4,5), substr(col1,6,7))))
  rgb2 = strtoi(paste0('0x',c(substr(col2,2,3), substr(col2,4,5), substr(col2,6,7))))
  rgb = floor(rgb1 + where*(rgb2-rgb1))
  hex = paste0('#',paste0(format(as.hexmode(rgb),width=2),collapse=''))
  return (hex)
}
blend_hex = function(col1, col2, where) {
  stopifnot(length(col1)==length(col2))
  if (length(col1)!=length(where) & length(col1)==1) {
    col1 = rep(col1, length(where))
    col2 = rep(col2, length(where))
  }
  return_value = character(length(col1))
  for (i in 1:length(col1)) {
    return_value[i] = blend_hex_atomic(col1[i], col2[i], where[i])
  }
  return (return_value)
}

# clinicaltrials.gov dates come in two formats - this handles both cases
parse_ct_dates = function(x) {
  returnvalue = as.Date(rep(NA,length(x)))
  returnvalue[grepl(',',x)] = as.Date(x[grepl(',',x)], format='%B %d, %Y', optional=T)
  returnvalue[!grepl(',',x)] = as.Date(paste0(x[!grepl(',',x)], ', 1'), format='%B %Y, %d', optional=T)
  return (returnvalue)
}

resx=300

###############
#### DATASET
###############

cat(file=stderr(), 'done!\nReading and processing datasets..'); flush.console()

# auxiliary datasets
genes = read_tsv('other_data/gencode_genes.tsv', col_names='gene', col_types=cols())

# clinicaltrials.gov curated data
trials = read_tsv('curated_dataset/trials.tsv',guess_max=5000, col_types=cols())
interventions = read_tsv('curated_dataset/interventions.tsv',guess_max=10000, col_types=cols())
targets = read_tsv('curated_dataset/targets.tsv',guess_max=5000, col_types=cols())
timatch = read_tsv('xml_parsing_output/trial_intervention_match.tsv',guess_max=5000,quote='',comment='', col_types=cols())


# drug data
imaging_agents = read_tsv('other_data/imaging_agents.tsv', col_types=cols()) # approved for imaging, not therapy
overrides = read_tsv('other_data/overrides.tsv', col_types=cols()) # manual overrides of questionable drugbank target assignments
drugbank = read_tsv('other_data/drugbank.tsv',guess_max=5000, col_types=cols())
drugbank = separate_rows(drugbank, target, sep=',') # separate out the comma-delimited target column
drugbank$status[drugbank$drug %in% imaging_agents$drug] = 'imaging'
drugbank$status[paste0(drugbank$drug, '-', drugbank$target) %in% paste0(overrides$drug, '-', overrides$target)] = 'override'
atc_summary = read_tsv('other_data/atc_summary.tsv', col_types=cols())

# sources for disease-specific FDA approvals
# https://www.apdaparkinson.org/what-is-parkinsons/treatment-medication/medication/
# https://www.alz.org/media/documents/fda-approved-treatments-alzheimers-ts.pdf
# https://www.theaftd.org/what-is-ftd/ftd-and-als-ftd-als/#:~:text=At%20present%2C%20there%20are%20no,both%20physical%20and%20cognitive%20abilities.
# https://www.als.org/navigating-als/living-with-als/fda-approved-drugs
# https://hdsa.org/hd-research/therapies-in-pipeline/
approved = read_tsv('other_data/approved.tsv', col_types=cols())
approved$fdastatus = 'approved'

interventions %>%
  filter(intervention_class=='molecular') %>%
  select(intervention_name, intervention_type, target_gene) %>%
  separate_rows(target_gene, sep=',') %>%
  mutate(target_gene=replace_na(target_gene, 'unknown')) %>%
  filter(!(target_gene %in% c('unknown','none'))) -> ig_match

interventions %>%
  transmute(intervention_name,
            n_trials,
            example_trial,
            intervention_type, # as annotated in ClinicalTrials.gov - used for joins
            intervention_class, # our determination - used for analysis
            target_gene,
            regulatory_status,
            best_name
  ) -> ivns

# QC to check for duplicates with different orderings of combination therapies
# ivns %>%
#   filter(intervention_class=='molecular') %>%
#   inner_join(tim, by=c('intervention_name','intervention_type')) %>%
#   inner_join(trls, by='nct') %>%
#   select(best_name) %>% pull() %>% unique() -> drug_names
# tibble(drug_names, name=drug_names) %>%
#   separate_rows(name, sep=',') %>%
#   ungroup() %>%
#   arrange(drug_names, name) %>%
#   group_by(drug_names) %>%
#   summarize(.groups='keep',
#             alpha_combined = toString(name)) -> recombined
# recombined %>%
#   group_by(alpha_combined) %>%
#   summarize(.groups='keep', n=n()) %>%
#   arrange(desc(n))

# QC to make sure we never disagree with ourselves about what the target is
# interventions %>% filter(incl_count > 0 & !(target_gene %in% c('none','unknown'))) %>% group_by(best_name, target_gene) %>% summarize(.groups='keep', n=n()) %>% arrange(desc(n)) -> targ_assign
# targ_assign %>% group_by(best_name) %>% summarize(.groups='keep', n_versions=n()) %>% arrange(desc(n_versions)) -> targ_multip

timatch %>%
  filter(nct %in% trials$nct[trials$dataset_inclusion=='i - include']) %>%
  rename(intervention_name_original = intervention_name) %>%
  left_join(interventions, by=c('intervention_name_original','intervention_type')) %>%
  group_by(nct, intervention_name, intervention_type, intervention_class) %>%
  slice(1) %>%
  ungroup() %>%
  select(nct, intervention_name, intervention_type, intervention_class) -> tim

# identify additional potential SOC arms
trials %>%
  select(nct, disease_area, year) %>%
  inner_join(tim, by='nct') %>%
  inner_join(interventions, by=c('intervention_name','intervention_type','intervention_class')) %>%
  left_join(approved, by=c('best_name'='drug', 'disease_area')) %>%
  mutate(already_approved = !is.na(approval_year) & year >= approval_year) %>%
  select(nct, disease_area, year, approval_year, intervention_name, intervention_type, intervention_class, best_name, already_approved) -> tim_temp1

tim_temp1 %>%
  group_by(nct) %>%
  summarize(.groups='keep', 
            any_potentially_novel = any(!already_approved & !(intervention_class %in% c('placebo','none')))) -> pnovel

tim_temp1 %>%
  inner_join(pnovel, by='nct') %>%
  mutate(probable_soc_arm = already_approved & any_potentially_novel) %>% # this drug is approved, and trial contains a separate, potentially novel intervention arm
  mutate(intervention_class = ifelse(probable_soc_arm, 'placebo', intervention_class)) -> tim_temp2

tim_temp2 %>%
  select(nct, intervention_type, intervention_name, probable_soc_arm) %>%
  inner_join(tim, by=c('nct','intervention_name','intervention_type')) %>%
  mutate(intervention_class = ifelse(probable_soc_arm, 'placebo', intervention_class)) -> tim2

tim2 %>%
  select(nct, intervention_name, intervention_type, intervention_class) -> tim

tim %>%
  group_by(nct) %>%
  summarize(.groups='keep',
            has_placebo = any(intervention_class %in% 'placebo'),
            iclass = case_when(any(intervention_class %in% 'behavioral') ~ 'behavioral',
                               any(intervention_class %in% 'procedure') ~ 'procedure',
                               any(intervention_class %in% 'device') ~ 'device',
                               any(intervention_class %in% 'molecular') ~ 'molecular',
                               TRUE ~ 'none')) -> trl_ivn_info

trials$has_placebo = trl_ivn_info$has_placebo[match(trials$nct, trl_ivn_info$nct)]
trials$iclass = trl_ivn_info$iclass[match(trials$nct, trl_ivn_info$nct)]

trials$n_incl = 0
trials$n_excl = 0

hyphened = grepl(' - ',trials$inclusion_text) | grepl(' - ',trials$exclusion_text)
trials$n_incl[hyphened] = str_count(trials$inclusion_text[hyphened], ' - ')
trials$n_excl[hyphened] = str_count(trials$exclusion_text[hyphened], ' - ')
numbered = !hyphened & (grepl('1\\.',trials$inclusion_text) | grepl('1\\.',trials$exclusion_text))
trials$n_incl[numbered] = str_count(trials$inclusion_text[numbered], '[0-9]+\\. ')
trials$n_excl[numbered] = str_count(trials$exclusion_text[numbered], '[0-9]+\\. ')
trials$n_incl_excl = trials$n_incl + trials$n_excl
trials$n_incl_excl[!is.na(trials$curated_ie)] = trials$curated_ie[!is.na(trials$curated_ie)] # override if curated
#  sum(!hyphened & !numbered & is.na(trials$curated_ie) & trials$dataset_inclusion=='i - include' & trials$year >= 2000) # 0

write(paste('Inclusion/exclusion criteria count curated: ',sum(!is.na(trials$curated_ie)),'\n',sep=''),text_stats_path,append=T)

trials %>% 
  transmute(nct, 
            dataset_inclusion = case_when(year < 2000 ~ 't - exclude, pre-2000', 
                                          iclass %in% 'none' ~ 'n - exclude, no intervention',
                                          TRUE ~ dataset_inclusion), 
            exclude_comment,
            year,
            disease_area = ifelse(disease_area %in% c('other','unspecified'), 'other/unspecified', disease_area),
            incl = dataset_inclusion=='i - include',
         phs = case_when(
           phase %in% c('1','1/2','Early 1') ~ 1,
           phase %in% c('2','2/3') ~ 2,
           phase %in% c('3') ~ 3,
           phase %in% c('4') ~ 4,
           TRUE ~ -1
         ),
         industry = case_when(
           sponsor_type == 'industry' ~ TRUE,
           TRUE ~ FALSE
         ),
         advanced_exclusion,
         advanced_exclusion_comments,
         score_test,
         score_least = as.numeric(score_least_advanced),
         score_most = case_when(
           score_most_advanced == '6b' ~ 6.2,
           TRUE ~ suppressWarnings(as.numeric(score_most_advanced))
         ),
         status = case_when(
           overall_status %in% c('Terminated','Withdrawn','Suspended') ~ 'Stopped',
           overall_status %in% c('Completed') ~ 'Completed',
           overall_status %in% c('Not yet recruiting', 'Recruiting', 'Active, not recruiting', 'Enrolling by invitation') ~ 'Active',
           overall_status %in% c('Unknown status') ~ 'Unknown',
           TRUE ~ 'Unknown'
         ),
         start_date_d = parse_ct_dates(start_date),
         primary_completion_d = parse_ct_dates(primary_completion_date),
         study_completion_d = parse_ct_dates(completion_date),
         completion_year = as.integer(format(study_completion_d, '%Y')),
         primary_duration = as.integer(difftime(primary_completion_d, start_date_d, units='days')),
         study_duration = as.integer(difftime(study_completion_d, start_date_d, units='days')),
         earliest = as.integer(substr(earliest_disease_stage, 1, 1)),
         latest = as.integer(substr(latest_disease_stage, 1, 1)),
         enrollment = ifelse(enrollment_type=='Actual',enrollment,as.numeric(NA)),
         enrollment_type,
         n_incl,
         n_excl,
         n_incl_excl,
         has_placebo,
         iclass
  ) -> trls_all

trls_all %>%
  filter(dataset_inclusion=='i - include') -> trls

write(paste('Trials considered: ',nrow(trls_all),'\n',sep=''),text_stats_path,append=T)
write(paste('Trials curated: ',nrow(trls),'\n',sep=''),text_stats_path,append=T)

write(paste('All interventions: ',nrow(ivns),'\n',sep=''),text_stats_path,append=T)
ivn_count = length(unique(tim$intervention_name[tim$nct %in% trls$nct[trls$dataset_inclusion=='i - include']]))
write(paste('Interventions curated in included trials: ',ivn_count,'\n',sep=''),text_stats_path,append=T)

ivns %>%
  inner_join(tim, by=c('intervention_name','intervention_type')) %>%
  inner_join(trls, by='nct') %>%
  group_by(intervention_name) %>%
  summarize(.groups='keep', n_nct=length(unique(nct))) -> nct_per_ivn
only1 = mean(nct_per_ivn$n_nct==1)
meanotherwise = mean(nct_per_ivn$n_nct[nct_per_ivn$n_nct > 1])
sdotherwise = sd(nct_per_ivn$n_nct[nct_per_ivn$n_nct > 1])
write(paste('Interventions: ',percent(only1),' had just 1 trial, rest had (mean±sd)  ',round(meanotherwise,1),'±',round(sdotherwise,1),' trials.','\n',sep=''),text_stats_path,append=T)

interventions %>% 
  filter(!is.na(target_gene_source)) %>%
  group_by(target_gene_source) %>% 
  summarize(.groups='keep', n_=n()) %>% 
  ungroup() %>% 
  mutate(target_gene_source=ifelse(n_ < 10,'other',target_gene_source)) %>% 
  group_by(target_gene_source) %>% 
  summarize(.groups='keep', ntot=sum(n_)) %>% 
  arrange(desc(ntot)) %>%
  mutate(text = paste(target_gene_source, ntot, sep=': ')) -> ivn_target_sources

write(paste('Intervention target sources: ',paste(ivn_target_sources$text,collapse='; '),'\n',sep=''),text_stats_path,append=T)


###############
#### METADATA
###############

cat(file=stderr(), 'done!\nLoading metadata..'); flush.console()

diseases = tibble(disease_area=c('alzheimer','parkinson','ftd/als','huntington','other/unspecified'),
                  disp = c('AD','PD','FTD/ALS','HD','other/NS'),
                  axis_order=1:5,
                  color=c('#880088','#FF6600','#008888','#66FF00','#A9A9A9'))

phases = tibble(phslet = c('i','ii','iii','o'),
                disp = c('I','II','III','other/NS'),
                axis_order=1:4,
                color = c('#A8DDB5','#7BCCC4','#4EB3D3','#DDDD44'))

sponsors = tibble(sponsor = c('industry','non'),
                 disp = c('industry','other'),
                 axis_order = 1:2,
                 color=c('#CC3A78','#004FBA'))

iclasses = tibble(iclass = c('molecular','device','behavioral','procedure'),
                  disp = c('drug', 'device', 'behavioral', 'procedure'),
                  axis_order = c(1,3,2,4),
                  color = c('#F0027F','#386CB0','#7FC97F','#FDC086'))

controls = tibble(placebo = c('yes','no'),
                  disp = c('yes','no'),
                  axis_order = 1:2,
                  color = c('#33A02C','#CDB7B5'))

cat(file=stderr(), 'done!\nCreating Figure 1..'); flush.console()

###############
#### FIGURE 1
###############

trls_all %>%
  group_by(dataset_inclusion, exclude_comment) %>%
  summarize(.groups='keep', n=n()) %>%
  mutate(exclusion_subcat = case_when(is.na(exclude_comment) ~ 'other',
                                      n >= 10 ~ exclude_comment,
                                      TRUE ~ 'other')) %>%
  group_by(dataset_inclusion, exclusion_subcat) %>%
  summarize(.groups='keep', ntot=sum(n)) %>%
  select(dataset_inclusion, exclusion_subcat, n=ntot) %>%
  arrange(grepl('^i',dataset_inclusion), desc(dataset_inclusion), exclusion_subcat=='other', desc(n)) -> include_exclude

include_exclude = include_exclude[!is.na(include_exclude$dataset_inclusion),]

incl_count = include_exclude$n[!is.na(include_exclude$dataset_inclusion) & include_exclude$dataset_inclusion=='i - include']
pre2000_count = sum(include_exclude$n[!is.na(include_exclude$dataset_inclusion) & include_exclude$dataset_inclusion=='t - exclude, pre-2000'])
ni_list = paste0('- ', 
                 include_exclude$exclusion_subcat[include_exclude$dataset_inclusion=='n - exclude, no intervention'],
                 ' (',
                 include_exclude$n[include_exclude$dataset_inclusion=='n - exclude, no intervention'],
                 ')',
                 collapse='\\n')
hv_count = include_exclude$n[!is.na(include_exclude$dataset_inclusion) & include_exclude$dataset_inclusion=='h - exclude, healthy volunteers']
o_list = paste0('- ', 
                 include_exclude$exclusion_subcat[include_exclude$dataset_inclusion=='e - exclude, other'],
                 ' (',
                 include_exclude$n[include_exclude$dataset_inclusion=='e - exclude, other'],
                 ')',
                 collapse='\\n')

grViz(paste0("digraph flowchart {
             # node definitions with substituted label text
             node [fontname = Helvetica, shape = rectangle]        
             search [label = '@@1']
             exnhub [label = '@@2' shape = point]
             exn [label = '@@3']
             exhhub [label = '@@4' shape = point]
             exh [label = '@@5']
             exohub [label = '@@6' shape = point]
             exo [label = '@@7']
             incl [label = '@@8']
             ivns [label = '@@9']
             exthub [label = '@@10' shape = point]
             ext [label = '@@11']
             # edge definitions with the node IDs
             edge [arrowhead=none] search -> exthub
             edge [arrowhead=none] exthub -> ext
             edge [arrowhead=none] exthub -> exnhub
             edge [arrowhead=none] exnhub -> exn
             edge [arrowhead=none] exnhub -> exhhub
             edge [arrowhead=none] exhhub -> exh
             edge [arrowhead=none] exhhub -> exohub
             edge [arrowhead=none] exohub -> exo
             edge [arrowhead=none] exohub -> incl
             edge [arrowhead=vee] incl -> ivns
             }
             
             [1]: 'Search results (",formatC(sum(include_exclude$n),big.mark=','),")'
             [2]: ''
             [3]: 'Excluded - no intervention\\n",ni_list,"'
             [4]: ''
             [5]: 'Excluded - healthy volunteers (",hv_count,")'
             [6]: ''
             [7]: 'Excluded - other\\n",o_list,"'
             [8]: 'Trial included for manual curation (",formatC(incl_count,big.mark=','),")'
             [9]: 'Unique interventions (",formatC(ivn_count,big.mark=','),")'
             [10]: ''
             [11]: 'Excluded - pre-2000 (",pre2000_count,")'
             ")) %>% export_svg %>% charToRaw %>% rsvg_png("display_items/figure-1.png",width=3.5*resx)






########################
#### BEGIN FIGURE 2 ####
########################

cat(file=stderr(), 'done!\nCreating Figure 2..'); flush.console()

resx=600
png('display_items/figure-2.png',width=6.5*resx,height=5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4,4,4,4, 4 ,
                         5,5,6,6,9,9,10,10,
                         5,5,7,7,9,9,10,10,
                         5,5,8,8,9,9,10,10),nrow=4,byrow=T)
layout(layout_matrix, widths=c(1.65, 1, 1, 1, 1, 1.5, 1.5, 1), heights=c(1,.33,.33,.33))

panel = 1



trls %>%
  filter(incl) %>%
  mutate(vrbl=case_when(is.na(disease_area) | disease_area %in% c('other','unspecified') ~ 'other/unspecified', TRUE ~ disease_area)) %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  inner_join(diseases, by=c('vrbl'='disease_area')) %>%
  arrange(axis_order)  -> univariate_by_disease


trls %>%
  filter(incl) %>%
  mutate(vrbl=case_when(industry ~ 'industry', TRUE ~ 'non')) %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  arrange(vrbl) %>%
  inner_join(sponsors, by=c('vrbl'='sponsor')) %>%
  arrange(axis_order)  -> univariate_by_sector

trls %>%
  filter(incl) %>%
  mutate(vrbl=case_when(phs == 1 ~ 'i', phs == 2 ~ 'ii', phs == 3 ~ 'iii', TRUE ~ 'o')) %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  inner_join(phases, by=c('vrbl'='phslet')) %>%
  arrange(axis_order) -> univariate_by_phase

trls %>%
  filter(incl) %>%
  mutate(vrbl=iclass) %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  inner_join(iclasses, by=c('vrbl'='iclass')) %>%
  arrange(axis_order) -> univariate_by_iclass

trls %>%
  filter(incl) %>%
  mutate(vrbl=case_when(has_placebo ~ 'yes', TRUE ~ 'no')) %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  inner_join(controls, by=c('vrbl'='placebo')) %>%
  arrange(axis_order)  ->  univariate_by_placebo

trls %>%
  filter(incl) %>%
  mutate(vrbl='TOTAL') %>%
  group_by(vrbl) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/n()) %>%
  ungroup() %>%
  mutate(disp='total', axis_order=1, color='#000000')  -> univariate_by_total

rbind(cbind(cat='disease',univariate_by_disease),
      cbind(cat='sponsor',univariate_by_sector),
      cbind(cat='intervention',univariate_by_iclass),
      cbind(cat='phase',univariate_by_phase),
      cbind(cat='placebo',univariate_by_placebo),
      cbind(cat='total',univariate_by_total)) %>%
  group_by(cat) %>%
  mutate(proportion_count = n_total/sum(n_total),
         proportion_py = py/sum(py)) -> univariate_table

write_tsv(univariate_table, 'qc/univariate_data.tsv')

write(paste('Industry sponsored: N=',
            univariate_table$n_total[univariate_table$vrbl=='industry'],
            ' (',percent(univariate_table$proportion_count[univariate_table$vrbl=='industry']),'), proportion PY=',
            percent(univariate_table$proportion_py[univariate_table$vrbl=='industry']),
            '\n',sep=''),text_stats_path,append=T)
write(paste('Drug trials: N=',
            univariate_table$n_total[univariate_table$vrbl=='molecular'],
            ' (',percent(univariate_table$proportion_count[univariate_table$vrbl=='molecular']),'), proportion PY=',
            percent(univariate_table$proportion_py[univariate_table$vrbl=='molecular']),
            '\n',sep=''),text_stats_path,append=T)

fig2 = univariate_table[1:(nrow(univariate_table)-1),]
fig2$y = nrow(fig2):1

ymeta = fig2 %>% group_by(cat) %>% summarize(.groups='keep', miny=min(y), maxy=max(y), midy=mean(y))
ybreaks = sort(unique(c(ymeta$miny - 0.5, ymeta$maxy + 0.5)))

par(mar=c(3,0,3,0))
xlims = c(0,1)
ylims = range(fig2$y) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4, line=-0.1, adj=1, las=2, at=fig2$y, text=fig2$disp, cex=0.6)
overhang = 0.4
for (i in 1:nrow(ymeta)) {
  axis(side=4, line=-3.85, at=c(ymeta$miny[i]-overhang, ymeta$maxy[i]+overhang), tck=-0.025, labels=NA)
  mtext(side=4, line=-3.95, at=ymeta$midy[i], text=ymeta$cat[i], cex=0.6, las=2, adj=1)
}

par(mar=c(3,0,3,0.5))
xlims = c(0,2100)
ylims = range(fig2$y) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:40*100, tck=-0.025, labels=NA)
axis(side=1, at=0:4*1000, tck=-0.05, labels=NA)
axis(side=1, at=0:2*1000, labels=c('0','1K','2K'), lwd.ticks=0, lwd=0, line=-1, cex.axis=.45)
mtext(side=1, line=1, text='N trials', cex=0.6)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
barwidth=0.4
rect(xleft=rep(0,nrow(fig2)), xright=fig2$n_total, ybottom=fig2$y-barwidth, ytop=fig2$y+barwidth, col=fig2$color, border=NA)
par(xpd=T)
abline(h=ybreaks, lwd=0.125)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.5); panel = panel + 1

par(mar=c(3,0,3,0.5))
xlims = c(0,5e5)
ylims = range(fig2$y) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:20*5e4, tck=-0.025, labels=NA)
axis(side=1, at=0:10*1e5, labels=NA, tck=-0.05)
axis(side=1, at=c(0,2e5,5e5), labels=c('0','200K','500K'), lwd.ticks=0, lwd=0, line=-1, cex.axis=.45)
mtext(side=1, line=1, text='patient-years', cex=0.6)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
barwidth=0.4
rect(xleft=rep(0,nrow(fig2)), xright=fig2$py, ybottom=fig2$y-barwidth, ytop=fig2$y+barwidth, col=fig2$color, border=NA) 
par(xpd=T)
abline(h=ybreaks, lwd=0.125)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.5); panel = panel + 1

trls %>%
  filter(incl) %>%
  mutate(vrbl1=case_when(industry ~ 'industry', TRUE ~ 'non')) %>%
  mutate(vrbl2=iclass) %>%
  group_by(vrbl1, vrbl2) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/sum(!is.na(enrollment) & !is.na(study_duration))) %>%
  inner_join(sponsors, by=c('vrbl1'='sponsor')) %>%
  inner_join(iclasses, by=c('vrbl2'='iclass'), suffix=c('_x','_y')) %>%
  rename(x=axis_order_x, y=axis_order_y) %>%
  arrange(desc(n_total)) -> 
  xtab_sector_iclass

trls %>%
  filter(incl) %>%
  mutate(vrbl1=case_when(phs == 1 ~ 'i', phs == 2 ~ 'ii', phs == 3 ~ 'iii', TRUE ~ 'o')) %>%
  mutate(vrbl2=iclass) %>%
  group_by(vrbl1, vrbl2) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/sum(!is.na(enrollment) & !is.na(study_duration))) %>%
  inner_join(phases, by=c('vrbl1'='phslet')) %>%
  inner_join(iclasses, by=c('vrbl2'='iclass'), suffix=c('_x','_y')) %>%
  rename(x=axis_order_x, y=axis_order_y) %>%
  arrange(desc(n_total)) -> 
  xtab_phs_iclass

trls %>%
  filter(incl) %>%
  mutate(vrbl1=case_when(has_placebo ~ 'yes', TRUE ~ 'no')) %>%
  mutate(vrbl2=iclass) %>%
  group_by(vrbl1, vrbl2) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/sum(!is.na(enrollment) & !is.na(study_duration))) %>%
  inner_join(controls, by=c('vrbl1'='placebo')) %>%
  inner_join(iclasses, by=c('vrbl2'='iclass'), suffix=c('_x','_y')) %>%
  rename(x=axis_order_x, y=axis_order_y) %>%
  arrange(desc(n_total))  -> 
  xtab_placebo_iclass


trls %>%
  filter(incl) %>%
  mutate(vrbl1=disease_area) %>%
  mutate(vrbl2=iclass) %>%
  group_by(vrbl1, vrbl2) %>%
  summarize(.groups='keep',
            n_total = n(),
            n_completed = sum(status=='Completed'),
            mean_enrollment = mean(enrollment, na.rm=T),
            mean_duration = mean(study_duration, na.rm=T)/365.25,
            py = sum(enrollment * study_duration, na.rm=T)/365.25,
            py_compl = sum(enrollment[status=='Completed'] * study_duration[status=='Completed']/365.25, na.rm=T),
            py_per = sum(enrollment * study_duration, na.rm=T)/365.25/sum(!is.na(enrollment) & !is.na(study_duration))) %>%
  inner_join(diseases, by=c('vrbl1'='disease_area')) %>%
  inner_join(iclasses, by=c('vrbl2'='iclass'), suffix=c('_x','_y')) %>%
  rename(x=axis_order_x, y=axis_order_y) %>%
  arrange(desc(n_total))  -> 
  xtab_disease_iclass

rbind(cbind(xtab_sector_iclass, xcat='sponsor', xcatorder=4),
      cbind(xtab_phs_iclass, xcat='phase', xcatorder=2),
      cbind(xtab_placebo_iclass, xcat='placebo', xcatorder=3),
      cbind(xtab_disease_iclass, xcat='disease', xcatorder=1)) %>% 
  ungroup() %>%
  mutate(ycat='intervention') %>% 
  mutate(y_combined = max(y) + 1 - y) %>%
  group_by(xcat, ycat) %>%
  mutate(proportion_count = n_total/sum(n_total),
         proportion_py = py/sum(py)) %>%
  ungroup() %>% 
  group_by(xcat, xcatorder, x, vrbl1) %>% 
  arrange(xcatorder, x) %>% 
  mutate(x_combined = cur_group_id()) -> xtab             

write(paste('Industry sponsored drug trials: N=',
            xtab$n_total[xtab$vrbl1=='industry' & xtab$vrbl2=='molecular'],
            ' (',percent(xtab$proportion_count[xtab$vrbl1=='industry' & xtab$vrbl2=='molecular']),'), proportion PY=',
            percent(xtab$proportion_py[xtab$vrbl1=='industry' & xtab$vrbl2=='molecular']),
            '\n',sep=''),text_stats_path,append=T)

cex_factor = 4.25
xtab$color = blend_hex('#AAAA22','#0001CD',pmin(1,xtab$py_per/1000))
xlims = range(xtab$x_combined) + c(-0.5, 0.5)
ylims = range(xtab$y) + c(-0.5, 0.5)
par(mar=c(4,4,4,1))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
points(xtab$x_combined, xtab$y_combined, pch=19, 
       cex=cex_factor*sqrt(xtab$n_total/1000),
       col=alpha(xtab$color, ci_alpha + (1-ci_alpha)*pmin(1,xtab$py_per/1000)))
mtext(side=3, at=xtab$x_combined, text=gsub('/','/\n',xtab$disp_x), cex=0.45)
mtext(side=2, at=xtab$y_combined, text=xtab$disp_y, cex=0.45)
xtab %>% 
  group_by(xcat) %>% 
  summarize(.groups='keep',
            minx = min(x_combined),
            maxx = max(x_combined),
            midx = mean(x_combined)) -> xmeta
overhang = 0.45
for (i in 1:nrow(xmeta)) {
  axis(side=3, line=1.75, at=unlist(xmeta[i,c('minx','maxx')])+c(-1,1)*overhang, tck=0.05, labels=NA)
  mtext(side=3, line=2.0, at=xmeta$midx[i], text=xmeta$xcat[i], cex=0.6)
}
xtab %>% 
  group_by(ycat) %>% 
  summarize(.groups='keep',
            miny = min(y),
            maxy = max(y),
            midy = mean(y)) -> ymeta
overhang = 0.4
for (i in 1:nrow(ymeta)) {
  axis(side=2, line=1.75, at=unlist(ymeta[i,c('miny','maxy')])+c(-1,1)*overhang, tck=0.05, labels=NA)
  mtext(side=2, line=2, at=ymeta$midy[i], text=ymeta$ycat[i], cex=0.7)
}

par(xpd=T)
legpy = c(100,300,1000)
legend(x=min(xlims),y=min(ylims),legend=formatC(legpy,big.mark=','),pch=19,pt.cex=3,
       col=alpha(blend_hex('#FFFF88','#0001CD',legpy/1000), ci_alpha + (1-ci_alpha)*legpy/1000),
       title='patient-years/trial',title.col='#000000',
       horiz=T,bty='n')
legnt = c(100,300,1000)
legend(x=mean(xlims),y=min(ylims),legend=formatC(legnt,big.mark=','),pch=19,
       pt.cex=cex_factor*sqrt(legnt/1000),
       col='#888888',title='N trials',title.col='#000000',
       horiz=T,bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

trls %>%
  mutate(indmol = industry & iclass=='molecular') %>%
  mutate(py = enrollment * study_duration / 365.25) %>%
  filter(!is.na(py) & py > 0) %>%
  select(indmol, nct, py) %>%
  arrange(-py) %>%
  group_by(indmol) %>%
  mutate(pyrank = rank(-py),
         pycum = cumsum(py)/sum(py)) -> pyranks

trls$pyrank = pyranks$pyrank[match(trls$nct, pyranks$nct)]
trls$pycum = pyranks$pycum[match(trls$nct, pyranks$nct)]
trls %>%
  mutate(indmol = industry & iclass=='molecular') %>%
  group_by(indmol) %>%
  summarize(.groups='keep', 
            n_nct = n(),
            py = sum(enrollment * study_duration/365.25, na.rm=T),
            pyper = mean(enrollment * study_duration/365.25, na.rm=T),
            sdpyper = sd(enrollment * study_duration/365.25, na.rm=T),
            l20p = mean(enrollment <= 20, na.rm=T),
            m33 = min(pyrank[pycum > 0.33],na.rm=T),
            placebo_no = sum(!has_placebo),
            placebo_yes = sum(has_placebo),
            completed_no = sum(!(status %in% 'Completed')),
            completed_yes = sum((status %in% 'Completed')),
            specphase_no = sum(phs %in% -1),
            specphase_yes = sum(phs %in% 1:3)
            ) %>%
  ungroup() %>%
  mutate(p_nct = n_nct/sum(n_nct),
         p_py = py / sum(py))-> ctables

write(paste('Top 1 trial PY: industry-drug: ',
            percent(pyranks$pycum[pyranks$indmol & pyranks$pyrank==1]),
            ' non:',
            percent(pyranks$pycum[!pyranks$indmol & pyranks$pyrank==1]),
            '\n',sep=''),text_stats_path,append=T)
write(paste('Top 33% of trials: industry-drug: ',
            ctables$m33[ctables$indmol],
            ' trials; non: ',
            ctables$m33[!ctables$indmol],
            ' trials.\n',sep=''),text_stats_path,append=T)
write(paste('≤20 PY trials: industry-drug: ',
            percent(ctables$l20p[ctables$indmol]),
            '; non: ',
            percent(ctables$l20p[!ctables$indmol]),
            '.\n',sep=''),text_stats_path,append=T)


fisher_placebo = fisher.test(ctables[,grepl('placebo',colnames(ctables))])
fisher_completed = fisher.test(ctables[,grepl('completed',colnames(ctables))])
fisher_specphase = fisher.test(ctables[,grepl('specphase',colnames(ctables))])

write(paste('Industry-drug vs. other - placebo OR: ',
            round(fisher_placebo$estimate,1),
            ', P = ',
            formatC(fisher_placebo$p.value,format='g',digits=1),
            '.\n',sep=''),text_stats_path,append=T)
write(paste('Industry-drug vs. other - completed OR: ',
            round(fisher_completed$estimate,1),
            ', P = ',
            formatC(fisher_completed$p.value,format='g',digits=1),
            '.\n',sep=''),text_stats_path,append=T)
write(paste('Industry-drug vs. other - specified phase OR: ',
            round(fisher_specphase$estimate,1),
            ', P = ',
            formatC(fisher_specphase$p.value,format='g',digits=1),
            '.\n',sep=''),text_stats_path,append=T)

par(mar=c(1,1,2,1))
im_col = '#CC3A78'
ot_col = '#33C587'
trls$py = trls$enrollment * trls$study_duration / 365.25
ot_py = rev(sort(trls$py[!(trls$industry & trls$iclass=='molecular') & !is.na(trls$py) & trls$py > 0]))
im_py = rev(sort(trls$py[(trls$industry & trls$iclass=='molecular') & !is.na(trls$py) & trls$py > 0]))
pie(c(im_py,ot_py), 
    col = alpha(c(rep(im_col,length(im_py)), rep(ot_col, length(ot_py))),c(.25,.5,.75)),
    labels=NA, border=NA)
par(xpd=T)
text(x=0,y=1,col=im_col,labels=paste0('industry drug trials'))
text(x=0,y=-1,col=ot_col,labels=paste0('all other'))
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.5); panel = panel + 1

yes_col = '#9A32CD'
no_col = '#808A87'
yes_y = 2
no_y = 1


par(mar=c(1,1,3,1))
xlims = c(-2000,2000)
ylims = c(0.5,2.5)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=-2:2*1000, labels=NA, tck=-0.05)
axis(side=1, at=-20:20*100, labels=NA, tck=-0.025)
axis(side=1, at=-2:2*1000, labels=c('2K','1K','0','1K','2K'), cex.axis=0.6, lwd=0, lwd.ticks=0, line=-1)
barwidth = 0.4
rect(xleft=0, xright=-ctables$placebo_no[!ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright=-ctables$placebo_yes[!ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
rect(xleft=0, xright=ctables$placebo_no[ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright=ctables$placebo_yes[ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
abline(v=0)
mtext(side=3, at=-1000, text='all other', cex=0.45)
mtext(side=3, at=1000, text='industry+drug', cex=0.45)
mtext(side=3, at=0, line=0.75, text='placebo-controlled', cex=0.6)
mtext(side=2, at=c(yes_y,no_y), line=0.25, las=2, text=c('yes','no'), col=c(yes_col, no_col), cex=0.45)
par(xpd=T)
text(x=600, y=mean(no_y,yes_y), pos=4, labels=paste0('OR=',formatC(fisher_placebo$estimate,format='f',digits=1)), cex=0.6)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 1.5); panel = panel + 1

par(mar=c(2,1,2,1))
xlims = c(-2000,2000)
ylims = c(0.5,2.5)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=-2:2*1000, labels=NA, tck=-0.05)
axis(side=1, at=-20:20*100, labels=NA, tck=-0.025)
axis(side=1, at=-2:2*1000, labels=c('2K','1K','0','1K','2K'), cex.axis=0.6, lwd=0, lwd.ticks=0, line=-1)
barwidth = 0.4
rect(xleft=0, xright=-ctables$completed_no[!ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright=-ctables$completed_yes[!ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
rect(xleft=0, xright= ctables$completed_no[ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright= ctables$completed_yes[ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
abline(v=0)
mtext(side=3, at=-1000, text='all other', cex=0.45)
mtext(side=3, at=1000, text='industry+drug', cex=0.45)
mtext(side=3, at=0, line=0.75, text='completed', cex=0.6)
mtext(side=2, at=c(yes_y,no_y), line=0.25, las=2, text=c('yes','no'), col=c(yes_col, no_col), cex=0.45)
par(xpd=T)
text(x=600, y=mean(no_y,yes_y), pos=4, labels=paste0('OR=',formatC(fisher_completed$estimate,format='f',digits=1)), cex=0.6)
par(xpd=F)


par(mar=c(3,1,1,1))
xlims = c(-2000,2000)
ylims = c(0.5,2.5)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=-2:2*1000, labels=NA, tck=-0.05)
axis(side=1, at=-20:20*100, labels=NA, tck=-0.025)
axis(side=1, at=-2:2*1000, labels=c('2K','1K','0','1K','2K'), cex.axis=0.6, lwd=0, lwd.ticks=0, line=-1)
barwidth = 0.4
rect(xleft=0, xright=-ctables$specphase_no[!ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright=-ctables$specphase_yes[!ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
rect(xleft=0, xright= ctables$specphase_no[ctables$indmol], ybottom=no_y-barwidth, ytop=no_y+barwidth, col=no_col, border=NA)
rect(xleft=0, xright= ctables$specphase_yes[ctables$indmol], ybottom=yes_y-barwidth, ytop=yes_y+barwidth, col=yes_col, border=NA)
abline(v=0)
mtext(side=3, at=-1000, text='all other', cex=0.45)
mtext(side=3, at=1000, text='industry+drug', cex=0.45)
mtext(side=3, at=0, line=0.75, text='specified phase', cex=0.6)
mtext(side=2, at=c(yes_y,no_y), line=0.25, las=2, text=c('yes','no'), col=c(yes_col, no_col), cex=0.45)
mtext(side=1, line=1.0, text='N trials', cex=0.6)
par(xpd=T)
text(x=600, y=mean(no_y,yes_y), pos=4, labels=paste0('OR=',formatC(fisher_specphase$estimate,format='f',digits=1)), cex=0.6)
par(xpd=F)


trls %>%
  mutate(indmol = industry & iclass=='molecular') %>%
  mutate(phs_fct=as.factor(case_when(phs == 1 ~ 'i', phs == 2 ~ 'ii', phs == 3 ~ 'iii', TRUE ~ 'o'))) %>%
  select(indmol, year, phs_fct, n_incl_excl) %>%
  inner_join(phases, by=c('phs_fct'='phslet')) -> ie_data
ie_data$color_im = ifelse(ie_data$indmol, im_col, ot_col)


ie_data %>%
  mutate(color = color_im) %>%
  group_by(indmol, year, color) %>%
  summarize(.groups='keep',
            n_nct=n(),
            mean=mean(n_incl_excl),
            l95 = mean(n_incl_excl) - 1.96*sd(n_incl_excl)/sqrt(n()),
            u95 = mean(n_incl_excl) + 1.96*sd(n_incl_excl)/sqrt(n())) -> ie_smry_by_indmol

ie_data %>%
  group_by(year, phs_fct, color) %>%
  summarize(.groups='keep',
            n_nct=n(),
            mean=mean(n_incl_excl),
            l95 = mean(n_incl_excl) - 1.96*sd(n_incl_excl)/sqrt(n()),
            u95 = mean(n_incl_excl) + 1.96*sd(n_incl_excl)/sqrt(n())) -> ie_smry_by_phs



par(mar=c(3,2,2,4))
xlims = c(2000,2020) + c(-0.5, 0.5)
ylims = c(0,22)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=1999:2021, labels=NA, tck=-0.025)
axis(side=1, at=seq(2000,2020,5), labels=NA, tck=-0.05)
axis(side=1, at=seq(2000,2020,5), lwd.ticks=0, lwd=0, line=-1, cex.axis=0.6)
mtext(side=1, line=1, text='year', cex=0.6)
axis(side=2, at=0:30, labels=NA, tck=-0.02)
axis(side=2, at=0:3*10, labels=NA, tck=-0.06)
axis(side=2, at=0:3*10, las=2, lwd=0, cex.axis=0.7, line=-0.5)
mtext(side=2, line=1, text='N incl/excl criteria', cex=0.6)
for (indmol_status in c(T,F)) {
  subs = subset(ie_data, indmol==indmol_status)
  #set.seed(1)
  #points(x=jitter(subs$year,amount=.25), y=jitter(subs$n_incl_excl,amount=.25), pch=20, cex=0.5, col=alpha(subs$color, ci_alpha))
  m = lm(n_incl_excl ~ year, data=subs)
  linear_p = summary(m)$coefficients['year','Pr(>|t|)']
  loess_model = loess(n_incl_excl ~ year, data=subs, span=1.5)
  x = seq(2000,2020,.1)
  y = predict(loess_model, x)
  points(x, y, type='l', col=subs$color_im[1], lwd=2)
  final_y = y[length(y)]
  mtext(side=4, line=0.25, at=final_y, col=subs$color_im[1], text=ifelse(indmol_status,'industry\ndrug\ntrials','all\nother'),las=2, cex=0.6)
  #abline(m, col=subs$color[1])
  smry_subs = subset(ie_smry_by_indmol, indmol==indmol_status)
  barwidth=0.4
  segments(x0=smry_subs$year-barwidth,x1=smry_subs$year+barwidth, y0=smry_subs$mean, col=alpha(smry_subs$color, ci_alpha), lwd=1.5)
  
  write(paste('Inclusion/exclusion criteria time trend -- ',
              ifelse(indmol_status,'Industry-drug','All other'), ': mean rose from ',
              paste0(formatC(predict(loess_model, c(2000,2020)),format='f', digits=1), collapse=' to '),
              ', P = ',
              formatC(linear_p,format='g',digits=2),
              '.\n',sep=''),text_stats_path,append=T)
}
#legend('bottomleft',c('industry drug trials','all other'),lwd=1.5,col=c(im_col,ot_col),text.col=c(im_col,ot_col),cex=0.6,bty='n')
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.5); panel = panel + 1

par(mar=c(3,2,2,4))
xlims = c(2000,2020) + c(-0.5, 0.5)
ylims = c(0,22)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=1999:2021, labels=NA, tck=-0.02)
axis(side=1, at=seq(2000,2020,5), labels=NA, tck=-0.06)
axis(side=1, at=seq(2000,2020,5), lwd.ticks=0, lwd=0, line=-1, cex.axis=0.6)
mtext(side=1, line=1, text='year', cex=0.6)
axis(side=2, at=0:30, labels=NA, tck=-0.025)
axis(side=2, at=0:3*10, labels=NA, tck=-0.05)
axis(side=2, at=0:3*10, las=2, lwd=0, cex.axis=0.7, line=-0.5)
mtext(side=2, line=1, text='N incl/excl criteria', cex=0.6)
for (phs in phases$phslet) {
  subs = subset(ie_data, phs_fct==phs)
  #set.seed(1)
  #points(x=jitter(subs$year,amount=.25), y=jitter(subs$n_incl_excl,amount=.25), pch=20, cex=0.5, col=alpha(subs$color, ci_alpha))
  m = lm(n_incl_excl ~ year, data=subs)
  linear_p = summary(m)$coefficients['year','Pr(>|t|)']
  loess_model = loess(n_incl_excl ~ year, data=subs, span=1.5)
  x = seq(2000,2020,.1)
  y = predict(loess_model, x)
  final_y = y[length(y)]
  if (phs=='ii') final_y = final_y+0.5
  mtext(side=4, line=0.25, at=final_y, col=subs$color[1], text=paste0(ifelse(phs=='o','','phase '),phases$disp[phases$phslet==phs]), las=2, cex=0.6)
  points(x, y, type='l', col=subs$color[1], lwd=2)
  #abline(m, col=subs$color[1])
  smry_subs = subset(ie_smry_by_phs, phs_fct==phs)
  barwidth=0.4
  segments(x0=smry_subs$year-barwidth,x1=smry_subs$year+barwidth, y0=smry_subs$mean, col=alpha(smry_subs$color, ci_alpha), lwd=1.5)
  
  
  write(paste('Inclusion/exclusion criteria time trend -- phase ',phs,': mean rose from ',
              paste0(formatC(predict(loess_model, c(min(x[!is.na(y)]),2020)),format='f', digits=1), collapse=' to '),
              ', P = ',
              formatC(linear_p,format='g',digits=2),
              '.\n',sep=''),text_stats_path,append=T)
}
#legend('bottomleft',phases$disp,lwd=1.5,col=phases$color,text.col=phases$color,title='phase',title.col='#000000',cex=0.6,bty='n')
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0, line = 0.5); panel = panel + 1

silence_message = dev.off()

######################
#### END FIGURE 2 ####
######################

















###############
#### FIGURE 3 #
###############

cat(file=stderr(), 'done!\nCreating Figure 3..'); flush.console()

png('display_items/figure-3.png',width=6.5*resx,height=5*resx,res=resx)

layout_matrix = matrix(c(1,2,2,3,3,
                         1,2,2,4,4,
                         5,5,6,6,7,
                         5,5,6,6,7), nrow=4, byrow=T)
layout(layout_matrix, widths=c(.6, .5, .5, .5, 1), heights=c(1,.9,.75, .75))

panel = 1

# legend panel
par(mar=c(3,0,3,0))
stages = tibble(stage=0:4,
                desc=c("at-risk","molecular","detectable","mild","diagnosed"),
                color=c('#BCED91','#61B329','#0099CC','#7D26CD','#E3170D'))
plot(NA, NA, xlim=c(0,1), ylim=range(stages$stage)+c(-0.5, 0.5), ann=F, axes=F)
points(x=rep(0.1,nrow(stages)),y=max(stages$stage)-stages$stage,pch=22,cex=2,bg=stages$color,lwd=1)
text(x=rep(0.15,nrow(stages)),y=max(stages$stage)-stages$stage,pos=4,labels=paste0(stages$stage,' - ',stages$desc),font=1)



# barplot panel
par(mar=c(3,3,3,1))
trls %>%
  filter(incl) %>%
  group_by(earliest) %>%
  summarize(n=n()) %>%
  filter(!is.na(earliest)) %>%
  mutate(p = n/sum(n)) -> stagecount
barwidth=0.8
xlims = range(stages$stage)+c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=c(0, 2000), ann=F, axes=F, xaxs='i', yaxs='i')
rect(xleft=stages$stage-barwidth/2, xright=stages$stage + barwidth/2, ybottom=0, ytop=stagecount$n, col=stages$color, border=NA)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=0.25, at=stages$stage, text=stages$stage, cex=.75)
axis(side=2, at=0:4*500, labels=formatC(0:4*500,big.mark=','), las=2)
mtext(side=1, line=1.5, text='earliest disease stage', cex=0.6)
mtext(side=2, line=3.5, text='N trials', cex=0.6)
segments(x0=-0.3, x1=2.3, y0=250)
text(x=1, y=250, pos=3, labels=paste0(percent(sum(stagecount$n[stagecount$earliest < 3])/sum(stagecount$n)),' of trials'))
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

write_tsv(stagecount, 'qc/disease_stage_count.tsv')

write(paste('Disease stages: ',paste0('stage ',stagecount$earliest,' N = ',stagecount$n,' (',percent(stagecount$p),')', collapse='; '), sep=''),text_stats_path,append=T)

# forest of preventive difference-in-kind
trls %>%
  filter(incl) %>%
  mutate(preventive = earliest < 3) %>%
  mutate(preventdisp = ifelse(preventive, 'stage 0-2', 'stage 3-4')) %>%
  mutate(pspec = phs %in% 1:3,
         behav = iclass=='behavioral',
         proced = iclass=='procedure',
         drug = iclass=='molecular',
         device = iclass=='device',
         complete = status=='Completed') -> trls_prev_prep

forest_rows = tibble(vrbl=c('industry','drug','behav','device','proced','pspec','complete','has_placebo'),
                     disp=c('industry-sponsored','drug','behavioral','device','procedure','specified phase','completed','placebo-controlled'))

if(exists('prev_forest')) {
  rm(prev_forest)
}
for (vrbl in forest_rows$vrbl) {
  ctable = table(trls_prev_prep[,c('preventive',vrbl)])
  fisher_obj = fisher.test(ctable)
  forest_row = tibble(vrbl=vrbl, or=fisher_obj$estimate, l95=fisher_obj$conf.int[1], u95=fisher_obj$conf.int[2], p=fisher_obj$p.value)
  if (exists('prev_forest')) {
    prev_forest = rbind(prev_forest, forest_row)
  } else {
    prev_forest = forest_row
  }
}
prev_forest$disp = forest_rows$disp[match(prev_forest$vrbl, forest_rows$vrbl)]
prev_forest$y = nrow(prev_forest):1
xlims = c(0,3.25)
ylims = range(prev_forest$y) + c(-0.5, 0.5)

write_tsv(prev_forest, 'qc/preventive_trials_forest.tsv')

prev_forest %>%
  mutate(text=paste0(disp,' OR=',round(or,2),' P=',formatC(p,format='fg',digits=2))) -> prevtext

write(paste('Preventive trial differences: ',paste0(prevtext$text, collapse='; '), sep=''),text_stats_path,append=T)


par(mar=c(2,6,3,1))
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=0:10/2, labels=NA, tck=-0.03)
axis(side=1, at=0:5, labels=NA, tck=-0.06)
axis(side=1, at=0:5, labels=formatC(0:5,format='f',digits=1), lwd=0, line=-0.5)
mtext(side=1, at=1, line=1.5, text='odds ratio', cex=0.8)
mtext(side=1, at=c(0.5,1.5), line=0.35, text=c('depleted','enriched'), cex=0.6)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
abline(v=1, lty=3)
mtext(side=2, line=0.25, at=prev_forest$y, text=prev_forest$disp, las=2, cex=0.6)
points(prev_forest$or, prev_forest$y, pch=19)
segments(x0=prev_forest$l95, x1=prev_forest$u95, y0=prev_forest$y, lwd=2)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1


trls %>%
  filter(incl) %>%
  mutate(preventive = earliest < 3) %>%
  group_by(preventive) %>%
  summarize(.groups='keep',
            n_dur = sum(!is.na(study_duration)),
            n_enr = sum(!is.na(enrollment)),
            dur = mean(study_duration/365.25, na.rm=T),
            enr = mean(enrollment, na.rm=T)) -> prev_duration_enr

duration_ks_obj = suppressWarnings(ks.test(trls$study_duration[trls$earliest < 3]/365.25, trls$study_duration[trls$earliest >= 3]/365.25))
enrollment_ks_obj = suppressWarnings(ks.test(trls$enrollment[trls$earliest < 3], trls$enrollment[trls$earliest >= 3]))

write(paste0('Preventive trials duration: ',formatC(prev_duration_enr$dur[prev_duration_enr$preventive],format='g',digits=3),
      ' vs. ',formatC(prev_duration_enr$dur[!prev_duration_enr$preventive],format='g',digits=3),' years, P = ',
      formatC(duration_ks_obj$p.value, format='g', digits=2),', N = ',
      prev_duration_enr$n_dur[prev_duration_enr$preventive],' vs ',prev_duration_enr$n_dur[!prev_duration_enr$preventive],
      '\n'),text_stats_path,append=T)

write(paste0('Preventive trials duration: ',formatC(prev_duration_enr$enr[prev_duration_enr$preventive],format='g',digits=3),
             ' vs. ',formatC(prev_duration_enr$enr[!prev_duration_enr$preventive],format='g',digits=3),' years, P = ',
             formatC(enrollment_ks_obj$p.value, format='g', digits=2),', N = ',
             prev_duration_enr$n_enr[prev_duration_enr$preventive],' vs ',prev_duration_enr$n_enr[!prev_duration_enr$preventive],
             '\n'),text_stats_path,append=T)

durmax = 5
enrmax = 500
trls %>%
  filter(incl) %>%
  mutate(preventive = earliest < 3) %>%
  mutate(preventdisp = ifelse(preventive, 'stage 0-2', 'stage 3-4')) %>%
  mutate(y = ifelse(preventive, 2, 1)) %>%
  mutate(dur = pmin(durmax, study_duration/365.25),
         enr = pmin(enrmax, enrollment)) -> de_scatter

trls %>%
  mutate(preventive = earliest < 3) %>%
  mutate(preventdisp = ifelse(preventive, 'stage 0-2', 'stage 3-4')) %>%
  group_by(preventive, preventdisp) %>%
  summarize(.groups='keep',
            dur_mean = mean(study_duration/365.25, na.rm=T),
            dur_l95  = mean(study_duration/365.25, na.rm=T) - 1.96*sd(study_duration/365.25, na.rm=T)/sqrt(sum(!is.na(study_duration))),
            dur_u95  = mean(study_duration/365.25, na.rm=T) + 1.96*sd(study_duration/365.25, na.rm=T)/sqrt(sum(!is.na(study_duration))),
            enr_mean = mean(enrollment, na.rm=T),
            enr_l95  = mean(enrollment, na.rm=T) - 1.96*sd(enrollment, na.rm=T)/sqrt(sum(!is.na(enrollment))),
            enr_u95  = mean(enrollment, na.rm=T) + 1.96*sd(enrollment, na.rm=T)/sqrt(sum(!is.na(enrollment)))
  ) -> de_smry

s02_col = blend_hex(stages$color[1],stages$color[3],.5)
s34_col = '#A9A9A9' #blend_hex(stages$color[4],stages$color[5],.5)
de_smry$color = ifelse(de_smry$preventive, s02_col, s34_col)

par(mar=c(2,4,2,6))
xlims = c(0,durmax*1.05)
ylims = c(0,enrmax*1.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:10, labels=NA, tck=-0.03)
axis(side=1, at=c(0,2.5,5), labels=NA, tck=-0.06)
axis(side=1, at=c(0,2.5,5), labels=c('0','2.5','5+'), lwd=0, line=-0.5)
mtext(side=1, line=1.25, text='duration (y)', cex=0.6)
axis(side=2, at=0:10*100, labels=NA, tck=-0.025)
axis(side=2, at=c(0,250,500), labels=NA, tck=-0.05)
axis(side=2, at=c(0,250,500), labels=c('0','250','500+'), lwd=0, line=-0.5, las=2)
mtext(side=2, line=2.5, text='enrollment', cex=0.6)
points(de_scatter$dur[!de_scatter$preventive], de_scatter$enr[!de_scatter$preventive], pch=20, col=alpha(s34_col, .1))
points(de_scatter$dur[ de_scatter$preventive], de_scatter$enr[ de_scatter$preventive], pch=20, col=alpha(s02_col, .5))
segments(x0=de_smry$dur_l95, x1=de_smry$dur_u95, y0=de_smry$enr_mean, lwd=3, col=de_smry$color)
segments(x0=de_smry$dur_mean, y0=de_smry$enr_l95, y1=de_smry$enr_u95, lwd=3, col=de_smry$color)
par(xpd=T)
legend(x=max(xlims), y=max(ylims), legend=c('stages 0-2', 'stages 3-4'), pch=20, col=c(s02_col,s34_col), text.col=c(s02_col,s34_col), bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1


expand.grid(year=2000:2020, earliest=0:4) %>% 
  inner_join(stages, by=c('earliest'='stage')) %>%
  as_tibble() -> stage_years

stage_years %>%
  left_join(trls, by=c('year' = 'year', 'earliest' = 'earliest')) %>%
  filter(incl) %>%
  group_by(year, earliest, color) %>%
  summarize(.groups='keep',
            n=n()) -> stage_by_year_upto

scale_2020 = 365.25/as.numeric(as.Date('2020-03-31') - as.Date('2020-01-01'))
stage_by_year_upto$n[stage_by_year_upto$year==2020] = stage_by_year_upto$n[stage_by_year_upto$year==2020] * scale_2020

# ordinal logit for whether stage is changing by year
trls %>%
  mutate(stage = factor(earliest, levels=0:4, ordered=T)) -> polr_data
ordered_logit_model = polr(stage ~ year, data=polr_data, Hess=T)
year_t_statistic = coef(summary(ordered_logit_model))['year','t value']
p_value = pt(year_t_statistic, df=nrow(coef(summary(ordered_logit_model))))*2


write(paste0('Ordinal logit model: [disease] stage ~ year: P = ',formatC(p_value,format='g',digits=1),'\n'),text_stats_path,append=T)

# absolute proportions by beginning and end of period
trls %>%
  mutate(non4 = earliest <4) %>%
  filter(year > 2016) %>%
  group_by(non4) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  mutate(p = n/sum(n)) -> lastfour
trls %>%
  mutate(non4 = earliest <4) %>%
  filter(year < 2004) %>%
  group_by(non4) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  mutate(p = n/sum(n)) -> firstfour

write(paste0('Trials at disease stage < 4: ',percent(firstfour$p[firstfour$non4]),' in 2000-2003, ',percent(lastfour$p[firstfour$non4]),' in 2017-2020, ','\n'),text_stats_path,append=T)

trls %>%
  mutate(score_test=replace_na(score_test,'none')) %>%
  group_by(score_test) %>%
  summarize(.groups='keep', n_trials=n()) %>%
  ungroup() %>%
  mutate(proportion = n_trials/sum(n_trials)) -> count_by_score_test

any_test = sum(count_by_score_test$n_trials) - count_by_score_test$n_trials[count_by_score_test$score_test=='none']

write(paste0('Score/test: N=',any_test,' (',percent(any_test/sum(count_by_score_test$n_trials)),') trials used any.\n'),text_stats_path,append=T)
write(paste0('Score/test: N=',count_by_score_test$n_trials[count_by_score_test$score_test=='MMSE'],' (',percent(count_by_score_test$proportion[count_by_score_test$score_test=='MMSE']),') trials used MMSE.\n'),text_stats_path,append=T)
write(paste0('Score/test: N=',count_by_score_test$n_trials[count_by_score_test$score_test=='Hoehn and Yahr'],' (',percent(count_by_score_test$proportion[count_by_score_test$score_test=='Hoehn and Yahr']),') trials used Hoehn and Yahr.\n'),text_stats_path,append=T)





par(mar=c(4,4,3,1))
xlims = c(2000, 2020)
ylims = c(0, ceiling(max(stage_by_year_upto$n, na.rm=T)/10)*10)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(2000,2020,1), tck=-0.025, labels=NA)
axis(side=1, at=seq(2000,2020,5), tck=-0.05)
mtext(side=1, line=2, text='year', cex=0.8)
axis(side=2, at=0:10*50, tck=-0.03, labels=NA)
axis(side=2, at=0:5*100, tck=-0.06, labels=NA)
axis(side=2, at=0:5*100, las=2, lwd=0, line=0, cex.axis=0.8)
mtext(side=2, line=2.5, text='trials/year', cex=0.8)
for (s in 4:0) {
  subs = subset(stage_by_year_upto, earliest==s)
  polygon(x=c(subs$year, rev(subs$year)), y=c(rep(0,nrow(subs)), rev(subs$n)), col=subs$color, border=NA)
}
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

par(mar=c(4,4,3,3))
trls %>%
  filter(incl) %>%
  filter(score_test=='MMSE') %>%
  mutate(score_most = ifelse(is.na(score_most),0,score_most),
         score_least = ifelse(is.na(score_least),30,score_least)) -> mmse_trls
mmse_least = lm(score_least ~ year, data=mmse_trls)
mmse_most = lm(score_most ~ year, data=mmse_trls)

write(paste0('MMSE: least impaired ',paste(formatC(predict.lm(mmse_least,data.frame(year=c(2000,2020))),format='f',digits=1), collapse=' to '),
             ' P = ',formatC(summary(mmse_least)$coefficients['year','Pr(>|t|)'], format='g', digits=2),
             ', most impaired ',paste(formatC(predict.lm(mmse_most,data.frame(year=c(2000,2020))),format='f',digits=1), collapse=' to '),
             ' P = ',formatC(summary(mmse_most)$coefficients['year','Pr(>|t|)'], format='g', digits=2),'\n'),text_stats_path,append=T)



set.seed(1)
mmse_trls$x = jitter(mmse_trls$year, amount=0.5)
barwidth = 0.25
barcol = '#810262'
baralpha = 0.1
fitcol = '#008837'

xlims = c(2000, 2020)
ylims = c(0, 30)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(2000,2020,1), tck=-0.025, labels=NA)
axis(side=1, at=seq(2000,2020,5), tck=-0.05)
mtext(side=1, line=2, text='year', cex=0.8)
axis(side=2, at=0:30, labels=NA, tck=-0.025)
axis(side=2, at=0:6*5, labels=NA, tck=-0.05)
axis(side=2, at=0:6*5, las=2, lwd=0, lwd.ticks=0, line=-0.5)
mtext(side=2, line=2.5, text='MMSE')
par(xpd=T)
arrows(y0=29, y1=1, x0=1996, code=2, angle=45, length=0.025)
text(x=1996, y=29, pos=2, srt=90, labels='normal', cex=0.6)
text(x=1996, y=10, pos=2, srt=90, labels='more impaired', cex=0.6)
par(xpd=F)
rect(xleft=mmse_trls$x-barwidth/2, xright=mmse_trls$x+barwidth/2, ybottom=mmse_trls$score_most, ytop=mmse_trls$score_least, col=alpha(barcol, baralpha), border=NA)
abline(mmse_least, col=fitcol, lwd=2)
abline(mmse_most, col=fitcol, lwd=2)
max_pos = mmse_least$coefficients['(Intercept)'] + 2020*mmse_least$coefficients['year']
min_pos = mmse_most$coefficients['(Intercept)'] + 2020*mmse_most$coefficients['year']
mtext(side=4, at=max_pos, col=fitcol, las=2, line=0.25, text='max')
mtext(side=4, at=min_pos, col=fitcol, las=2, line=0.25, text='min')
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1
# text(x=2010, 
#      y=predict.lm(mmse_least,data.frame(year=2010)),
#      pos=3, 
#      srt=mmse_least$coefficients[2]*45, 
#      labels=paste(formatC(mmse_least$coefficients[2],format='fg',digits=1),'/y, P=',formatC(summary(mmse_least)$coefficients[2,4],format='fg',digits=2)),
#      col=fitcol)

trls %>%
  filter(incl) %>%
  filter(score_test=='Hoehn and Yahr') %>%
  mutate(score_most = ifelse(is.na(score_most),5,score_most),
         score_least = ifelse(is.na(score_least),0,score_least)) -> hy_trls


hy_least = lm(score_least ~ year, data=hy_trls)
hy_most = lm(score_most ~ year, data=hy_trls)

write(paste0('H&Y: least impaired ',paste(formatC(predict.lm(hy_least,data.frame(year=c(2000,2020))),format='f',digits=1), collapse=' to '),
             ' P = ',formatC(summary(hy_least)$coefficients['year','Pr(>|t|)'], format='g', digits=2),
             ', most impaired ',paste(formatC(predict.lm(hy_most,data.frame(year=c(2000,2020))),format='f',digits=1), collapse=' to '),
             ' P = ',formatC(summary(hy_most)$coefficients['year','Pr(>|t|)'], format='g', digits=2),'\n'),text_stats_path,append=T)




set.seed(1)
hy_trls$x = jitter(hy_trls$year, amount=0.5)
barwidth = 0.25
barcol = '#810262'
baralpha = 0.1
fitcol = '#008837'

xlims = c(2000, 2020)
ylims = c(0, 5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(2000,2020,1), tck=-0.025, labels=NA)
axis(side=1, at=seq(2000,2020,5), tck=-0.05)
mtext(side=1, line=2, text='year', cex=0.8)
axis(side=2, at=seq(0,5,.5), labels=NA, tck=-0.025)
axis(side=2, at=0:5, labels=NA, tck=-0.05)
axis(side=2, at=0:5, las=2, lwd=0, lwd.ticks=0, line=-0.5)
mtext(side=2, line=2.5, text='Hoehn and Yahr')
par(xpd=T)
arrows(y0=0.15, y1=4.85, x0=1996, code=2, angle=45, length=0.025)
text(x=1996, y=0.75, pos=2, srt=90, labels='normal', cex=0.6)
text(x=1996, y=4.9, pos=2, srt=90, labels='more impaired', cex=0.6)
par(xpd=F)
rect(xleft=hy_trls$x-barwidth/2, xright=hy_trls$x+barwidth/2, ybottom=hy_trls$score_most, ytop=hy_trls$score_least, col=alpha(barcol, baralpha), border=NA)
abline(hy_least, col=fitcol, lwd=2)
abline(hy_most, col=fitcol, lwd=2)
min_pos = hy_least$coefficients['(Intercept)'] + 2020*hy_least$coefficients['year']
max_pos = hy_most$coefficients['(Intercept)']  + 2020*hy_most$coefficients['year']
mtext(side=4, at=max_pos, col=fitcol, las=2, line=0.25, text='max')
mtext(side=4, at=min_pos, col=fitcol, las=2, line=0.25, text='min')
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

silence_message = dev.off()











###############
#### FIGURE 4 # 
###############

cat(file=stderr(), 'done!\nCreating Figure 4..'); flush.console()

# eulerr uses grid and so does not respect multi-panel layouts in base R
# need to go through the analog hole to force the euler plot into this layout.
png('display_items/figure-4c.png',width=4.5*resx,height=4.5*resx,res=resx)
genes$drugbank = genes$gene %in% drugbank$target[drugbank$status %in% 'approved']
genes$developed = genes$gene %in% ig_match$target_gene
genes$supported = genes$gene %in% targets$gene
euler_obj = euler(genes[,c('drugbank','supported','developed')])
plot(euler_obj, quantities=T,
     labels=c('approved','supported','developed'),
     pos=2, 
     fill=c('#A9A9A959','#EEAD0E59','#D0A9AA59'))
developed_drugbank = table(genes[,c('developed','drugbank')])
fisher_developed_approved = fisher.test(developed_drugbank)
developed_supported = table(genes[,c('developed','supported')])
fisher_developed_supported = fisher.test(developed_supported)
silence_message = dev.off()
figure_4c = image_convert(image_read('display_items/figure-4c.png'),'png')
silence_message = file.remove('display_items/figure-4c.png')



write(paste0('Enrichment for approved targets: OR=',round(fisher_developed_approved$estimate,1),
             ' P=',formatC(fisher_developed_approved$p.value,format='f',digits=2),'\n'),text_stats_path,append=T)
write(paste0('Enrichment for supported targets: OR=',round(fisher_developed_supported$estimate,1),
             ' P=',formatC(fisher_developed_supported$p.value,format='f',digits=2),'\n'),text_stats_path,append=T)


png('display_items/figure-4.png',width=6.5*resx,height=6.5*resx,res=resx)


layout_matrix = matrix(c(1,1,1,1,2,2,
                         3,3,4,4,4,5), nrow=2, byrow=T)
layout(layout_matrix, heights=c(1,1), widths=c(1,1,1,1,.8,1.5))

panel = 1


write(paste0('Drug trials: N=',sum(trls$iclass=='molecular'),'\n'),text_stats_path,append=T)

ivns %>%
  filter(intervention_class=='molecular') %>%
  inner_join(tim, by=c('intervention_name','intervention_type')) %>%
  inner_join(trls, by='nct') %>%
  select(best_name) %>% pull() %>% unique -> drug_names
  
write(paste0('Unique molecular entities: N=',length(drug_names),'\n'),text_stats_path,append=T)


tirc_meta = tibble(priority = c(6, 7, 1, 3, 4, 5, 2),
                   shortname = c('pre-approval','post-approval','novel, supported','repurposing','combination','metoo','novel, unsup'),
                   yorder = c(7,6,2,5,3,4,1),
                   color = c('#00E5EE','#CDCDCD','#EEAD0E','#72587F','#CD919E','#D0A9AA','#9D8851'),
                   tranche = c('approved','approved','novel','repurposing','repurposing','repurposing','novel'),
                   disp = c('pre-approval','post-approval','novel genetically\nsupported targrets','repurposed drugs','combinations','repurposed targets','novel other hypotheses'))



trls %>%
  filter(incl) %>%
  filter(iclass=='molecular') %>%
  mutate(patient_years = enrollment * study_duration/365.25) %>%
  inner_join(tim, by=c('nct')) %>%
  inner_join(ivns, by=c('intervention_name','intervention_type','intervention_class')) %>%
  filter(intervention_class=='molecular') %>%
  select(-target_gene) %>%
  left_join(approved, by=c('disease_area','best_name'='drug')) %>%
  left_join(ig_match, by=c('intervention_name'='intervention_name')) %>%
  left_join(targets, by=c('disease_area','target_gene'='gene')) %>%
  select(disease_area, nct, year, intervention_name, drug=best_name, regulatory_status, approval_year, patient_years, gen_assoc=source, target_gene) -> tir

ig_match$target_status = ifelse(ig_match$target_gene %in% drugbank$target[drugbank$status %in% 'approved'], 'approved', 'non')
tir$target_status = ifelse(tir$intervention_name %in% ig_match$intervention_name[ig_match$target_status=='approved'], 'approved', 'non')




tir %>%
  mutate(shortname = case_when(!is.na(approval_year) & year < approval_year  ~ 'pre-approval',
                               !is.na(approval_year) & year >= approval_year  ~ 'post-approval',
                               !(regulatory_status %in% 'approved') & target_status=='non' & !is.na(gen_assoc) ~ 'novel, supported',
                               regulatory_status %in% 'approved' & is.na(approval_year) ~ 'repurposing',
                               regulatory_status %in% 'combination' ~ 'combination',
                               !(regulatory_status %in% 'approved') & target_status=='approved' ~ 'metoo',
                               !(regulatory_status %in% 'approved') & target_status=='non' & is.na(gen_assoc) ~ 'novel, unsup',
                               TRUE ~ 'other')) %>%
  inner_join(select(tirc_meta, priority, shortname), by=c('shortname')) -> tirc

tirc %>%
  group_by(nct) %>%
  arrange(priority) %>%
  slice(1) %>%
  ungroup()  %>%
  rename(classification=shortname) -> tirc1

write_tsv(tirc1, 'qc/drug_trial_classifications.tsv', na='')

#  QC of post-approval trials:
tirc1 %>%
  filter(classification=='post-approval') %>%
  inner_join(select(trials, nct, brief_title),by='nct') %>%
  select(nct, disease_area, year, drug, approval_year, brief_title) -> post_approval_qc
write_tsv(post_approval_qc, 'qc/post_approval_trials.tsv')


tirc1 %>% 
  filter(classification=='repurposing') %>% 
  group_by(drug) %>% 
  summarize(.groups='keep', n=n()) %>% 
  arrange(desc(n)) -> repurposing_count

write_tsv(repurposing_count, 'qc/repurposing_trial_count.tsv', na='')

tirc1 %>%
  group_by(priority, classification) %>%
  summarize(.groups='keep',
            n=n(),
            py = sum(patient_years, na.rm=T)) %>%
  ungroup() %>%
  mutate(pn = n/sum(n),
         ppy = py/sum(py)) -> tirc_smry

write_tsv(tirc_smry, 'qc/drug_trial_classifications_summary.tsv')


write(paste('Pre-approval and post-approval: total proportion PY = ',percent(sum(tirc_smry$ppy[tirc_smry$classification %in% c('pre-approval','post-approval')])),'\n',sep=''),text_stats_path,append=T)

write(paste('Pre-approval: proportion PY = ',percent(sum(tirc_smry$ppy[tirc_smry$classification %in% c('pre-approval')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Post-approval: proportion PY = ',percent(sum(tirc_smry$ppy[tirc_smry$classification %in% c('post-approval')])),'\n',sep=''),text_stats_path,append=T)

write(paste('Pre-approval vs. post-approval: N trials = ',
            paste0(c(tirc_smry$n[tirc_smry$classification=='pre-approval'], tirc_smry$n[tirc_smry$classification=='post-approval']),collapse=' vs '),
            ' (ratio: ',round(tirc_smry$n[tirc_smry$classification=='post-approval']/tirc_smry$n[tirc_smry$classification=='pre-approval'],1),
            ')',' PY = ',paste0(formatC(c(tirc_smry$py[tirc_smry$classification=='pre-approval'], tirc_smry$py[tirc_smry$classification=='post-approval']),big.mark=',',format='f',digits=0),collapse=' vs '),
            ' (ratio: ',round(tirc_smry$py[tirc_smry$classification=='post-approval']/tirc_smry$py[tirc_smry$classification=='pre-approval'],1),'\n',sep=''),text_stats_path,append=T)

write(paste('Donepezil N trials = ',sum(tirc1$drug=='donepezil'),'\n',sep=''),text_stats_path,append=T)

write(paste('Novel hypotheses: total proportion PY = ',percent(sum(tirc_smry$ppy[tirc_smry$classification %in% c('novel, supported','novel, unsup')])),'\n',sep=''),text_stats_path,append=T)
write(paste('Novel hypotheses: fraction lacking genetic support = ',
            paste0(c(tirc_smry$n[tirc_smry$classification=='novel, unsup'], sum(tirc_smry$n[tirc_smry$classification %in% c('novel, supported','novel, unsup')])),collapse='/')
            ,'\n',sep=''),text_stats_path,append=T)

write(paste('Genetically supported vs. unsupported novel hypotheses: N trials = ',
            paste0(c(tirc_smry$n[tirc_smry$classification=='novel, supported'], tirc_smry$n[tirc_smry$classification=='novel, unsup']),collapse=' vs '),
            ' (ratio: ',round(tirc_smry$n[tirc_smry$classification=='novel, unsup']/tirc_smry$n[tirc_smry$classification=='novel, supported'],1),
            ')',' PY = ',paste0(formatC(c(tirc_smry$py[tirc_smry$classification=='novel, supported'], tirc_smry$py[tirc_smry$classification=='novel, unsup']),big.mark=',',format='f',digits=0),collapse=' vs '),
            ' (ratio: ',round(tirc_smry$py[tirc_smry$classification=='novel, unsup']/tirc_smry$py[tirc_smry$classification=='novel, supported'],1),'\n',sep=''),text_stats_path,append=T)


tirc1 %>% 
  filter(classification=='novel, unsup') %>% 
  group_by(disease_area, target_gene) %>% 
  summarize(.groups='keep', n=n()) %>% 
  arrange(desc(n)) -> unsup_count


write_tsv(unsup_count, 'qc/unsupported_targets_count.tsv', na='')


tirc1 %>%
  filter(classification=='novel, unsup' & is.na(target_gene)) %>%
  select(drug) %>% pull() %>% unique() %>% length() -> unique_notarget_unsup_drugs


write(paste('Novel unsupported targets: ',sum(unsup_count$n[is.na(unsup_count$target_gene)]),'/',sum(unsup_count$n),
            ' (',percent(sum(unsup_count$n[is.na(unsup_count$target_gene)])/sum(unsup_count$n)),')',
            ' no target assigned (',unique_notarget_unsup_drugs,' unique drugs); ',unsup_count$n[unsup_count$disease_area=='alzheimer' & unsup_count$target_gene %in% 'MAPT'],'  MAPT; ',
            unsup_count$n[unsup_count$disease_area=='alzheimer' & unsup_count$target_gene %in% 'BACE1'],'  BACE1','\n',sep=''),text_stats_path,append=T)



tirc_smry %>%
  inner_join(tirc_meta, by=c('priority','classification'='shortname')) %>%
  arrange(desc(yorder)) %>%
  mutate(pyper=py/n,
         ppyper = pyper/sum(pyper),
         ppyper_cum = cumsum(ppyper),
         ybottom = max(yorder) * (1-ppyper_cum),
         ytop = max(yorder) * (1-ppyper_cum + ppyper)) -> tirc_smry_o

tirc_smry_o %>%
  group_by(tranche) %>%
  summarize(.groups='keep', ymid = (min(ybottom) + max(ytop))/2, ymax=max(ytop), ymin=min(ybottom)) -> tirc_ymeta

par(mar=c(4,5,3,2))
ylims = c(0.0, 7.0)
xlims = c(0, 1000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:16*50, tck=-0.025, labels=NA)
axis(side=1, at=0:8*100, tck=-0.05)
mtext(side=1, line=2.0, text='N trials', cex=0.6)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
rect(xleft=rep(0,nrow(tirc_smry_o)), xright=tirc_smry_o$n, ybottom=tirc_smry_o$ybottom, ytop=tirc_smry_o$ytop, border=NA, col=tirc_smry_o$color)
mtext(side=2, at=tirc_ymeta$ymid, line=0.25, text=tirc_ymeta$tranche, las=2, cex=0.6)
par(xpd=T)
ybreaks = unique(c(tirc_ymeta$ymax, tirc_ymeta$ymin))
segments(x0=rep(-100,length(ybreaks)), x1=rep(800, length(ybreaks)), y0=ybreaks, y1=ybreaks, lwd=0.125)
par(xpd=F)
text(x=tirc_smry_o$n, y=(tirc_smry_o$ybottom + tirc_smry_o$ytop)/2, pos=4, labels=tirc_smry_o$disp, cex=0.8)
text(x=tirc_smry_o$n/2, y=(tirc_smry_o$ybottom + tirc_smry_o$ytop)/2, labels=percent(tirc_smry_o$ppy,digits=0), font=2, col="#FFFFFF", cex=0.8)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1


tirc %>%
  filter(!is.na(approval_year)) %>%
  group_by(drug, disease_area, year, approval_year) %>%
  summarize(.groups='keep', n = length(unique(nct))) -> approved_n_by_year

tirc %>%
  filter(!is.na(approval_year)) %>%
  group_by(drug, disease_area, year, approval_year) %>%
  summarize(.groups='keep', n_nct = length(unique(nct))) %>%
  inner_join(diseases, by=c('disease_area')) %>%
  mutate(color=ifelse(year < approval_year, tirc_meta$color[tirc_meta$shortname=='pre-approval'], tirc_meta$color[tirc_meta$shortname=='post-approval'])) %>%
  group_by(-axis_order, desc(drug)) %>%
  mutate(y = cur_group_id()) -> abyy


abyy %>%
  mutate(timing=ifelse(year < approval_year, 'pre', 'post')) %>%
  group_by(drug, approval_year, timing) %>%
  summarize(.groups='keep', n_trials = sum(n_nct)) -> abyy_summary

write_tsv(abyy_summary, 'qc/trials_pre_and_post_approval.tsv')

abyy %>%
  group_by(y, drug, approval_year, disease_area) %>%
  summarize(.groups='keep', n=n()) -> ymeta

ymeta %>%
  inner_join(diseases, by=c('disease_area')) %>%
  group_by(disease_area, disp) %>%
  summarize(.groups='keep',
            midy=mean(y),
            miny=min(y),
            maxy=max(y)) -> yaxes

par(mar=c(4,7,3,1))
xlims = c(1999.5, 2020.5)
ylims = range(abyy$y) + c(-0.5, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=seq(2000,2020,5), tck=-0.05)
axis(side=1, at=2000:2020, labels=NA, tck=-0.025)
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
mtext(side=2, line=0.25, at=ymeta$y+0.5, text=ymeta$drug, las=2, cex=0.6)
overhang = 0.4
for (i in 1:nrow(yaxes)) {
  axis(side=2, line=6.5, tck=0.05, at=unlist(yaxes[i,c('miny','maxy')])+c(-1,1)*overhang+0.5, labels=NA)
  mtext(side=2, line=7, at=yaxes$midy[i]+0.5, text=yaxes$disp[i], las=2, cex=0.7)
}
abline(h=ymeta$y, lwd=0.125)
barwidth=0.4
rect(xleft=abyy$year-barwidth, xright=abyy$year+barwidth, ybottom=abyy$y, ytop=abyy$y + abyy$n_nct / max(abyy$n_nct), col=abyy$color, border=NA)
points(x=pmax(ymeta$approval_year,1999.5), y=ymeta$y+0.2, pch=17, cex=0.8)
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

par(mar=c(0,0,0,0))
plot(as.raster(figure_4c))
mtext(LETTERS[panel], side=3, cex=1.5, adj=0.2, line =-1); panel = panel + 1


write(paste('Repurposing categories: ',percent(sum(tirc_smry$ppy[tirc_smry$classification %in% c('repurposing','metoo','combination')])),'\n'),text_stats_path,append=T)

write(paste('Targets tested: ',length(unique(genes$gene[genes$developed])),'\n'),text_stats_path,append=T)
write(paste('Repurposed targets: ',length(unique(genes$gene[genes$developed & genes$drugbank])),'\n'),text_stats_path,append=T)

drugbank %>%
  filter(status=='approved' & !is.na(atc)) %>%
  select(drug, atc) %>%
  group_by(drug, atc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(atc1 = substr(atc,1,1)) %>%
  left_join(atc_summary, by=c('atc1'='level1')) %>%
  mutate(atc_name = replace_na(disease_area,'other')) %>%
  select(drug, atc, atc1, atc_name) -> db_atc

approved %>%
  select(approved_dz=disease_area, drug) %>%
  right_join(tirc1, by='drug') %>%
  mutate(approved_dz = replace_na(approved_dz, 'other')) %>%
  left_join(db_atc, by='drug') %>%
  mutate(atc_name = case_when(drug %in% 'insulin' ~ 'metabolic',
                              is.na(atc_name) ~ 'other',
                              atc_name %in% c('respiratory','skeletomuscular','dermatological','hematologic') ~ 'other',
                              approved_dz != 'other' ~ 'other neurodegenerative',
                              atc_name == 'neurology' ~ 'other neurological',
                              TRUE ~ atc_name)) %>%
  filter(classification=='repurposing') -> repurposing 

repurposing %>% filter(atc_name=='other neurodegenerative') %>% group_by(drug, disease_area) %>% summarize(.groups='keep', n_nct=n()) -> between_nd_pivots

write_tsv(between_nd_pivots, 'qc/pivots_between_neurodegenerative_diseases.tsv')

write(paste('Pivots from one neurodegenerative disease to another: ',sum(between_nd_pivots$n_nct),'\n'),text_stats_path,append=T)

repurposing %>%
  group_by(atc_name) %>%
  summarize(.groups='keep', n_nct = n()) -> repur_atc

write(paste('Other repurposing sources: neurology N=',repur_atc$n_nct[repur_atc$atc_name=='other neurological'],
            ' metabolic N=',repur_atc$n_nct[repur_atc$atc_name=='metabolic'],
            ' cardiovascular N=',repur_atc$n_nct[repur_atc$atc_name=='cardiovascular'],'\n'),text_stats_path,append=T)

write(paste('All non-ND repurposing sources: N=',sum(repur_atc$n_nct[!(repur_atc$atc_name %in% 'other neurodegenerative')]),
            '\n'),text_stats_path,append=T)

repurposing %>%
  group_by(drug) %>%
  summarize(.groups='keep', n=n()) %>%
  arrange(desc(n)) -> top_repur_drugs

write(paste('Top repurposed drug: ',top_repur_drugs$drug[1],', N=',top_repur_drugs$n[1],
            '\n'),text_stats_path,append=T)


ig_match %>%
  inner_join(tim, by=c('intervention_name','intervention_type')) %>%
  filter(intervention_class=='molecular') %>%
  inner_join(trls, by=c('nct')) %>%
  left_join(targets, by=c('target_gene'='gene', 'disease_area')) %>%
  mutate(patient_years = enrollment * study_duration/365.25) %>%
  filter(!is.na(source)) %>%
  inner_join(diseases, by='disease_area') %>%
  inner_join(select(interventions,intervention_name,intervention_type,best_name), by=c('intervention_name','intervention_type')) %>%
  filter(!target_gene %in% c('PSEN2','APH1B')) %>%
  mutate(target_gene = ifelse(target_gene=='PSEN1','PSEN1*',target_gene)) %>%
  group_by(target_gene, disease_area, axis_order, assoc_year, color) %>%
  summarize(.groups='keep', 
            fih_year = min(year),
            drugs = toString(unique(best_name)),
            n_trls = length(unique(nct)),
            n_i = length(unique(nct[phs==1])),
            n_ii = length(unique(nct[phs==2])),
            n_iii = length(unique(nct[phs==3])),
            n_o = length(unique(nct[phs %in% c(-1, 4)])),
            py = sum(patient_years, na.rm=T),
            py_i = sum(patient_years[phs==1], na.rm=T),
            py_ii = sum(patient_years[phs==2], na.rm=T),
            py_iii = sum(patient_years[phs==3], na.rm=T),
            py_o = sum(patient_years[phs %in% c(-1, 4)], na.rm=T)) %>%
  ungroup() %>%
  mutate(n_prop = n_trls/sum(n_trls),
         py_prop = py/sum(py)) %>%
  arrange(axis_order, target_gene) -> t_adv
t_adv$y = nrow(t_adv):1

write_tsv(t_adv, 'qc/target_advancement.tsv')


write(paste('Gene-disease links identified: ',nrow(targets),'\n'),text_stats_path,append=T)

write(paste('Novel supported hypotheses: ',nrow(t_adv),'\n'),text_stats_path,append=T)

# a lot of these trials don't have final enrollment and/or duration
# numbers so many cells have pyper==0, so this plot is nonsensical:
# t_adv %>%
#   mutate(pyper = py/n_trls,
#          ppyper = pyper/sum(pyper),
#          ppyper_cum = cumsum(ppyper),
#          ybottom = max(y)*(1-ppyper_cum),
#          ytop = max(y)*(1-ppyper_cum + ppyper)) -> t_adv

t_adv %>%
  filter(fih_year > assoc_year) -> prospective 

write(paste('Prospective trials based on genetics: mean ',round(mean(prospective$fih_year - prospective$assoc_year)),' years, min ',round(min(prospective$fih_year - prospective$assoc_year)),'\n'),text_stats_path,append=T)

t_adv$pos = 4
t_adv$pos[t_adv$target_gene %in% c('PSEN1','BST1','GCH1','APOE')] = 3
t_adv$pos[t_adv$target_gene %in% c('C9orf72','LRRK2')] = 1
t_adv$pos[t_adv$target_gene %in% c('APH1B','GBA')] = 2
par(mar=c(3,5,2,4))
lims = c(1990, 2020)
plot(NA, NA, xlim=lims, ylim=lims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=1990:2020, tck=-0.025, labels=NA)
axis(side=1, at=seq(1990,2020,10), labels=NA, tck=-0.05)
axis(side=1, at=seq(1990,2020,10), lwd=0, line=-0.25)
mtext(side=1, line=1.75, text='genetic association', cex=0.8)
axis(side=2, at=1990:2020, tck=-0.025, labels=NA)
axis(side=2, at=seq(1990,2020,10), tck=-0.05, labels=NA)
axis(side=2, at=seq(1990,2020,10), lwd=0, line=-0.25, las=2)
mtext(side=2, line=2.5, text='first trial', cex=0.8)
points(t_adv$assoc_year, t_adv$fih_year, pch=19, col=t_adv$color)
par(xpd=T)
text(x=t_adv$assoc_year, y=t_adv$fih_year, pos=t_adv$pos, font=3, labels=t_adv$target_gene, col=t_adv$color)
par(xpd=F)
abline(a=0, b=1, lwd=0.25)
legend('bottomright', diseases$disp[1:4], col=diseases$color[1:4], text.col=diseases$color[1:4], pch=19, bty='n')
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

tirc1 %>%
  filter(classification=='novel, supported') %>%
  group_by(target_gene) %>%
  summarize(.groups='keep',
            n_trials = n(),
            sum_py = sum(patient_years, na.rm=T)) %>%
  ungroup() %>%
  mutate(proportion_trials = n_trials/sum(n_trials),
         proportion_py = sum_py / sum(sum_py)) -> genetic_targets_effort

write_tsv(genetic_targets_effort, 'qc/genetic_targets_effort.tsv')

write(paste('APP: ',percent(genetic_targets_effort$proportion_trials[genetic_targets_effort$target_gene=='APP']),
            ' of trials, ',percent(genetic_targets_effort$proportion_py[genetic_targets_effort$target_gene=='APP']),
            ' of PY','\n',sep=''),text_stats_path,append=T)

bg_color = '#000000'
par(mar=c(3,5,2,1))
xlims = c(0, 120)
ylims = range(t_adv$y) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:50*10, labels=NA, tck=-0.025)
axis(side=1, at=0:10*50, labels=NA, tck=-0.05)
axis(side=1, at=0:5*100, lwd=0, line=-0.25)
mtext(side=1, line=2, text='N trials',cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, at=t_adv$y, text=t_adv$target_gene, font=3, las=2, cex=0.6)
barwidth = 0.4
rect(xleft=rep(0,nrow(t_adv)), xright=t_adv$n_trls, ybottom=t_adv$y-barwidth, ytop=t_adv$y+barwidth, col=t_adv$color, border=NA)
t_adv %>%
  inner_join(diseases, by='disease_area') %>%
  group_by(disp) %>%
  summarize(.groups='keep', miny=min(y), maxy=max(y), midy=mean(y)) -> ymeta
overhang = 0.4
for (i in 1:nrow(ymeta)) {
  axis(side=2, line=3.5, at=c(ymeta$miny[i]-overhang, ymeta$maxy[i]+overhang), tck=0.025, labels=NA)
  mtext(side=2, line=4, at=ymeta$midy[i], text=ymeta$disp[i], las=2, cex=0.7)
}
mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.2, line = 0.5); panel = panel + 1

trls[trls$earliest %in% 0:2 & trls$nct %in% tirc1$nct[tirc1$classification=='novel, supported'],] -> prev_gsup_trials

write(paste('Preventive, genetically supported trials: ',paste(prev_gsup_trials$nct, collapse=', '),'\n'),text_stats_path,append=T)

silence_message = dev.off()


cat(file=stderr(), paste0('done!\nAll tasks completed in ',round(as.numeric(Sys.time() - start_time),1),' seconds.\n')); flush.console()


