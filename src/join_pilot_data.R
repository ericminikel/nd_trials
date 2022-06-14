options(stringsAsFactors = F)
if (interactive()) {
  setwd('~/d/sci/src/neurology_prevention/')
}

library(sqldf)

trials = read.table('output/trials.tsv',sep='\t',header=T,quote='',comment.char='')
interventions_grouped = read.table('output/interventions_grouped.tsv',sep='\t',header=T,quote='',comment.char='')

pilot_trials = read.table('data_received/pilot_trials_2020-04-16.tsv',sep='\t',header=T,quote='',comment.char='')
pilot_interventions_grouped = read.table('data_received/pilot_interventions_grouped_2020-04-16.tsv',sep='\t',header=T,quote='',comment.char='')

curated_columns = c('dataset_inclusion','exclude_comment','disease_area','sponsor_type','earliest_disease_stage','latest_disease_stage','advanced_exclusion','advanced_exclusion_comments','comments')
trial_match_indices = match(trials$nct, pilot_trials$nct)
for (column in curated_columns) {
  trials[,column] = pilot_trials[trial_match_indices,column]
}

# double check if any got lost
#sum(!(pilot_trials$nct %in% trials$nct))
#pilot_trials[!(pilot_trials$nct %in% trials$nct),] # just 1 - because it had not launchd as of 3/31/2020. ok to lose this one

ivn_curated_columns = c('intervention_class','target_gene','target_gene_comments')
ivn_match_indices = match(interventions_grouped$intervention_name, tolower(pilot_interventions_grouped$intervention_name))
for (column in ivn_curated_columns) {
  interventions_grouped[,column] = pilot_interventions_grouped[ivn_match_indices,column]
}

write.table(trials,'output/trials.tsv',sep='\t',col.names=T,row.names=F,quote=F,na='')
write.table(interventions_grouped,'output/interventions_grouped.tsv',sep='\t',col.names=T,row.names=F,quote=F,na='')


