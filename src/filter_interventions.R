options(stringsAsFactors = F)
if (interactive()) {
  setwd('~/d/sci/src/neurology_prevention/')
}

library(sqldf)

drugbank = read_tsv('data/drugbank.tsv', guess_max=5000)

trials = read.table('data/curated/trials.tsv',sep='\t',header=T,quote='',comment.char='')
interventions_grouped = read.table('output/interventions_grouped.tsv',sep='\t',header=T,quote='',comment.char='')
tim = read.table('output/trial_intervention_match.tsv',sep='\t',header=T,quote='',comment.char='')
tim$intervention_name = tolower(tim$intervention_name)
table(trials$dataset_inclusion)
mean(tim$intervention_name %in% interventions_grouped$intervention_name)
mean(interventions_grouped$intervention_name %in% tim$intervention_name)
mean(tim$nct %in% trials$nct)
mean(trials$nct %in% tim$nct)

interventions_grouped$db_target = drugbank$target[match(interventions_grouped$intervention_name, drugbank$drug)]
interventions_grouped$db_status = drugbank$status[match(interventions_grouped$intervention_name, drugbank$drug)]

incl_count = sqldf("
select   i.intervention_name, sum(case when t.dataset_inclusion in ('i - include','') then 1 else 0 end) n_incl
from     trials t, tim ti, interventions_grouped i
where    t.nct = ti.nct and ti.intervention_name = i.intervention_name
group by 1
order by 1
;")

mean(interventions_grouped$intervention_name %in% incl_count$intervention_name)
sum(duplicated(interventions_grouped$intervention_name)) # ah some are dups b/c one trial called them a device and the other a drug, e.g.

interventions_grouped$incl_count = incl_count$n_incl[match(interventions_grouped$intervention_name, incl_count$intervention_name)]

write.table(interventions_grouped,'output/interventions_grouped_with_count.tsv',sep='\t',col.names=T,row.names=F,quote=F,na='')
