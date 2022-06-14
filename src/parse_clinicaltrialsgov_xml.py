#!/usr/bin/python

import sys
import zipfile
import random
from xml.etree import cElementTree

srch_ncts = []
with open('data/2020-04-16-0956-search-result.tsv',mode='rb') as f:
	header = f.readline()
	colnames = header.strip().split('\t')
	for line in f.readlines():
		cells = line.strip().split('\t')
		nct = cells[colnames.index('NCT Number')]
		srch_ncts.append(nct)

# # select a random set of 25 trials from the search
# random.seed(1)
# pilot_random_ncts = random.sample(srch_ncts,25)

# # and include the poscons meredith already looked at
# # note NCT02008357 = A4
# pilot_poscon_ncts = ['NCT00007189','NCT00592995','NCT01998841','NCT02008357','NCT01760005','NCT02519036','NCT02215616','NCT02197130','NCT02006472','NCT01795859','NCT00712426','NCT03344601','NCT02535884','NCT03854019','NCT03842969','NCT03761849','NCT03787758','NCT03575676','NCT03515213','NCT03764215','NCT03342053','NCT03225833','NCT03225846','NCT02453061','NCT02509793','NCT02481674','NCT02336633','NCT00514774']

# random.seed(1)
# pilot_ncts = pilot_random_ncts + pilot_poscon_ncts
# pilot_ncts = random.sample(pilot_ncts, k=len(pilot_ncts))

# # double check none missing:
# missing_ncts = list(set(pilot_ncts) - set(srch_ncts))
# print missing_ncts
# # good, none missing

# put in random order
random.seed(1)
ncts = random.sample(srch_ncts, k=len(srch_ncts))

trl = open('output/trials.tsv',mode='wb')
trl.write('\t'.join(['nct', 'gs_nct', 'dataset_inclusion', 'exclude_comment', 'year', 'brief_title', 'official_title', 'disease_area', 'phase', 'healthy_volunteers', 'intervention_list', 'lead_sponsor_agency', 'lead_sponsor_agency_class', 'sponsor_type', 'inclusion_text', 'exclusion_text', 'earliest_disease_stage', 'latest_disease_stage', 'advanced_exclusion', 'advanced_exclusion_comments', 'score_test', 'score_least_advanced', 'score_most_advanced', 'enrollment', 'enrollment_type', 'overall_status', 'comments', 'start_date', 'primary_completion_date', 'completion_date', 'results_first_posted_date'])+'\n')
#ivn = open('output/interventions.tsv',mode='wb')
ivn_grp = open('output/interventions_grouped.tsv',mode='wb')
timatch = open('output/trial_intervention_match.tsv',mode='wb')
timatch.write('\t'.join(['nct', 'intervention_name','intervention_type'])+'\n')

all_interventions = []

# shortlist of ncts used for implementation/debugging
# ncts = ['NCT02590276','NCT03981536','NCT04140136','NCT01565395','NCT03160898','NCT03103815','NCT02714036','NCT03241784','NCT03427086','NCT04172792','NCT03068754','NCT03482050','NCT02868580','NCT02437110','NCT03168711','NCT03293394','NCT03843710','NCT03892382','NCT01492686','NCT01884571']

in_path = 'bigdata/AllPublicXML_2020-04-16.zip'

z = zipfile.ZipFile(in_path) 
for nct in ncts:
	subpath = nct[:7]+'xxxx/'+nct+'.xml'
	f = z.open(subpath)
	link = 'https://clinicaltrials.gov/ct2/show/' + nct
	gs_nct = '=HYPERLINK("'+link+'","'+nct+'")'
	iterparser = iter(cElementTree.iterparse(f, events=("start", "end")))	
	interventions = []
	brief_title = ''
	official_title = ''
	lead_sponsor_agency = ''
	lead_sponsor_agency_class = ''
	year = ''
	enrollment = ''
	enrollment_type = ''
	overall_status = ''
	inclusion_text = ''
	exclusion_text = ''
	healthy_volunteers = ''
	phase = ''
	intervention_list = ''
	start_date = ''
	primary_completion_date = ''
	completion_date = ''
	results_first_posted_date = ''
	for event, elem in iterparser:
		if elem.tag == 'brief_title' and event=='end':
			brief_title = elem.text
		elif elem.tag == 'official_title' and event=='end':
			official_title = elem.text
		elif elem.tag == 'lead_sponsor' and event=='end':
			lead_sponsor_node = elem
			lead_sponsor_agency_node = elem.find('agency')
			lead_sponsor_agency = lead_sponsor_agency_node.text
			lead_sponsor_agency_class_node = elem.find('agency_class')
			lead_sponsor_agency_class = lead_sponsor_agency_class_node.text
		elif elem.tag == 'start_date' and event=='end':
			start_date_node = elem
			year = start_date_node.text[-4:]
			start_date = start_date_node.text
		elif elem.tag == 'enrollment' and event=='end':
			enrollment_node = elem
			enrollment = elem.text
			if 'type' in elem.attrib:
				enrollment_type = elem.attrib['type']
			else:
				enrollment_type = 'Unspecified'
		elif elem.tag == 'overall_status' and event=='end':
			overall_status_node = elem
			overall_status = overall_status_node.text
		elif elem.tag == 'eligibility' and event=='end':
			eligibility_node = elem
			criteria_node = eligibility_node.find('criteria')
			textblock_node = criteria_node.find('textblock')
			eligibility_text = ' '.join(textblock_node.text.strip().split())
			# try to split out inclusion and exclusion where possible
			inclusion_start = eligibility_text.lower().find('inclusion')
			exclusion_start = eligibility_text.lower().find('exclusion')
			inclusion_text = eligibility_text[inclusion_start:exclusion_start]
			exclusion_text = eligibility_text[exclusion_start:]
			healthy_volunteers_node = eligibility_node.find('healthy_volunteers')
			if healthy_volunteers_node is not None:
				healthy_volunteers = healthy_volunteers_node.text
				if healthy_volunteers == 'Accepts Healthy Volunteers':
					healthy_volunteers = 'Yes'
			else:
				healthy_volunteers = '?'
		elif elem.tag == 'intervention' and event=='end':
			intervention_node = elem
			intervention_name_node = intervention_node.find('intervention_name')
			intervention_type_node = intervention_node.find('intervention_type')
			interventions.append([intervention_name_node.text, intervention_type_node.text])
			timatch.write('\t'.join([nct,intervention_name_node.text.lower(),intervention_type_node.text.lower()]).encode('utf-8')+'\n')
		elif elem.tag == 'phase' and event=='end':
			phase_node = elem
			phase = phase_node.text.replace('Phase ','')
		elif elem.tag == 'primary_completion_date' and event=='end':
			primary_completion_date_node = elem
			primary_completion_date = primary_completion_date_node.text
		elif elem.tag == 'completion_date' and event=='end':
			completion_date_node = elem
			completion_date = completion_date_node.text
		elif elem.tag == 'results_first_posted' and event=='end':
			results_first_posted_node = elem
			results_first_posted_date = results_first_posted_node.text
		elif elem.tag == 'clinical_study' and event=='end':
			intervention_list = '; '.join([x[0]+' ('+x[1]+')' for x in interventions])
			all_interventions += [[nct] + intervention for intervention in interventions]
			# fields for curation
			exclude_comment = ''
			sponsor_type = ''
			disease_area = ''
			earliest_disease_stage = ''
			latest_disease_stage = ''
			dataset_inclusion = ''
			advanced_exclusion = ''
			advanced_exclusion_comments = ''
			score_test = ''
			score_least_advanced = ''
			score_most_advanced = ''
			target_gene = ''
			target_comments = ''
			target_genetic_support = ''
			target_genetic_support_comments = ''
			comments = ''
			trl.write('\t'.join([nct, gs_nct, dataset_inclusion, exclude_comment, year, brief_title, official_title, disease_area, phase, healthy_volunteers, intervention_list, lead_sponsor_agency, lead_sponsor_agency_class, sponsor_type, inclusion_text, exclusion_text, earliest_disease_stage, latest_disease_stage, advanced_exclusion, advanced_exclusion_comments, score_test, score_least_advanced, score_most_advanced, enrollment, enrollment_type, overall_status, comments, start_date, primary_completion_date, completion_date ,results_first_posted_date]).encode('utf-8')+'\n')
			#for intervention in interventions:
			#	ivn.write('\t'.join([gs_nct,intervention[0],intervention[1],target_gene,target_comments,target_genetic_support,target_genetic_support_comments]).encode('utf-8')+'\n') # nct, name, type
	# if elem.tag == namespace + 'uniprot':
	# 	root = elem
	# 	continue
	# if elem.tag == namespace + 'entry' and event == 'end': # wait until that node is done parsing
	# 	gene = ''
	# 	comments = ''
	# 	entry = elem
	# 	entries_read += 1

# group interventions for curation
ivn_dict = {} # [name, type] key to [example_nct, count_nct]
for intervention in all_interventions:
	name_and_type = tuple(map(lambda x:x.lower(), intervention[1:]))
	if name_and_type in ivn_dict.keys():
		ivn_dict[name_and_type][0] += 1
	else:
		ivn_dict[name_and_type] = [1, intervention[0]]

ivn_rows = []
for ivn_key in ivn_dict.iterkeys():
	ivn_row = list(ivn_key) + ivn_dict[ivn_key]
	ivn_row[2] = str(ivn_row[2])
	ivn_row[3] = '=HYPERLINK("https://clinicaltrials.gov/ct2/show/'+ivn_row[3]+'","'+ivn_row[3]+'")'
	ivn_row_text = '\t'.join(ivn_row)
	ivn_rows.append(ivn_row_text)
ivn_rows.sort()

ivn_grp.write('\t'.join(['intervention_name', 'intervention_type', 'n_trials', 'example_trial', 'intervention_class', 'target_gene', 'target_gene_comments'])+'\n')

for ivn_row_text in ivn_rows:
	ivn_grp.write(ivn_row_text.encode('utf-8')+'\t\t\t\n')

