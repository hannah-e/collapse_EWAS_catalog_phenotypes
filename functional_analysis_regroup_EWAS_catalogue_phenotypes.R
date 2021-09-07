############################################################################################################################################################
#
# 
# Use data from EWAS catalog to regroup phenotypes into a smaller number of categories
# 
#
###########################################################################################################################################################

#load data and libraries
library(readr)

EWAS_catalog <- as.data.frame(read_delim("http://ewascatalog.org/static//docs/ewascatalog-results.txt.gz", delim="\t"))


#### reclassify columns in the object "EWAS_catalog" #######################################################################################################


EWAS_catalog$phenotype<-EWAS_catalog$StudyID

temp<-grepl("age|aging|fetal_vs_adult_liver", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"age"

temp<-grepl("tissue", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"tissue"

temp<-grepl("smok|cotinine", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"smoking"

temp<-grepl("alcohol", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"alcohol"

temp<-grepl("sex", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"sex"

temp<-grepl("ancestry|ethnicity", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"ancestry"

temp<-grepl("cancer|carcinoma|adenoma|melanoma", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"cancer"

temp<-grepl("rheumatoid|ulcerative_colitis|lupus|sjogrens|crohn|inflammatory_bowel_disease|atopy|graves|psoriasis|multiple_sclerosis|cow_milk_allergy|26292806", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"autoimmune"
#26292806 study captures atopy

temp<-grepl("hiv|human_immunodeficiency_virus", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"infection"

temp<-grepl("blood_pressure|triglycerides|hdl|highdensity_lipoprotein|insulin|glucose|diabetes|chronic_kidney_disease|lipemia|creactive|ige|obesity|hepatic_fat|cholesterol|lipoprotein|ldl|vldl|statin_use|ischaemic_stroke|hypertension|c-reactive_protein|homair|resistin|hba1c|atrial_fibrillation|adiponectin|myocardial_infarction|leptin|liver_fat|coronary_heart_disease|idl|proinsulin|phospholipids|lp_a", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"cardiometabolic"

temp<-grepl("perinatal|birth_weight|birthweight|maternal_underweight|plasma_folate|prenatal|pregnancy|preterm_birth|season_of_birth|breastfeeding|31230546|33396735|utero|fetal_intolerance_of_labor|Starling-PS_maternal_serum|gestational_weight_gain|parity|fetal_brain_development", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"perinatal"
#31230546 study captures hypertension in pregnancy, preeclampsia
#33396735 study captures bottle, breast and mixed feeding behaviours

temp<-grepl("copd|fev1|fvc|chronic_obstructive_pulmonary_disease|asthma|lung_function|28621160|26681806", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"cardiovascular"
#28621160 study captures fev1 and fvc
#26681806 study captures cardiovascular biomarker GDF-15

temp<-grepl("bmi|body_mass_index|waist_circumference|arm_circumference|head_circumference|hip_circumference|fat|mass|weight|height|skinfold|waist|bone_mineral_density|Battram-T_head|Battram-T_hip|Battram-T_leg|Battram-T_pelvis|Battram-T_ribs|Battram-T_spine|Battram-T_total_body|Battram-T_trunk", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"anthropometric"

temp<-grepl("dementia|schizophrenia|palsy|alzheimer|depressive_disorder|depressive_symptoms|attention_deficit_hyperactivity_disorder|wellbeing|amyloid_plaques|depression|cognitive|personality_disorder|tic_disorders|aggressive_behaviour|infant_attention|cortical|stress|anxiety|neurobehavioural_scale|seizures|conduct_problems|parkinson|social_communication_deficits|hippocampus_volume|thalamus_volume|antidepressant_use|response_to_antidepressants|apolipoprotein|apoe", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"neurological"

temp<-grepl("socioeconomic_position|maternal_education|educational_attainment|31062658", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"SEP"
#31062658 captures various markers of maternal SEP

temp<-grepl("child_abuse|exposure_to_community_and_domestic_violence|exposure_to_intimate_partner_violence|bullying|perceived_discrimination|victimization", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"social_adversity"

temp<-grepl("healthy_eating| mediterranean_diet| diet_quality|pufa_intake|tea_consumption|fruit_consumption|juice_consumption|30101351|folate_intake|serum_copper|selenium|cadmium|arsenic|ppdde|substance_use|illicit_drug_use|esterogen_exposure|altitude| noise_pollution|noise_polution| sunlight_duration", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"diet_environment"
#30101351 study captures one-carbon metabolism nutritents

temp<-grepl("atmospheric_iron|atmospheric_nickel|atmospheric_vanadium|nitrogen_dioxide_exposure|proximity_to_major_roadways| particulate_matter| air_pollution", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"air_pollutants"

temp<-grepl("24014485|33413638|trimethylaminenoxide", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"metabolites"
#24014485 study captures study of blood serum metabolites
#33413638 NMR of predominantly lipids

temp<-grepl("31282290", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"miRNA"
#31282290 study captures EWAS of miRNA expression

#classify everything else as "other":
temp<-grepl("age|tissue|smoking|alcohol|sex|ancestry|cancer|autoimmune|infection|cardiometabolic|perinatal|cardiovascular|anthropometric|neurological|SEP|social_adversity|diet_environment|air_pollutants|metabolites|miRNA", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==F]<-"other"


#### check classifications ################################################################################################################################

table(EWAS_catalog$phenotype)

#### Create table of unique CpGs per category to use in enrichment analysis ###############################################################################

#cut to required p-value threshold
EWAS_catalog<-EWAS_catalog[EWAS_catalog$P<(0.05/484781),] #EDIT TO REDUCE TO REQUIRED P-THRESHOLD

#create output table
EWAS_catalog_collapsed_freq_table<-table(EWAS_catalog$phenotype)

#get number of unique CpGs in each category
for (i in names(EWAS_catalog_collapsed_freq_table)){

temp<-EWAS_catalog[EWAS_catalog$phenotype==i,]
temp<-length(temp$CpG[duplicated(temp$CpG)==F])
EWAS_catalog_collapsed_freq_table[names(EWAS_catalog_collapsed_freq_table)==i]<-temp
}

rm(temp)


#### end ################################################################################################################################



