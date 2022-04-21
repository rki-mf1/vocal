## To Run: 
## Rscript --vanilla vocal/Script_VOCAL_DESH.autopilot.V2.R -d /scratch/kongkitimanonk/Vocal/vocal-test-preproduction --date 2021-11-09 -l 2 -o /scratch/kongkitimanonk/Vocal/vocal-test-preproduction
##
library(optparse)
library(digest)
library(tidyverse, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(glue, warn.conflicts = FALSE)
library(igraph, warn.conflicts = FALSE)

options(error=traceback)
option_list <- list(
    make_option(c("-f", "--file_variants_table"), default="input/",
                help="The path to the file the variants annotation (output from Mutations2Function.py)"),
    make_option(c("-v", "--anno_vocal_path"), default=file.path(getwd(),"data"),
                help="directory path where both ECDC_assigned_variants.csv and escape_data_bloom_lab.csv are located"),

    make_option(c("-o", "--outdir"), default="~/Vocal/vocal-test-preproduction/",
                help="Output directory [default %default]")
    )

args <- parse_args(OptionParser(option_list = option_list))

theme_set(theme_minimal())

theme_update(axis.text = element_text(size = 14),
             axis.title = element_text(size = 14),
             strip.text = element_text(size = 14),
             legend.position = "bottom", 
             legend.title = element_blank())

##### Functions

write_csv_wcomment = function(data, fout, comments = c()){
  comment_lines = str_c("#", comment, sep = " ")
  write(comment_lines, file = fout, append = FALSE)
  write_csv(data, fout, append = TRUE)  
}


#####
##### Setting up some variables

indel_max_size = 5 #larger indel will not be considered
indel_bad_qc_size = 10

Hamming_dist_clustering = 1 #For clustering the alerts detected by Vocal

#file_variant_table = "variants_with_phenotype_sc2-global-2021-04-21.tsv"
file_variant_table = args$file_variants_table


out_path = str_c(args$outdir)
if (!dir.exists(out_path)){
  dir.create(out_path)
}

vocal_path = args$anno_vocal_path

file_ECDCvariants_csv = file.path(vocal_path, 
                                  "ECDC_assigned_variants.csv")
file_BLOOM_mutation_csv = file.path(vocal_path, 
                                    'escape_data_bloom_lab.csv')

antibody_escape_score_raw = read_csv(file_BLOOM_mutation_csv,
                                 col_types = "cccicciccdcddddcic")

antibody_escape_score_summary = antibody_escape_score_raw %>% 
  filter(eliciting_virus == "SARS-CoV-2") %>% 
  group_by(wildtype, site, mutation) %>% 
  summarize(n_measure = n(),
            across(ends_with("_escape"), 
                   .fns = list(avg = mean),
                   .names = "{.fn}.{.col}")) %>% 
  ungroup() %>% 
  mutate(aa_pattern = glue("{wildtype}{site}{mutation}"))

###Info on the variants of concern

ECDC_variants = read_csv(file_ECDCvariants_csv) #,  show_col_types = FALSE)

VariantsOfConcern = ECDC_variants %>% filter(Status == "VOC") %>% 
  pull(`Pango lineage`)

VariantsOfInterest = ECDC_variants %>% filter(Status == "VOI") %>% 
  pull(`Pango lineage`) ## %>% grep("B.1.617",., invert = TRUE, value = TRUE)##Test for pres


##Not used for the moment
VariantsMonitoring = ECDC_variants %>% filter(Status == "Monitoring") %>% 
  pull(`Pango lineage`)

####
#### Reading up the results
####

##Vocal file
var_pheno = read_tsv(file_variant_table,
                     col_type = "ccccicicc") %>% 
                    mutate( type = factor(type, 
                                  levels = c("LineageDefiningMutation", 
                                              "MutationOfConcern", "RegionOfInterest", 
                                              "NotAnnotated", "VariantType"))) %>%
                    complete(type = type)

##metadata files

###
### Computing the table of phenotypes in wide format
### Each mutation in each sample is now on one line
### and filtering the variant that do not make sense
###  to be merged after with the metadata
###

var_pheno_wide_filter = var_pheno %>% 
  filter(variant_size <= indel_max_size) %>% 
  mutate(Test = 1L) %>% 
  pivot_wider(id_cols = ID:variant_size, 
              names_from = "type", 
              values_from = "Test",
              values_fn = max, 
              values_fill = 0) 

### 
### Our different types of variants that can be scored
### 
scores_combination = expand_grid(LineageDefiningMutation = c(TRUE,FALSE), 
                                 MutationOfConcern = c(TRUE,FALSE), 
                                 RegionOfInterest = c(TRUE,FALSE),
                                 NotAnnotated = c(TRUE,FALSE), 
                                 VariantType = c("VariantOfConcern", 
                                                 "VariantOfInterest", 
                                                 "Other")) %>% 
          mutate(s_pm = case_when(
                  VariantType != "Other" & NotAnnotated & !RegionOfInterest ~ 1,
                  VariantType == "Other" & NotAnnotated & !RegionOfInterest ~ 1,
                  VariantType == "Other" & !NotAnnotated & !MutationOfConcern & !RegionOfInterest ~ 1,
                  TRUE ~ 0
                ), 
                s_moc = case_when(
                  VariantType != "Other" & !NotAnnotated & MutationOfConcern & !LineageDefiningMutation ~ 1,
                  VariantType == "Other" & !NotAnnotated & MutationOfConcern ~ 1,
                  TRUE ~ 0
                ) ,
                 s_roi = case_when(
                   NotAnnotated & RegionOfInterest ~ 1,
                   VariantType == "Other" & !NotAnnotated & !MutationOfConcern & RegionOfInterest ~ 1,
                   TRUE ~ 0
                 )) %>% mutate(across(LineageDefiningMutation:NotAnnotated, as.integer))

var_pheno_wide_filter_with_score_bloom = inner_join(var_pheno_wide_filter, 
                                                 scores_combination) %>% 
  left_join(antibody_escape_score_summary, by = c("aa_pattern"))



##TEmp for 
var_pheno_summary_simple = var_pheno_wide_filter_with_score_bloom %>% 
  group_by(ID, variant_type) %>% 
  summarise(ListMutationsSelected = str_c(aa_pattern[s_roi > 0 | s_moc > 0 | s_pm >0], collapse = ","),
            nMutationsTotal = n(),
            nLineageDefining = sum(LineageDefiningMutation),
            across(starts_with("s_"), sum),
            across(ends_with("_escape"), ~ sum(., na.rm =TRUE),
                   .names = "sum_of_{.col}")) %>% 
  ungroup()


alerts_colors = c("red" = "#FF2400", "pink" = "deeppink", "orange" = "orange", 
           "lila" = "orchid", "grey" = "slategrey")
alert_codes = factor(names(alerts_colors),
                     ordered = TRUE
                )

# ###The list of all candidates with some alert level
# isolates_select_4_alert= var_pheno_score_summary %>% 
#   filter(s_moc>=1 | s_roi >= 3 | s_pm >= 5) %>% 
#   select(ID) %>% distinct()


var_pheno_summary_wide = var_pheno_summary_simple %>% 
#  inner_join(isolates_select_4_alert) %>% 
  pivot_wider(names_from = "variant_type", 
              values_from = c("ListMutationsSelected", "nMutationsTotal",
                              "nLineageDefining", 
                              "s_moc","s_roi","s_pm",
                              ends_with("_escape")))

var_pheno_summary_wide_with_alert = var_pheno_summary_wide %>% 
  mutate(VariantType = "Other", ##FIX change to your Variant Type
         LINEAGE = "None", ##FIX change to your LINEAGE Type
          s_moc_roi_tot = s_moc_M + s_moc_D + s_roi_M + s_roi_D,
          alert_level = 
           case_when(
             VariantType == "Other" & 
               (s_moc_M >= 3 | s_moc_D >= 2) & 
                    s_pm_M >= 0 ~ "red",
             VariantType == "Other" & 
               (s_moc_M >= 2 | s_moc_D >= 1) & 
               s_pm_M >= 4 ~ "red",
             VariantType == "Other" & 
               (s_moc_M >= 2 | s_moc_D >= 1) & 
               s_pm_M >= 2 ~ "orange",
             VariantType %in% c("VariantOfConcern", 
                                "VariantOfInterest") & 
               (s_moc_M >= 1 | s_moc_D >= 1) ~ "pink",
             VariantType == "Other" & 
               (s_moc_M >= 1 | s_moc_D >= 1) &
               (s_roi_M >= 1 | s_roi_D >= 1) ~ "lila",
             TRUE ~ "grey",
           ))

var_pheno_summary_wide_with_alert %>% group_by(alert_level) %>% count()
write_tsv(var_pheno_summary_wide_with_alert, 
                    file = file.path(out_path, "var_pheno_summary_wide_with_alert.tsv"))
error_log_outputFile <-file.path(out_path, "R_error.output.txt")

### Some simple functions for clustering
###
tryCatch( {
hamming <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}

hamming_tibble = function(tab){
  mat = as.matrix(tab[,-1])
  hamming(mat)
}

##Get the cluster from a table tab where first column is the ID and
##the other ones are the mutation patterns
get_sequence_clusters = function(tab, max_dist = Hamming_dist_clustering){
  matdists = hamming_tibble(tab)
  g = graph.adjacency(matdists <= max_dist, mode = "undirected")
  clus = components(g)
  tibble(ID = tab$ID, 
         cluster_ID = clus$membership,
         cluster_size = clus$csize[clus$membership],)
}
          },
          error = function(e){
           cat(paste("error at clustering ", as.character(e)), file=error_log_outputFile,append=TRUE,sep="\n")
          }
        )

###
### Let's prepare the table with the mutations 
### for all the alert levels.
### We create a nested dataframe for each alert level
###

tryCatch( {
mutations_per_alert_level = var_pheno_summary_wide_with_alert  %>% 
              ungroup() %>% filter(alert_level != "grey") %>%
              select(ID, alert_level) %>% 
              distinct() %>% ##WARNING THIS IS A HOTFIX
              inner_join(var_pheno_wide_filter %>% 
                           ungroup() %>% select(ID, aa_pattern)) %>% 
              mutate(value = 1L) %>%
              group_by(alert_level) %>% 
              nest() %>% 
              mutate(data = map(data, 
                                  ~ pivot_wider(., 
                                                names_from = aa_pattern, 
                                                values_from = value, 
                                                values_fill = 0)) )

alert_level_groups_with_clusters = mutations_per_alert_level %>% 
  filter(alert_level != "grey") %>% 
  mutate(clusters = map(data, get_sequence_clusters))

alerts_with_clusters_ID = alert_level_groups_with_clusters%>% select(-data) %>% 
  unnest( c(alert_level, clusters)) %>% 
  rename(cluster_ID_in_alert_level = cluster_ID)

vocal_list_samples_with_alert = var_pheno_summary_wide_with_alert %>% 
  inner_join(alerts_with_clusters_ID) %>% 
  arrange(desc(alert_level), desc(s_moc_roi_tot), desc(cluster_size))

vocal_common_mutations_in_clusters = vocal_list_samples_with_alert %>% 
    inner_join(var_pheno_wide_filter) %>% 
  group_by(alert_level, cluster_ID_in_alert_level) %>% 
  summarize(FreqMutThr = round( length( ID %>% unique() ) * 0.3, d = 0 ) + 1,
            aa_pattern_common = fct_lump_min(aa_pattern, FreqMutThr) #,
            #aa_pattern_uncommon = fct_lump_min(aa_pattern, -FreqMutThr)
            ) %>% 
  summarize(
            ListFrequentMutations_gt30perc = str_c(aa_pattern_common %>% 
                                             unique() %>% 
                                             grep("Other", ., invert = TRUE, 
                                                  fixed = TRUE, value=TRUE), 
                                          collapse = ",") ,
            # ListRareMutations_le30perc = str_c(aa_pattern_uncommon %>% 
            #                              unique() %>% 
            #                              grep("Other", ., invert = TRUE, 
            #                                    fixed = TRUE, value=TRUE), 
            #                               collapse = ","),
            ) %>% 
  ungroup() %>% arrange(desc(alert_level), cluster_ID_in_alert_level)  

vocal_list_clusters_properties = vocal_list_samples_with_alert %>% 
  group_by(alert_level, cluster_ID_in_alert_level, cluster_size) %>% 
  summarize(
            across(starts_with("s_"), ~ round(mean(.), d= 1), .names = "{.col}.avg"),
            across(contains("_escape_"), ~ round(mean(.), d= 1), .names = "{.col}"),
            Lineages = str_c(LINEAGE %>% unique(), collapse=  ","))

vocal_list_clusters_properties_with_mutations = vocal_common_mutations_in_clusters %>% 
  inner_join(vocal_list_clusters_properties) %>% 
  arrange(desc(alert_level), desc(s_moc_roi_tot.avg)) %>% 
  select(!c(ends_with("escape_D"), ends_with("escape_T"), ends_with("escape_I"), 
            ends_with("_T.avg"))) %>% 
  relocate(c(cluster_size, Lineages), .after = alert_level) %>% 
  relocate(cluster_ID_in_alert_level, .after = last_col()) %>% 
  rename(n_samples = cluster_size) %>% 
  mutate(across(where(is.numeric), ~ replace_na(.x , 0)))
  
write_tsv(vocal_list_clusters_properties_with_mutations, 
          file = file.path(out_path, 
                           "vocal_alerts_clusters_summaries_all.tsv"))

          },
          error = function(e){
           cat(paste("error at Let's prepare the table with the mutations ", as.character(e)), file=error_log_outputFile,append=TRUE,sep="\n")
          }
        )


tryCatch( {
          vocal_samples_out = vocal_list_samples_with_alert %>%
            mutate(across(where(is.numeric), ~ replace_na(.x , 0)),
                  across(where(is.character), ~ replace_na(.x , "")),
            ) %>%
            select(!c(ends_with("escape_D"), ends_with("escape_T"), ends_with("escape_I"), 
                      ends_with("_T.avg"), ends_with("MutationsSelected_T"))) %>% 
            relocate(ID, alert_level, SAMPLING_DATE, s_moc_roi_tot, s_moc_M, s_moc_D,
                    s_moc_I, s_roi_M,	s_roi_D,	s_roi_I,	
                    s_pm_M, s_pm_D, s_pm_I, starts_with("ListMutationsSelected"),
                    LINEAGE, VariantType, cluster_size, cluster_ID_in_alert_level) %>% 
            arrange(desc(alert_level), desc(s_moc_roi_tot), 
                    desc(s_moc_D), desc(s_moc_M), 
                    desc(s_roi_D), desc(s_roi_M))
            
            write_tsv(vocal_samples_out, 
                      file = file.path(out_path, "vocal_alerts_samples_all.tsv"))
          },
          error = function(e){
           cat(paste("error at vocal_samples_out", as.character(e)), file=error_log_outputFile,append=TRUE,sep="\n")
          }
        )
tryCatch( 
          {
          all_samples_with_alert_out = var_pheno_summary_wide_with_alert %>% 
            mutate(across(where(is.numeric), ~ replace_na(.x , 0)),
                  across(where(is.character), ~ replace_na(.x , "")),
            ) %>%
            select(!c(ends_with("escape_D"), ends_with("escape_T"), ends_with("escape_I"), 
                      ends_with("_T.avg"), ends_with("MutationsSelected_T"))) %>% 
            relocate(ID, alert_level, SAMPLING_DATE, 
                    s_moc_roi_tot,
                    s_moc_M, s_moc_D,
                    s_moc_I, s_roi_M,	s_roi_D,	s_roi_I,	
                    s_pm_M, s_pm_D, s_pm_I, starts_with("ListMutationsSelected"),
                    LINEAGE, LINEAGE.LATEST, VariantType) %>% 
            arrange(desc(alert_level), desc(s_moc_roi_tot),
                    desc(s_moc_D), desc(s_moc_M), 
                    desc(s_roi_D), desc(s_roi_M), 
                    desc(SAMPLING_DATE))
            
          write_tsv(all_samples_with_alert_out, 
                    file = file.path(out_path, "vocal_all_samples_with_alert_code.tsv"))
          },
          error = function(e){
            cat(paste("error at all_samples_with_alert_out",as.character(e)), file=error_log_outputFile, append=TRUE,sep="\n")
          }
        )

