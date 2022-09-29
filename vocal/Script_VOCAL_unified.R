##
## VOCAL unified script V.1
## This script takes a table of variant annotation
## and metadata table and generates a Vocal report
##
## The focus of this script is on genericity so all informations that would
## be specific to DESH database have been taken out
## To Run, from the vocal root directory:
## Rscript --vanilla vocal/Script_VOCAL_unified.R -f variants_with_phenotype_sc2-global.tsv -a report.meta.tsv
##
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr)) 
suppressPackageStartupMessages(library(readr))  
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(logger, warn.conflicts = FALSE))

# Suppress summarise() info
options(dplyr.summarise.inform = FALSE)

# options(error = traceback)

option_list <- list(
  make_option(
    c("-f", "--file_variant_table"),
    help = "file corresponding to the variants"
  ),
  make_option(
    c("-a", "--file_annotations"),
    default = "",
    help = "file containing metadata on the samples (optional)"
  ),
  make_option(
    c("-v", "--anno_vocal_path"),
    default = file.path(getwd(), "data"),
    help = "directory path where Vocal database is stored (files concerned: ECDC_assigned_variants.csv and escape_data_bloom_lab.csv and filiation_information) [default %default]"
  ),
  make_option(c("-o", "--outdir"), default = "results/",
              help = "Output directory [default %default]"),
  make_option(c("--id_column"), default = "ID",
              help = "Column name for the sample ID (this argument will be used if file_annotations is used) [default %default]"),
  make_option(c("--lineage_column"), default = "LINEAGE",
              help = "Column name reporting LINEAGE information in the metadata (this argument will be used if file_annotations is used) [default %default]"),
  make_option(c("--date_column"), default = "SAMPLING_DATE",
              help = "Column name reporting sampling date information (this argument will be used if file_annotations is used) [default %default]"),
  make_option(c("--geoloc_column"), default = "PRIMARY_DIAGNOSTIC_LAB_PLZ",
              help = "Column name for geolocation (this argument will be used if file_annotations is used) [default %default]"),
  make_option(c("--verbose"), default = FALSE, action = "store_true",
              help = "Write more output to the console")
)

parser = OptionParser(option_list = option_list)
exit <- function() { invokeRestart("abort") }  
args <- parse_args(parser)
# test if there is at least one argument: if not, return an error
args_tester <- function(args){
  if (length(args)==0) {
    print_help(parser)
    stop("At least one argument must be supplied.")
    exit()
  }else if (is.null(args$file_variant_table)){
    print_help(parser)
    stop("Require --file_variant_table")
    exit()
  }else{
    # print("Good to Go!")
  }
  return (args)
}

args = args_tester(args)

# Setup logging
if (args$verbose) {
  log_threshold(TRACE)
} else {
  log_threshold(WARN)
}

theme_set(theme_minimal())

theme_update(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  strip.text = element_text(size = 14),
  legend.position = "bottom",
  legend.title = element_blank()
)

## Define Functions

write_csv_wcomment = function(data, fout, comments = c()) {
  comment_lines = str_c("#", comment, sep = " ")
  write(comment_lines, file = fout, append = FALSE)
  write_csv(data, fout, append = TRUE)
}

## Setting up some variables and parameters
VARIANT_TYPE_LEVELS = c(
  "LineageDefiningMutation",
  "MutationOfConcern",
  "RegionOfInterest",
  "NotAnnotated",
  "VariantType"
)
VOC_KEY = "VariantOfConcern"
VOI_KEY = "VariantOfInterest"
OTHER_KEY = "Other"


ID_COL = args$id_column
DATE_COL = args$date_column
LINEAGE_COL = args$lineage_column
GEOLOC_COL = args$geoloc_column

MUT_TYPE_COL = "variant_type"
VARIANT_CLASS_COL = "VariantType"

indel_max_size = 5 # larger indel will not be considered
indel_bad_qc_size = 10

Hamming_dist_clustering = 1 # For clustering the alerts detected by Vocal

## This should all be deduced from the information in the metadata file
## Should go out as well
lookback_duration = NULL
## File paths and dates
earliestdate = "2020-01-01"

# file_variant_table = "variants_with_phenotype_sc2-global-2021-04-21.tsv"
file_variant_table = args$file_variant_table
file_annotations = args$file_annotations


out_path = str_c(args$outdir)
if (!dir.exists(out_path)) {
  dir.create(out_path)
}
######## VOCAL parameter files ########
## Preparing the information data

vocal_path = args$anno_vocal_path

file_ECDCvariants_csv = file.path(vocal_path,
                                  "ECDC_assigned_variants.csv")
file_BLOOM_mutation_csv = file.path(vocal_path,
                                    'escape_data_bloom_lab.csv')

file_variant_filiations_csv = file.path(vocal_path, "lineage.all.tsv")

antibody_escape_score_raw = read_csv(file_BLOOM_mutation_csv,
                                     col_types = "cccicciccdcddddcic")

## Prepare the BLOOM score
antibody_escape_score_summary = suppressMessages( antibody_escape_score_raw %>%
  filter(eliciting_virus == "SARS-CoV-2") %>%
  group_by(wildtype, site, mutation) %>%
  summarize(n_measure = n(),
            across(
              ends_with("_escape"),
              .fns = list(avg = mean),
              .names = "{.fn}.{.col}"
            )) %>%
  ungroup() %>%
  mutate(aa_pattern = glue("{wildtype}{site}{mutation}"))
)

## Gathering Information about the set of Variants of Concern and Variants of Interest
## Make the parsing of the files more robust (especially column names)
ECDC_variants = read_csv(file_ECDCvariants_csv, , col_types = cols()) #,  show_col_types = FALSE)
variants_filiation = read_tsv(file_variant_filiations_csv, , col_types = cols())

VariantsOfConcern = suppressMessages(ECDC_variants %>% filter(Status == "VOC") %>%
  inner_join(variants_filiation, by = c(`Pango lineage` = "lineage")) %>%
  separate_rows(sublineage, sep = ",") %>%
  pivot_longer(
    cols = c(`Pango lineage`, sublineage),
    names_to = "lineage_type",
    values_to = "absolute.lineage"
  ) %>%
  filter(absolute.lineage %>% str_starts("[A-Z]")) %>% 
  pull(absolute.lineage) %>% unique()
)

VariantsOfInterest = suppressMessages(ECDC_variants %>% filter(Status == "VOI") %>%
  inner_join(variants_filiation, by = c(`Pango lineage` = "lineage")) %>%
  separate_rows(sublineage, sep = ",") %>%
  pivot_longer(
    cols = c(`Pango lineage`, sublineage),
    names_to = "lineage_type",
    values_to = "absolute.lineage"
  ) %>%
  filter(absolute.lineage %>% str_starts("[A-Z]")) %>% 
  pull(absolute.lineage) %>% unique()
)

############ Reading up the results ################

## VOCAL file

var_pheno = read_tsv(file_variant_table,
                     col_type = "ccccicicc") %>%
  mutate(type = factor(type,
                       levels = VARIANT_TYPE_LEVELS)) %>%
  complete(type = type)

## Metadata files

if (file.exists(file_annotations)) {
  metadata = read_tsv(file_annotations, col_types = cols()) #, show_col_types = FALSE)
  if (!(LINEAGE_COL %in% colnames(metadata))) {
    stop( "No column with LINEAGE information in the metadata table, this column is required, exiting")
  }
  if (!(DATE_COL %in% colnames(metadata))) {
    stop(
      "No column with time information in the metadata table. this column is required, exiting"
    )
  }
  if (!(GEOLOC_COL %in% colnames(metadata))) {
    warning(
      "No column with geolocalisation information in the metadata table some infos will be missing"
    )
  }
  
  # metadata = rename(metadata, "lineage" = LINEAGE_COL)
  all_reports = suppressMessages(metadata %>%
    mutate(across(any_of(DATE_COL),  ~ date(.x))) %>%
    filter(across(any_of(DATE_COL), ~ . >= date(earliestdate))) %>%
    mutate(
      VariantType = case_when(
        LINEAGE %in% VariantsOfConcern ~ VOC_KEY,
        # lineage
        LINEAGE %in% VariantsOfInterest ~ VOI_KEY,
        # lineage  , have to handle this by change column at the begining
        TRUE ~ OTHER_KEY
      )
    )
  )
  
} else{
  warning("No meta information is given \n")
  
  # metadata <-  data.frame()
  metadata = var_pheno %>% distinct(ID)
  metadata[LINEAGE_COL] <- 'Other'
  metadata[DATE_COL] <- NA
  
  all_reports = metadata %>%
    mutate(
      VariantType = case_when(
        LINEAGE %in% VariantsOfConcern ~ VOC_KEY,
        # lineage
        LINEAGE %in% VariantsOfInterest ~ VOI_KEY,
        # lineage  , have to handle this by change column at the begining
        TRUE ~ OTHER_KEY
      )
    )
  
}

###
### Computing the table of phenotypes in wide format
### Each mutation in each sample is now on one line
### and filtering the variant that do not make sense
###  to be merged after with the metadata
### This results in a table with a binary mask of the
###

#all_reports = rename(all_reports, "ID" = ID_COL)
#ID_COL = "ID"
var_pheno_wide_filter = var_pheno %>%
  filter(variant_size <= indel_max_size) %>%
  mutate(Test = 1L) %>%
  pivot_wider(
    id_cols = ID:variant_size,
    names_from = "type",
    values_from = "Test",
    values_fn = max,
    values_fill = 0
  ) %>%
  inner_join(all_reports, by = 'ID')

######### Compute score #########
## Adding the mask for scoring,
## Our different types of variants that can be scored

score_mutation <- function(pheno_table) {
  ##
  if (!(all(VARIANT_TYPE_LEVELS %in% colnames(pheno_table)))) {
    simpleError("Table of phenotypes is missing some columns for scoring mutations")
  }
  
  pheno_table_infomask = pheno_table %>%
    mutate(
      infomask = str_c(
        LineageDefiningMutation,
        MutationOfConcern,
        RegionOfInterest,
        NotAnnotated,
        VariantType,
        sep = "."
      )
    )
  
  scores_combination = expand_grid(
    LineageDefiningMutation = c(TRUE, FALSE),
    MutationOfConcern = c(TRUE, FALSE),
    RegionOfInterest = c(TRUE, FALSE),
    NotAnnotated = c(TRUE, FALSE),
    VariantType = c(VOC_KEY,
                    VOI_KEY,
                    OTHER_KEY)
  ) %>%
    mutate(
      s_pm = case_when(
        VariantType != OTHER_KEY & NotAnnotated & !RegionOfInterest ~ 1,
        VariantType == OTHER_KEY &
          NotAnnotated & !RegionOfInterest ~ 1,
        VariantType == OTHER_KEY &
          !NotAnnotated &
          !MutationOfConcern & !RegionOfInterest ~ 1,
        TRUE ~ 0
      ),
      s_moc = case_when(
        VariantType != OTHER_KEY &
          !NotAnnotated &
          MutationOfConcern & !LineageDefiningMutation ~ 1,
        VariantType == OTHER_KEY &
          !NotAnnotated & MutationOfConcern ~ 1,
        TRUE ~ 0
      ) ,
      s_roi = case_when(
        NotAnnotated & RegionOfInterest ~ 1,
        VariantType == OTHER_KEY &
          !NotAnnotated & !MutationOfConcern & RegionOfInterest ~ 1,
        TRUE ~ 0
      )
    ) %>% mutate(across(LineageDefiningMutation:NotAnnotated, as.integer))
  pheno_table_with_score = suppressMessages(inner_join(pheno_table_infomask,
                                      scores_combination))
  return(pheno_table_with_score)
}

var_pheno_wide_filter_with_score = score_mutation(var_pheno_wide_filter)

var_pheno_score_summary = suppressMessages(var_pheno_wide_filter_with_score %>%
  left_join(antibody_escape_score_summary, by = c("aa_pattern")) %>%
  group_by(across(any_of(
    c(
      ID_COL,
      DATE_COL,
      GEOLOC_COL,
      MUT_TYPE_COL,
      VARIANT_CLASS_COL,
      LINEAGE_COL
    )
  )))  %>%
  summarise(
    ListMutationsSelected = str_c(aa_pattern[s_roi > 0 |
                                               s_moc > 0 |
                                               s_pm > 0],
                                  collapse = ","),
    nMutationsTotal = n(),
    nLineageDefining = sum(LineageDefiningMutation),
    across(starts_with("s_"), sum),
    across(ends_with("_escape"), ~ sum(., na.rm = TRUE),
           .names = "sum_of_{.col}")
  ) %>%
  ungroup()
)
### If a latest annotation for the
### lineage is provided, we add a column
if (FALSE) {
  # file.exists(file_annot_latest)
  print(file_annot_latest)
  desh_latest = read_tsv(file_annot_latest) %>%
    select(all_of(ID_COL, LINEAGE_COL)) %>%
    rename(LINEAGE.LATEST, LINEAGE_COL)
  
  var_pheno_score_summary = suppressMessages(inner_join(var_pheno_score_summary, desh_latest))
}
#write_tsv(var_pheno_score_summary, "./variants_desh_with_score.tsv")
#write_tsv(var_pheno_score_summary_with_latest_lineage, "./variants_desh_with_score_with_latest_lineage.tsv")


######### ALERT #########
alerts_colors = c(
  "red" = "#FF2400",
  "pink" = "deeppink",
  "orange" = "orange",
  "lila" = "orchid",
  "grey" = "slategrey"
)
alert_codes = factor(names(alerts_colors), ordered = TRUE)

var_pheno_summary_wide = var_pheno_score_summary %>%
  pivot_wider(
    names_from = all_of(MUT_TYPE_COL),
    values_from = c(
      "ListMutationsSelected",
      "nMutationsTotal",
      "nLineageDefining",
      "s_moc",
      "s_roi",
      "s_pm",
      ends_with("_escape")
    ),
    values_fill = list(
      ListMutationsSelected = "",
      nMutationsTotal = 0 ,
      nLineageDefining = 0,
      s_moc = 0 ,
      s_roi = 0,
      s_pm = 0
    )
  )

## This fonction returns the alert levels for VOCAL
## Version 1 with hard set thresholds
compute_alert_levels_v1 <- function(pheno_table_wide) {
  ## WARNING: all variable names are hard coded in this function

  if(! "s_roi_I" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_roi_I"] = 0}
  if(! "s_moc_I" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_moc_I"] = 0}
  if(! "s_moc_M" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_moc_M"] = 0}
  if(! "s_roi_M" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_roi_M"] = 0}
  if(! "s_moc_D" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_moc_D"] = 0}
  if(! "s_roi_D" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_roi_D"] = 0}
  if(! "s_pm_M" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_pm_M"] = 0}
    if(! "s_pm_D" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_pm_D"] = 0}
    if(! "s_pm_I" %in% colnames(pheno_table_wide))
  {pheno_table_wide["s_pm_I"] = 0}

    if(! "nMutationsTotal_M" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nMutationsTotal_M"] = 0}
    if(! "nMutationsTotal_D" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nMutationsTotal_D"] = 0}
    if(! "nMutationsTotal_I" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nMutationsTotal_I"] = 0}
    if(! "nLineageDefining_M" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nLineageDefining_M"] = 0}
    if(! "nLineageDefining_D" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nLineageDefining_D"] = 0}
    if(! "nLineageDefining_I" %in% colnames(pheno_table_wide))
  {pheno_table_wide["nLineageDefining_I"] = 0}

  pheno_table_wide_with_alert = pheno_table_wide %>%
    mutate(
      s_moc_roi_tot = (s_moc_M + s_moc_D +
                         s_roi_M + s_roi_D +
                         s_moc_I + s_roi_I),
      alert_level =
        case_when(
          VariantType %in% c(VOC_KEY, VOI_KEY) &
            (s_moc_M >= 1 | s_moc_D >= 1) ~ "pink",
          VariantType == OTHER_KEY &
            (s_moc_M >= 3 | s_moc_D >= 2) &
            s_pm_M >= 0 ~ "red",
          VariantType == OTHER_KEY &
            (s_moc_M >= 2 | s_moc_D >= 1) &
            s_pm_M >= 4 ~ "red",
          VariantType == OTHER_KEY &
            (s_moc_M >= 2 | s_moc_D >= 1) &
            s_pm_M >= 2 ~ "orange",
          VariantType == OTHER_KEY &
            (s_moc_M >= 1 | s_moc_D >= 1) &
            (s_roi_M >= 1 | s_roi_D >= 1) ~ "lila",
          TRUE ~ "grey",
        ),
        impact =
          case_when(
            alert_level ==  "pink" | alert_level ==  "red"
              ~ "HIGH",
            alert_level ==  "orange"
              ~ "MODERATE",
            alert_level ==  "lila" | alert_level ==  "grey"
              ~ "LOW",
          ),
    )
  return(pheno_table_wide_with_alert)
}

var_pheno_summary_wide_with_alert = compute_alert_levels_v1(var_pheno_summary_wide)

prediction_overview = suppressMessages(var_pheno_summary_wide_with_alert %>% group_by(alert_level) %>% count())
log_info("Prediction Results:")
if (log_threshold() >= TRACE) {
print(prediction_overview)
}

write.table(
  prediction_overview,
  file = file.path(out_path, "prediction-overview.txt"),
  sep = "\t"
)

### Some simple functions for clustering
hamming <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}

hamming_tibble = function(tab) {
  mat = as.matrix(tab[, -1])
  hamming(mat)
}

## Get the cluster from a table tab where first column is the ID and
## the other ones are the mutation patterns
get_sequence_clusters = function(tab, max_dist = Hamming_dist_clustering) {
  matdists = hamming_tibble(tab)
  g = graph.adjacency(matdists <= max_dist, mode = "undirected")
  clus = components(g)
  tibble(
    ID = tab$ID,
    cluster_ID = clus$membership,
    cluster_size = clus$csize[clus$membership],
  )
}

###
### Let's prepare the table with the mutations
### for all the alert levels.
### We create a nested dataframe for each alert level
###

mutations_per_alert_level = suppressMessages(var_pheno_summary_wide_with_alert  %>%
  ungroup() %>%
  filter(alert_level != "grey") %>%
  select(ID, alert_level) %>%
  distinct() %>% ## WARNING THIS IS A HOTFIX
  inner_join(var_pheno_wide_filter %>%
               ungroup() %>% select(ID, aa_pattern)) %>%
  mutate(value = 1L) %>%
  group_by(alert_level) %>%
  nest() %>%
  mutate(data = map(
    data,
    ~ pivot_wider(
      .,
      names_from = aa_pattern,
      values_from = value,
      values_fill = 0
    )
  ))
  )

alert_level_groups_with_clusters = mutations_per_alert_level %>%
  mutate(clusters = map(data, get_sequence_clusters))
print("Test")
if(nrow(alert_level_groups_with_clusters) != 0){
alerts_with_clusters_ID = alert_level_groups_with_clusters %>% select(-data) %>%
  unnest(c(alert_level, clusters)) %>%
  rename(cluster_ID_in_alert_level = cluster_ID)
}else{
# assign empty
alerts_with_clusters_ID <- data.frame(alert_level=character(), ID=character(), cluster_ID_in_alert_level=character(),
                         cluster_size=numeric())
}


vocal_list_samples_with_alert = suppressMessages(var_pheno_summary_wide_with_alert %>%
  left_join(alerts_with_clusters_ID) %>%
  arrange(desc(alert_level),
          desc(s_moc_roi_tot),
          desc(cluster_size),
          desc(DATE_COL))
  )
########### Output goes Here ###########
log_debug("Write Results at: {out_path}")

vocal_common_mutations_in_clusters = suppressMessages( vocal_list_samples_with_alert %>%
  filter(alert_level != "grey") %>%
  inner_join(var_pheno_wide_filter) %>%
  group_by(alert_level, cluster_ID_in_alert_level) %>%
  summarize(
    FreqMutThr = round(length(ID %>% unique()) * 0.3, d = 0) + 1,
    aa_pattern_common = fct_lump_min(aa_pattern, FreqMutThr) #,
    #aa_pattern_uncommon = fct_lump_min(aa_pattern, -FreqMutThr)
  ) %>%
  summarize(
    ListFrequentMutations_gt30perc = str_c(
      aa_pattern_common %>%
        unique() %>%
        grep(
          OTHER_KEY,
          .,
          invert = TRUE,
          fixed = TRUE,
          value = TRUE
        ),
      collapse = ","
    ) ,
    # ListRareMutations_le30perc = str_c(aa_pattern_uncommon %>%
    #                              unique() %>%
    #                              grep(OTHER_KEY, ., invert = TRUE,
    #                                    fixed = TRUE, value=TRUE),
    #                               collapse = ","),
  ) %>%
  ungroup() %>% arrange(desc(alert_level), cluster_ID_in_alert_level))

vocal_list_clusters_properties = suppressMessages(vocal_list_samples_with_alert %>%
    filter(!is.na(cluster_ID_in_alert_level)) %>%
    group_by(alert_level, cluster_ID_in_alert_level, cluster_size) %>%
    summarize(
      last_seen_isolate = max(.data[[DATE_COL]]),
      first_seen_isolate = min(.data[[DATE_COL]]),
      date_range = difftime(last_seen_isolate,
                            first_seen_isolate,
                            units = "days"),
      across(starts_with("s_"), ~ round(mean(.), d = 1), .names = "{.col}.avg"),
      across(contains("_escape_"), ~ round(mean(.), d = 1), .names = "{.col}"),
      Lineages = str_c(.data[[LINEAGE_COL]] %>% unique(), collapse =  ",")
    )
  )

vocal_list_clusters_properties_with_mutations = suppressMessages(vocal_common_mutations_in_clusters %>%
  inner_join(vocal_list_clusters_properties) %>%
  arrange(desc(alert_level),
          desc(s_moc_roi_tot.avg),
          desc(last_seen_isolate)) %>%
  select(!c(
    ends_with("escape_D"),
    ends_with("escape_T"),
    ends_with("escape_I"),
    ends_with("_T.avg")
  )) %>%
  relocate(c(cluster_size, last_seen_isolate, date_range, Lineages),
           .after = alert_level) %>%
  relocate(cluster_ID_in_alert_level, .after = last_col()) %>%
  rename(n_samples = cluster_size) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x , 0)))
)
write_csv(
  vocal_list_clusters_properties_with_mutations,
  file = file.path(out_path,
                   "vocal-alerts-clusters-summaries-all.csv")
)

error_log_outputFile <- file.path(out_path, "R-error-output.txt")

tryCatch({
  vocal_samples_out = vocal_list_samples_with_alert %>%
    mutate(across(where(is.numeric), ~ replace_na(.x , 0)),
           across(where(is.character), ~ replace_na(.x , "")),) %>%
    select(!c(
      ends_with("escape_D"),
      ends_with("escape_T"),
      ends_with("escape_I"),
      ends_with("_T.avg"),
      ends_with("MutationsSelected_T")
    )) %>%
    relocate(
      c(
        all_of(ID_COL),
        alert_level,
        all_of(DATE_COL),
        s_moc_roi_tot,
        s_moc_M,
        s_moc_D,
        s_moc_I,
        s_roi_M,
        s_roi_D,
        s_roi_I,
        s_pm_M,
        s_pm_D,
        s_pm_I,
        starts_with("ListMutationsSelected"),
        all_of(LINEAGE_COL),
        any_of("LINEAGE.LATEST"),
        all_of(VARIANT_CLASS_COL),
        cluster_size,
        cluster_ID_in_alert_level
      )
    ) %>%
    arrange(across(
      c(
        alert_level,
        s_moc_roi_tot,
        s_moc_D,
        s_moc_M,
        s_roi_D,
        s_roi_M,
        any_of(DATE_COL)
      ),
      desc
    ))
  
  write_csv(vocal_samples_out,
            file = file.path(out_path, "vocal-alerts-samples-all.csv"))
},
error = function(e) {
  cat(
    paste("error at vocal-samples-out", as.character(e)),
    file = error_log_outputFile,
    append = TRUE,
    sep = "\n"
  )
  stop("error at vocal-samples-out", as.character(e))
})
log_info("** Success **")
