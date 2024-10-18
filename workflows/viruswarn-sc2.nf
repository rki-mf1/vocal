/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// help message
if (params.help) { exit 0, helpMSG() }

// Parameters sanity checking
Set valid_params = ['cores', 'max_cores', 'memory', 'help',
                    'fasta', 'metadata', 'year', 'psl', 'strict', 
                    'output', 'preprocess_dir', 'vocal_dir', 
                    'annot_dir', 'report_dir', 'runinfo_dir',
                    'publish_dir_mode', 'conda_cache_dir',
                    'cloudProcess', 'cloud-process']

def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VOCAL_SUB } from '../subworkflows/local/sub_vocal'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRUSWARN_SC2 {

ref_nt = Channel.fromPath( file("test/ref.fna", checkIfExists: true) )

input_fasta = Channel.fromPath( file("${params.fasta}", checkIfExists: true) )

if (params.year == 2022) {
    log.info"INFO: VirusWarn-SC2 uses mutation, lineage and VOC/VOI/VUM information from November 2022"
    mutation_table =  Channel.fromPath( file("data/2022-11/table_cov2_mutations_annotation.tsv", checkIfExists: true) )
    ecdc =  Channel.fromPath( file("data/2022-11/ECDC_assigned_variants.csv", checkIfExists: true) )
    lineages =  Channel.fromPath( file("data/2022-11/lineage.all.tsv", checkIfExists: true) )
} else if (params.year == 2021) {
    log.info"INFO: VirusWarn-SC2 uses mutation, lineage and VOC/VOI/VUM information from September 2021"
    mutation_table =  Channel.fromPath( file("data/2021-09/table_cov2_mutations_annotation.tsv", checkIfExists: true) )
    ecdc =  Channel.fromPath( file("data/2021-09/ECDC_assigned_variants.csv", checkIfExists: true) )
    lineages =  Channel.fromPath( file("data/2021-09/lineage.all.tsv", checkIfExists: true) )
} else {
    exit 1,
    "ERROR: $params.year is an invalid input for the parameter year!"
}

if (params.metadata != '') {
    metadata = Channel.fromPath( file("${params.metadata}", checkIfExists: true) )
} else {
    log.warn"WARNING! No metadata file was given. This can lead to problems, like not correctly identifying pink alerts!"
    metadata = params.metadata
}

bloom =  Channel.fromPath( file("data/escape_data_bloom_lab.csv", checkIfExists: true) )

vocal_version = Channel.fromPath( file(".version", checkIfExists: true) )
db_version = Channel.fromPath( file("data/.db_version", checkIfExists: true) )

email = Channel.fromPath( file("templates/email.html") )
email_sum = Channel.fromPath( file("templates/email.sum.html") )

VOCAL_SUB ( 
    ref_nt, input_fasta, mutation_table, metadata,
    ecdc, bloom, lineages, 
    vocal_version, db_version, email, email_sum 
)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\u001B[31m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

    Workflow: VirusWarn-SC2

    ${c_yellow}Usage examples:${c_reset}
    nextflow run main.nf -profile conda,local --fasta 'test/sample-test.fasta'

    ${c_yellow}Input options:${c_reset}
    ${c_green} --fasta ${c_reset}           REQUIRED! Path to the input fasta file.
                        [ default: $params.fasta ]
    ${c_green} --metadata ${c_reset}        Path to the metadata file.
                        [ default: $params.metadata ]
    ${c_green} --year ${c_reset}            Specify the year from which the information should 
                        be used for the ranking.
                        [ default: $params.year ]
    ${c_green} --psl ${c_reset}             Run process with pblat alignment.
                        [ default: $params.psl ]
    ${c_green} --strict ${c_reset}          Run process with strict alert levels (without orange).
                        [ default: $params.strict ]

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]
    ${c_yellow}Output options:${c_reset}
    --output                 Name of the result folder [default: $params.output]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: $params.publish_dir_mode]
                        ${c_dim}With 'symlink' results are lost when removing the work directory.${c_reset}
    ${c_yellow}Caching:${c_reset}
    --conda_cache_dir        Location for storing the conda environments [default: $params.conda_cache_dir]
    
    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executors${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executor${c_reset} (choose one):
        local
        slurm
    
    ${c_blue}Engines${c_reset} (choose one):
        conda
        mamba
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/