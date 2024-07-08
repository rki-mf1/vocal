process REPORT {
    label 'report'

    publishDir "${params.output}/${params.report_dir}", mode: params.publish_dir_mode

    input:
        path variants_with_phenotypes
        path ecdc
        path bloom
        path lineages
        path vocal_version
        path db_version
        val email
        val email_sum
        val metadata

    output:
        path "vocal-alerts-samples-all.csv",                emit: alerts_samples
        path "vocal-alerts-clusters-summaries-all.csv",     emit: alerts_clusters
        path "vocal-report.html",                           emit: report

    script:
    """
    echo "Step 3: Detect and alert emerging variants"

    echo "Preparing csv for report..."
    
    Script_VOCAL_unified.R \
        -f ${variants_with_phenotypes} \
        -s "vocal-alerts-samples-all.csv" \
        -c "vocal-alerts-clusters-summaries-all.csv" \
        -a "${metadata}" \
        --ecdc ${ecdc} \
        --bloom ${bloom} \
        --lineages ${lineages}
    
    echo "Building HTML report..."

    Reporter.py  \
        -s "vocal-alerts-samples-all.csv" \
        -c "vocal-alerts-clusters-summaries-all.csv" \
        -v ${vocal_version} \
        -d ${db_version} \
        -e ${ecdc} \
        -m ${email} \
        -n ${email_sum} \
        -o "vocal-report.html"
    """

    stub:
    """
    touch vocal-alerts-samples-all.csv 
    touch vocal-alerts-clusters-summaries-all.csv 
    touch vocal-report.html 
    """
}