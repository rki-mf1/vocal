process PREPROCESS {
    label 'preprocess'

    publishDir "${params.output}/${params.preprocess_dir}", mode: params.publish_dir_mode

    input:
        path metadata

    output:
        path "metadata.tsv",   emit: metadata

    script:
    """
    echo "Step 0: Preprocessing metadata file"

    file=${metadata}
    required_columns=("ID" "LINEAGE")

    header=\$(head -n 1 "\$file")

    IFS=\$'\t' read -r -a header_columns <<< "\$header"

    missing_columns=()
    for col in "\${required_columns[@]}"; do
        if [[ ! " \${header_columns[@]} " =~ " \${col} " ]]; then
            missing_columns+=("\$col")
        fi
    done

    expected_header="accession\tdescription\tlab\tsource\tcollection\ttechnology\tplatform\tchemistry\tmaterial\tct\tsoftware\tsoftware_version\tgisaid\tena\tzip\tdate\tsubmission_date\tlineage\tseqhash\tdna_profile\taa_profile\tfs_profile" 

    if [ \${#missing_columns[@]} -gt 0 ]; then
        if [[ "\$header" == "\$expected_header" ]]; then
            awk -v FS='\t' -v OFS='\t' \
                'NR==1 {print "ID", "PRIMARY_DIAGNOSTIC_LAB_PLZ", "SAMPLING_DATE", "LINEAGE"} NR>1 {print \$1, \$15, \$16, \$18}' \
                ${metadata} > metadata.tsv
            echo "The metadata tsv from covSonar was adjusted to VirusWarn-SC2 format."
        else
            echo "ERROR: The columns \${missing_columns[*]} are missing in the metadata!"
            echo "Please check if the required columns are existing. You may need to rename them."
            exit 1
        fi
    else
        mv \$file metadata.tsv
        echo "All required columns of the metadata are present."
    fi
    """

    stub:
    """
    touch metadata.tsv
    """
}