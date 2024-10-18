include { PREPROCESS }      from '../../modules/local/preprocess'
include { VOCAL }           from '../../modules/local/vocal'
include { PSL }             from '../../modules/local/psl'
include { ANNOTATION }      from '../../modules/local/annotation'
include { REPORT }          from '../../modules/local/report'

workflow VOCAL_SUB {
    take:
        ref_nt
        input_fasta
        mutation_table
        metadata
        ecdc
        bloom
        lineages
        vocal_version
        db_version
        email
        email_sum

    main:
        if (metadata != ''){
            PREPROCESS ( metadata )
        }

        if (params.psl) {
            PSL ( ref_nt, input_fasta )

            ANNOTATION ( PSL.out.variant_table, mutation_table )
        } else {
            VOCAL ( input_fasta )

            ANNOTATION ( VOCAL.out.variant_table, mutation_table )
        }
        if (metadata != ''){
            REPORT ( 
                ANNOTATION.out.variants_with_phenotypes, ecdc, bloom, lineages, 
                vocal_version, db_version, email, email_sum, PREPROCESS.out.metadata
            )
        } else {
            REPORT ( 
                ANNOTATION.out.variants_with_phenotypes, ecdc, bloom, lineages, 
                vocal_version, db_version, email, email_sum, metadata
            )
        }

    emit:
        report = REPORT.out.report

}