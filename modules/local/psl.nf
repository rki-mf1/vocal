process PSL {
    label 'psl'

    publishDir "${params.output}/${params.vocal_dir}", mode: params.publish_dir_mode

    input:
        path ref_nt
        path input_fasta

    output:
        path "output.psl",          emit: output_psl
        path "variant_table.tsv",   emit: variant_table

    script:
    """
    echo "Step 1: Generate a PSL file with alingments and annotate mutations in the proteins"

    echo "Aligning sequences..."

    pblat ${ref_nt} \
        ${input_fasta} \
        -threads=4 \
        "output.psl"

    echo "Annotation by VOCAL..."

    vocal.py \
        -i ${input_fasta} \
        --PSL "output.psl" \
        -o "variant_table.tsv"
    """

    stub:
    """
    touch output.psl
    touch variant_table.tsv
    """
}