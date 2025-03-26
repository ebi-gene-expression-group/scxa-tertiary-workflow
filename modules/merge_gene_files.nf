/*
 * mergeGeneFiles: Merges gene file with genemeta on column 1, and keeps column1 and 4
 */
process MERGEGENEFILES {
    input:
      path gene
      path filtered_genemeta

    output:
      path 'merged_genemeta.tsv'

    script:
    def args    = task.ext.args ?: ""
    """
        # Sort both files by the first column for join compatibility
        sort -k1,1 "$gene" > sorted_gene.txt
        sort -k1,1 "$filtered_genemeta" > sorted_genemeta.txt
        
        # Perform a left join to keep all data from gene file
        join -a 1 -t \$'\t' -o 0,1.2,2.2 sorted_gene.txt sorted_genemeta.txt | cut -f1,3 > merged_genemeta.tsv
    """
    stub:
    """
      touch merged_genemeta.tsv
    """
}