process COLUMN_REARRANGE_1 {
    input:
      path genemeta
      val col

    output:
      path 'filtered_genemeta.txt'

    script:
    def args    = task.ext.args ?: ""
    """
      # Find the column number of the specified gene_id column name
      col_num=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col\$" | cut -d: -f1)
  
      # If column is found, extract it; otherwise, raise an error
      if [[ -z "\$col_num" ]]; then
          echo "Error: Column '$col' not found in $genemeta" >&2
          exit 1
      fi
  
      # Extract the gene_id column (without the header)
      tail -n +2 "$genemeta" | cut -f\$col_num > filtered_genemeta.txt
    """

    stub:
    """
      touch filtered_genemeta.txt
    """
}