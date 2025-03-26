process COLUMN_REARRANGE_2 {
    
    input:
      path genemeta
      val col1
      val col2

    output:
      path 'filtered_genemeta_2.txt'

    script:
    def args    = task.ext.args ?: ""
    """
      # Find the column number of the specified gene_id column name
      col_num_1=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col1\$" | cut -d: -f1)
      col_num_2=\$(head -n1 "$genemeta" | tr '\\t' '\\n' | grep -n "^$col2\$" | cut -d: -f1)
  
      # If either column is not found, raise an error
      if [[ -z "\$col_num_1" || -z "\$col_num_2" ]]; then
          echo "Error: Column '$col1' or '$col2' not found in $genemeta" >&2
          exit 1
      fi
  
      # Extract the gene_id column (without the header)
      tail -n +2 "$genemeta" | cut -f\$col_num_1,\$col_num_2 > filtered_genemeta_2.txt
    """
    stub:
    """
      touch filtered_genemeta_2.txt
    """
}