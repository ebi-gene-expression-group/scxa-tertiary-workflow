if [ -z "$SCXA_WORKFLOW_ROOT" ]; then
    echo "Variable SCXA_WORKFLOW_ROOT is not defined or empty. Please load SC env."
    echo "Exiting..."
    exit 1;
fi

EXP_ID=$1

echo "Creating ${pwd}/${EXP_ID} directory"
mkdir ${EXP_ID}
cd ${EXP_ID}

echo "Copying data to ${pwd}/${EXP_ID}"
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/${EXP_ID}.cell_metadata.tsv cell_metadata.tsv
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/genes.tsv.gz . && gunzip genes.tsv.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/matrix.mtx.gz . && gunzip matrix.mtx.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/barcodes.tsv.gz . && gunzip barcodes.tsv.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/reference/gene_annotation.txt genes_metadata.tsv
echo "Copying data for ${EXP_ID} finished"
