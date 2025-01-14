#!/usr/bin/env bash

# This is EMBL-EBI specific script to fetch data from workflow root and put it in a place for downstream workflow to use 

if [ -z "$SCXA_WORKFLOW_ROOT" ]; then
    echo "Variable SCXA_WORKFLOW_ROOT is not defined or empty. Please load SC env."
    echo "Exiting..."
    exit 1;
fi

if [ -z "$1" ]; then
    echo "Experiment ID is not provided. Please provide EXP ID"
    echo "bash data_prep.sh <EXP-ID> [output path]"
    echo "Exiting..."
    exit 1;
fi

EXP_ID=$1

outdir="$(pwd)"

if [ "$2" ]; then
    outdir=$2
fi


echo "Creating ${outdir}/${EXP_ID} directory"
mkdir -p ${outdir}/${EXP_ID}
cd ${outdir}/${EXP_ID}

echo "Copying data to ${outdir}/${EXP_ID}"

cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/${EXP_ID}.cell_metadata.tsv cell_metadata.tsv
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/genes.tsv.gz . && gunzip -f genes.tsv.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/matrix.mtx.gz . && gunzip -f matrix.mtx.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/filtered_normalised/barcodes.tsv.gz . && gunzip -f barcodes.tsv.gz
cp ${SCXA_WORKFLOW_ROOT}/results/${EXP_ID}/*/bundle/reference/gene_annotation.txt genes_metadata.tsv

echo "Copying data for ${EXP_ID} finished"
