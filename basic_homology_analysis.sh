#!/bin/bash
# Basic homology analysis using NCBI E-utilities

echo "Starting homology analysis for syn3A proteins..."

# Extract first protein sequence for testing
echo "Extracting first protein sequence..."
head -20 syn3A_proteins.fasta > first_protein.fasta

# Get the sequence without header
SEQUENCE=$(grep -v "^>" first_protein.fasta | tr -d '\n')
PROTEIN_ID=$(grep "^>" first_protein.fasta | head -1 | cut -d' ' -f1 | sed 's/>//')

echo "Analyzing protein: $PROTEIN_ID"
echo "Sequence length: ${#SEQUENCE} amino acids"

# Create a temporary file for BLAST submission
cat > blast_params.txt << EOF
CMD=Put
PROGRAM=blastp
DATABASE=nr
QUERY=$SEQUENCE
FORMAT_TYPE=XML
HITLIST_SIZE=50
EXPECT=0.01
EOF

echo "Submitting BLAST search..."
RID_RESPONSE=$(curl -s -X POST -d @blast_params.txt "https://blast.ncbi.nlm.nih.gov/Blast.cgi")

# Extract RID
RID=$(echo "$RID_RESPONSE" | grep -o "RID = [A-Z0-9]*" | cut -d' ' -f3)

if [ -z "$RID" ]; then
    echo "Failed to submit BLAST search"
    exit 1
fi

echo "BLAST job submitted with RID: $RID"

# Wait for results
echo "Waiting for BLAST results..."
for i in {1..30}; do
    echo "Checking results (attempt $i/30)..."
    
    RESULT=$(curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=$RID")
    
    if echo "$RESULT" | grep -q "Status=WAITING"; then
        echo "Still waiting..."
        sleep 15
    elif echo "$RESULT" | grep -q "Status=FAILED"; then
        echo "BLAST search failed"
        exit 1
    elif echo "$RESULT" | grep -q "Status=UNKNOWN"; then
        echo "BLAST RID expired"
        exit 1
    elif [ ${#RESULT} -gt 1000 ]; then
        echo "BLAST results received!"
        echo "$RESULT" > blast_result_${PROTEIN_ID}.xml
        
        # Count mycoplasma hits
        MYCOPLASMA_HITS=$(echo "$RESULT" | grep -i mycoplasma | wc -l)
        echo "Found $MYCOPLASMA_HITS potential mycoplasma hits"
        
        # Count total hits
        TOTAL_HITS=$(echo "$RESULT" | grep -o "<Hit>" | wc -l)
        echo "Total BLAST hits: $TOTAL_HITS"
        
        break
    else
        sleep 15
    fi
done

echo "Analysis complete for protein $PROTEIN_ID"

# Clean up
rm -f blast_params.txt

echo "Results saved to blast_result_${PROTEIN_ID}.xml"