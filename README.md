# General structure

#Source Code Overview
![module diagram](Structure.png "Source Code Overview")
 Exon hash

|key|description|
|---|-----------|
|number|the exon number. This starts from 1|
|start|the start position in the DNA sequence|
|end|the end position in the DNA sequence|

### get_overall_codon_frequencies
returns a hash of codon to counts for all of genes in chromosome 16
### calculate_custom_restriction_enzyme_sites
input: this takes two inputs

1.  gene_id
2.  enzyme cutting pattern. It should be a string contain bases, and also a caret '^' within the string.  e.g. "G^GATCC"

Returns an array of positions where this enzyme will cut.

## example usage
See `run.sh` for sample code usage

