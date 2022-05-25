# HapCheck

Inspect final marker haplotypes. Identify which haplotypes might be redundant/chimeric and are created from fragments of other two marker haplotypes. Comparisons are done for each pair.

# Input 
* `--dir` - directory
* `--haplotypes` - FASTA with all haplotype sequences. Sequences needs to be named `marker-id` f. ex. trap-1 all throughout the file.
* `--marker` - name of the marker
* `--indices` - comma separates ids for marker haplotypes

# Output
* `tmp.fasta` - FASTA with full sequences in marker haplotypes
* `tmp.snps` - FASTA with segregating SNPs in marker haplotypes


# Set-up
## Clone repository
```
git clone https://github.com/juuzia/HapCheck.git
```
## Create conda environment
```
cd HapCheck
conda env create -n HapCheck --file env.yml
conda activate HapCheck
```

## Example usage
```
python HapCheck/get_haplotype_combinations.py --dir . --haplotypes <haplotype_file>.fasta --marker trap --indices 2,20,22,23,29,41,44,46
```