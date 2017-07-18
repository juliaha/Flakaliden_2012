## Preprocessing of 16S and ITS1 rRNA data

### Activate QIIME (1.9.0) and create a working directory with the raw data (read R1/R2 and index I1/I2 files) and mapping files
### launch Python (2.7.8) and download additional scripts from https://github.com/druvus/16S-demultiplexing

### Validate the mapping file:
```
validate_mapping_file.py -m map.txt -o check_id_map_output
```

### Create separate mapping files for the two index reads with one of the custom scripts:
```
fix_mappingfile.py map.txt
```

### Do some initial quality filtering on read quality and barcode errors using default settings: 
```
split_libraries_fastq.py -o mapping_1 -i R1.fastq -b I1.fastq --rev_comp_mapping_barcodes -m map_1.txt --store_demultiplexed_fastq

split_libraries_fastq.py -o mapping_2 -i R2.fastq -b I1.fastq --rev_comp_mapping_barcodes -m map _1.txt --store_demultiplexed_fastq
```

### Sort the read files to have the same order with a custom script:
```
ln -s mapping_1/seqs.fastq only_R1.fastq 

ln -s mapping_2/seqs.fastq only_R2.fastq

syncsort_fq only_R1.fastq only_R2.fastq
```

### Activate Cutadapt (1.4.1), works with Python (2.7.3) and remove the overhanging sequencing primer sequence from reads:
#### ITS1
```
cutadapt -a ATGCTGCGTTCTTCATCGATGC only_R1.fastq_paired.fq -o only_R1_clean.fq 

cutadapt -a CTCTTGGTCATTTAGAGGAAGTAA only_R2.fastq_paired.fq -o only_R2_clean.fq
```

#### 16S
```
cutadapt -a CCGGACTACHVGGGTWTCTAAT only_R1.fastq_paired.fq -o only_R1_clean.fq

cutadapt -a GTGTGCCAGCMGCCGCGGTAA only_R2.fastq_paired.fq -o only_R2_clean.fq
```

### 16S – Merge reads using FLASH (1.2.9) (300/500 cycle kit –r 151/251)
```
flash -m 20 -M 100 -r 151 -f 253 -o extendedFrags.fastq only_R1_clean.fq  only_R2_clean.fq

syncsort_fq_readindex extendedFrags.fastq I1.fastq I2.fastq
```

### For ITS1 continue only with R1 reads

#### Bring the index read files and the reads in the same order:
```
syncsort_fq_readindex only_R1_clean.fq I1.fastq I2.fastq
```

### 16S and ITS1
#### Demultiplex both read files (R1; R2) 
#### First using the reverse index file (I1) and the corresponding mapping file of that index:
```
split_libraries_fastq.py -o mapping -i only_R1_clean.fq_synced.fq -b I1.fastq_synced.fq  --rev_comp_mapping_barcodes -m map_1.txt  --store_demultiplexed_fastq --phred_offset 33
```

### Generate separate files for each reverse barcode:
```
split_fastq.py mapping/seqs.fastq I2.fastq_synced.fq 
```

### Split every file by the forward index:
```
split_libraries_fastq.py -o mapping_${barcode} –i mapping/seqs_${barcode}.fastq -b mapping/seqs_${barcode}_barcode.fastq  --rev_comp_mapping_barcodes -m map_2.txt --store_demultiplexed_fastq --rev_comp_barcode --phred_offset 33
```

### Rename and move the output files:
```
mv mapping_${barcode}/seqs.fna mapping_all/seqs_${barcode}.fna
```

### Fix the header of each read to include the sample information:
```
fix_header.py map_.txt mapping_all/seqs_${barcode}.fna
```

### Combine all the files:
```
cat corrected_* > corrected_tmp.fna 
```

### Add consecutive numbers to all sequences:
```
fix_ID.py corrected_tmp.fna
```

### The preprocessing was performed separately for each sample type (needles, roots, soil) for 16S and ITS1. 
### At this point all non experiment related samples (mock community, extraction and amplification controls) are removed and the demultiplexed files get merged for taxonomic assignment and clustering (ITS_all.fna/16S_all.fna). 

#### Only ITS1 - Non-fungal and chimeric sequences are identified using ITSx (Version 1.0) and filtered from the file:
```
ITSx -i ITS_all.fna -o ITSx_all -p ITSx/HMMs/ -t f --cpu 20 --preserve T -E 1 --allow_single_domain 1e-5,0

filter_fasta.py -f ITS_all.fna -o ITS_all_nochimera.fna -s all_chimera.txt -n

filter_fasta.py -f ITS_all_nochimera.fna -o ITS_all_nochimera_2.fna -s all_no_detect.txt –n
```

#### ITSx displays in detail the composition of each read as of SSU,ITS1,5.8S,ITS2 region and provides separate files containing reads with the complete ITS1 regions only. These were first filtered from the rest of the sequences, as they are already cleaned of the conservative SSU region, which might complicate taxonomy assignment or cluster formation:
```
grep "^>" all_ITS1.fasta | sed -e 's/>//g;s/ /\t/g' | cut -f 1 > all_ITS1_header.txt

filter_fasta.py -f ITS_all_nochimera_2.fna -o ITS_all_nochimera_3.fna -s all_ITS1_header.txt -n 
```

#### Using ITS1F-ITS2 primers and reads of the ~length of 250bp, reads should contain parts of the conservative SSU, partial or full ITS1 region and a partial 5.8S region. The length of the conservative SSU was in all reads 46bp and was removed from all reads, before merging with the clean ITS1 region reads from the previous step:
```
awk '/^#/ {next} /^>/ { print $0 } /^[^>]/ { print substr($0, 47, length($0) - 46) }' ITS_all_nochimera_3.fna > ITS_all_nochimera_trim.fna

cat split_all/all_ITS1.fasta ITS_all_nochimera_trim.fna > ITS_all_ITS1.fna
```

## 16s and ITS1 - clustering with VSEARCH (Version 1.10.2) 
### First dereplication of reads:
```
vsearch --derep_fulllength ITS_all_ITS1.fna --output derep_ITS1.fa --sizeout &
```

### Size sorting (removes singletons):
```
vsearch --sortbysize derep_ITS1.fa --output sorted_ITS1.fa --minsize 2 --sizein --relabel OUT
```

### Cluster at 95% (ITS1)/ 97% (16S):
```
vsearch --cluster_size sorted_ITS1.fa --id 0.95 --strand both --centroids ITS1_centroids.fa --relabel OTU_ --sizeout --uc ITS1_centroids.uc
```

### Perform another chimera checking (de novo):
```
vsearch --uchime_denovo ITS1_centroids.fa --nonchimera ITS1_v_nochimera.fa --chimeras ITS1_v_chimera.fa
```

### Populate clusters (--id 0.95 ITS1/ 0.97 16S):
```
vsearch --usearch_global ITS_all_ITS1.fna --db ITS1_v_nochimera.fa --id 0.95 --strand both --uc ITS1_v_.uc
```

### Transform into OTU table (tab delimited) using script from here https://github.com/leffj/helper-code-for-uparse/blob/master/create_otu_table_from_uc_file.py:
```
create_otu_table_from_uc_file.py -i ITS1_v_.uc -o ITS1_v_.txt
```

### Taxonomy assignment using BLAST (2.2.26) and QIIME and UNITE (Version 2016_01_31)/ SILVA (Version 119) database files:
```
biom convert --table-type="OTU table" -i ITS1_v_.txt -o ITS1_v_.biom --to-json 
```

#### ITS1:
```
parallel_assign_taxonomy_blast.py -i ITS1_v_nochimera.fa -o v_nochimera_taxonomy -T --jobs_to_start 10 --reference_seqs_fp qiime_databases/unite/v9/sh_refs_qiime_ver7_dynamic_31.01.2016.fasta --id_to_taxonomy_fp qiime_databases/unite/v9/sh_taxonomy_qiime_ver7_dynamic_31.01.2016.txt 
```
##### No. samples: 327 (including three mock community samples)
##### No. observations: 5166
##### Total read count: 17565367

#### 16S:
```
assign_taxonomy.py -i 16S_v_nochimera.fa -o rdp_assigned_taxonomy -m rdp --rdp_max_memory 20000 --reference_seqs_fp qiime_databases/silva/Silva_119/Silva119_release/rep_set/97/Silva_119_rep_set97.fna --id_to_taxonomy_fp qiime_databases/silva/Silva_119/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt 

parallel_align_seqs_pynast.py -i 16S_v_nochimera.fa --template_fp qiime_databases/silva/Silva_119/Silva119_release/core_alignment/core_Silva119_alignment.fna --min_percent_id 0.60 -o pynast_aligned_seqs -T --jobs_to_start 4

filter_alignment.py -o pynast_aligned_seqs/ -i pynast_aligned_seqs/16S_v_nochimera_aligned.fasta 

make_phylogeny.py -i pynast_aligned_seqs/16S_v_nochimera_aligned_pfiltered.fasta -o pynast_aligned_seqs/rep_set.tre
```

##### No. samples: 333 (including each three mock, negative water control, ecoli culture samples)
##### No. observations: 22755
##### Total read count: 22673362



### ITS1 and 16S - Remove every otu with counts less than or with 10 reads in every out using the script filter_otus_per_sample.py (https://gist.github.com/adamrp/7591573):
```
filter_otus_per_sample.py -i ITS1_v_.biom -o otu_v_10.biom -n 11
```
#### ITS1:
```
biom add-metadata -i otu_v_10.biom --observation-metadata-fp v_nochimera_taxonomy/ITS1_v_nochimera_tax_assignments.txt -o otu_table_v_10_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy 
```
##### No. observations: 2630
##### Total read count: 17263282

#### 16S:
```
biom add-metadata -i 16S_v_10.biom --observation-metadata-fp rdp_assigned_taxonomy/16S_v_nochimera_tax_assignments.txt -o 16S_otu_table_v_10_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy

filter_otus_from_otu_table.py -i 16S_otu_table_v_10_tax.biom -o 16S_otu_pynast.biom -e pynast_aligned_seqs/16S_v_nochimera_failures.fasta
```

##### No. observations: 7739
##### Total read count: 20446149

#### ITS1 - Remove reads with an assignments other than the kingdom fungi and reads without taxonomic assignment:
```
filter_taxa_from_otu_table.py -i otu_v_10.biom -o otu_v_10_2.biom -n k__Protista 

filter_taxa_from_otu_table.py -i otu_v_10_2.biom -o otu_v_10_3.biom -n "No blast hit" 
```
##### No. observations: 1867
##### Total read count: 15276438

#### 16S - Remove plant organellar and Archaea sequences and reads without taxonomic assignment:
```
filter_taxa_from_otu_table.py -i 16S_otu_pynast.biom -o 16S_all.biom -n D_2__Chloroplast,D_4__mitochondria,D_0__Archaea &

filter_taxa_from_otu_table.py -i 16S_all.biom -o 16S_all.biom -n "Unclassified" &
```

##### No. observations: 6704
##### Total read count: 16809660


### ITS1 and 16S - Split into different sample types (needle, root or soil), remove singletons after:
```
filter_samples_from_otu_table.py -i otu_v_10_3.biom -o needle_v_10.biom -m all_tissue_map.txt -s 'SampleType:needle' 

filter_otus_from_otu_table.py -i needle_v_10.biom -o needle_v_10_1.biom -n 1
```

##### Needle ITS1 
##### No. samples: 108
##### No. observations: 1191
##### Total read count: 3036165

```
filter_samples_from_otu_table.py -i otu_v_10_3.biom -o soil_v_10.biom -m all_tissue_map.txt -s 'SampleType:soil'
```

##### after –n 1:
##### Soil ITS1 
##### No. samples: 108
##### No. observations: 844
##### Total read count: 5968902

```
filter_samples_from_otu_table.py -i otu_v_10_3.biom -o root_v_10.biom -m all_tissue_map.txt -s 'SampleType:root'
```
##### after –n 1:
##### Root ITS1 
##### No. samples: 108
##### No. observations: 691
##### Total read count: 6112784

##### 16S Soil:
##### No. observations: 5331
##### Total read count: 7040923

##### 16S Needle:
##### No. observations: 1009
##### Total read count: 1202870

##### 16S Root:
##### No. observations: 3891
##### Total read count: 7946020

### Restrict to 0.005% abundance per sample type:
#### For ITS1 needles (at least 151 reads of 3036165)
#### For ITS1 soil (at least 298 reads of 5968902)
#### For ITS1 roots (at least 305 reads of 6112784)
```
filter_otus_from_otu_table.py -i #tissue#_v_10_1.biom -o #tissue#_v_10_000005.biom --min_count_fraction 0.00005
```

##### Soil ITS1 
##### No. observations: 360
##### Total read count: 5930724 

##### Needle ITS1 
##### No. observations: 407
##### Total read count: 3001175

##### Root ITS1 
##### No.  observations: 311
##### Total read count: 6081910

##### 16S soil:
##### No. observations: 1579
##### Total read count: 6740639

##### 16S Needle:
##### No. observations: 408
##### Total read count: 1188901

##### 16S Root:
##### No. observations: 1240
##### Total read count: 7704727

### Combine technical replicates:
```
collapse_samples.py -b #tissue#.biom -m map.txt --output_biom_fp #tissue#bioreps.biom --output_mapping_fp bioreps_map.txt --collapse_fields SampleType,Treatment,Date,SoilPlot &
```

### Rarefy ITS1/16S (39000/16500):
```
single_rarefaction.py -i #tissue#bioreps.biom -o #tissue#_obs_39000.biom -d 39000 &
```
### Continue with analysis in R
