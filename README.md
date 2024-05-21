# Nanopore demux-count pipeline

A Nextflow pipeline for the processing of reads from Oxford Nanopore Technology. The pipeline locates the specified forward and reverse primer sequences and trims the read N bases upstream of the forward primer and N bases downstream of the reverse primer, where N is the barcode length.

If demultiplexing by barcodes is required, the pipeline can split the data into separate fastq files by sample. Sequences of interest can also be counted if a fasta file is provided of guide sequences.

The read structure is typically:

`[sequence][fwd_index][fwd_primer][sequence_of_interest][rev_primer][rev_index][sequence]`

The pipeline uses the following tools:

- [splitcode](https://github.com/pachterlab/splitcode) for demultiplexing
- [minimap2](https://github.com/lh3/minimap2) for alignment

## How to install (WEHI only)

The easiest way to run the pipeline is to use the [Seqera Platform](https://seqera.services.biocommons.org.au/) service provided to WEHI by Australian Biocommons. You can find more information about Seqera Platform (formerly Nextflow Tower) on WEHI's [Research Computing page](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Nextflow-Tower.aspx). See the [Configuration](https://github.com/WEHIGenomicsRnD/nf-qc-pipe#tower-configuration) section for more info.

## How to install (general users)

Clone the repository as follows:

```bash
git clone https://github.com/WEHIGenomicsRnD/nanopore-overhang.git
cd nanopore-overhang-process
```

## How to run (via command line)

If you want to run the pipeline via the command line, you can do this via:

```bash
nextflow run main.nf \
    --input_dir $path_to_fastqs \
    --outdir $output_directory \
    --fwd_primer $fwd_primer \
    --rev_primer $rev_primer \
    --primer_mismatches 3 \
    --barcode_length 12 \
    --idx_5p_mismatch 1 \
    --idx_3p_mismatch 1 \
    --index_template_file $index_file \
    --guides_fasta $guides_fasta \
    --use_db false
```

If you are running on WEHI's Milton HPC:

```bash
module load nextflow
nextflow run main.nf \
    --input_dir $path_to_fastqs \
    --outdir $output_directory \
    --fwd_primer $fwd_primer \
    --rev_primer $rev_primer \
    --primer_mismatches 3 \
    --barcode_length 12 \
    --idx_5p_mismatch 1 \
    --idx_3p_mismatch 1 \
    --index_template_file $index_file \
    --guides_fasta $guides_fasta \
    --use_db true \
    -profile log,milton
```

Note that this uses the Genomics database for the primer lookup. To do this you will need to either run through Seqera or have set your `DB` secret key to access the database.

If you are running against a large number of highly similar reference sequence, you may need to set the `--minimap_f` parameter. For exmaple, if you have 10s of thousands of sequences, set this to something like `1000`.

## Configuration

Here are the parameters you will need to set:

- `--input_dir`: full path to folder containing your read files. Files must be FASTQ format (.fq or .fastq extension). Optionally, they may be gzipped.
- `--outdir`: directory path where output is written.
- `--fwd_primer`: forward primer sequence.
- `--rev_primer`: reverse primer sequence.
- `--primer_mismatches`: how many mismatches are allowed in the primer sequences. Calculated as the levehnstein edit distance using [edlib](https://github.com/Martinsos/edlib). You may want to set this higher for longer primer sequnces.
- `--barcode_length`: how many bases to trim to the left and right of the primer sequences. If your barcode includes spacers make sure to take that into account (i.e., non-informative bases between the index and primer). Set this to 0 if you do not have barcodes.
- `--idx_5p_mismatches`: mismatches allowed in the 5' index.
- `--idx_3p_mismatches`: mismatches allowed in the 3' index.
- `--index_template_file`: if demultiplexing, use this index file to specify or lookup indexes (see below for format).
- `--guides_fasta`: (optional) fasta file contains guide sequences to count.
- `--use_db`: boolean value, default: false. Whether or not to look up indexes in the Genomics database.
- `--lenient_counts`: boolean value, default: false. If true, reads do not have to span the whole guide sequence to be counted (they will be counted as a partial map).
- `--count_only`: only perform counting of fastq files input via `--input_dir`, i.e., skip primer trimming and demultiplexing.

### Index template file format

If you are using the Genomics database for index lookup, your index file should look like this:

```
index_name
Fwd_01
Fwd_02
Rev_01
Rev_02
```

This will fetch the index names from the database. If your indexes are custom ones, or you don not want to use the database, use the following file format:

```
index_name,direction,sequence
Fwd_01,F,TAGATCGC
Fwd_02,F,CTCTCTAT
Rev_01,R,TCGCCTTA
Rev_02,R,CTAGTACG
```

Note that both sequences must match the forward direction. We do not perform any reverse complementing of the reverse sequence.
