# Nanopore Overhang Preprocess

A Nextflow pipeline for the preprocessing of reads using the overhang protocol run on Oxford Nanopore Technology. The pipeline locates the specified forward and reverse primer sequences and trims the read N bases upstream of the forward primer and N bases downstream of the reverse primer. The 5' and 3' ends will then contain index information ready for demultiplexing.

The read structure is typically `[sequence][fwd_index][fwd_primer][sequence_of_interest][rev_primer][rev_index][sequence]`. The pipeline trims off the 5' and 3' superfluous 'sequence'. 

## How to install

The easiest way to run the pipeline is to use the [Nextflow Tower](https://tower.services.biocommons.org.au) service provided to us by Australian Biocommons. You can find more information about Tower on WEHI's [Research Computing page](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Nextflow-Tower.aspx). See the [Tower Configuration](https://github.com/WEHIGenomicsRnD/nf-qc-pipe#tower-configuration) section for more info.

If you want to install manually, you can clone the repository as follows:

```bash
git clone https://github.com/WEHIGenomicsRnD/nanopore-overhang-process.git
cd nanopore-overhang-process
```

## How to run (Nextflow pipeline via command line)

If you want to run the pipeline via the command line, you can do this via:

```bash
module load nextflow
nextflow run main.nf --input_dir $path_to_fastqs \
    --outdir $output_directory \
    --fwd_primer AGCTGTCGTAGTGGTAGTC \
    --rev_primer GTCTGATCGTGCTGCTGAT \
    --mismatches 3 \
    --barcode_length 12 \
    -profile log,milton
```

## Configuration

Here are the parameters you will need to set:

- `--inputdir`: full path to folder containing your read files. Files must be FASTQ format (.fq or .fastq extension). Optionally, they may be gzipped.
- `--output_dir`: directory path where output is written.
- `--fwd_primer`: forward primer sequence.
- `--rev_primer`: reverse primer sequence.
- `--mismatches`: how many mismatches are allowed in the primer sequences. Calculated as the levehnstein edit distance using [edlib](https://github.com/Martinsos/edlib).
- `--barcode_length`: how many bases to trim to the left and right of the primer sequences. If your barcode includes spacers make sure to take that into account (i.e., non-informative bases between the index and primer).
- `--log`: enables logging, which writes the IDs of failed reads to a log file and specifies why each read could not be processed.
