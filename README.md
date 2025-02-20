# Debugging University of Guleph Demuxing Pipeline

## SUMMARY
KAC and Noura AlMansouri have downloaded the University of Guleph ONT demuxxer. The raw data can be found [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.41ns1rnp1#readme)

## BIOINFORMATIC TOOLS MISSING FROM THE OMEN LAPTOP
| Bioinformatic tool    | Purpose |
| -------- | ------- |
| seqtk  | "stuff here"   |
| cutadapt  | Demuxing and filtering   |
| R  | "stuff here"   |

## Directory structure:
to keep us organised, directory structure recommended by Hebert et al, 2024 is added

/home/Guleph_barcoding
├── DATA_INPUT/
│   ├── parameters_XXXXX.xlsx
│   ├── taxonomy_XXXXX.xlsx (optional)
│   └── Raw_Reads_XXXXX.fastq.gz
├── PRIMERS/
│   └── primerDB.fasta
├── REFS/
│   └── BOLD_UniqueBINS_sintax.fasta
└── SCRIPTS/
    ├── ONT.sh
    └── SUBSCRIPTS/
        ├── get_parameters.R
        ├── primary_contig_consensus.R
        ├── final_contig_consensus.R
        ├── autotrim_taxfilter_dominants.R
        └── generate_histograms.R
