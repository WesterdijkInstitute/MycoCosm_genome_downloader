# MycoCosm genome downloader

Helps to download genomes (fasta + gff) files from MycoCosm.


# Requirements

* `ete3` and `libxml2`. Alternatively, you can use the provided `yml` file with conda (`conda env create --file mycocosmdownloader.yml`).
* [`curl`](https://curl.se/) (should already be installed in Mac or Linux systems)
* A JGI account


# Usage (tested 2022-03)

## Update taxonomy database

As part of the output of this program, a "taxonomy" file is created, which includes the lineage for each species (this lineage is also used to create the folder structure for the outpu). For this, ete3 is used. Use the following command to update the database:
```
python mycocosm_genome_downloader.py --update
```

## Obtain data files

Obtain the [list of genomes](https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=asc), and the list of files:
```
python mycocosm_genome_downloader.py --getgenomelist --getxml
```

The list of genomes is a comma-separated value file which enlists all current projects ("**portals**") in MycoCosm (`MycoCosm_Genome_list.csv`, Note: it uses ISO-8859-15 encoding); The list of files is stored as `MycoCosm_data.xml` (Note: this is ~40Mb file which may take several tries to download. After downloading, it is re-formatted to be human-readable).

These files will be stored in the folder specified by the `--outputfolder` argument (default: 'output'). You will need your credentials to retrieve the file list.


## Download data

Once you have the genome and file list, use them to download the data with
```
python mycocosm_genome_downloader.py --csv [path to MycoCosm_Genome_list.csv] --xml [path to MycoCosm_data.xml]
```

This will create a folder structure that follows the fungal taxonomy tree as MycoCosm's [main page](https://mycocosm.jgi.doe.gov/mycocosm/home)

Optional: use parameter `--simulate` to only create the structure, but don't download any file.

Note: as of 2022-01, the downloaded files sum about 33 GiB


## Output

Each 'leaf' folder within the folder structure represents a genome from MycoCosm and contains the fasta + gff files (compressed as .gz). Additionally, three files are found inside the base output folder:
* `JGI_taxonomy.tsv`. Columns: `Short name` (Portal), `Accession` (same as short name in this case), `TaxId` (from NCBI), `Name` (project name), `Path` (relative path from output folder to Portal), `Assembly file`, `GFF file`, `lineage` (comma-separated)
* `JGI_download_list_[date].txt`: All the files that were downloaded
* `List_gene_gff_filenames.txt`: All gff files for each portal


### Optional: re-use data

If you already have downloaded data, you can simply use your local copy of the files instead of downloading them. This makes it easier to update with new genomes. First, use option `--getprevious`, which will scan a base folder with previous results and create a file called `previously_downloaded_files.tsv`. Use this file with option `--previous` to skip files already downloaded

## Known limitations, notes

* The incorrect gff might have been chosen
* Some old gff files can't be processed correctly by other tools (e.g. antiSMASH)
* Can't handle dikaryon genomes ("primary/secondary alleles")

## More details

### Target download folders

These are the folders that are scanned for files, as they appear while in each Portal's download section:

* Assemblies: "Assembled scaffolds (unmasked)"
* Gene annotations: "Filtered Models ('best')"

### Choosing the right annotation file

There are usually a few files to choose from in the "Filtered Models ('best')" folder. The criteria to choose a single file is implemented in the following order:
* In some cases, the choosing was done manually and incorporated to the `hardcoded_gff_files.tsv` file. If the Portal has an entry here, use the indicated file. This is a work in progress!
* Skip some specific files: "Aciri1_meta_GeneCatalog_genes_20111216.gff.gz", "Exoaq1_GeneCatalog_20160901.gff3.gz", "Exoaq1_GeneCatalog_20160828.gff3.gz", "Fonpe1_GeneCatalog_20160901.gff3.gz", "Copmic2_FM1_removed_alleles.gff.gz"
* Skip files with `proteins` or `secondary_alleles` in their name
* Skip files with the following extensions: `gtf.gz`, `tgz` and any other that is *not* `gz`
* Prefer `gff3` files
* If we have two `gff3` (or two `gff`) files, prefer the most recent


### Ignored assembly files / portals

* Files with `MitoAssembly`, `MitoScaffolds`, `PrimaryAssemblyScaffolds` or `SecondaryAssemblyScaffolds` in their name
* Excluded assemblies (hardcoded). Ignore these filenames as they're not related to assemblies, are old versions or
 are assemblies of meta-samples: "1034997.Tuber_borchii_Tbo3840.standard.main.scaffolds.fasta.gz", "Spofi1.draft.mito.scaffolds.fasta.gz", "Patat1.draft.mito.scaffolds.fasta.gz", "PleosPC9_1_Assembly_scaffolds.fasta.gz", "Neuhi1_PlasmidAssemblyScaffolds.fasta.gz", "CocheC5_1_assembly_scaffolds.fasta.gz", "Alternaria_brassicicola_masked_assembly.fasta.gz", "Aciri1_meta_AssemblyScaffolds.fasta.gz", "Rhoto_IFO0880_2_AssemblyScaffolds.fasta.gz"
* Ignored portals (hardcoded). Metaprojects or old versions: "Rhoto_IFO0880_2", "Aciri1_meta", "Pospl1"
