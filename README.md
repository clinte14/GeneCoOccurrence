<h1 align="center">GeneCoOccurrence</h1>

<p align="center">
  <img width="256" height="256" src="https://user-images.githubusercontent.com/35710809/167226809-ea5ec455-674c-4111-a1b0-acb84f29b3ee.png">
</p>
<p align="center"><b>
A Tool for Calculating Bacterial Gene Co-occurrence Without Phylogenetic Inference
</p></b>

[Description](#description) <br>
[Installation](#installation) <br>
[Basic Usage](#basic-usage) <br>
[Workflow](#workflow) <br>
[Tutorial With Sample Data](#tutorial-with-sample-data)

## Description
**Jan 25th 2024:** A link to bioRxiv software announcement paper will be provided soon.

 **An example of GeneCoOccurrence usage** is found in this [Nature Microbiology](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9830645) study of unique genomic islands in the current pandemic strain of *Vibrio cholerae* .
 
 **GeneCoOccurrence enables investigation of the possible functional pairing of co-occurring bacterial genes on genetic elements that lack informative phylogeny.** This uniquely allows investigation of horizontally acquired regions of DNA in genomic islands, mobile genetic elements, etc. The software works by calculating the frequency at which a pair of genes are jointly present within individual genomes across a given set of sequenced bacteria. The output includes co-occurrence scores, heatmaps, and graphical networks to provide context to the co-occurrence of gene pairs.

 *See [Tutorial With Sample Data](#tutorial-with-sample-data) below to see sample input and step-by-step examples.*
## Installation
We recommend creating an isolated environment using [Conda](https://www.anaconda.com/download), followed by installation using the Python package manager pip.
1. Create Conda enviornment named 'gco' (i.e. GeneCoOccurrence) 
```
conda create --name gco python=3.11
```
2. Enter the 'gco' isolated enviornment
```
conda activate gco
```

3. Install GeneCoOccurrence using pip with our package at testpypi
```
python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ genecooccurrence
```

## Basic Usage

```
gco -i blast_output.xml -o /home/name/my_project -c protein_id_to_common_names.csv
```

* `-i <Input File>` *required*

  BLAST results of your genes of interest (saved as a single .xml file, which is an available output format for BLAST on NCBI website and stand-alone command line tool).
  **OR**
Any binary presence/absence matrix in comma-separated values (.csv) file. It should be utf-8 formatted, comma separated, and have a '.csv' suffix.
* `-o <Output Directory>` *optional* 
  *Optional: If not provided defaults to current working directory*
* `-c <Common Name>` *optional* 
  *Optional:  A .csv file that allows conversion of BLAST query protein ID's to common gene names of your choosing.* 
* `-c <Common Name>` *optional* 
  *Optional:  A .csv file that allows conversion of BLAST query protein ID's to common gene names of your choosing.* 

## Workflow
![image](https://github.com/clinte14/GeneCoOccurrence/assets/35710809/98b89fc7-4e34-4efc-befb-c476e5ceec60)
**Workflow of GeneCoOccurrence**. Output folders in grey, input file in blue, intermediate output files in red, final output files in green. **A.** User input is either a presence/absence matrix OR BLAST results. **B.** A presence/absence matrix is generated if BLAST results were chosen as input. **C.** The co-occurrence for all GOIs *i* to *j* is summed and fed into a Pearson Correlation followed by a partial correlation correction which results in co-occurrence score. **D.** Output includes a co-occurrence heatmap of all genes *i* to *j*, maximum related subnetwork visuals, and a co-occurrence table.

## Tutorial with Sample Data

[1. Using BLAST Results](#1-using-blast-results)
*You may wish to use BLAST results of genes of interest to calculate gene pair co-occurrence, most likely to infer functional relatedness.*

[2. Using A Custom Presence/Absence Matrix](#2-using-a-custom-presenceabsence-matrix)
*You may wish to use a custom presence/absence matrix, for example to calculate the co-occurrence of antibiotic resistance genes within metagenomic data, etc.*

### 1. Using BLAST Results
1. First [install](#installation) GeneCoOccurrence as described above.
2. Right-click this [link](https://github.com/clinte14/GeneCoOccurrence/blob/master/tutorial.zip) and choose 'Save As' to save a zipped directory containing the example data and output.
3. Navigate to and unzip the previous directory, which should have this structure:
```
└── tutorial
    ├── manure_metagenomic_study
    │  ├── unenriched_gene
    │  └── unenriched_gene_matrix.csv
    ├── README.txt
    └── vibrio_VSPI_study
        ├── prot_ID_to_common_name.csv
        ├── VSPI
        └── VSPI_genes_BLAST_results.xml
```
4. Enter the *vibrio_VSPI_study* directory, which contains the two input files: *VSPI_genes_BLAST_results.xml* and *prot_ID_to_common_name.csv*
5. Run the command:
```
gco -i ./VSPI_genes_BLAST_results.xml -o ./VSPI -c ./protein_id_to_common_names.csv
```
6. Note this command will overwrite the existing output directory *VSPI*. See [workflow](#workflow) above to see the output files available in directory *VSPI*.

**Helpful Information for Generating Your Own Compatible BLAST Data**
It is important that BLAST matches have enough sequence similarity to the GOI to reasonably infer homology. For this reason, the BLAST program (*e.g.,* BLASTp or BLASTn) and E-value cutoff value should be carefully chosen given the context of your organism(s) and possible MGEs of interest. Our study into functionally related gene pairs in [*Vibrio cholerae*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9830645/) used BLASTp with an E-value cutoff of 10^-4 in the NCBI protein non-redundant (nr) database limited to taxid:2 (bacteria). Additionally, we note that the online BLAST tool provided by NCBI limits the total number of total hits returned. This means that hits meeting the search criteria could be arbitrarily dropped, reducing the input available to GeneCoOccurrence to make predictions. In these cases, users should consider using [localized BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) searches.  **The BLAST search should have an output format of ‘5’ (single-file .xml).**


**Helpful Information for Creating a Common Name .csv**
We offer the ability to convert less helpful BLAST query gene IDs (_e.g. BAF33440_) to common gene names (_e.g. KfrC_). These common gene names are subsequently used in intermediate and final output files (including visualizations). To enable this feature, use the ‘-c’ flag (_e.g._ -c _gene_name_conversions.csv_) pointing to a comma-separated .csv file. This file should be formatted with query gene IDs in the first column and corresponding common gene names on an equivalent row in the second column.


### 2. Using a Custom Presence/Absence Matrix
1. First [install](#installation) GeneCoOccurrence as described above.
2. Right-click this [link](https://github.com/clinte14/GeneCoOccurrence/blob/master/tutorial.zip) and choose 'Save As' to save a zipped directory containing the example data and output.
3. Navigate to and unzip the previous directory, which should have this structure:
```
└── tutorial
    ├── manure_metagenomic_study
    │  ├── unenriched_gene
    │  └── unenriched_gene_matrix.csv
    ├── README.txt
    └── vibrio_VSPI_study
        ├── prot_ID_to_common_name.csv
        ├── VSPI
        └── VSPI_genes_BLAST_results.xml
```
4. Enter the *manure_metagenomic_study* directory, which contains the single input file: *unenriched_gene_matrix.csv*
5. Run the command:
```
gco -i ./unenriched_gene_matrix.csv -o ./unenriched_gene
```
6. Note this command will overwrite the existing output directory *unenriched_gene*. See [workflow](#workflow) above to see the output files available in directory *unenriched_gene*.

## Authors
Copyright (C) 2018-2024  Clinton A. Elg
Please contact via Github



## License

This project is licensed under the GNU General Public License, Version 2. Please see the <a href="https://github.com/clinte14/GeneCoOccurrence/blob/master/LICENSE.md" title="LICENSE.md">LICENSE.md</a> link for details.


## Acknowledgments
<a href="https://www.flaticon.com/free-icons/genetic" title="genetic icons">Genetic icons created by Freepik - Flaticon</a>
<!--
Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
THIS README from https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc
-->
