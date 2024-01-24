
<h1 align="center">GeneCoOccurrence</h1>

<p align="center">
  <img width="256" height="256" src="https://user-images.githubusercontent.com/35710809/167226809-ea5ec455-674c-4111-a1b0-acb84f29b3ee.png">
</p>
<p align="center"><b>
A Tool for Exploring Relationships Within Small Sets of Unclassified Microbial Genes
</p></b>


## Description

GeneCoOccurrence provides an accessible tool for biologist to explore interactions and functional relationships between small sets of unclassified genes. It works by calculating the co-occurrence of homologues genes in individual strains from across the sequenced bacterial domain. **May 6th, 2022: Expect a biorxiv draft soon!**


## Getting Started

### Installing
For now, you must install dependencies (see below) and run via Python3 interpreter in a BASH shell (Linux/Unix/MacOS/WSL2). **May 6th, 2022: Docker image coming soon.**

### Basic Usage

<pre>
main.py -i blast_output.xml -o /home/name/my_project
</pre>
* `-i <Input File>`<br>
  *Required*<br>
  BLAST results of your genes of interest (saved as a single .xml file, which is an available output format for BLAST on NCBI website and stand-alone command line tool).

* `-o <Output Directory>`<br>
  *Optional: If not provided defaults to current working directory*<br>
  Desired output directory

### Helpful Flags / Options
<pre>
main.py -i blast_output.xml -o /home/name/my_project <i>-c prot_ID_to_common_name.csv</i>
</pre>
* `-c <Common Name>`<br>
  *Optional: See details below for default (a) or secondary default (b)*<br>
 A .csv file that allows conversion of BLAST query protein ID's to common gene names of your choosing. This ensures all visuals created by GCO use names you chose rather than *(a)* names stripped from the BLAST output (if available), or *(b)* reference numbers from original BLAST query (last resort).<br>
 
  * The .csv should be comma seperated with the format of `reference_number,common_name`, *i.e.:*
     ```
    BAF33440.1,kfrC
    BAF33441.1,kfrB
    BAF33442.1,kfrA
    BAF33443.1,korB
    BAF33444.1,incC1
    BAF33445.1,incC2
     ```

### Workflow
![image](https://github.com/clinte14/GeneCoOccurrence/assets/35710809/98b89fc7-4e34-4efc-befb-c476e5ceec60)
**Workflow of GeneCoOccurrence**. Output folders in grey, input file in blue, intermediate output files in red, final output files in green. **A.** User input is either a presence/absence matrix OR BLAST results. **B.** A presence/absence matrix is generated if BLAST results were chosen as input. **C.** The co-occurrence for all GOIs *i* to *j* is summed and fed into a Pearson Correlation followed by a partial correlation correction which results in co-occurrence score. **D.** Output includes a co-occurrence heatmap of all genes *i* to *j*, maximum related subnetwork visuals, and a co-occurrence table.

## Tutorial with Sample Data
<!-- Any advise for common problems or issues.
```
command to run if program contains helper info
``` -->


## Authors
Copyright (C) 2018-2022  Clinton A. Elg
Please contact via Github

## Version History
Currently pre-production / alpha version

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
