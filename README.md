# kappa

### Description
This repository contains the files to typset my doctoral thesis (or kappa) titled **"Novel methods for dose-response meta-analysis"** (Karolinska Institutet, 2018). It is based on the tex code written by [Andrea Discacciati](https://github.com/anddis/phd-thesis).

### Structure
The thesis is written combining [`R`](https://www.r-project.org/) and [LaTex](http://www.latex-project.org/) by generating the tex files with [`knitr`](https://yihui.name/knitr/).

The thesis *Rnw* has been deviding into multiple *Rnw* input files. The root/master file is `kappa.Rnw`.
The main folders are:
- chapters: different *Rnw*s chapters of the kappa
- additional: useful *tex* files, included *bib* library
- style: styles for latex sources
- packrat: for handling `R` package dependencies and enhancing reproducibility
- figure: pictures and time consuming figures
- data
- code

### Author
Alessio Crippa,  Karolinska Institutet, Stockholm, Sweden
