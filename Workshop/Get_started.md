# Haney Lab RNA-Seq Analysis Workshop

## Downloads
Please download:
* [R >= 3.5.0](https://cran.r-project.org)
* [RStudio](https://www.rstudio.com/products/rstudio/download/) (Desktop version)

## R Packages
Open RStudio and try to get familiar with the layout. By default, on the left side is a console. This is where you can type
code. Type into the console
```R
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```
