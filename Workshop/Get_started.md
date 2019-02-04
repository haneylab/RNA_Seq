# Haney Lab RNA-Seq Analysis Workshop

## Downloads
Please download:
* [R >= 3.5.0](https://cran.r-project.org)
* [RStudio](https://www.rstudio.com/products/rstudio/download/) (Desktop version)
* Files – from OwnCloud in the `Data_Andrew/` directory

It is helpful if you download the workshop_files into a new folder dedicated to this workshop.  

## R Packages
Open RStudio and try to get familiar with the layout. By default, on the left side is a console. This is where you can type
code. Type into the console:
```R
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("DESeq2", "pheatmap", "ashr"))
```
