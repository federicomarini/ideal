*This information is also contained in the `ideal` package vignette.*

`ideal` is a Bioconductor package containing a Shiny application for
interactively analyzing RNA-seq expression data, by interactive exploration of the 
results of a Differential Expression analysis.




## Getting started

*[ideal](http://bioconductor.org/packages/ideal)* is an R package distributed as part of the [Bioconductor](http://bioconductor.org) project. To install the package, start R and enter:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ideal")
```

If you prefer, you can install and use the development version, which can be retrieved via Github (https://github.com/federicomarini/ideal). To do so, use


```r
library("devtools")
install_github("federicomarini/ideal")
```

Once *ideal* is installed, it can be loaded by the following command.


```r
library("ideal")
```



## Introduction

*[ideal](http://bioconductor.org/packages/ideal)* is a Bioconductor package containing a Shiny application for analyzing RNA-Seq data in the context of differential expression. This enables an interactive and at the same time analysis, keeping the functionality accessible, and yet providing a comprehensive selection of graphs and tables to mine the dataset at hand.

*[ideal](http://bioconductor.org/packages/ideal)* is an R package which fully leverages the infrastructure of the Bioconductor project in order to deliver an interactive yet reproducible analysis for the detection of differentially expressed genes in RNA-Seq datasets. Graphs, tables, and interactive HTML reports can be readily exported and shared across collaborators. The dynamic user interface displays a broad level of content and information, subdivided by thematic tasks. All in all, it aims to enforce a proper analysis, by reaching out both life scientists and experienced bioinformaticians, and also fosters the communication between the two sides, offering robust statistical methods and high standard of accessible documentation.

It is structured in a similar way to the *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)*, also designed as an interactive companion tool for RNA-seq analysis focused rather on the exploratory data analysis e.g. using principal components analysis as a main tool.

The interactive/reactive design of the app, with a dynamically generated user interface makes it easy and immediate to apply the gold standard methods (in the current implementation, based on *[DESeq2](http://bioconductor.org/packages/DESeq2)*) in a way that is information-rich and accessible also to the bench biologist, while also providing additional insight also for the experienced data analyst. Reproducibility is supported via state saving and automated report generation.

### Citation info

If you use *[ideal](http://bioconductor.org/packages/ideal)* for your analysis, please cite it as here below:


```r
citation("ideal")
```

```

To cite package 'ideal' in publications use:

  Federico Marini (2017). ideal: Interactive Differential Expression
  AnaLysis. R package version 0.9.0.
  https://github.com/federicomarini/ideal

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {ideal: Interactive Differential Expression AnaLysis},
    author = {Federico Marini},
    year = {2017},
    note = {R package version 0.9.0},
    url = {https://github.com/federicomarini/ideal},
  }
```



## Using the application

There are different ways to use `ideal` for interactive differential expression analysis.

### Launching `ideal` locally

First load the library


```r
library("ideal")
```

and then launch the app with the `ideal` function. This takes the following essential parameters as input:

- `dds_obj` - a `DESeqDataSet` object. If not provided, then a `countmatrix` and a `expdesign` need to be provided. If none of the above is provided, it is possible to upload the data during the execution of the Shiny App
- `res_obj` -  a `DESeqResults` object. If not provided, it can be computed during the execution of the application
- `annotation_obj` - a `data.frame` object, with row.names as gene identifiers (e.g. ENSEMBL ids) and a column, `gene_name`, containing e.g. HGNC-based gene
symbols. If not provided, it can be constructed during the execution via the `org.eg.XX.db` packages
- `countmatrix` - a count matrix, with genes as rows and samples as columns. If not provided, it is possible to upload the data during the execution of the Shiny App
- `expdesign` -a `data.frame` containing the info on the experimental covariates of each sample. If not provided, it is possible to upload the data during the execution of the Shiny App

Different modalities are supported to launch the application:

- `ideal(dds_obj = dds, res_obj = res, annotation_obj = anno)`, where the objects are precomputed in the current session and provided as parameters
- `ideal(dds_obj = dds)`, as in the command above, but where the result object is assembled at runtime 
- `ideal(countmatrix = countmatrix, expdesign = expdesign)`, where instead of passing the defined `DESeqDataSet` object, its components are given, namely the count matrix (e.g. generated after a run of featureCounts or HTSeq-count) and a data frame with the experimental covariates. The design formula can be constructed interactively at runtime
- `ideal()`, where the count matrix and experimental design can simply be uploaded at runtime, where all the derived objects can be extracted and computed live. These files have to be formatted as tabular text files, and a function in the package tries to guess the separator, based on heuristics of occurrencies per line of commonly used characters

### Accessing the public instance of `ideal` 

To use *[ideal](http://bioconductor.org/packages/ideal)* without installing any additional software, you can access the public instance of the Shiny Server made available at the Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI) in Mainz.

This resource is accessible at this address: 

http://shiny.imbei.uni-mainz.de:3838/ideal

### Deploying to a Shiny Server

A deployment-oriented version of the package is available at https://github.com/federicomarini/ideal_serveredition. This repository contains also detailed instruction to setup the running instance of a Shiny Server, where `ideal` can be run without further installation for the end-users.

Please note that you still need `ideal` to be installed there once during the setup phase - for this operation, you might require root administrator permissions. 


<!-- Following the structure of \texttt{pcaExplorer}, the deployment as a standalone web application is easily achieved also for \texttt{ideal}, with the only requirement of a running Shiny Server. Handling of sensitive human patient data can thus occur in a secured way, e.g. with restricted access to an internal instance of the application. This will also avoid the need of external software installation for the users, which could moreover seamlessly access higher performing computing infrastructures. -->

<!-- This implementation of \texttt{ideal} is made available at github.com/federicomarini/ideal\_serveredition  -->



## Getting to know the user interace and the functionality

The user interface is dynamically displayed according to the provided and computed objects, with tabs that are actively usable only once the required input is effectively available.

Moreover, for some relevant UI widgets, the user can receive additional information by hovering over with the mouse, with the functionality powered by the *[shinyBS](http://cran.fhcrc.org/web/packages/shinyBS/index.html)* package.

For the user which is either new with the app UI/functionality, or not extensively familiar with the topic of differential expression, it is possible to obtain a small *guided tour* of the App by clicking on the respective help buttons, marked in the app like this.

<button id="btn" type="button" class="btn btn-default action-button" style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4">
  <i class="fa fa-info"></i>
  Click me for a quick tour
</button>


These trigger the start of a step-by-step guide and feature introduction, powered by the *[rintrojs](http://cran.fhcrc.org/web/packages/rintrojs/index.html)* package.

### The controls sidebar

Some of the input controls which affect different tabs are located in the sidebar, while others are as well in the individual tabs of the app. By changing one or more of the input parameters, the user can get a fine control on what is computed and displayed.

#### App settings

- **Group/color by** - Select the group of samples to stratify the analysis for plotting. Can also assume multiple values.
- **Select the gene(s) of interest - ids** - Select a subset of genes for deeper analysis. If an annotation object is provided, the user can handily select the genes e.g. based on their HGNC symbol
- **False Discovery Rate** - Set as default to 0.05, it is the FDR value for the Benjamini-Hochberg procedure for adjusting p-values in the multiple testing comparison scenario

#### Plot export settings   

**Width** and **Height** for the figures to export are input here in cm.

#### Quick viewer

This displays a list of the underlying objects with which basically all of the analysis can be performed. A green tick icon appears close to each when the respective component is either provided or calculated. For obtaining the best analysis experience in `ideal`, it is recommended to provide all of them.

#### First steps help

Clicking on this button activated the `intro.js` based tour for getting to know the components and the structure of the app. Dedicated step-by-step procedures are also available in each individual tab.

### The task menu

The task menu, accessible by clicking on the cog icon in the upper right part of the application, provides two functionalities:

- `Exit ideal & save` will close the application and store the content of the `input` and `values` reactive objects in a list of two elements in the `ideal_env` environment, respectively called `ideal_inputs_YYYYMMDD_HHMMSS` and `ideal_values_YYYYMMDD_HHMMSS`
- `Save State as .RData` will similarly store `LiveInputs` and `r_data` in a binary file named `idealState_YYYYMMDD_HHMMSS.Rdata`, without closing the application 


           
## The main app panels

The *[ideal](http://bioconductor.org/packages/ideal)* app is a one-paged dashboard, structured in different panels, where each of them is focused on a different aspect of the data exploration. 

On top of the panels, three `valueBox` objects serve as guiding elements for having an overview of the data at hand: how many genes and samples are in the data, how many entries are in the annotation object, and how many genes were found to be differentially expressed in the results. Whenever each of the underlying objects is available, the background color turns from red to green.

For the main analysis, the available panels are described in the following subsections.

### Welcome!

The landing page for the app is also where you might likely be reading this text (otherwise in the package vignette).



### Data Setup


The Data Setup panel is where you can upload or inspect the required inputs for running the app. This builds on the primary idea used by *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* and extends it with the following aspects:

- the panel structure appears dynamically in three consecutive mandatory steps, marked with color from red to yellow to green, with optional steps in light blue.
- the optional step of retrieving the annotation on the fly relieves the user from the task of composing the `data.frame` in advance, and is based on the widely adopted `org.XX.eg.db` Bioconductor packages.
- when the objects are already passed as parameters, or computed, a brief overview/summary for them is displayed
- to tighten the concert operations between similar tools with different scope (as *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* and *[ideal](http://bioconductor.org/packages/ideal)* are), the information flow can move from the data exploration to decisions taken at the moment of testing


A diagnostic mean-dispersion plot is also provided in a collapsible element at the bottom of the panel, shown when the `DESeqDataSet` is generated and the `DESeq` command from the `DESeq2` package has been applied.




### Counts Overview


As in pcaExplorer, interactive tables for the raw, normalized or (r)log-transformed counts are shown in this tab. The user can also generate a sample-to-sample correlation scatter plot with the selected data.

Additionally, `ideal` has an option to include a filter step at the gene level by removing genes with low absolute or averages low values. After this, it might be possible to have to re-run the analysis in step 3 from the Data Setup panel. 

### Extract Results


This tab is an interface for generating the summary tables after testing for DE. It is usually based on the Wald test, as implemented in DESeq2, but when the factor of interest is assuming more than two levels, the user can also perform an ANOVA-like test across the groups with the likelihood ratio test. Options for enabling/disabling automated independent filtering, adding the additional column of unshrunken log2 fold change values (instead of the moderated estimates used by default), as well as using the Independent Hypothesis Weighting (IHW) framework, are provided.


The False Discovery Rate (FDR) can be set from the sidebar panel, and a couple of diagnostic plots, such as the histogram of raw p-values and the distribution of log2fc, are shown below the interactive enhanced version of the table - with clickable elements to link to ENSEMBL database and NCBI website.

### Summary Plots



In this tab an interactive MA plot for the contrast selected in the Extract Results tab is displayed. Clicking on a single gene in the zoomed plot (enabled by brushing in the main plot), it is possible to obtain a boxplot for its expression values, flanked by an overview of information accessed live from the Entrez database. Alternatively, a volcano plot of -log10(p-value) versus log fold change can provide a slightly different perspective. The subset of selected genes are also here presented in static and interactive heatmaps, with the underlying data accessible from the collapsible box element.


### Gene Finder


The functionality in the Gene Finder builds upon the one provided by `pcaExplorer`, and allows to query up to four genes in the same view, which can here be selected from a dropdown input list which supports autocompletion. 

A combined summary table (with both normalized counts and results statistics) is located below an MA plot where the selected genes are marked and annotated on the plot. To avoid repeating this manually, the user can also quickly upload a list of genes as text file (one gene identifier per line), such as members of gene families (e.g. all cytokines, all immunoglobulines, ...) or defined by common function (e.g. all housekeeping genes, or others based on any annotation).


### Functional Analysis


The Functional Analysis tab takes the user from the simple lists of DE genes to insight on the affected biological pathways, with three approaches based on the Gene Ontology (GO) databases. This panel of ideal has a slim interface to 

- `limma::goana` for the quick yet standard implementation
- `topGO`, particularly valuable for pruning terms which are topologically less meaningful than their specific nodes
- `goseq`, which accounts for the specific length bias intrinsic in RNA-Seq assays (longer genes have higher chances of being called DE).

*[ideal](http://bioconductor.org/packages/ideal)* allows the user to work simultaneously with more gene lists, two of which can be uploaded in a custom way (e.g. list of gene families, or extracted from other existing publications). 

The interaction among these lists can be visually represented in Venn diagrams, as well as with the appealing alternative from the UpSetR package, where all combination of sets are explicitly shown. 

Each of the methods for GO enrichment delivers its own interactive `DT`-based table, which can then be explored interactively with the display of a heatmap for all the (DE) genes annotated to a particular term, picking the normalized transformed values for comparing robustly the expression values. This is simply triggered by clicking any of the rows for the results tables. Another useful feature is provided by the clickable link to the AmiGO database on each of the GO term identifiers.


### Report Editor


The Report Editor tab works in the same way of `pcaExplorer`, with the scope of providing an interface to full computational reproducibility of the analyses.

General `Markdown options` and `Editor options` are available, and the text editor, based on the `shinyAce` package, contains a comprehensive template report, that can be edited to the best convenience of the user.

The code contained in the template report fetches the latest state of the reactive values in the ongoing session, and its output is a comprehensive HTML file that can be expanded, edited, previewed in the tab itself, downloaded, and shared with a few mouse clicks.

### About


The About tab contains the output of `sessionInfo`, plus general information on *[ideal](http://bioconductor.org/packages/ideal)*, including the link to the Github development version. If requested, the modular structure of the app can be easily expanded, and many new operations on the same set of input data and derived results can be embedded in the same framework. 




## Running `ideal` on an exemplary data set

We can run *[ideal](http://bioconductor.org/packages/ideal)* for demonstration purpose on published datasets that are available as SummarizedExperiment in an experiment Bioconductor packages.

We will use the *[airway](http://bioconductor.org/packages/airway)* dataset, which can be installed with this command

<!-- [][]or use pasilla[][] -->



```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("airway")
```

This package provides a `RangedSummarizedExperiment` object of read counts in genes for an RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone. More details such as gene models and count quantifications can be found in the *[airway](http://bioconductor.org/packages/airway)* package vignette. 

To run *[ideal](http://bioconductor.org/packages/ideal)* on this dataset, the following commands are required. First, prepare the objects to be passed as parameters of *[ideal](http://bioconductor.org/packages/ideal)*


```r
library(airway)
library(DESeq2)

data(airway)

dds_airway <- DESeqDataSet(airway,design= ~ cell + dex)
dds_airway
```

```
class: DESeqDataSet 
dim: 64102 8 
metadata(2): '' version
assays(1): counts
rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
rowData names(0):
colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
colData names(9): SampleName cell ... Sample BioSample
```

```r
# run deseq on it
dds_airway <- DESeq(dds_airway)
# extract the results
res_airway <- results(dds_airway, contrast = c("dex","trt","untrt"),alpha = 0.05)
```

Then launch the app itself


```r
ideal(dds_obj = dds_airway)
# or also providing the results object
ideal(dds_obj = dds_airway,res_obj = res_airway)
```

The `annotation` for this dataset can be built manually by exploiting the *[org.Hs.eg.db](http://bioconductor.org/packages/org.Hs.eg.db)* package


```r
library(org.Hs.eg.db)
genenames_airway <- mapIds(org.Hs.eg.db,keys = rownames(dds_airway),column = "SYMBOL",keytype="ENSEMBL")
annotation_airway <- data.frame(gene_id = rownames(dds_airway),
                                gene_name = genenames_airway,
                                row.names = rownames(dds_airway),
                                stringsAsFactors = FALSE)
head(annotation_airway)                                
```

```
                        gene_id gene_name
ENSG00000000003 ENSG00000000003    TSPAN6
ENSG00000000005 ENSG00000000005      TNMD
ENSG00000000419 ENSG00000000419      DPM1
ENSG00000000457 ENSG00000000457     SCYL3
ENSG00000000460 ENSG00000000460  C1orf112
ENSG00000000938 ENSG00000000938       FGR
```

or alternatively, can be handily created at runtime in the optional step.

Then again, the app can be launched with 


```r
ideal(dds_obj = dds_airway,
      annotation_obj = annotation_airway)
```

If desired, alternatives can be used. See the well written annotation workflow available at the Bioconductor site (https://bioconductor.org/help/workflows/annotation/annotation/).

<!-- # Running `ideal` on synthetic datasets -->

## Functions exported by the package for standalone usage

The functions exported by the *[ideal](http://bioconductor.org/packages/ideal)* package can be also used in a standalone scenario, provided the required objects are in the working environment. They are listed here for an overview, but please refer to the documentation for additional details. Where possible, for each function a code snippet will be provided for its typical usage.

### `deseqresult2DEgenes` and `deseqresult2tbl`

`deseqresult2DEgenes` and `deseqresult2tbl` generate a tidy table with the results of DESeq2, sorted by the values in the `padj` column.


```r
summary(res_airway)
```

```

out of 33469 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 2206, 6.6% 
LFC < 0 (down)   : 1801, 5.4% 
outliers [1]     : 0, 0% 
low counts [2]   : 16058, 48% 
(mean count < 6)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

```r
res_airway
```

```
log2 fold change (MAP): dex trt vs untrt 
Wald test p-value: dex trt vs untrt 
DataFrame with 64102 rows and 6 columns
                 baseMean log2FoldChange      lfcSE       stat       pvalue
                <numeric>      <numeric>  <numeric>  <numeric>    <numeric>
ENSG00000000003 708.60217    -0.37415246 0.09884435 -3.7852692 0.0001535423
ENSG00000000005   0.00000             NA         NA         NA           NA
ENSG00000000419 520.29790     0.20206175 0.10974241  1.8412367 0.0655868795
ENSG00000000457 237.16304     0.03616686 0.13834540  0.2614244 0.7937652416
ENSG00000000460  57.93263    -0.08445399 0.24990709 -0.3379415 0.7354072415
...                   ...            ...        ...        ...          ...
LRG_94                  0             NA         NA         NA           NA
LRG_96                  0             NA         NA         NA           NA
LRG_97                  0             NA         NA         NA           NA
LRG_98                  0             NA         NA         NA           NA
LRG_99                  0             NA         NA         NA           NA
                       padj
                  <numeric>
ENSG00000000003 0.001245144
ENSG00000000005          NA
ENSG00000000419 0.191760396
ENSG00000000457 0.909231300
ENSG00000000460 0.878782741
...                     ...
LRG_94                   NA
LRG_96                   NA
LRG_97                   NA
LRG_98                   NA
LRG_99                   NA
```

```r
head(deseqresult2tbl(res_airway))
```

```
               id   baseMean log2FoldChange     lfcSE     stat        pvalue
1 ENSG00000152583   997.4398       4.313962 0.1721373 25.06116 1.319237e-138
2 ENSG00000165995   495.0929       3.186823 0.1281565 24.86664 1.708565e-136
3 ENSG00000101347 12703.3871       3.618734 0.1489434 24.29604 2.158637e-130
4 ENSG00000120129  3409.0294       2.871488 0.1182491 24.28338 2.937247e-130
5 ENSG00000189221  2341.7673       3.230395 0.1366745 23.63569 1.656535e-123
6 ENSG00000211445 12285.6151       3.553360 0.1579821 22.49217 4.952260e-112
           padj
1 2.296923e-134
2 1.487391e-132
3 1.252801e-126
4 1.278510e-126
5 5.768386e-120
6 1.437063e-108
```

In particular, `deseqresult2DEgenes` only includes genes detected as DE


```r
head(deseqresult2DEgenes(res_airway,FDR = 0.05))
```

```
               id   baseMean log2FoldChange     lfcSE     stat        pvalue
1 ENSG00000152583   997.4398       4.313962 0.1721373 25.06116 1.319237e-138
2 ENSG00000165995   495.0929       3.186823 0.1281565 24.86664 1.708565e-136
3 ENSG00000101347 12703.3871       3.618734 0.1489434 24.29604 2.158637e-130
4 ENSG00000120129  3409.0294       2.871488 0.1182491 24.28338 2.937247e-130
5 ENSG00000189221  2341.7673       3.230395 0.1366745 23.63569 1.656535e-123
6 ENSG00000211445 12285.6151       3.553360 0.1579821 22.49217 4.952260e-112
           padj
1 2.296923e-134
2 1.487391e-132
3 1.252801e-126
4 1.278510e-126
5 5.768386e-120
6 1.437063e-108
```

```r
# the output in the first lines is the same, but
dim(res_airway)
```

```
[1] 64102     6
```

```r
dim(deseqresult2DEgenes(res_airway))
```

```
[1] 4007    7
```

This tables can be enhanced with clickable links to the ENSEMBL and NCBI gene databases by the following commands


```r
myde <- deseqresult2DEgenes(res_airway,FDR = 0.05)
myde$symbol <- mapIds(org.Hs.eg.db,keys = as.character(myde$id),column = "SYMBOL",keytype="ENSEMBL")
```

```
'select()' returned 1:many mapping between keys and columns
```

```r
myde_enhanced <- myde
myde_enhanced$id <- ideal:::createLinkENS(myde_enhanced$id,species = "Homo_sapiens")
myde_enhanced$symbol <- ideal:::createLinkGeneSymbol(myde_enhanced$symbol)
DT::datatable(myde_enhanced[1:100,], escape = FALSE)
```




### `ggplotCounts`

`ggplotCounts` extends the functionality of the `plotCounts` function of *[DESeq2](http://bioconductor.org/packages/DESeq2)*, and plots the normalized counts of a single gene as a boxplot, with jittered points superimposed.


```r
ggplotCounts(dds = dds_airway,
             gene = "ENSG00000103196", # CRISPLD2 in the original publication
             intgroup = "dex")
```


If an `annotation_obj` is provided, their gene name can also be included in the title.


```r
ggplotCounts(dds = dds_airway,
             gene = "ENSG00000103196", # CRISPLD2 in the original publication
             intgroup = "dex",
             annotation_obj = annotation_airway)
```


When used in the context of the app, it is possible to seamless search for genes also by their gene name, making exploration even more immediate.

### `goseqTable`

`goseqTable` is a wrapper to extract the functional GO terms enriched in  in a list of (DE) genes, based on the algorithm and the implementation in the *[goseq](http://bioconductor.org/packages/goseq)* package.

Its counterpart, based on the *[topGO](http://bioconductor.org/packages/topGO)* package, can be found in the *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* package.

<!-- maybe eval=FALSE here? -->


```r
res_subset <- deseqresult2DEgenes(res_airway)[1:100,]
myde <- res_subset$id
myassayed <- rownames(res_airway)
mygo <- goseqTable(de.genes = myde,
                   assayed.genes = myassayed,
                   genome = "hg38",
                   id = "ensGene",
                   testCats = "GO:BP",
                   addGeneToTerms = FALSE)
```

```
Can't find hg38/ensGene length data in genLenDataBase...
```

```
Found the annotaion package, TxDb.Hsapiens.UCSC.hg38.knownGene
```

```
Trying to get the gene lengths from it.
```

```
Loading required package: GenomicFeatures
```

```

Attaching package: 'GenomicFeatures'
```

```
The following object is masked from 'package:topGO':

    genes
```

```
Fetching GO annotations...
```

```
For 44913 genes, we could not find any categories. These genes will be excluded.
```

```
To force their use, please run with use_genes_without_cat=TRUE (see documentation).
```

```
This was the default behavior for version 1.15.1 and earlier.
```

```
Calculating the p-values...
```

```
'select()' returned 1:1 mapping between keys and columns
```

```r
head(mygo)
```

```
       category over_represented_pvalue under_represented_pvalue numDEInCat
8485 GO:0048856            6.883358e-06                0.9999974         50
5122 GO:0032879            1.861302e-05                0.9999936         29
1123 GO:0003300            2.620204e-05                0.9999987          5
3180 GO:0014897            2.996555e-05                0.9999984          5
3179 GO:0014896            3.519238e-05                0.9999981          5
2663 GO:0010243            3.982549e-05                0.9999909         14
     numInCat                                term ontology      p.adj
8485     5487    anatomical structure development       BP 0.09903272
5122     2488          regulation of localization       BP 0.09903272
1123       69          cardiac muscle hypertrophy       BP 0.09903272
3180       71         striated muscle hypertrophy       BP 0.09903272
3179       73                  muscle hypertrophy       BP 0.09903272
2663      777 response to organonitrogen compound       BP 0.09903272
```

As for the results, this table can be enhanced by adding the links for each category to the AmiGO database

<!-- maybe eval=FALSE here? -->


```r
mygo_enhanced <- mygo
mygo_enhanced$category <- ideal:::createLinkGO(mygo_enhanced$category)
DT::datatable(mygo_enhanced,escape=FALSE)
```



### `plot_ma`

The MA plot provided by *[ideal](http://bioconductor.org/packages/ideal)* displays the gene-wise log2-fold changes (logFCs) versus the average expression value. As a main input parameter, a `DESeqResults` object is required. Control on the appearance of the plot can be applied by selecting the False Discovery Rate (`FDR`), the point transparency (`point_alpha`), adding horizontal lines at particular logFC values (`hlines`). The user can also decide to add rug plots in the margins (setting `add_rug` to `TRUE`).

To facilitate the inspection of a particular gene or gene set, `intgenes` can assume the value of a vector of genes (either the IDs or the gene symbols if `symbol` column is provided in `res_obj`. Labels can be added via `labels_intgenes`, while classical labels/title can be also edited as preferred (see `plot_ma` for all details). 


```r
plot_ma(res_airway, FDR = 0.05, hlines = 1, title ="Adding horizontal lines")
```



```r
plot_ma(res_airway, FDR = 0.1,
        intgenes = c("ENSG00000103196",  # CRISPLD2
                     "ENSG00000120129",  # DUSP1
                     "ENSG00000163884",  # KLF15
                     "ENSG00000179094"), # PER1
        title = "Providing a shortlist of genes"
       )
```


```r
res_airway2 <- res_airway
res_airway2$symbol <- mapIds(org.Hs.eg.db,keys = rownames(res_airway2),column = "SYMBOL",keytype="ENSEMBL")
```

```
'select()' returned 1:many mapping between keys and columns
```

```r
plot_ma(res_airway2, FDR = 0.05,
        intgenes = c("CRISPLD2",  # CRISPLD2
                     "DUSP1",  # DUSP1
                     "KLF15",  # KLF15
                     "PER1"), # PER1
        annotation_obj = annotation_airway,
        hlines = 2,
        add_rug = FALSE,
        title = "Putting gene names..."
       )
```


### `plot_volcano`

The volcano plot gives a different flavor for the gene overview, displaying log2-fold changes and log p-values

In a way similar to `plot_ma`, genes can be annotated with `intgenes`, and vertical lines can be added via `vlines`. `ylim_up` controls the y axis upper limit to visualize better the bulk of genes - keep in mind that very small p-values due to robust differences/large effect sizes can be "cut out".


```r
plot_volcano(res_airway)
```



```r
plot_volcano(res_airway2, FDR = 0.05,
        intgenes = c("CRISPLD2",  # CRISPLD2
                     "DUSP1",  # DUSP1
                     "KLF15",  # KLF15
                     "PER1"), # PER1
        title = "Selecting the handful of genes - and displaying the full range for the -log10(pvalue)",
        ylim_up = -log10(min(res_airway2$pvalue, na.rm =TRUE)))
```



### `sepguesser`   

`sepguesser` makes an educated guess on the separator character for the input text file (`file`). The separator list can be provided as a vector in `sep_list` (defaults to comma, tab, semicolon, and whitespace - which ideally could cover most of the cases). The heuristics is based on the number of occurrencies of each separator in each line.


```r
sepguesser(system.file("extdata/design_commas.txt",package = "ideal"))
```

```
[1] ","
```

```r
sepguesser(system.file("extdata/design_semicolons.txt",package = "ideal"))
```

```
[1] ";"
```

```r
sepguesser(system.file("extdata/design_spaces.txt",package = "ideal"))
```

```
[1] " "
```

```r
anyfile <- system.file("extdata/design_tabs.txt",package = "ideal") # we know it is going to be TAB
guessed_sep <- sepguesser(anyfile)
guessed_sep
```

```
[1] "\t"
```

```r
# to be used for reading in the same file, without having to specify the sep
read.delim(anyfile, header = TRUE, as.is = TRUE, 
           sep = guessed_sep, 
           quote = "", row.names = 1, check.names = FALSE)
```

```
           SampleName    cell   dex albut        Run avgLength Experiment
SRR1039508 GSM1275862  N61311 untrt untrt SRR1039508       126  SRX384345
SRR1039509 GSM1275863  N61311   trt untrt SRR1039509       126  SRX384346
SRR1039512 GSM1275866 N052611 untrt untrt SRR1039512       126  SRX384349
SRR1039513 GSM1275867 N052611   trt untrt SRR1039513        87  SRX384350
SRR1039516 GSM1275870 N080611 untrt untrt SRR1039516       120  SRX384353
SRR1039517 GSM1275871 N080611   trt untrt SRR1039517       126  SRX384354
SRR1039520 GSM1275874 N061011 untrt untrt SRR1039520       101  SRX384357
SRR1039521 GSM1275875 N061011   trt untrt SRR1039521        98  SRX384358
              Sample    BioSample sizeFactor
SRR1039508 SRS508568 SAMN02422669  1.0236476
SRR1039509 SRS508567 SAMN02422675  0.8961667
SRR1039512 SRS508571 SAMN02422678  1.1794861
SRR1039513 SRS508572 SAMN02422670  0.6700538
SRR1039516 SRS508575 SAMN02422682  1.1776714
SRR1039517 SRS508576 SAMN02422673  1.3990365
SRR1039520 SRS508579 SAMN02422683  0.9207787
SRR1039521 SRS508580 SAMN02422677  0.9445141
```



## Creating and sharing output objects

While running the app, the user can

- generate and save graphics
- create and export tables
- generate, preview, download/export an HTML report
- save the values of the `reactiveValues` in an environment, or in binary format

This functionality to retrieve and share the output is provided by action buttons that are placed close to each element of interest. 



## Further development

Additional functionality for the *[ideal](http://bioconductor.org/packages/ideal)* will be added in the future, as it is tightly related to a topic under current development research. 

Improvements, suggestions, bugs, issues and feedback of any type can be sent to marinif@uni-mainz.de.



## Session Info {.unnumbered}


```r
sessionInfo()
```

```
R version 3.3.0 (2016-05-03)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] TxDb.Hsapiens.UCSC.hg38.knownGene_3.4.0
 [2] GenomicFeatures_1.26.3                 
 [3] ideal_0.9.0                            
 [4] BiocStyle_2.3.30                       
 [5] covr_2.2.2                             
 [6] airway_0.108.0                         
 [7] org.Hs.eg.db_3.4.0                     
 [8] BiocParallel_1.8.1                     
 [9] shinyAce_0.2.1                         
[10] rmarkdown_1.3                          
[11] knitr_1.15.13                          
[12] rintrojs_0.1.2                         
[13] rentrez_1.0.4                          
[14] shinyBS_0.61                           
[15] shinydashboard_0.5.3                   
[16] GOstats_2.40.0                         
[17] Category_2.40.0                        
[18] Matrix_1.2-8                           
[19] limma_3.30.11                          
[20] dplyr_0.5.0                            
[21] plyr_1.8.4                             
[22] stringr_1.2.0                          
[23] goseq_1.26.0                           
[24] geneLenDataBase_1.10.0                 
[25] BiasedUrn_1.07                         
[26] UpSetR_1.3.2                           
[27] gplots_3.0.1                           
[28] IHW_1.2.0                              
[29] pcaExplorer_2.0.0                      
[30] pheatmap_1.0.8                         
[31] d3heatmap_0.6.1.1                      
[32] ggplot2_2.2.1                          
[33] DESeq2_1.14.1                          
[34] SummarizedExperiment_1.4.0             
[35] GenomicRanges_1.26.3                   
[36] GenomeInfoDb_1.10.3                    
[37] topGO_2.26.0                           
[38] SparseM_1.74                           
[39] GO.db_3.4.0                            
[40] AnnotationDbi_1.36.2                   
[41] IRanges_2.8.1                          
[42] S4Vectors_0.12.1                       
[43] Biobase_2.34.0                         
[44] graph_1.52.0                           
[45] BiocGenerics_0.20.0                    
[46] debrowser_0.99.0                       
[47] shiny_1.0.0                            
[48] DT_0.2                                 
[49] ggvis_0.4.3                            
[50] clusterProfiler_3.2.0                  
[51] DOSE_3.0.0                             

loaded via a namespace (and not attached):
 [1] backports_1.0.5          Hmisc_4.0-2             
 [3] fastmatch_1.1-0          NMF_0.20.6              
 [5] igraph_1.0.1             lazyeval_0.2.0.9000     
 [7] GSEABase_1.36.0          splines_3.3.0           
 [9] gridBase_0.4-7           lpsymphony_1.2.0        
[11] digest_0.6.12            BiocInstaller_1.24.0    
[13] foreach_1.4.3            htmltools_0.3.5         
[15] GOSemSim_2.0.4           rsconnect_0.7           
[17] gdata_2.17.0             magrittr_1.5            
[19] checkmate_1.8.2          memoise_1.0.0           
[21] cluster_2.0.5            doParallel_1.0.10       
[23] Biostrings_2.42.1        annotate_1.52.1         
[25] matrixStats_0.51.0       colorspace_1.3-2        
[27] ggrepel_0.6.5            RCurl_1.95-4.8          
[29] jsonlite_1.3             genefilter_1.56.0       
[31] survival_2.40-1          iterators_1.0.8         
[33] registry_0.3             gtable_0.2.0            
[35] zlibbioc_1.20.0          XVector_0.14.0          
[37] scales_0.4.1.9000        DBI_0.5-1               
[39] edgeR_3.16.5             rngtools_1.2.4          
[41] Rcpp_0.12.9.2            xtable_1.8-2            
[43] htmlTable_1.9            foreign_0.8-67          
[45] Formula_1.2-1            AnnotationForge_1.16.1  
[47] rex_1.1.1                htmlwidgets_0.8         
[49] threejs_0.2.2            fgsea_1.0.2             
[51] RColorBrewer_1.1-2       acepack_1.4.1           
[53] XML_3.98-1.5             nnet_7.3-12             
[55] locfit_1.5-9.1           labeling_0.3            
[57] reshape2_1.4.2           munsell_0.4.3           
[59] tools_3.3.0              RSQLite_1.1-2           
[61] devtools_1.12.0          fdrtool_1.2.15          
[63] evaluate_0.10            yaml_2.1.14             
[65] caTools_1.17.1           RBGL_1.50.0             
[67] nlme_3.1-131             mime_0.5                
[69] slam_0.1-37              DO.db_2.9               
[71] biomaRt_2.30.0           png_0.1-7               
[73] tibble_1.2               geneplotter_1.52.0      
[75] stringi_1.1.2            lattice_0.20-34         
[77] data.table_1.10.4        bitops_1.0-6            
[79] httpuv_1.3.3             rtracklayer_1.34.2      
[81] qvalue_2.6.0             R6_2.2.0                
[83] latticeExtra_0.6-28      bookdown_0.3.9          
[85] KernSmooth_2.23-15       gridExtra_2.2.1         
[87] codetools_0.2-15         gtools_3.5.0            
[89] assertthat_0.1           pkgmaker_0.22           
[91] rprojroot_1.2            withr_1.0.2             
[93] GenomicAlignments_1.10.0 Rsamtools_1.26.1        
[95] mgcv_1.8-17              grid_3.3.0              
[97] rpart_4.1-10             tidyr_0.6.1             
[99] base64enc_0.1-3         
```

```r
# devtools::session_info()
```


