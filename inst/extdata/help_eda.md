## Performing EDA before DE

For the purpose of Exploratory Data Analysis, fundamental in -omics data like in many low-throughput experimental datasets, we can recommend another package we developed, `pcaExplorer` (https://bioconductor.org/packages/pcaExplorer/, published at https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2879-1).
There you can have PCA plots, heatmaps of distances across samples, and much more (e.g. hints towards a functional interpretation of principal components).

By splitting the steps of Exploratory Data Analysis and Differential Expression analysis, while still keeping a common framework (with similar input objects and formats), we hope to raise the awareness of the importance of proper exploratory steps, before testing for differential expression - and at the same time, not increase the burden for users by having to learn _yet another software interface_ and its requirements. 

You can use `pcaExplorer` both as an interactive Shiny app, as well as with its functionality from the exported functions, inserted e.g. in fully fledged analysis reports.

If you want to try out `pcaExplorer`, you can do so by checking out http://shiny.imbei.uni-mainz.de:3838/pcaExplorer/ - or using it on your local machine (it is imported directed by `ideal`).
