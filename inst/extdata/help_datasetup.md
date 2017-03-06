## Data Setup

### Upload a dataset

Use the widgets below in Step 1 to upload a text tab/comma/semicolon/space delimited data files. You will need two pieces of information, namely

- the count matrix (as output e.g. from featureCounts or HTSeq-count)
- the experimental design matrix/data frame

In both cases, the first row is assumed to contain headers and the first column contains the feature names (i.e. gene IDs/names). 

You can check how the data look like in the collapsible elements.

### Select the DE design

Select the experimental factors that you want to account for when setting up your comparison of interest.

In the easiest case, this field will be one of the columns you provided in the experimental design matrix

### Optional steps

You can additionally

- create an annotation object - recommended, as it will enhance your experience throughout the app runs
- remove samples - if you think some have to be deemed as outliers (please consider using  [`pcaExplorer`](http://bioconductor.org/packages/pcaExplorer/) for this purpose)

### Run the DESeq2 pipeline

Once everything is set, you can just click on the run button. This operation might take a while depending on the dataset size, but will compute most of the components you will require for the further analyses.

You can also inspect a diagnostic mean-dispersion plot for the current dataset in the collapsible element.

Once you are done with this help box, you can close it by clicking on the `Help` text.
