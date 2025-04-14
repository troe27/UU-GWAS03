# Genome-wide association study
## Table of content

- [Introduction](#introduction)
- [Task 1: A simple GWAS using R](#task-1--a-simple-gwas-using-r)
    - [Reading in the data](#reading-in-the-data)
    - [Test run on the first locus](#test-run-on-the-first-locus)
        - [Step 1](#step-1)
        - [Step 2](#step-2)
        - [Step 3: Chi-square test](#step-3-chi-square-test)
- [Task 2: Whole-genome association study](#task-2-whole-genome-association-study)
    - [Step 1: Create a data frame with 4 columns](#step-1-create-a-data-frame-with-4-columns)
    - [Step 2: For-loop](#step-2-for-loop)
    - [Step 3: Manhattan plot](#step-3-manhattan-plot)
    - [Step 4: Q-Q plot](#step-4-q-q-plot)
- [Task 3: Simple GWAS using PLINK](#task-3-simple-gwas-using-plink)
    - [Step 1: Run association study by PLINK in the terminal](#step-1-run-association-study-by-plink-in-the-terminal)
    - [Step 2: Loading the results](#step-2-loading-the-results)
    - [Step 3: Manhattan plot](#step-3-manhattan-plot-1)
    - [Step 4: Q-Q plot](#step-4-q-q-plot-1)
- [Task 4: [Advanced] GWAS linear mixed model](#task-4-advanced-gwas-linear-mixed-model)
    - [Performing GWAS evaluation](#performing-gwas-evaluation)
    - [Making the Manhattan plot](#making-the-manhattan-plot)

## Introduction
Welcome to the lab practice. Today, we will work on the human data that passed our quality control yesterday. There will be three sections in today’s practice:

- A simple association study using R
- A simple association study using PLINK
- Optional: GWAS linear mixed model using R  

The data set we will be using today can be found [here]().

## Task 1 : A simple GWAS using R
Here’s an example of conducting a chi-square test in R. You will be conducting a GWAS test utilizing the chi-square method.

```R
### Your code here ###
```


#### Reading in the data
Let’s start the GWAS test. First, read the three data-files with the ```human_all``` prefix into ```R``` 
using the ```read_plink()```function in genio package.

```R
### Your code here ###
```


#### test run on the first locus

##### Step 1:
Generate a data frame with only two columns. Remove samples with missing phenotypes or genotypes.

- **Phenotype:** extract phenotype information from the fam data frame (1, 2 indicate this is a case/control study).

- **Genotype:** extract genotype information from the X matrix (the first row). Remove samples which missing genotype information.


```R
### Your code here ###
```


##### Step 2:
use the function ```table()``` to convert the data frame into a contingency table . (HINT: ```?table```)

_**Note:** the ```table()``` function will automatically remove rows with missing values._

```R
### Your code here ###
```

#### Step 3. Chi-square test

```R
### Your code here ###
```
  
## Task 2: Whole-genome association study  
Please repeat the same steps you took with the first marker for all the other markers.

#### Step 1: Create a data frame with 4 columns.
- **CHR**: Chromosome number  
- **SNP**: SNP marker ID  
- **BP**: Marker position  
- **P**: p-value (NA for all cells)  


```R
### Your code here ###
```

#### Step 2: For-loop
 A for loop can be utilized to conduct a chi-square test on all SNP markers. The p-values can be extracted from each iteration and added to the data frame created in Step 1. It’s essential to remember that R may provide a warning message if the expected value is small; you may ignore them.

```R
### Your code here ###
```


#### Step 3: Manhattan plot
 To create a Manhattan plot, utilize the manhattan() function within the qqman package. The data frame will already be in the default input format for this function if you use the same column names as shown in Step 1. Additionally, use the annotatePval = 0.01 option to annotate peaks with a p-value < 10-4.

Check [this webpage](https://r-graph-gallery.com/101_Manhattan_plot.html) for more information.

```R
### Your code here ###
```


#### Step 4. Q-Q plot.

 The [QQ plot](https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot) serves as a crucial tool for identifying issues in a GWAS. To use it, input the vector of p-values from your association result into the ```qq()``` function.

```R
### Your code here ###
```


## Task 3: simple GWAS using PLINK
There are more efficient ways to conduct association tests than using R. For instance, one can use well-developed software like PLINK. For more information, check the [PLINK website](https://www.cog-genomics.org/plink/1.9/assoc). 

#### Step 1: Run association study by PLINK in the terminal.

```bash
### Your code here ###
```

#### Step 2: loading the results
 Read association test results (```*.assoc```) table by into R with the ```read.table()``` function.

```R
### Your code here ###
```



#### Step 3. Manhattan plot.
 The result table generated from PLINK is the exact format for the plotting function.

```R
### Your code here ###
```


#### Step 4. Q-Q plot

```R
### Your code here ###
```



You may notice a difference in results from R and PLINK. To prevent any confusion, please reference the help manual for more information on the chi-square test (```chisq.test()```) and [Yate’s correction for continuity](https://en.wikipedia.org/wiki/Yates%27s_correction_for_continuity).

_Based on the Q-Q plot outcome, which result is superior?_

## Task 4: [Advanced] GWAS linear mixed model
You can choose from different R packages for fitting GWAS linear mixed model. Some of these packages are ```GWAStools```, ```rrBLUP```, and ```sommer```. Among these choices, my personal preference is the ```sommer``` package.

We will use a [rice dataset from Cornell University in 2011](https://www.nature.com/articles/ncomms1467). Raw data can be found on the [rice diversity website](http://www.ricediversity.org/data/sets/44kgwas/).

```R
### Your code here ###
```


Read PLINK binary format files (rice.bed, rice.bim, rice.fam). This dataset contains various trait measurements. You can select one from the file rice_phenotype.txt and incorporate phenotype information into the fam data frame.

Read binary files by the ```genio()``` function.

```R
### Your code here ###
```

Read the phenotype table (I will choose plant height in the answer version)

```R
### Your code here ###
```

Merge ```data$fam``` and ```peno```. The ```NSFTVID``` in the phenotype table is its individual ID.

```R
### Your code here ###
```

Start preparing the input data for the ```GWAS()``` function in the ```sommer``` package.

- The genotype matrix needs to be transposed (row: samples, column: markers).
- Genotype coding should be converted from (0, 1, 2) to (-1, 0, 1).
- Missing values in genotype coded as 0. (Today, we’re keeping it simple, but there are more effective methods for imputing missing values.)
- The row names of the genotype matrix should be identical to the ID in the phenotype table.
Make sure the sample id is factorized.

```R
### Your code here ###
```

Now we need to generate an additive relationship matrix using the ```A.mat()``` function.

```R
### Your code here ###
```

Then, we can fit the linear mixed model using the ```GWAS()``` function.

```R
### Your code here ###
```

#### Performing GWAS evaluation
Create a dataframe for plotting the function. You should have ```CHR```, ```SNP```, ```BP```, and ```P``` in four columns. The p-value could be calculated by ```10^(-as.numeric(fit$scores))).```

```R
### Your code here ###
```

#### Making the Manhattan plot.
To create a Manhattan plot, please keep in mind that if the p-value is 0, its log value becomes infinity, which cannot be plotted. Therefore, kindly remove the data point before proceeding with the plot.

```R
### Your code here ###
```
