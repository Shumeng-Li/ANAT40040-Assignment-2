# ANAT40040 - Assignment 2
# Gene Expression Analysis and Interpretation

## Explanation of different code parts

### 1. Import Required Libraries
Load the necessary R libraries for data processing, visualization, differential expression analysis, pathway enrichment, and survival analysis.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L1-L22

### 2. Set Working Directory, untar the folder and extract the files
Set the working directory and extracts the tar.gz file.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L24-L34

### 3. Read the necessary files
Load the RNA-seq, patient clinical data, and CNA (copy number aberration) data into R.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L36-L43

### 4. Match the RNA-seq patient ids with the CNA ids and the Patient Data ids.
Clean the column names of the RNA-seq data by extracting the first 12 characters and replacing ‘. ‘ to “-” to standardise patient IDs, followed by matching Patient_ID then removes duplicate entries, and ensures only rows with valid gene symbols are kept.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L45-L69

### 5. Create metadata using the CNA level of ERBB2+ 
Extract ERBB2-specific CNA data, convert it into long format, and label patient samples as ‘Amplified’ or ‘Not Amplified’  based on ERBB2 CNA levels (greater than 0 means amplified).
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L71-L99

### 6. Normalize data using DESeq2
Create a DESeq2 object for differential expression analysis and performs normalization.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L101-L111

### 7. Obtain Differentially Expressed Genes
Conduct differential expression analysis (DEA) between ERBB2 groups. Shrinks the log2 fold changes to reduce noise.
Generates MA plots (unshrunken & shrunken) showing the mean expression against log fold change for all genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L113-L130

### 8. Top 10 Differentially Expressed Genes Ranked by Fold Change
#### 8.1 Top 10 Differentially Expressed Genes Table
Generate a table listing the top 10 differentially expressed genes ranked by absolute log2 fold change, including statistical information such as base mean, log2 fold change, standard error and adjusted p-value.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L132-L158

#### 8.2 Top 10 Differentially Expressed Genes Boxplot
Generate boxplots of the normalised expression levels (log10 scale) of the first 10 differentially expressed genes in the ERBB2 Amplified and Not Amplified groups.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L160-L188

### 9. Perform Pathway Enrichment Analysis
Perform KEGG pathway enrichment analysis on the up- and down-regulated genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L190-L224

### 10. Principal Component Analysis (PCA) Plot
Perform PCA to reduce the dimensionality of the normalized expression data and visualize sample clustering.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L226-L247

### 11. Heatmap of DE Genes
Generate heatmaps for the top 50 and top 500 differentially expressed genes, respectively, showing their correlation and clustering.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/52bf176e9a79e5d6edae6418765454a9505bfcfa/Assignment2.R#L249-L292

### 12. Survival Analysis 
#### 12.1 Lasso Regularized Cox Regression with DE Genes
Fit a Lasso regularized Cox regression model to predict survival outcomes using differentially expressed genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L294-L371

#### 12.2 Kaplan-Meier Survival Curve
Plot the Kaplan-Meier survival curves for CNA level of ERBB2+ to visualize survival differences.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/4271a73fa1fcab774d82d8de2c6e9060ee39e113/Assignment2.R#L373-L393

