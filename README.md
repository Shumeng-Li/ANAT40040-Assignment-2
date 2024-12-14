# ANAT40040 - Assignment 2
# Gene Expression Analysis and Interpretation

## Explanation of different code parts

### 1. Import Required Libraries
Load the necessary R libraries for data processing, visualization, differential expression analysis, pathway enrichment, and survival analysis.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L1-L24

### 2. Set Working Directory, untar the folder and extract the files
Set the working directory and extracts the tar.gz file.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L26-L36

### 3. Read the necessary files
Load the RNA-seq, patient clinical data, and CNA (copy number aberration) data into R.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L38-L45

### 4. Match the RNA-seq patient ids with the CNA ids and the Patient Data ids.
Clean the column names of the RNA-seq data by extracting the first 12 characters and replacing ‘. ‘ to “-” to standardise patient IDs, followed by matching Patient_ID then removes duplicate entries, and ensures only rows with valid gene symbols are kept.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L47-L71

### 5. Create metadata using the CNA level of ERBB2+ 
Extract ERBB2-specific CNA data, convert it into long format, and label patient samples as ‘Amplified’ or ‘Not Amplified’  based on ERBB2 CNA levels (greater than 0 means amplified).
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L73-L101

### 6. Normalize data using DESeq2
Create a DESeq2 object for differential expression analysis and performs normalization.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L103-L113

### 7. Obtain Differentially Expressed Genes
Conduct differential expression analysis (DEA) between ERBB2 groups. Shrinks the log2 fold changes to reduce noise.
Generates MA plots (unshrunken & shrunken) showing the mean expression against log fold change for all genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L115-L132

### 8. Top 10 Differentially Expressed Genes Ranked by Fold Change

https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L134-L190

### 9. Perform Pathway Enrichment Analysis
Perform KEGG pathway enrichment analysis on the up- and down-regulated genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L192-L227

### 10. Principal Component Analysis (PCA) Plot
Perform PCA to reduce the dimensionality of the normalized expression data and visualize sample clustering.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L229-L250

### 11. Heatmap of DE Genes
Generate heatmaps for the top 10 and top 100 differentially expressed genes, respectively, showing their correlation and clustering.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L252-L295

### 12. Survival Analysis 
#### 12.1 Lasso Regularized Cox Regression with DE Genes
Fit a Lasso regularized Cox regression model to predict survival outcomes using differentially expressed genes.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L297-L366

#### 12.2 Kaplan-Meier Survival Curve
Plot the Kaplan-Meier survival curves for CNA level of ERBB2+ to visualize survival differences.
https://github.com/Shumeng-Li/ANAT40040-Assignment-2/blob/a5e5f0c381a443214108025dff38520c266084cb/Assignment2.R#L369-L383

