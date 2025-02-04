# Compact transcription factor cassettes generate functional, engraftable motor neurons by direct conversion
[<img src="https://img.shields.io/badge/code_license-MIT-green">](./LICENSE)
[<img src="https://img.shields.io/badge/text_license-CC--BY--4.0-green">](https://creativecommons.org/licenses/by/4.0/)

| Item                     | DOI                       |
| ------------------------ |:-------------------------:|
| Article                              | [![Dataset DOI](https://img.shields.io/badge/Article_DOI-TBD-blue)](TBD)      |
| Dataset                              | [![Dataset DOI](https://img.shields.io/badge/Dataset_DOI-10.5281/zenodo.14743950-blue)](https://doi.org/10.5281/zenodo.14743950)       |
| scRNAseq data I                       | [![Dataset GSE](https://img.shields.io/badge/Dataset_GSE-GSE287882-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287882)       |
| scRNAseq data II                      | [![Dataset GSE](https://img.shields.io/badge/Dataset_GSE-GSE287783-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287783)       |


## For python analysis

### Setting data directory
We use rushd to organge our data analysis. rushd requires you to specify a data directory. After you decide where to put the data locally, create a file called `datadir.txt` in the GitHub repo that contains the code. The `datadir.txt` file will contain only the **absolute path** to the folder in which you placed the data.

For example, if you placed the GitHub repo in `/Users/username/Documents/GitHub/article-engraftable-neurons/` and you place the data in `/Users/username/Documents/GitHub/article-engraftable-neurons/data/`, then you would place the `datadir.txt` file in `/Users/username/Documents/GitHub/article-engraftable-neurons/` and it should contain one line:

```
/Users/username/Documents/GitHub/article-engraftable-neurons/data/

```

### Instructions

1. Download local copy of data
2. Update datadir.txt to point to data folder
3. Use virtual environment with Python 3.9.6
4. Install packages
```
pip install -r requirements.txt
```

## For scRNA-seq (R) analysis

### Setting up the scRNA-seq data directory

#### Download data 
First download required GSE files (see table at top)

#### Rename data and put into folders
After downloading all required *barcodes.tsv.gz, *features.tsv.gz, and *matrix.mtx.gz, you will have to rename them and place them into subfolders for the `Read10(data.dir=)` function. If you place it into a folder called `scRNAseq-data` it should look like this:

```
scRNAseq-data
├── 6FDDRR-14dpi-1-iMN
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── 6FDDRR-14dpi-2-iMN
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── EmbMN-1
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── EmbMN-2
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── LNI-DDRR-4dpi
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── LNI-DDRR-iMN
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── LNI-DDRR-nosupplements
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── MEFs
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

#### Add datadir.txt file to point to data
In the r directory where the .R files are, add a `datadir.txt` file that points to where the folders with the scRNA-seq data is. For example, if you placed the `scRNAseq-data` folder on your Desktop, you would put a `datadir.txt` file in
`/Users/username/Documents/GitHub/article-engraftable-neurons/r/` that contains the following text

````
/Users/username/Desktop/scRNAseq-data/

```

This allows the .R scripts to change the working directory to the one containing the data.

### Analysis
Run `2023.10.16_embMNfilter.R` first to isolate the median motor column (MMC) motor neurons. It will produce a .Rds file with just the MMC cells.

Then run `2023.10.17_scRNA-analysis.R`
