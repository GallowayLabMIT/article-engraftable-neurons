# Compact transcription factor cassettes generate functional, engraftable motor neurons by direct conversion
[<img src="https://img.shields.io/badge/code_license-MIT-green">](./LICENSE)
[<img src="https://img.shields.io/badge/text_license-CC--BY--4.0-green">](https://creativecommons.org/licenses/by/4.0/)

| Item                     | DOI                       |
| ------------------------ |:-------------------------:|
| Article                              | [![Dataset DOI](https://img.shields.io/badge/Article_DOI-TBD-blue)](TBD)      |
| Dataset                              | [![Dataset DOI](https://img.shields.io/badge/Dataset_DOI-10.5281/zenodo.14743950-blue)](https://doi.org/10.5281/zenodo.14743950)       |


## Setting data directory
We use rushd to organge our data analysis. rushd requires you to specify a data directory. After you decide where to put the data locally, create a file called datadir.txt in the GitHub repo that contains the code. The datadir.txt file will contain only the **absolute path** to the folder in which you placed the data.

For example, if you placed the GitHub repo in /Users/username/Documents/GitHub/article-engraftable-neurons/ and you place the data in /Users/username/Documents/GitHub/article-engraftable-neurons/data/, then you would place the datadir.txt file in /Users/username/Documents/GitHub/article-engraftable-neurons/ and it should contain one line:

```
/Users/username/Documents/GitHub/article-engraftable-neurons/data/
```

## Instructions

1. Download local copy of data
2. Update datadir.txt to point to data folder
3. Use virtual environment with Python 3.9.6
4. Install packages
```
pip install -r requirements.txt
```
