# sGLMM (Sparse Graph-structrued Lasso Mixed model)

Implementation of sGLMM in this paper:

    ''Ye, W. Liu, X. Wang, H. & Xing, EP. A Sparse Graph-structrued Lasso Mixed model for genetic association with confounding correction''

## Introduction

sGLMM is a state-of-art model that could extract the genetic variation efficiently from a confounded dataset.

#File Structure:

* **models/** -   main method for sGLMM
* **utility/** -  other helper files 
* **runsGLMM.py**  -  main entry point of using sGLMM to work with your own data

## How to Use

### An Example Command

```
python runsGLMM.py --snum=50 -n data/mice.plink
```
this command will run the program and specify the the number of relevant SNPs that the model selects to 50. 

Options:

    -h, --help          show this help message and exit
    -t FILETYPE         choices of input file type (only csv or plink)
    -n FILENAME         name of the input file
    --lambda=LMBD       The weight of the penalizer. If neither lambda or snum
                        is given, cross validation will be run.
    --snum=SNUM         the number of targeted variables the model selects. If
                        neither lambda or snum is given, cross validation will
                        be run.
    --threshold=THRESHOLD
                        The threshold to mask the weak genotype relatedness
    -q                  Run in quiet mode
    -m                  Run without missing genotype imputation

#### Data Support

* sGLMM currently supports CSV and binary PLINK files.
* Extensions to other data format can be easily implemented through `FileReader` in `utility/dataLoadear`. Feel free to contact us for the support of other data format.

## Installation
First, you need to download this repository by using :

```
git clone https://github.com/YeWenting/sGLMM.git
```

Then, before you run the program, you need to install following package in your Python, that is :

- Numpy
- Scipy
- Pysnptool

You can install the dependency by using:

```
pip install -r requirement.txt
```

## Contact
[Wenting Ye](mailto:wenting_ye@bupt.edu.cn)