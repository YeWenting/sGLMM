# sGLMM (Sparse Graph-structrued Lasso Mixed model)

Implementation of sGLMM in this paper:

    ''Ye, W. Liu, X. Wang, H. & Xing, EP. A Sparse Graph-structrued Lasso Mixed model for genetic association with confounding correction''

## Introduction

sGLMM is a state-of-art model that could extract the genetic variation efficiently from a confounded dataset.

#ff# File Structure:

* models/ -  main method for LRVA
* utility/ other helper files
* lrva.py main entry point of using LRVA to work with your own datafff

## An Example Command:

```
python lrva.py -n data/mice.plink
```
#### Data Support
* LRVA currently supports CSV and binary PLINK files.
* Extensions to other data format can be easily implemented through `FileReader` in `utility/dataLoadear`. Feel free to contact us for the support of other data format.

## Installation
You will need to have numpy, scipy and pysnptool installed on your current system.
You can install LRVA using pip by doing the following

```
   pip install git+https://github.com/HaohanWang/LRVA
```
This should make the module pylmm available as well as the two scripts pylmmGWAS.py and pylmmKinship.py.

You can also clone the repository and do a manual install.
```
   git clone https://github.com/HaohanWang/LRVA
   python setup.py install
```

## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
