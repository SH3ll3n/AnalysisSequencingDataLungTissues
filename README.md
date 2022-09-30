# Analyze sequencing data collected from healthy lung tissues and tumoral lung tissues

# What is it
This analyses of collected data is an educational university project 
to understand which effect do air pollutants have on the investigated tissue (lung),
how do healthy and treated tissues respond to pollutant exposure and which pollutant 
agent impacts the tissue the most?

# How to install
For installation you need to use R language and install the libraries on your PC.

# How to use
To use, you need two csv files, I call "count_matrix_G1.csv" and "coldata.csv". Inside is like:

count_matrix_G1.csv contains IDs of the genes:
```ENSG00000223972.5,0,0,0,0,0,0
ENSG00000227232.5,18,124,29,24,98,60
ENSG00000278267.1,0,0,0,0,0,2
ENSG00000243485.5,0,0,0,0,0,0
ENSG00000284332.1,0,0,0,0,0,0
ENSG00000237613.2,0,0,0,0,0,0
```

coldata.csv contains samples, conditions, gender...:
```,"Sample","Condition","Gender","Age","ApoE","Group","Smell_test"
5,"T109","control","male",64,"34","B",5
6,"T110","control","female",72,"33","C",12
```

Then you simply run the code on your laptop.
