# BIBOT

## Overview

BIBOT stand for bibliography bot, it uses natural language processing (NLP) approaches to parse the content of abstracts of large number of publications; NLP is an emerging field of machine learning, which aims at capturing the meaning of sentences and texts written in a natural language (English), such as scientific articles. BIBOT is written in python 2.7 language, and is designed to interact with the Medline database through the NCBI API, in order to retrieve abstracts, articles and meta-data about theses articles (year of publication, author list, journal of publication, language, conflict of interest statement, etc.).

## Dependencies
* Biopython
* unidecode
* nltk
* itertools
* os
* time
* shutil
* datetime
* glob
* getopt
* sys

## Installation
Most of the python module used by BIBOT are native for the 2.7 python version, just type the following lines to install the other depedencies:

Command line:
```Bash
pip install biopython
pip install bioservices
pip install nltk
pip install unidecode
```
Python console:
```python
import nltk
nltk.download('punkt')
nltk.download('averaged_perceptron_tagger')
nltk.download('maxent_ne_chunker')
nltk.download('words')
```

You can then use BIBOT as a standard python script.

## Usage

BIBOT v1 take 2 mandatory parameters and one optionnal parameters

* -a --action (mandatory)
* -r --rterm (mandatory)
* -c --conf (optionnal)

### action

action is the global action to perform, only 2 options available for now!
* run
* debug

run for a classic run of BIBOT, debug to display a few informations about the initialisation of variables.

### rterms
The list of terms you want to use to screen the NCBI database. Each combination of at least two keywords is used to generate a query, which leads to (2^n)-n-1 generated queries for a run, where n is the number of keywords provided to the program. rterms can be a semi-col delimited list of terms or the name of a txt file containing a semi-col delimited list of terms.

### conf
conf is an optionnal parameter, the name of a configuration file. The configuration file should be a semi-col delimited file and should wontain the following lines:

* min year
* authorized languages
* validation keywords

several validation keywords lines cab be add, only articles matching a term in each of the validation keywords list will be selected.

conf file exemple:
```Bash
min year;2015
authorized languages;eng
validation keywords;autoimmunity,SLE,RA
validation keywords;machine learning, big data, artificial intelligence
```

if no configuration file is provide, default values are used.

### Usage exemple

```Bash
python bibot.py -a run -r "machine learning;SjS;big data" -c myconf.csv
```

When the run is complete, bibot store the selected articles in the abstract subolder and the corresponding meta data in the meta subfolder. each articles is designed by it's pmid.
