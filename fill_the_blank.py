import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import shutil
import datetime
import sys
import glob
import re
import nltk
import operator
import numpy as np
import collections
from Bio import Entrez
from Bio.Entrez import efetch, read


##-----------##
## FUNCTIONS ##########################################################################################
##-----------##
##
## => Manifest :
##    - describe_articles_type(abstract_folder)
##    - get_date_from_meta_save(meta_file)
##    - plot_pulbications_years(meta_data_folder):
##    - plot_word_evolution(item_list, run_folder):
##    - get_all_subjects_from_abstract(abstract_folder):
##    - get_number_of_articles_from_log_file(log_file):
##    - get_number_of_keywords_from_log_file(log_file):
##    - get_articles_filter_stat(log_file):
##    - get_number_of_articles_for_year(run_folder, year):
##    - get_country_publication_stat(run_folder):
##    - get_article_domain(abstract):
##    - retrieve_techniques(display_details):
##    - generate_techniques_figures(category_to_techniques_to_count):
##    - get_cohorte_size(abstract_file):
##    - get_cohorte_size_for_category():
##    - generate_cohorte_size_figures(run_folder):
##
##

def describe_articles_type(abstract_folder):
	##
	## 
	## Count the number of articles talking about
	## several diseases, return a dictionnary.
	## currently work on:
	##		- SjS
	##		- SLE
	##		- RA
	##

	print "[INFO] => Running"

	abstract_to_disease = {}
	abstract_to_disease["SjS"] = 0
	abstract_to_disease["SLE"] = 0
	abstract_to_disease["RA"] = 0 
	abstract_to_disease["Other"] = 0
	abstract_list = glob.glob(str(abstract_folder)+"/*.txt")
	for abstract_file in abstract_list:

		talking_about_SjS = False
		talking_about_RA = False
		talking_about_SLE = False

		abstract = open(abstract_file, "r")

		for line in abstract:

			## SjS
			m = re.findall(r"([S,s]j.gren)", line)
			if(len(m)>0):
				talking_about_SjS = True
			m = re.findall(r"(SjS)", line)
			if(len(m)>0):
				talking_about_SjS = True

			## SLE
			m = re.findall(r"([L,l]upus)", line)
			if(len(m)>0):
				talking_about_SLE = True
			m = re.findall(r"(SLE)", line)
			if(len(m)>0):
				talking_about_SLE = True

			## RA
			m = re.findall(r"([A,a])rthrisis", line)
			if(len(m)>0):
				talking_about_RA = True
			m = re.findall(r"(RA)", line)
			if(len(m)>0):
				talking_about_RA = True

		if(talking_about_SjS):
			abstract_to_disease["SjS"] += 1

			print "[SJS] => "+str(abstract_file)


		if(talking_about_SLE):
			abstract_to_disease["SLE"] += 1


		if(talking_about_RA):
			abstract_to_disease["RA"] += 1
		else:
			abstract_to_disease["Other"] += 1			
		abstract.close()

	return abstract_to_disease



def get_date_from_meta_save(meta_file):
	##
	## Get the date of an article using the
	## meta data file created on local device,
	## no connection needed to NCBI server
	##
	## -> return the year of publication
	##

	## Retrieve the year of publication
	## from the meta data file.
	year = "NA"
	meta_data = open(meta_file, "r")
	for line in meta_data:
		
		line = line.replace("\n", "")
		if(line[0] == ">"):

			line_in_array = line.split(";")
			if(line_in_array[0] == ">Date"):
				date_in_array = line_in_array[1].split("/")
				year = date_in_array[2]

	meta_data.close()

	## return only the year of publication
	return year



def plot_pulbications_years(meta_data_folder):
	##
	## Retrieve the year of publications of all
	## articles from the meta_data_folder and
	## plot the histogramm of publications over
	## the years
	##

	## create the structure
	year_to_count = {}
	for meta_file in glob.glob(meta_data_folder+"/*.csv"):
		year = get_date_from_meta_save(meta_file)
		
		if(int(year) < 2018):

			if(year not in year_to_count.keys()):
				year_to_count[year] = 1
			else:
				year_to_count[year] += 1

	
	## add for publi, to remove
	for key in year_to_count.keys():
		print "[DATA] => "+str(key)+ " : " +str(year_to_count[key])

	## plot graphe
	plt.bar(year_to_count.keys(), year_to_count.values(), color='b', align='center', width=0.3)
	plt.savefig("images/years_publications_evolution.png")
	plt.close()




def plot_word_evolution(item_list, run_folder):
	##
	## -> Plot word frequency evolution over
	## the last decade
	##

	for item in item_list:

		## get year to pmid
		year_to_pmid = {}
		for meta_file in glob.glob(run_folder+"/meta/*.csv"):
			pmid = meta_file.split("/")
			pmid = pmid[-1]
			pmid = pmid.split(".")
			pmid = pmid[0]
			year = get_date_from_meta_save(meta_file)
			if(year not in year_to_pmid.keys()):
				year_to_pmid[year] = []
				year_to_pmid[year].append(pmid)
			else:
				year_to_pmid[year].append(pmid)

		## get apparition count in articles
		year_to_count = {}
		for year in year_to_pmid.keys():
			list_of_pmis_to_check = year_to_pmid[year]
			year_to_count[year] = 0
			for pmid in list_of_pmis_to_check:
				abstract_file = run_folder+"/"+"abstract/"+str(pmid)+"_abstract.txt"
				try:
					abstract_data = open(abstract_file, "r")
					for line in abstract_data:	
						m = re.findall(r"("+str(item)+")", line)		
						if(m is not None):
							if(len(m) > 0):
								year_to_count[year] += 1
					abstract_data.close()
				except:
					## do nothing
					tardis = 1

		## get apparition frequencies
		year_to_frequency = {}
		for year in year_to_pmid.keys():
			frequency = float(year_to_count[year])/float(len(year_to_pmid[year]))
			year_to_frequency[year] = frequency

		## save the figure
		y_vector = []
		x_vector = sorted(year_to_frequency.keys())

		for year in x_vector:
			y_vector.append(year_to_frequency[year])
		fig = plt.figure()

		plt.plot(x_vector, y_vector)
		plt.savefig("images/"+str(item)+"_frequency.png")
		plt.close()

		## get apparition count
		y_vector = []
		x_vector = sorted(year_to_count.keys())
		for year in x_vector:
			y_vector.append(year_to_count[year])
		fig = plt.figure()

		plt.plot(x_vector, y_vector)
		plt.savefig("images/"+str(item)+"_count.png")
		plt.close()






def get_all_subjects_from_abstract(abstract_folder):
	##
	## The idea is to find all key word from
	## each abstract in a run folder and
	## return a list of key word.
	##
	## Word are selected by the number of apparition in
	## all abstracts, if found in more than 10 % of the articles,
	## the word is selected
	##

	## get list of all abstract files
	abstract_list = glob.glob(abstract_folder+"/*.txt")
	name_list = []
	name_to_count = {}

	treshold = 0.1

	abstract_parsed = 0
	for abstract in abstract_list:
		names_found_in_abstract = []
		abstract_data = open(abstract, "r")
		for line in abstract_data:
			try:
				tokens = nltk.word_tokenize(line.encode('utf8'))
				tagged = nltk.pos_tag(tokens)
				entities = nltk.chunk.ne_chunk(tagged)
				abstract_parsed += 1
			except:
				## Something went wrong
				entities = []

			for item in entities:
				try:
					if(item[1] in ["NN", "NNS", "NNP"]):
						if(item[0] not in names_found_in_abstract):
							names_found_in_abstract.append(item[0])
				except:
					## Somethig went wrong
					choucroute = True

		for name in names_found_in_abstract:
			if(name not in name_to_count.keys()):
				name_to_count[name] = 1
			else:
				name_to_count[name] += 1
		abstract_data.close()

	for name in name_to_count.keys():
		if(float(name_to_count[name]/float(abstract_parsed)) > treshold and name not in name_list):
			name_list.append(name)	

	return name_list


def get_number_of_articles_from_log_file(log_file):
	##
	## Get the total number of articles found
	## by the run by looking to the run's log
	## file.
	##
	## return an int or "NA" if nothing is found. 
	##

	number_of_articles = "NA"
	log_data = open(log_file, "r")
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(line_in_array[0] == "Total_number_of_articles"):
			try:
				number_of_articles = int(line_in_array[1])
			except:
				number_of_articles = "NA"

	log_data.close

	return number_of_articles


def get_number_of_keywords_from_log_file(log_file):
	##
	## Get the total number of keywords used to
	## generate the combination of request by
	## the run.
	##
	## return an int or "NA" if nothing is found. 
	##

	number_of_keywords = "NA"
	log_data = open(log_file, "r")
	cmpt = 0
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(cmpt == 0):
			try:
				number_of_keywords = len(line_in_array)
			except:
				number_of_keywords = "NA"
		else:
			break

		cmpt += 1

	log_data.close()

	return number_of_keywords


def get_articles_filter_stat(log_file):
	##
	## Get the count of articles that passed the first
	## filter, the second filter, and both.
	## return a tuple of int
	##

	first_filter_pass_cmpt = 0
	second_filter_pass_cmpt = 0
	pass_both_filter_cmpt = 0

	log_data = open(log_file, "r")
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(line[0] == ">"):
			
			filter_1_info = line_in_array[1].split("=")
			if(filter_1_info[1] == "PASSED"):
				first_filter_pass_cmpt += 1

			filter_2_info = line_in_array[2].split("=")
			if(filter_2_info[1] == "PASSED"):
				second_filter_pass_cmpt += 1


			if(filter_1_info[1] == "PASSED" and filter_2_info[1] == "PASSED"):
				pass_both_filter_cmpt += 1

	log_data.close()

	return (first_filter_pass_cmpt, second_filter_pass_cmpt, pass_both_filter_cmpt)



def get_number_of_articles_for_year(run_folder, year):
	##
	## get the number of articles published in a
	## specific year, parse data from meta files
	## in the given run folder
	##

	article_cmpt = 0

	for meta_file in glob.glob(run_folder+"/meta/*.csv"):
		meta_data = open(meta_file, "r")
		for line in meta_data:
			line = line.replace("\n", "")
			line_in_array = line.split(";")
			if(line_in_array[0] == ">Date"):
				date = line_in_array[1]
				date_in_array = date.split("/")
				year_fetched = date_in_array[-1]
				if(int(year_fetched) == year):
					article_cmpt += 1
		meta_data.close()


	## add for publi, to remove
	print "[DATA] => "+str(year)+" : "+str(article_cmpt)

	return article_cmpt


def get_country_publication_stat(run_folder):
	##
	## Get the publications stats for country.
	## get the list of pmid retrieved from the
	## meta folder and connect to the NCBI to fecth
	## publications informations, parse it to get the
	## country of publication.
	## 
	## return a dictionnary
	##

	## init structure
	country_to_count = {}

	## get list of PMID to process
	meta_file_list = glob.glob(run_folder+"/meta/*.csv")
	for meta_file in meta_file_list:
		meta_file_in_array = meta_file.split("/")
		file_name = meta_file_in_array[-1]
		file_name_in_array = file_name.split(".")
		pmid = file_name_in_array[0]

		## get country publication
		try:
			handle = efetch(db='pubmed', id=pmid, retmode='xml', )
			informations = read(handle)
			stuff = informations[u'PubmedArticle'][0]
			country = stuff[u'MedlineCitation'][u'MedlineJournalInfo'][u'Country']
			print country # to delete
		except:
			country = "NA"

		## fill dictionnary
		if(country not in country_to_count.keys()):
			country_to_count[country] = 1
		else:
			country_to_count[country] += 1

	return country_to_count




def get_article_domain(abstract):
	##
	## IN PROGRESS
	##
	## Fit an article on a pre existing
	## class based on the content of its abstract
	##
	##
	## Class:
	##
	##	- Diagnostic : predict the diagnostic of a patient
	##
	##	- Therapeutic : test a treatment
	##
	##	- Modelisation : try to understand the disease
	##
	## 	- Unclassified : not sure what we are talking about
	##

	## Control structure
	diagnostic_keywords = ["classification", "Classification", "criteria", "diagnostic", "diagnosis", "prevalence", "epidemiological"]
	therapeutic_keywords = ["therapeutic", "therapy", "treatment", "treatments", "rituximab"]
	modelisation_keywords = ["model", "models", "modelisation", "components", "dynamics", "composition", "pathway", "regulatory", "regulates", "mechanistic", "mechanism", "mechanisms"]

	diagnostic_score = 0
	therapeutic_score = 0
	modelisation_score = 0

	article_theme = "Unclassified"

	## Looking for keyword in the abstract with nltk
	names_found_in_abstract = []
	abstract_data = open(abstract, "r")
	for line in abstract_data:
		try:
			tokens = nltk.word_tokenize(line.encode('utf8'))
			tagged = nltk.pos_tag(tokens)
			entities = nltk.chunk.ne_chunk(tagged)
		except:
			## Something went wrong
			entities = []
		for item in entities:
			try:
				if(item[1] in ["NN", "NNS", "NNP"]):
					if(item[0] not in names_found_in_abstract):
						names_found_in_abstract.append(item[0])
			except:
				## Somethig went wrong
				choucroute = True

	## compute score
	for name in names_found_in_abstract:
		if(name in diagnostic_keywords):
			diagnostic_score += 1
		elif(name in therapeutic_keywords):
			therapeutic_score += 1
		elif(name in modelisation_keywords):
			modelisation_score += 1

	abstract_data.close()
	
	## looking for regular expression
	## split abstract into sentences, look for combination
	## of keywords in each sentences.
	## Get the list of words
	abstract_data = open(abstract, "r")
	text = ""
	for line in abstract_data:
		text +=line
	abstract_data.close()
	sentences = text.split(". ")
	words = text.replace(".", "")
	words = words.replace(",", "")
	words = words.replace(":", "")
	words = words.replace(";", "")
	words = words.split(" ")


	## Looking for co-occurences
	for sentence in sentences:
		
		## modelisation 
		if(("mechanisms" in sentence or "mechanism" in sentence) and ("understood" in sentence or "understand" in sentence)):
			modelisation_score += 1
		if(("summarize" in sentence) and ("mechanism" in sentence or "mechanisms" in sentence or "pathway" in sentence or "pathways" in sentence)):
			modelisation_score += 1
		if(("present" in sentence and "study" in sentence) and "investigated" in sentence):
			modelisation_score += 1
		if("investigated" in sentence and "interactions" in sentence):
			modelisation_score += 1

		## Diagnostic
		if("identification" in sentence and ("phenotype" in sentence or "phenotypes" in sentence or "subphenotype" in sentence or "subphenotypes" in sentence)):
			diagnostic_score += 1
		if("severity" in sentence and "marker" in sentence):
			diagnostic_score += 1
		if("differentiate" in sentence and "patients" in sentence):
			diagnostic_score += 1
		if(("biomarker" in sentence or "biomarkers" in sentence) and ("disease" in sentence or "diseases" in sentence)):
			diagnostic_score += 1

		## Therapeutic
		if("impact" in sentence and "course of disease" in sentence):
			therapeutic_score += 1

	## Checking words
	for word in words:
		if(word in diagnostic_keywords):
			diagnostic_score += 1
		elif(word in therapeutic_keywords):
			therapeutic_score += 1
		elif(word in modelisation_keywords):
			modelisation_score += 1

	## compute score
	scores = {"Diagnostic":diagnostic_score, "Therapeutic":therapeutic_score, "Modelisation":modelisation_score}
	sorted_scores = sorted(scores.items(), key=operator.itemgetter(1))

	if(sorted_scores[-1][1] > 0):
		article_theme = sorted_scores[-1][0]

	return article_theme



def retrieve_techniques(run_folder, display_details):
	##
	## Retrieve techniques used in articles from parsing the abstracts
	## run_folder is the name of the run_folder
	## display_details is a boolean:
	##    - True: display details
	##    - False: do nothing 
	##
	##
	## Read abstracts in the articles_classification subfolders
	##    - Diagnostic
	##    - Therapeutic
	##    - Modelisation
	##    - Unclassified
	## 
	## Parse each abstract and look for a list of statistical analaysis appraoches:
	##    - t-test
	##    - logistic regression
	##    - cox regression
	##    - multivariate regression
	##    - univariate regression
	##    - linear regression
	##    - Poisson regression
	##    - regression tree
	##    - regression analysis
	##    - clustering
	##    - hierarchical clustering
	##    - decision tree
	##    - random forest
	##    - artificial neural network
	##    - support vector machine
	##    - machine learning
	##    - Wilcoxon test
	##    - PCA
	##    - chi-squared
	##    - multivariate analysis
	##    - univariate analysis
	##    - discriminant analysis
	##    - bioinformatics analysis
	##    - GLM (generalized linear model)
	##    - mixed models
	##    - normalization
	##    - k mean
	##
	## Return a tuple of dictionnary ({category : technique : count}, {year: techniques : count}
	##


	## initiate data structure
	category_to_techniques = {"Diagnostic":{}, "Therapeutic":{}, "Modelisation":{}, "Unclassified":{}}
	date_to_techniques_to_count = {}
	## Parse the abstracts
	for category in ["Diagnostic", "Therapeutic", "Modelisation", "Unclassified"]:
		for file in glob.glob("articles_classification/"+str(category)+"/*"):

			## get pmid
			pmid = file.split("/")
			pmid = pmid[-1]
			pmid = pmid.split("_")
			pmid = pmid[0]

			## get date
			date = get_date_from_meta_save(str(run_folder)+"/meta/"+str(pmid)+".csv")

			## initialise dict
			if(date not in date_to_techniques_to_count.keys()):
				date_to_techniques_to_count[date] = {}

			abstract = open(file, "r")
			abstract_in_line = ""
			for line in abstract:
				abstract_in_line += line
			abstract.close()

			abstract_in_array = abstract_in_line.split(" ")

			cmpt = 0
			for word in abstract_in_array:

				## t-test
				if(word == "t-test"):
					if("t-test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["t-test"] = 1
					else:
						category_to_techniques[str(category)]["t-test"] += 1

					## get information by year
					if("t-test" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["t-test"] = 1
					else:
						date_to_techniques_to_count[date]["t-test"] += 1

				elif(word == "t" and cmpt + 1 < len(abstract_in_array)):
					if(abstract_in_array[cmpt+1] in ["test", "tests"]):
						if("t-test" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["t-test"] = 1
						else:
							category_to_techniques[str(category)]["t-test"] += 1

						## get information by year
						if("t-test" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["t-test"] = 1
						else:
							date_to_techniques_to_count[date]["t-test"] += 1


				## regression
				##    - logistic
				##    - cox
				##    - multivariate
				##    - univariate
				##    - linear
				##    - tree
				##    - Poisson
				if(word in ["regression", "Regression"] and cmpt-1 >= 0):

					## logistic
					if(abstract_in_array[cmpt-1] in ["logistic", "Logistic", "logistical", "Logistical"]):
						if("logistic regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["logistic regression"] = 1
						else:
							category_to_techniques[str(category)]["logistic regression"] += 1

						## get information by year
						if("logistic regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["logistic regression"] = 1
						else:
							date_to_techniques_to_count[date]["logistic regression"] += 1

					## multivariate
					elif(abstract_in_array[cmpt-1] in ["multivariate", "Multivariate", "multiple", "Multiple"]):
						if("multivariate regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["multivariate regression"] = 1
						else:
							category_to_techniques[str(category)]["multivariate regression"] += 1

						## get information by year
						if("multivariate regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["multivariate regression"] = 1
						else:
							date_to_techniques_to_count[date]["multivariate regression"] += 1
					
					## univariate
					elif(abstract_in_array[cmpt-1] in ["simple", "Simple", "univariate", "Univariate"]):
						if("univariate regression" not in category_to_techniques[str(category)]):
							category_to_techniques[str(category)]["univariate regression"] = 1
						else:
							category_to_techniques[str(category)]["univariate regression"] += 1

						## get information by year
						if("univariate regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["univariate regression"] = 1
						else:
							date_to_techniques_to_count[date]["univariate regression"] += 1

					## linear
					elif(abstract_in_array[cmpt-1] in ["Linear", "linear"]):
						if("linear regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["linear regression"] = 1
						else:
							category_to_techniques[str(category)]["linear regression"] += 1

						## get information by year
						if("linear regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["linear regression"] = 1
						else:
							date_to_techniques_to_count[date]["linear regression"] += 1

					## cox
					elif(abstract_in_array[cmpt-1] in ["Cox", "cox", "Cox's", "cox's"]):
						if("Cox regression" not in category_to_techniques["Diagnostic"].keys()):
							category_to_techniques[str(category)]["Cox regression"] = 1
						else:
							category_to_techniques[str(category)]["Cox regression"] += 1

						## get information by year
						if("Cox regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["Cox regression"] = 1
						else:
							date_to_techniques_to_count[date]["Cox regression"] += 1

					## proportion hazard regression model
					elif(abstract_in_array[cmpt-1] in ["hazard", "hazards"]):
						if(abstract_in_array[cmpt-2] in ["proportion", "proportional"] and abstract_in_array[cmpt+1] in ["analysis", "model"]):
							if("proportion hazard regression model" not in category_to_techniques[str(category)].keys()):
								category_to_techniques[str(category)]["proportion hazard regression model"] = 1
							else:
								category_to_techniques[str(category)]["proportion hazard regression model"] += 1
					
							## get information by year
							if("proportion hazard regression model" not in date_to_techniques_to_count[date].keys()):
								date_to_techniques_to_count[date]["proportion hazard regression model"] = 1
							else:
								date_to_techniques_to_count[date]["proportion hazard regression model"] += 1

					## regression tree
					elif(abstract_in_array[cmpt+1] == "tree"):
						if("regression tree" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["regression tree"] = 1
						else:
							category_to_techniques[str(category)]["regression tree"] += 1

						## get information by year
						if("regression tree" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["regression tree"] = 1
						else:
							date_to_techniques_to_count[date]["regression tree"] += 1

					
					## Poisson regression
					elif(abstract_in_array[cmpt-1] in ["Poisson", "poisson"]):
						if("Poisson regression" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["Poisson regression"] = 1
						else:
							category_to_techniques[str(category)]["Poisson regression"] += 1

						## get information by year
						if("Poisson regression" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["Poisson regression"] = 1
						else:
							date_to_techniques_to_count[date]["Poisson regression"] += 1

					elif(abstract_in_array[cmpt+1] in ["analysis", "test", "tests"]):
						if("regression analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["regression analysis"] = 1
						else:
							category_to_techniques[str(category)]["regression analysis"] += 1

						## get information by year
						if("regression analysis" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["regression analysis"] = 1
						else:
							date_to_techniques_to_count[date]["regression analysis"] += 1

				## Cox's hazard model
				elif(word in ["Cox's", "cox's", "cox", "Cox"] and abstract_in_array[cmpt+1] in ["proportional", "hazards", "hazard", "model", "models"]):
					if("proportion hazard regression model" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["proportion hazard regression model"] = 1
					else:
						category_to_techniques[str(category)]["proportion hazard regression model"] += 1

					## get information by year
					if("proportion hazard regression model" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["proportion hazard regression model"] = 1
					else:
						date_to_techniques_to_count[date]["proportion hazard regression model"] += 1

				## clustering
				elif(word in ["clustering", "Clustering"]):
					if(abstract_in_array[cmpt-1] in ["Hierarchical", "hierarchical"]):
						if("hierarchical clustering" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["hierarchical clustering"] = 1
						else:
							category_to_techniques[str(category)]["hierarchical clustering"] += 1

						## get information by year
						if("hierarchical clustering" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["hierarchical clustering"] = 1
						else:
							date_to_techniques_to_count[date]["hierarchical clustering"] += 1

					else:
						if("clustering" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["clustering"] = 1
						else:
							category_to_techniques[str(category)]["clustering"] += 1

						## get information by year
						if("clustering" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["clustering"] = 1
						else:
							date_to_techniques_to_count[date]["clustering"] += 1

				## tree
				elif(word in ["tree", "Tree"]):
					if(abstract_in_array[cmpt-1] in ["decision", "Decision"]):
						if("decision tree" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["decision tree"] = 1
						else:
							category_to_techniques[str(category)]["decision tree"] += 1

						## get information by year
						if("decision tree" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["decision tree"] = 1
						else:
							date_to_techniques_to_count[date]["decision tree"] += 1

				## random forest
				elif(word in ["forest", "Forest"]):
					if(abstract_in_array[cmpt-1] in ["random", "Random"]):
						if("random forest" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["random forest"] = 1
						else:
							category_to_techniques[str(category)]["random forest"] += 1

						## get information by year
						if("random forest" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["random forest"] = 1
						else:
							date_to_techniques_to_count[date]["random forest"] += 1

				## artificial neural network
				elif(word in ["neural", "Neural"]):
					if(abstract_in_array[cmpt+1] in ["network", "networks"] and (abstract_in_array[cmpt-1] in ["artificial", "Artificial"] or abstract_in_array[cmpt+2] in ["algorithm", "algorithms"])):
						if("artificial neural network" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["artificial neural network"] = 1
						else:
							category_to_techniques[str(category)]["artificial neural network"] += 1

						## get information by year
						if("artificial neural network" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["artificial neural network"] = 1
						else:
							date_to_techniques_to_count[date]["artificial neural network"] += 1

				## machine learning
				## support vector machine
				elif(word in ["machine", "Machine"]):
					if(abstract_in_array[cmpt+1] in ["learning", "Learning"]):
						if("machine learning" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["machine learning"] = 1
						else:
							category_to_techniques[str(category)]["machine learning"] += 1

						## get information by year
						if("machine learning" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["machine learning"] = 1
						else:
							date_to_techniques_to_count[date]["machine learning"] += 1

					elif(abstract_in_array[cmpt-1] == "vector" and abstract_in_array[cmpt-2] in ["support", "Support"]):
						if("support vector machine" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["support vector machine"] = 1
						else:
							category_to_techniques[str(category)]["support vector machine"] += 1

						## get information by year
						if("support vector machine" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["support vector machine"] = 1
						else:
							date_to_techniques_to_count[date]["support vector machine"] += 1


				elif(word in ["machine-learning", "Machine-learning"]):
					if("machine learning" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["machine learning"] = 1
					else:
						category_to_techniques[str(category)]["machine learning"] += 1

					## get information by year
					if("machine learning" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["machine learning"] = 1
					else:
						date_to_techniques_to_count[date]["machine learning"] += 1

				## Wilcoxon test
				elif(word in ["Wilcoxon", "wilcoxon"]):
					if("Wilcoxon test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["Wilcoxon test"] = 1
					else:
						category_to_techniques[str(category)]["Wilcoxon test"] += 1

					## get information by year
					if("Wilcoxon test" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["Wilcoxon test"] = 1
					else:
						date_to_techniques_to_count[date]["Wilcoxon test"] += 1

				## PCA
				elif(word in ["Principal", "principal"]):
					if(abstract_in_array[cmpt+1] == "component" and abstract_in_array[cmpt+2] == "analysis"):
						if("PCA" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["PCA"] = 1
						else:
							category_to_techniques[str(category)]["PCA"] += 1

						## get information by year
						if("PCA" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["PCA"] = 1
						else:
							date_to_techniques_to_count[date]["PCA"] += 1

				## chi-squared
				elif(word in ["chi-squared", "Chi-squared", "Chi-sqare", "ch-square"]):
					if("chi-squared" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["chi-squared"] = 1
					else:
						category_to_techniques[str(category)]["chi-squared"] += 1

					## get information by year
					if("chi-squared" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["chi-squared"] = 1
					else:
						date_to_techniques_to_count[date]["chi-squared"] += 1

				## Fisher's test
				elif(word in ["Fisher's", "Fisher"]):
					if("Fisher test" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["Fisher test"] = 1
					else:
						category_to_techniques[str(category)]["Fisher test"] += 1

					## get information by year
					if("Fisher test" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["Fisher test"] = 1
					else:
						date_to_techniques_to_count[date]["Fisher test"] += 1

				## analysis
				## - multivariate
				## - univariate
				## - discriminant
				## - bioinformatics
				elif(word in ["analysis", "Analysis"]):

					if(abstract_in_array[cmpt-1] in ["multivariate", "Multivariate"]):
						if("multivariate analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["multivariate analysis"] = 1
						else:
							category_to_techniques[str(category)]["multivariate analysis"] += 1

						## get information by year
						if("multivariate analysis" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["multivariate analysis"] = 1
						else:
							date_to_techniques_to_count[date]["multivariate analysis"] += 1

					elif(abstract_in_array[cmpt-1] in ["univariate", "Univariate"]):
						if("univariate analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["univariate analysis"] = 1
						else:
							category_to_techniques[str(category)]["univariate analysis"] += 1

						## get information by year
						if("univariate analysis" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["univariate analysis"] = 1
						else:
							date_to_techniques_to_count[date]["univariate analysis"] += 1

					elif(abstract_in_array[cmpt-1] in ["discriminant", "Discriminant"]):
						if("discriminant analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["discriminant analysis"] = 1
						else:
							category_to_techniques[str(category)]["discriminant analysis"] += 1

						## get information by year
						if("discriminant analysis" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["discriminant analysis"] = 1
						else:
							date_to_techniques_to_count[date]["discriminant analysis"] += 1

					elif(abstract_in_array[cmpt-1] in ["bioinformatics", "Bioinformatics", "bioinformatic", "Bioinformatic"]):
						if("bioinformatics analysis" not in category_to_techniques[str(category)].keys()):
							category_to_techniques[str(category)]["bioinformatics analysis"] = 1
						else:
							category_to_techniques[str(category)]["bioinformatics analysis"] += 1

						## get information by year
						if("bioinformatics analysis" not in date_to_techniques_to_count[date].keys()):
							date_to_techniques_to_count[date]["bioinformatics analysis"] = 1
						else:
							date_to_techniques_to_count[date]["bioinformatics analysis"] += 1

				## GLM (generalized linear model)
				elif(word in ["GLM", "(GLM)", "(GLM, ", "GLM)"]):
					if("GLM" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["GLM"] = 1
					else:
						category_to_techniques[str(category)]["GLM"] += 1

					## get information by year
					if("GLM" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["GLM"] = 1
					else:
						date_to_techniques_to_count[date]["GLM"] += 1

				## mixed models
				elif(word in ["models", "model"] and abstract_in_array[cmpt-1] in ["mixed", "Mixed"]):
					if("mixed models" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["mixed models"] = 1
					else:
						category_to_techniques[str(category)]["mixed models"] += 1

					## get information by year
					if("mixed models" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["mixed models"] = 1
					else:
						date_to_techniques_to_count[date]["mixed models"] += 1


				## normalization
				elif(word in ["Normalization", "normalization"]):
					if("normalization" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["normalization"] = 1
					else:
						category_to_techniques[str(category)]["normalization"] += 1

					## get information by year
					if("normalization" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["normalization"] = 1
					else:
						date_to_techniques_to_count[date]["normalization"] += 1

				## k mean
				elif(word in ["kmean", "Kmean", "K-mean", "K-mean-clustering"]):
					if("k mean" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["k mean"] = 1
					else:
						category_to_techniques[str(category)]["k mean"] += 1

					## get information by year
					if("k mean" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["k mean"] = 1
					else:
						date_to_techniques_to_count[date]["k mean"] += 1

				elif(word in ["k", "K"] and abstract_in_array[cmpt+1] in ["mean", "means", "Mean", "Means"]):
					if("k mean" not in category_to_techniques[str(category)].keys()):
						category_to_techniques[str(category)]["k mean"] = 1
					else:
						category_to_techniques[str(category)]["k mean"] += 1

					## get information by year
					if("k mean" not in date_to_techniques_to_count[date].keys()):
						date_to_techniques_to_count[date]["k mean"] = 1
					else:
						date_to_techniques_to_count[date]["k mean"] += 1

				cmpt += 1

	## Display details of the returned structure
	if(display_details):
		total = 119 + 604 +726 +19
		detected = 0
		for key in category_to_techniques.keys():
			print "=> " +str(key) +" <="
			for tech in category_to_techniques[key].keys():
				print "    -> "+str(tech) +" : " +str(category_to_techniques[key][tech])
				detected += category_to_techniques[key][tech]
		print "=> "+str(detected) +" / " + str(total) +" || [ "+str(float(float(detected)/float(total))*100) +" % ]"


	## return dictionnary
	return (category_to_techniques, date_to_techniques_to_count)


def generate_techniques_figures(category_to_techniques_to_count):
	##
	## Generate 8 figures from techniques count data
	## category_to_techniques_to_count is a dictionnary obtain by the
	## retrieve_techniques() function.
	##
	## save generated figures in images folder.
	##

	## Define category of analysis
	regression_analysis = ["logistic regression", "cox regression", "multivariate regression",
	"univariate regression", "linear regression", "Poisson regression", "regression tree", "regression analysis"]
	machine_learning_analysis = ["random forest", "artificial neural network", "support vector machine", "machine learning",
	"GLM", "mixed models", "k mean", "clustering", "hierarchical clustering", "decision tree"]
	other_analysis = ["t-test", "Wilcoxon test", "PCA", "chi-squared", "multivariate analysis", "univariate analysis", "discriminant analysis",
	"bioinformatic analysis", "normalization"]

	## Detailed histogramm
	Total_number_of_techniques = 0
	total_regression_techniques = 0
	total_machine_learning_techniques = 0
	total_other_techniques = 0
	global_techniques_data = {}
	category_to_proportion = {}
	for category in category_to_techniques_to_count.keys():

		category_to_proportion[category] = {}
		regression_techniques = 0
		machine_learning_techniques = 0
		other_techniques = 0
		number_of_techniques = 0
		for tech in category_to_techniques_to_count[category].keys():
			number_of_techniques += category_to_techniques_to_count[category][tech]
			Total_number_of_techniques += category_to_techniques_to_count[category][tech]

			if(tech not in global_techniques_data.keys()):
				global_techniques_data[tech] = category_to_techniques_to_count[category][tech]
			else:
				global_techniques_data[tech] += category_to_techniques_to_count[category][tech]

			if(tech in regression_analysis):
				regression_techniques += category_to_techniques_to_count[category][tech]
				total_regression_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in machine_learning_analysis):
				machine_learning_techniques += category_to_techniques_to_count[category][tech]
				total_machine_learning_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in other_analysis):
				other_techniques += category_to_techniques_to_count[category][tech]
				total_other_techniques += category_to_techniques_to_count[category][tech]

		techniques_to_count = category_to_techniques_to_count[category]
		for tech in techniques_to_count.keys():
			techniques_to_count[tech] = float(techniques_to_count[tech]) / float(number_of_techniques) *100


		## Detailed Histogram
		plt.bar(techniques_to_count.keys(), techniques_to_count.values(), color='b', align='center', width=0.3)
		plt.xticks(rotation=45)
		plt.savefig("images/techniques_histogram_"+str(category)+".png")
		plt.close()

		## Pie chart
		regression_techniques = float(regression_techniques) / float(number_of_techniques) * 100
		machine_learning_techniques = float(machine_learning_techniques) / float(number_of_techniques) * 100
		other_techniques = float(other_techniques) / float(number_of_techniques) * 100
		x_vector = [regression_techniques, machine_learning_techniques, other_techniques]
		labels = ["regression techniques", "machine learning techniques", "other techniques"]
		plt.pie(x_vector, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
		plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
		plt.savefig("images/techniques_pie_"+str(category)+".png")
		plt.close()

		category_to_proportion[category]["regression"] = regression_techniques
		category_to_proportion[category]["machineLearning"] = machine_learning_techniques
		category_to_proportion[category]["other"] = other_techniques


	## Global histogram
	plt.bar(global_techniques_data.keys(), global_techniques_data.values(), color='b', align='center', width=0.3)
	plt.xticks(rotation=45)
	plt.savefig("images/techniques_histogram_global.png")
	plt.close()

	## Global pie chart
	total_regression_techniques = float(total_regression_techniques) / float(Total_number_of_techniques) * 100
	total_machine_learning_techniques = float(total_machine_learning_techniques) / float(Total_number_of_techniques) * 100
	total_other_techniques = float(total_other_techniques) / float(Total_number_of_techniques) * 100
	x_vector = [total_regression_techniques, total_machine_learning_techniques, total_other_techniques]
	labels = ["total regression techniques", "total machine learning techniques", "total other techniques"]
	plt.pie(x_vector, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
	plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	plt.savefig("images/techniques_pie_global.png")
	plt.close()

	category_to_proportion["all"] = {}
	category_to_proportion["all"]["regression"] = total_regression_techniques
	category_to_proportion["all"]["machineLearning"] = total_machine_learning_techniques
	category_to_proportion["all"]["other"] = total_other_techniques

	return category_to_proportion


def get_cohorte_size(abstract_file):
	##
	## IN PROGRESS
	##
	## -> Try to retrieve the number of patients/case
	## in the dataset presented in the abstract
	## 
	## TODO : add more regex for size detection
	## TODO : add to bibotlite.py
	## TODO : get MAJ from trashlib
	##

	## useful_data_structure


	## read abstract
	abstract = open(abstract_file, "r")

	cohorte_size = -1
	for line in abstract:


		## try to catch a written number before the apparition 
		## of the term in catch list
		catch_terms = ["patients", "cases", "observations", "subjects"]
		line_in_array = line.split(" ")
		cmpt_in_line = 0
		for elt in line_in_array:	
			if(elt in catch_terms):
				try:
					fetched_cohorte_size = w2n.word_to_num(line_in_array[cmpt_in_line-1])
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
						nothing_retrieved = False
				except:
					nothing_retrieved = True

				if(nothing_retrieved):
					try:
						fetched_cohorte_size = w2n.word_to_num(line_in_array[cmpt_in_line-2])
						if(fetched_cohorte_size > cohorte_size):
							cohorte_size = fetched_cohorte_size
							nothing_retrieved = False
					except:
						nothing_retrieved = True

				if(nothing_retrieved):
					try:
						fetched_cohorte_size = w2n.word_to_num(line_in_array[cmpt_in_line-3])
						if(fetched_cohorte_size > cohorte_size):
							cohorte_size = fetched_cohorte_size
							nothing_retrieved = False
					except:
						nothing_retrieved = True
			cmpt_in_line += 1



		## Try to catch a number before the apparition of the 
		## term in catch terms list. keep the biggest number
		## found in the abstract.
		catch_terms = ["patients", "cases", "observations", "subjects"]
		for item in catch_terms:
			m = re.findall(r"([0-9]+[\.|,]?[0-9]+) "+str(item), line)		
			if(m is not None):
				for match in m:
					match = match.replace(",", "")
					try:
						fetched_cohorte_size = int(match)
						if(fetched_cohorte_size > cohorte_size):
							cohorte_size = fetched_cohorte_size
					except:
						## do nothing
						tardis = 1


		## Try to catch a number after the apparition of the 
		## term "n = ". keep the biggest number found in the
		## abstract.
		m = re.findall(r"n = ([0-9]+[\.|,]?[0-9]+)", line)		
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1

		## Try to catch a number in a classic expression
		m = re.findall(r"[I,i]n total, ([0-9]+[\.|,]?[0-9]+) subjects", line)
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1

		"""
		## Try to catch a number in a classic expression
		m = re.findall(r"([0-9]+[\.|,]?[0-9]+)", line)
		if(m is not None):
			for match in m:
				match = match.replace(",", "")
				try:
					fetched_cohorte_size = int(match)
					if(fetched_cohorte_size > cohorte_size):
						cohorte_size = fetched_cohorte_size
				except:
					## do nothing
					tardis = 1
		"""


	## close abstract
	abstract.close()

	## check if the cohorte_size variable
	## is still set to -1, set it to NA if
	## True
	if(cohorte_size == -1):
		cohorte_size = "NA"

	## return the detected size of the cohorte
	return cohorte_size



def get_cohorte_size_for_category():
	##
	## Get the count of patients parsed from abstract for each 
	## each category.
	##
	## return a dictionnary
	##

	## retrieve data
	category_list = ["Diagnostic", "Therapeutic", "Unclassified", "Modelisation"]
	category_to_size = {}
	category_to_size["Retrieved"] = 0
	for category in category_list:
		abstract_list = glob.glob("articles_classification/"+str(category)+"/*")
		category_to_size[category] = 0
		for abstract in abstract_list:
			size = get_cohorte_size(abstract)
			if(size != "NA"):
				category_to_size[category] += size
				category_to_size["Retrieved"] += 1

				if(int(size) > 1000):
					print "[WARNINGS] => flag size "+str(size)+" in "+str(abstract)

	for key in category_to_size.keys():
		print "[DATA] => "+str(key) +" : " +str(category_to_size[key])
			
	return category_to_size



def generate_cohorte_size_figures(run_folder):
	##
	## Generate 3 figures about the evolution of 
	## the sizes of cohorte, save the image in
	## the images folder:
	##    - category_sizes_histogram.png
	##    - count_cohorte_size_evolution.png
	##    - median_cohorte_size_evolution.png
	##
	## run_folder is the name of the run_folder, 
	## should contain at least:
	##    - abstract directory
	##    - meta directory
	##
	##

	## Figure 1
	## Histogram with category
	category_to_size = get_cohorte_size_for_category()
	total_number_of_articles = 0
	for value in category_to_size.values():
		total_number_of_articles += value
	for key in category_to_size.keys():
		category_to_size[key] = float(category_to_size[key]) / float(total_number_of_articles) * 100
	plt.bar(category_to_size.keys(), category_to_size.values(), color='b', align='center', width=0.3)
	plt.xticks(rotation=45)
	plt.savefig("images/category_sizes_histogram.png")
	plt.close()

	## Figure 2 and 3
	## Evolution curve
	year_to_pmid = {}
	for meta_file in glob.glob(run_folder+"/meta/*.csv"):
		pmid = meta_file.split("/")
		pmid = pmid[-1]
		pmid = pmid.split(".")
		pmid = pmid[0]
		year = get_date_from_meta_save(meta_file)
		if(year not in year_to_pmid.keys()):
			year_to_pmid[year] = []
			year_to_pmid[year].append(pmid)
		else:
			year_to_pmid[year].append(pmid)

	## -> Global count evolution
	year_to_patient = {}
	for year in year_to_pmid.keys():
		if(int(year) < 2018):
			number_of_patients = 0
			for pmid in year_to_pmid[year]:
				abstract_file = run_folder+"/abstract/"+str(pmid)+"_abstract.txt"
				patients_retrieved = get_cohorte_size(abstract_file)
				if(patients_retrieved != "NA"):
					number_of_patients += int(patients_retrieved)
			year_to_patient[year] = number_of_patients
	ordered_data = collections.OrderedDict(sorted(year_to_patient.items()))
	x = []
	for elt in ordered_data.keys():
		x.append(int(elt))
	y = ordered_data.values()
	m,b = np.polyfit(x, y, 1)
	print "[INFO] => coeff : "+str(m)+" * x + "+str(b)
	fit = np.polyfit(x, y, 1)
	fit_fn = np.poly1d(fit) 
	plt.plot(x, y, 'b-*', label="cohorte sizes")
	plt.plot(x, fit_fn(x), '--k', label="linear regression")
	plt.legend()
	plt.ylim(ymin=0)
	plt.savefig("images/count_cohorte_size_evolution.png")
	plt.close()

	## -> Median evolution
	year_to_median = {}
	for year in year_to_pmid.keys():
		if(int(year) < 2018):
			cohort_size = []
			for pmid in year_to_pmid[year]:
				abstract_file = run_folder+"/abstract/"+str(pmid)+"_abstract.txt"
				patients_retrieved = get_cohorte_size(abstract_file)
				if(patients_retrieved != "NA"):
					cohort_size.append(patients_retrieved)
			try:
				cohorte_median = np.median(cohort_size)
			except:
				cohorte_median = 0
			year_to_median[year] = cohorte_median
	ordered_data = collections.OrderedDict(sorted(year_to_median.items()))
	plt.plot(ordered_data.keys(), ordered_data.values())
	plt.savefig("images/median_cohorte_size_evolution.png")
	plt.close()



def get_techniques_proportions(category_to_techniques_to_count):
	##
	##
	## Generate a dictionnary with proportion data for the
	## differents techniques retrieve in abstracts
	## category_to_techniques_to_count is a dictionnary obtain by the
	## retrieve_techniques() function.
	##
	## return a dictionnary
	##

	## Define category of analysis
	regression_analysis = ["logistic regression", "cox regression", "multivariate regression",
	"univariate regression", "linear regression", "Poisson regression", "regression tree", "regression analysis"]
	machine_learning_analysis = ["random forest", "artificial neural network", "support vector machine", "machine learning",
	"GLM", "mixed models", "k mean", "clustering", "hierarchical clustering", "decision tree"]
	other_analysis = ["t-test", "Wilcoxon test", "PCA", "chi-squared", "multivariate analysis", "univariate analysis", "discriminant analysis",
	"bioinformatic analysis", "normalization"]

	## Detailed histogramm
	Total_number_of_techniques = 0
	total_regression_techniques = 0
	total_machine_learning_techniques = 0
	total_other_techniques = 0
	global_techniques_data = {}
	category_to_proportion = {}
	for category in category_to_techniques_to_count.keys():

		category_to_proportion[category] = {}

		regression_techniques = 0
		machine_learning_techniques = 0
		other_techniques = 0
		number_of_techniques = 0
		for tech in category_to_techniques_to_count[category].keys():
			number_of_techniques += category_to_techniques_to_count[category][tech]
			Total_number_of_techniques += category_to_techniques_to_count[category][tech]

			if(tech not in global_techniques_data.keys()):
				global_techniques_data[tech] = category_to_techniques_to_count[category][tech]
			else:
				global_techniques_data[tech] += category_to_techniques_to_count[category][tech]

			if(tech in regression_analysis):
				regression_techniques += category_to_techniques_to_count[category][tech]
				total_regression_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in machine_learning_analysis):
				machine_learning_techniques += category_to_techniques_to_count[category][tech]
				total_machine_learning_techniques += category_to_techniques_to_count[category][tech]
			elif(tech in other_analysis):
				other_techniques += category_to_techniques_to_count[category][tech]
				total_other_techniques += category_to_techniques_to_count[category][tech]

		techniques_to_count = category_to_techniques_to_count[category]
		for tech in techniques_to_count.keys():
			techniques_to_count[tech] = float(techniques_to_count[tech]) / float(number_of_techniques) *100


		## proportion by category of techniques
		regression_techniques = float(regression_techniques) / float(number_of_techniques) * 100
		machine_learning_techniques = float(machine_learning_techniques) / float(number_of_techniques) * 100
		other_techniques = float(other_techniques) / float(number_of_techniques) * 100
		category_to_proportion[category]["regression"] = regression_techniques
		category_to_proportion[category]["machineLearning"] = machine_learning_techniques
		category_to_proportion[category]["other"] = other_techniques

	## global proportions
	total_regression_techniques = float(total_regression_techniques) / float(Total_number_of_techniques) * 100
	total_machine_learning_techniques = float(total_machine_learning_techniques) / float(Total_number_of_techniques) * 100
	total_other_techniques = float(total_other_techniques) / float(Total_number_of_techniques) * 100	
	category_to_proportion["all"] = {}
	category_to_proportion["all"]["regression"] = total_regression_techniques
	category_to_proportion["all"]["machineLearning"] = total_machine_learning_techniques
	category_to_proportion["all"]["other"] = total_other_techniques

	return category_to_proportion


def generate_tech_evolution_figures(year_to_tech):
	##
	## generate graph evolution of techniques (by count in abstract)
	## over the last decade
	##
	## year_to_tech is a dictionnary return by 
	## the retrieve_techniques functions
	##
	## save figure as images/count_techniques_evolution.png
	##
	##

	## Define category of analysis
	regression_analysis = ["logistic regression", "cox regression", "multivariate regression",
	"univariate regression", "linear regression", "Poisson regression", "regression tree", "regression analysis"]
	machine_learning_analysis = ["random forest", "artificial neural network", "support vector machine", "machine learning",
	"GLM", "mixed models", "k mean", "clustering", "hierarchical clustering", "decision tree"]
	other_analysis = ["t-test", "Wilcoxon test", "PCA", "chi-squared", "multivariate analysis", "univariate analysis", "discriminant analysis",
	"bioinformatic analysis", "normalization"]

	year_to_count = {}

	for year in year_to_tech.keys():
		year_to_count[year] = {}
				
		regression_techniques = 0
		machine_learning_techniques = 0
		other_techniques = 0
		
		for tech in year_to_tech[year].keys():
			
			if(tech in regression_analysis):
				regression_techniques += year_to_tech[year][tech]
			elif(tech in machine_learning_analysis):
				machine_learning_techniques += year_to_tech[year][tech]
			elif(tech in other_analysis):
				other_techniques += year_to_tech[year][tech]

		year_to_count[year]["regression"] = regression_techniques
		year_to_count[year]["machine_learning"] = machine_learning_techniques
		year_to_count[year]["other"] = other_techniques


	## create the figure
	y_vector = []
	for year in year_to_count.keys():
		if(int(year) < 2018):
			y_vector.append(year)
	y_vector = sorted(y_vector)
	r_vector = []
	for year in y_vector:
		r_vector.append(year_to_count[year]["regression"])
	m_vector = []
	for year in y_vector:
		m_vector.append(year_to_count[year]["machine_learning"])
	o_vector = []
	for year in y_vector:
		o_vector.append(year_to_count[year]["other"])

	plt.plot(y_vector, r_vector, 'r--', label="regression")
	plt.plot(y_vector, m_vector, 'b-*', label="machine learning")
	plt.plot(y_vector, o_vector, 'g-^', label="other")
	plt.legend()
	plt.savefig("images/count_techniques_evolution.png")
	plt.close()


def get_cohorte_size_for_year(year, run_folder):
	##
	## Get the cohorte size (total count of patients and
	## median size) for a given year
	##
	## year is an int, the year of interest,
	## run_folder is the run_folder
	##
	## return a tuple (total_count, median)
	##


	total_size = 0
	size_median = 0
	sizes_list = []

	category_list = ["Diagnostic", "Therapeutic", "Unclassified", "Modelisation"]
	for category in category_list:
		abstract_list = glob.glob("articles_classification/"+str(category)+"/*")
		for abstract in abstract_list:
			abstract_in_array = abstract.split("/")
			abstract_in_array = abstract_in_array[-1]
			pmid = abstract_in_array.split("_")
			pmid = pmid[0]
			meta_file = run_folder+"/meta/"+str(pmid)+".csv"
			year_retrieved = get_date_from_meta_save(meta_file)

			if(int(year) == int(year_retrieved)):
				size = get_cohorte_size(abstract)
				if(size != "NA"):
					total_size += int(size)
					sizes_list.append(int(size))

	if(len(sizes_list) > 0):
		size_median = np.median(sizes_list)

	return (total_size, size_median)



############
## SCRIPT #############################################################################################
############


##--------------------------------------##
## Get the Run Folder from command line ##
##--------------------------------------##
if(os.path.isdir(str(sys.argv[1])) and os.path.isdir(str(sys.argv[1])+"/meta") and os.path.isdir(str(sys.argv[1])+"/abstract")):
	run_folder = str(sys.argv[1])
	log_file = run_folder+"/bibotlite.log"
else:
	print "[ERROR] => Can't use the folder "+str(sys.argv[1])+" as run folder"
	sys.exit()


##-----------------------------##
## Write the FIXEME table file ##
##-----------------------------##

## variables initialisation
## Automatic retrieve
total_number_of_articles = "NA"						# OK
number_of_input_keywords = "NA"						# OK
articles_pass_the_first_filter = "NA"				# OK
articles_pass_the_last_filter = "NA"				# OK
articles_pass_both_filter = "NA"					# OK
machine_learning_SjS = "NA"							# OK
machine_learning_Other = "NA" 						# OK
machine_learning_SLE = "NA" 						# OK
machine_learning_RA = "NA" 							# OK
article_final_selection = "NA"						# OK
articles_pass_the_first_filter_proportion = "NA"	# OK
articles_pass_the_last_filter_proportion = "NA"		# OK
drop_1 = "NA"										# OK
drop_2 = "NA"										# OK
number_of_publication_from_first_contry = "NA"		# OK
number_of_publication_from_second_contry = "NA"		# OK
first_country = "NA"								# OK
second_country = "NA"								# OK
number_of_other_publications = "NA"					# OK
other_country_list = "NA"							# OK
articles_2008 = "NA"								# OK
articles_2017 = "NA"								# OK
parsed_abstract_for_size = "NA"						# OK
size_2008 = "NA"									# OK
size_2017 = "NA"									# OK
diagnostic_techniques_list = "NA"					# OK
diagnostic_number_of_articles = "NA"				# OK
diagnostic_number_of_patients = "NA"				# OK
diagnostic_sjogren = "NA"							# OK
diagnostic_other = "NA"								# OK
therapeutic_techniques_list = "NA"					# OK
therapeutic_number_of_articles = "NA"				# OK
therapeutic_number_of_patients = "NA"				# OK
therapeutic_sjogren = "NA"							# OK
therapeutic_other = "NA"							# OK
modelisation_techniques_list = "NA"					# OK
modelisation_number_of_articles = "NA"				# OK
modelisation_number_of_patients = "NA"				# OK
modelisation_sjogren = "NA"							# OK
modelisation_other = "NA"							# OK
unclassified_techniques_list = "NA"					# OK
unclassified_number_of_articles = "NA"				# OK
unclassified_number_of_patients = "NA"				# OK
unclassified_sjogren = "NA"							# OK
unclassified_other = "NA"							# OK
diagnostic_sjogren_proportion = "NA"				# OK
therapeutic_sjogren_proportion = "NA"				# OK
modelisation_sjogren_proportion = "NA"				# OK
unclassified_sjogren_proportion = "NA"				# OK
diagnostic_machine_learning_proportion = "NA"		# OK
diagnostic_regression_proportion = "NA"				# OK
diagnostic_other_proportion = "NA"					# OK
therapeutic_machine_learning_proportion = "NA"		# OK
therapeutic_regression_proportion = "NA"			# OK
therapeutic_other_proportion = "NA"					# OK
modelisation_machine_learning_proportion = "NA"		# OK
modelisation_regression_proportion = "NA"			# OK
modelisation_other_proportion = "NA"				# OK
unclassified_machine_learning_proportion = "NA"		# OK
unclassified_regression_proportion = "NA"			# OK
unclassified_other_proportion = "NA"				# OK
global_machine_learning_proportion = "NA"			# OK
global_regression_proportion = "NA"					# OK
global_other_proportion = "NA"						# OK



## Manual retrieve
term1 = "NA"										# TODO
term2 = "NA"										# TODO
term3 = "NA"										# TODO
term4 = "NA"										# TODO
term5 = "NA"										# TODO
term6 = "NA"										# TODO
article_manually_retrieve = "NA"					# TODO


## save existing manifest file if it already exist
if(os.path.isfile("draft_FIXEME_table.csv")):
	now = datetime.datetime.now()
	tag = str(now.minute)+"_"+str(now.hour)+"_"+str(now.day)+"_"+str(now.month)
	shutil.copy("draft_FIXEME_table.csv", "draft_FIXEME_table_"+tag+".csv",)


## Get variables values
disease_to_articles = describe_articles_type(run_folder+"/abstract")
machine_learning_SjS = disease_to_articles["SjS"]
machine_learning_Other = disease_to_articles["Other"]
machine_learning_SLE = disease_to_articles["SLE"]
machine_learning_RA = disease_to_articles["RA"]
total_number_of_articles = get_number_of_articles_from_log_file(log_file)
number_of_input_keywords = get_number_of_keywords_from_log_file(log_file)
article_stats = get_articles_filter_stat(log_file)
articles_pass_the_first_filter = article_stats[0]
articles_pass_the_last_filter = article_stats[1]
articles_pass_both_filter = article_stats[2]
articles_pass_the_first_filter_proportion = float(float(articles_pass_the_first_filter)/float(total_number_of_articles))*100
articles_pass_the_last_filter_proportion = float(float(articles_pass_the_last_filter)/float(total_number_of_articles))*100
drop_1 = total_number_of_articles - articles_pass_the_first_filter
drop_2 = articles_pass_the_first_filter - articles_pass_both_filter
articles_2008 = get_number_of_articles_for_year(run_folder, 2008)
articles_2017 = get_number_of_articles_for_year(run_folder, 2017)
if(article_manually_retrieve == "NA"):
	article_final_selection = articles_pass_both_filter
else:
	article_final_selection = articles_pass_both_filter + article_manually_retrieve
country_stats = get_country_publication_stat(run_folder)
print country_stats
first_country = max(country_stats, key=country_stats.get)
number_of_publication_from_first_contry = country_stats[first_country]
del country_stats[first_country]
if(len(country_stats.keys()) > 0):
	second_country = max(country_stats, key=country_stats.get)
	number_of_publication_from_second_contry = country_stats[second_country]
	del country_stats[second_country]
	number_of_other_publications = 0
	for val in country_stats.values():
		number_of_other_publications += val
	other_country_list = ""
	sorted_stats = sorted(country_stats.items(), key=operator.itemgetter(1))
	for key in sorted_stats[::-1]:
		other_country_list += str(key[0]) +", "
	other_country_list = other_country_list[:-2]


class_to_count = {"Diagnostic":0, "Therapeutic":0, "Modelisation":0, "Unclassified":0}
for abstract in glob.glob(str(run_folder)+"/abstract/*_abstract.txt"):
	theme = get_article_domain(abstract)
	class_to_count[theme] += 1
	pmid = abstract.split("/")
	pmid = pmid[-1]
	shutil.copy(abstract, "articles_classification/"+theme+"/"+pmid)
unclassified_number_of_articles = class_to_count["Unclassified"]
modelisation_number_of_articles = class_to_count["Modelisation"]
diagnostic_number_of_articles = class_to_count["Diagnostic"]
therapeutic_number_of_articles = class_to_count["Therapeutic"]


for theme in class_to_count.keys():
	disease_to_articles = describe_articles_type("articles_classification/"+theme)

	print "=> " + str(theme)
	print disease_to_articles

	sjs_patients = disease_to_articles["SjS"]
	if(theme == "Diagnostic"):
		diagnostic_sjogren = sjs_patients
		diagnostic_other = diagnostic_number_of_articles - sjs_patients
		diagnostic_sjogren_proportion = float(diagnostic_sjogren) / float(diagnostic_number_of_articles) * 100
	elif(theme == "Therapeutic"):
		therapeutic_sjogren = sjs_patients
		therapeutic_other = therapeutic_number_of_articles - sjs_patients
		therapeutic_sjogren_proportion = float(therapeutic_sjogren) / float(therapeutic_number_of_articles) * 100
	elif(theme == "Modelisation"):
		modelisation_sjogren = sjs_patients
		modelisation_other = modelisation_number_of_articles - sjs_patients
		modelisation_sjogren_proportion = float(modelisation_sjogren) / float(modelisation_number_of_articles) * 100
	elif(theme == "Unclassified"):
		unclassified_sjogren = sjs_patients
		unclassified_other = unclassified_number_of_articles - sjs_patients
		unclassified_sjogren_proportion = float(unclassified_sjogren) / float(unclassified_number_of_articles) * 100

## Techniques
diagnostic_techniques_list = "\\begin{tabular}{c}"
modelisation_techniques_list = "\\begin{tabular}{c}"
therapeutic_techniques_list = "\\begin{tabular}{c}"
unclassified_techniques_list = "\\begin{tabular}{c}"

category_to_techniques_to_count, year_to_techniques_to_count = retrieve_techniques(run_folder, True)
for tech in category_to_techniques_to_count["Diagnostic"].keys():
	diagnostic_techniques_list += str(tech)+" \\\\"
for tech in category_to_techniques_to_count["Therapeutic"].keys():
	therapeutic_techniques_list += str(tech)+" \\\\"
for tech in category_to_techniques_to_count["Modelisation"].keys():
	modelisation_techniques_list += str(tech)+" \\\\"
for tech in category_to_techniques_to_count["Unclassified"].keys():
	unclassified_techniques_list += str(tech)+" \\\\"

diagnostic_techniques_list += "\\end{tabular}"
modelisation_techniques_list += "\\end{tabular}"
therapeutic_techniques_list += "\\end{tabular}"
unclassified_techniques_list += "\\end{tabular}"

## Techniques proportion
#category_to_proportion = get_techniques_proportions(category_to_techniques_to_count)
category_to_proportion = generate_techniques_figures(category_to_techniques_to_count)
diagnostic_machine_learning_proportion = category_to_proportion["Diagnostic"]["machineLearning"]
diagnostic_regression_proportion = category_to_proportion["Diagnostic"]["regression"]
diagnostic_other_proportion = category_to_proportion["Diagnostic"]["other"]
therapeutic_machine_learning_proportion = category_to_proportion["Therapeutic"]["machineLearning"]
therapeutic_regression_proportion = category_to_proportion["Therapeutic"]["regression"]
therapeutic_other_proportion = category_to_proportion["Therapeutic"]["other"]
modelisation_machine_learning_proportion = category_to_proportion["Modelisation"]["machineLearning"]
modelisation_regression_proportion = category_to_proportion["Modelisation"]["regression"]
modelisation_other_proportion = category_to_proportion["Modelisation"]["other"]
unclassified_machine_learning_proportion = category_to_proportion["Unclassified"]["machineLearning"]
unclassified_regression_proportion = category_to_proportion["Unclassified"]["regression"]
unclassified_other_proportion = category_to_proportion["Unclassified"]["other"]
global_machine_learning_proportion = category_to_proportion["all"]["machineLearning"]
global_regression_proportion = category_to_proportion["all"]["regression"]
global_other_proportion = category_to_proportion["all"]["other"]

## cohorte sizes
category_to_size = get_cohorte_size_for_category()
diagnostic_number_of_patients = category_to_size["Diagnostic"]
therapeutic_number_of_patients = category_to_size["Therapeutic"]
modelisation_number_of_patients = category_to_size["Modelisation"]
unclassified_number_of_patients = category_to_size["Unclassified"]
parsed_abstract_for_size = category_to_size["Retrieved"]
size_2008 = get_cohorte_size_for_year(2008, run_folder)
size_2017 = get_cohorte_size_for_year(2017, run_folder)
size_2008 = size_2008[1] # get the median
size_2017 = size_2017[1] # get the median




## write a new manifest file
manifest_file = open("draft_FIXEME_table.csv", "w")
manifest_file.write("FIXME<totalnumberofarticles>;"+str(total_number_of_articles)+"\n")
manifest_file.write("FIXME<numberofinputkeywords>;"+str(number_of_input_keywords)+"\n")
manifest_file.write("FIXME<articlespassthefirstfilter>;"+str(articles_pass_the_first_filter)+"\n")
manifest_file.write("FIXME<articlespasslastfilter>;"+str(articles_pass_the_last_filter)+"\n")
manifest_file.write("FIXME<machinelearningSjS>;"+str(machine_learning_SjS)+"\n")
manifest_file.write("FIXME<machinelearningOther>;"+str(machine_learning_Other)+"\n")
manifest_file.write("FIXME<machinelearningSLE>;"+str(machine_learning_SLE)+"\n")
manifest_file.write("FIXME<machinelearningRA>;"+str(machine_learning_RA)+"\n")
manifest_file.write("FIXME<articlemanuallyretrieve>;"+str(article_manually_retrieve)+"\n")
manifest_file.write("FIXME<articlespassbothfilter>;"+str(articles_pass_both_filter)+"\n")
manifest_file.write("FIXME<articlefinalselection>;"+str(article_final_selection)+"\n")
manifest_file.write("FIXME<articlespassthefirstfilterproportion>;"+str(articles_pass_the_first_filter_proportion)+"\n")
manifest_file.write("FIXME<articlespasslastfilterproportion>;"+str(articles_pass_the_last_filter_proportion)+"\n")
manifest_file.write("FIXME<drop1>;"+str(drop_1)+"\n")
manifest_file.write("FIXME<drop2>;"+str(drop_2)+"\n")
manifest_file.write("FIXME<numberofpublicationfromfirstcontry>;"+str(number_of_publication_from_first_contry)+"\n")
manifest_file.write("FIXME<numberofpublicationfromsecondcontry>;"+str(number_of_publication_from_second_contry)+"\n")
manifest_file.write("FIXME<firstcountry>;"+str(first_country)+"\n")
manifest_file.write("FIXME<secondcountry>;"+str(second_country)+"\n")
manifest_file.write("FIXME<numberofotherpublications>;"+str(number_of_other_publications)+"\n")
manifest_file.write("FIXME<othercountrylist>;"+str(other_country_list)+"\n")
manifest_file.write("FIXME<articles2008>;"+str(articles_2008)+"\n")
manifest_file.write("FIXME<articles2017>;"+str(articles_2017)+"\n")
manifest_file.write("FIXME<parsedabstractforsize>;"+str(parsed_abstract_for_size)+"\n")
manifest_file.write("FIXME<size2008>;"+str(size_2008)+"\n")
manifest_file.write("FIXME<size2017>;"+str(size_2017)+"\n")
manifest_file.write("FIXME<term1>;"+str(term1)+"\n")
manifest_file.write("FIXME<term2>;"+str(term2)+"\n")
manifest_file.write("FIXME<term3>;"+str(term3)+"\n")
manifest_file.write("FIXME<term4>;"+str(term4)+"\n")
manifest_file.write("FIXME<term5>;"+str(term5)+"\n")
manifest_file.write("FIXME<term6>;"+str(term6)+"\n")
manifest_file.write("FIXME<diagnosticTechniquesList>;"+str(diagnostic_techniques_list)+"\n")
manifest_file.write("FIXME<diagnosticNumberOfArticles>;"+str(diagnostic_number_of_articles)+"\n")
manifest_file.write("FIXME<diagnosticNumberOfPatients>;"+str(diagnostic_number_of_patients)+"\n")
manifest_file.write("FIXME<diagnosticSjogren>;"+str(diagnostic_sjogren)+"\n")
manifest_file.write("FIXME<diagnosticOther>;"+str(diagnostic_other)+"\n")
manifest_file.write("FIXME<therapeuticTechniquesList>;"+str(therapeutic_techniques_list)+"\n")
manifest_file.write("FIXME<therapeuticNumberOfArticles>;"+str(therapeutic_number_of_articles)+"\n")
manifest_file.write("FIXME<therapeuticNumberOfPatients>;"+str(therapeutic_number_of_patients)+"\n")
manifest_file.write("FIXME<therapeuticSjogren>;"+str(therapeutic_sjogren)+"\n")
manifest_file.write("FIXME<therapeuticOther>;"+str(therapeutic_other)+"\n")
manifest_file.write("FIXME<modelisationTechniquesList>;"+str(modelisation_techniques_list)+"\n")
manifest_file.write("FIXME<modelisationNumberOfArticles>;"+str(modelisation_number_of_articles)+"\n")
manifest_file.write("FIXME<modelisationNumberOfPatients>;"+str(modelisation_number_of_patients)+"\n")
manifest_file.write("FIXME<modelisationSjogren>;"+str(modelisation_sjogren)+"\n")
manifest_file.write("FIXME<modelisationOther>;"+str(modelisation_other)+"\n")
manifest_file.write("FIXME<unclassifiedTechniquesList>;"+str(unclassified_techniques_list)+"\n")
manifest_file.write("FIXME<unclassifiedNumberOfArticles>;"+str(unclassified_number_of_articles)+"\n")
manifest_file.write("FIXME<unclassifiedNumberOfPatients>;"+str(unclassified_number_of_patients)+"\n")
manifest_file.write("FIXME<unclassifiedSjogren>;"+str(unclassified_sjogren)+"\n")
manifest_file.write("FIXME<unclassifiedOther>;"+str(unclassified_other)+"\n")
manifest_file.write("FIXME<diagnosticSjSproportion>;"+str(diagnostic_sjogren_proportion)+"\n")
manifest_file.write("FIXME<therapeuticSjSproportion>;"+str(therapeutic_sjogren_proportion)+"\n")
manifest_file.write("FIXME<modelisationSjSproportion>;"+str(modelisation_sjogren_proportion)+"\n")
manifest_file.write("FIXME<unclassifiedSjSproportion>;"+str(unclassified_sjogren_proportion)+"\n")
manifest_file.write("FIXME<diagnosticMachineLearningProportion>;"+str(diagnostic_machine_learning_proportion)+"\n")
manifest_file.write("FIXME<diagnosticRegressionProportion>;"+str(diagnostic_regression_proportion)+"\n")
manifest_file.write("FIXME<diagnosticOtherProportion>;"+str(diagnostic_other_proportion)+"\n")
manifest_file.write("FIXME<therapeuticMachineLearningProportion>;"+str(therapeutic_machine_learning_proportion)+"\n")
manifest_file.write("FIXME<therapeuticRegressionProportion>;"+str(therapeutic_regression_proportion)+"\n")
manifest_file.write("FIXME<therapeuticOtherProportion>;"+str(therapeutic_other_proportion)+"\n")
manifest_file.write("FIXME<modelisationMachineLearningProportion>;"+str(modelisation_machine_learning_proportion)+"\n")
manifest_file.write("FIXME<modelisationRegressionProportion>;"+str(modelisation_regression_proportion)+"\n")
manifest_file.write("FIXME<modelisationOtherProportion>;"+str(modelisation_other_proportion)+"\n")
manifest_file.write("FIXME<unclassifiedMachineLearningProportion>;"+str(unclassified_machine_learning_proportion)+"\n")
manifest_file.write("FIXME<unclassifiedRegressionProportion>;"+str(unclassified_regression_proportion)+"\n")
manifest_file.write("FIXME<unclassifiedOtherProportion>;"+str(unclassified_other_proportion)+"\n")
manifest_file.write("FIXME<globalMachineLearningProportion>;"+str(global_machine_learning_proportion)+"\n")
manifest_file.write("FIXME<globalRegressionProportion>;"+str(global_regression_proportion)+"\n")
manifest_file.write("FIXME<globalOtherProportion>;"+str(global_other_proportion)+"\n")
manifest_file.close()



##--------------------##
## Create The Figures ##
##--------------------##

## Check if the year publication figure already exist,
## save it with a time_tag if yes.
if(os.path.isfile("images/years_publications_evolution.png")):
	now = datetime.datetime.now()
	tag = str(now.minute)+"_"+str(now.hour)+"_"+str(now.day)+"_"+str(now.month)
	shutil.copy("images/years_publications_evolution.png", "images/years_publications_evolution_"+tag+".png")

## Create new year publication figure
plot_pulbications_years(str(run_folder)+"/meta")

## create the techniques figure
generate_techniques_figures(category_to_techniques_to_count)
generate_tech_evolution_figures(year_to_techniques_to_count)

## Create database size figure
generate_cohorte_size_figures(run_folder)




## Create frquency plot
## TODO : find the optimal graphe among the generated graphe
## TODO : drop year 2019 and 2018
## => DROP THIS PART
"""
item_list = get_all_subjects_from_abstract(run_folder+"/abstract")
print item_list
plot_word_evolution(item_list, run_folder)
"""

##-----------------##
## Write The Paper ##
##-----------------##

## Open manifest of FIXEME varables
FIXEME_to_value = {}
manifest = open("draft_FIXEME_table.csv", "r")
for line in manifest:
	line = line.replace("\n", "")
	line_in_array = line.split(";")
	FIXEME_to_value[line_in_array[0]] = line_in_array[1]
manifest.close()

## Create a copy of tex file with the good values
original_tex_file = open("draft.tex", "r")
output_tex_file = open("draft_ready.tex", "w")
for line in original_tex_file:

	for key in FIXEME_to_value.keys():
		line = line.replace(key, FIXEME_to_value[key])

	output_tex_file.write(line)

output_tex_file.close()
original_tex_file.close()





