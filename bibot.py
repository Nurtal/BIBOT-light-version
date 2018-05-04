#!/usr/bin/python
#coding: utf8 
from __future__ import unicode_literals
##==================================================================##
##=======================>BIBOTLIGHT<===============================##
##==================================================================##
## -> Light version of BIBOT, focus on pubmed fetching and natural  ## 
## langage processing for the selected articles.                    ##
## Current version is 1.0.											##
## Because reviewing is not funny enough.                           ##
##==================================================================##


##-------------##
## IMPORTATION #######################################################
##-------------##
from Bio import Entrez
from Bio.Entrez import efetch, read
from unidecode import unidecode
import nltk
import itertools
import os
import time
import shutil
import datetime
import glob
import getopt
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


##------------##
## PARAMETERS ########################################################
##------------##
Entrez.email = 'murlock.raspberypi@gmail.com'




##-----------##
## FUNCTIONS #########################################################
##-----------##


def fetch_abstract(pmid):
	##
	## Return abstract of a given
	## article using pmid
	##
	## => Return None when pmid can't be return
	## (can happen when article is in chinese)
	##
	
	handle = efetch(db='pubmed', id=pmid, retmode='xml', )
	xml_data = read(handle)
	xml_data = xml_data['PubmedArticle'][0]

	try:
		article = xml_data['MedlineCitation']['Article']
		abstract = article['Abstract']['AbstractText'][0]
		return abstract
	except IndexError:
		return None
	except KeyError:
		return None


def get_ListOfArticles(term, max_count):
	##
	## return the list of pmid article conatining
	## the term.
	##

	h = Entrez.esearch(db='pubmed', retmax=max_count, term=term)
	result = Entrez.read(h)
	listOfArticles = result["IdList"]

	return listOfArticles;



def save_abstract(abstract, save_file):
	##
	## -> Save the abstract in a text file
	## convert the abstract to unicode.
	##

	## preprocess abstract
	#abstract_preprocess = unicode(abstract)
	abstract_preprocess = abstract.encode('utf8')

	## save abstract in file
	output = open(save_file, "w")
	output.write(abstract_preprocess)
	output.close()


def load_text(text_file):
	##
	## -> Create  and return nltk Text 
	## object from text_file
	##

	## -> Create the Text object from the input
	## text file
	text_file=open(text_file,'rU')
	raw=text_file.read()
	raw = raw.decode('utf8')
	tokens = nltk.word_tokenize(raw)
	text = nltk.Text(tokens)

	## Return the nltk text object
	return text


def evaluate_article(pmid):
	##
	## [IN PROGRESS]
	##
	## -> Test if the abstract is cool
	## -> return true or false
	##
	## TODO : write doc
	##

	##------------------------##
	## Parameters for filters ##
	##------------------------##
	

	## initialize parameters
	oldest_year_authorized = "NA"
	authorized_languages = []
	valid_article = False
	check_date = True
	check_language = True
	validation_check = {}
	validation_keywords = {}

	exclusion_check = {}
	exclusion_keywords = {}

	exclusion_keywords_found = False

	## test if config file exist
	if(os.path.isfile("config.conf")):
		config_data = open("config.conf", "r")
		validation_keywords_cmpt = 0
		exclusion_keywords_cmpt = 0
		for line in config_data:
			line = line.replace("\n", "")
			line_in_array = line.split(";")

			if(line_in_array[0] == "min year"):
				oldest_year_authorized = line_in_array[1]
			elif(line_in_array[0] == "authorized languages"):
				languages_list = line_in_array[1].split(",")
				for elt in languages_list:
					authorized_languages.append(unicode(elt))
			elif(line_in_array[0] == "validation keywords"):
				validation_keywords_cmpt += 1
				validation_check["keywords_"+str(validation_keywords_cmpt)] = False
				validation_keywords["keywords_"+str(validation_keywords_cmpt)] = []
				keywords_list = line_in_array[1].split(",")
				for elt in keywords_list:
					if(elt not in validation_keywords["keywords_"+str(validation_keywords_cmpt)]):
						validation_keywords["keywords_"+str(validation_keywords_cmpt)].append(str(elt))

			## Retrieve Exclusion list
			elif(line_in_array[0] == "exclusion keywords"):
				exclusion_keywords_found = True
				exclusion_keywords_cmpt += 1
				exclusion_check["exclusion_"+str(exclusion_keywords_cmpt)] = False
				exclusion_keywords["exclusion_"+str(exclusion_keywords_cmpt)] = []
				keywords_list = line_in_array[1].split(",")
				for elt in keywords_list:
					if(elt not in exclusion_keywords["exclusion_"+str(exclusion_keywords_cmpt)]):
						exclusion_keywords["exclusion_"+str(exclusion_keywords_cmpt)].append(str(elt))

		config_data.close()

	## default configuration
	else:
		oldest_year_authorized = 2008
		authorized_languages = [u'eng']
		validation_check["keywords_1"] = False
		validation_check["keywords_2"] = False
		validation_keywords["keywords_1"]= ["algorithm", "machine" "learning", "neural", "network", "statistic", "deep", "classification", "model"]
		validation_keywords["keywords_2"] = ["Sjogren" ,"sjogren", "lupus", "autoimmunity", "rhumatoid", "arthrisis", "RA", "SjS", "SLE"]
		validation_check["exclusion_1"] = False
		exclusion_keywords["exclusion_1"]= []

	if(not exclusion_keywords_found):
		validation_check["exclusion_1"] = False
		exclusion_keywords["exclusion_1"]= []

	##---------------##
	## The Easy Part ##
	##---------------##
	## get meta data on the articles
	try:
		handle = efetch(db='pubmed', id=pmid, retmode='xml', )
		informations = read(handle)
		stuff = informations[u'PubmedArticle'][0] 
		
		## get date from the history attribute, select
		## the date of acceptation.
		date = stuff[u'PubmedData']["History"][1]
		month = date[u'Month']
		day = date[u'Day']
		year = date[u'Year']

		## get the name of the review
		journal_name = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'MedlineJournalInfo'][u'MedlineTA']
		
		## get the keywords for the articles
		## the format is a bit strange, may have to be carreful
		## with this data (mix of strings and unicode elements)
		keywords_list = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'KeywordList']

		## Get the author's conflict of interest,
		## because we can.
		try:
			conflict_of_interest = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'CoiStatement']
		except:
			conflict_of_interest = "NA"

		## Get title of the article
		article_title = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'Article'][u'ArticleTitle']

		## Get language of the article
		article_language = informations[u'PubmedArticle'][0][u'MedlineCitation'][u'Article'][u'Language'][0]

		## Get country of publications
		country = stuff[u'MedlineCitation'][u'MedlineJournalInfo'][u'Country']
		
	except:
		return (False,False,False)

	##----------------##
	## The Smart Part ## 
	##----------------##
	## run further analysis on the abstract using nltk

	##
	## WORKING ON EXCLUSION LIST
	##

	## fetch the abstract and convert it to
	## a nltk text object.
	abstract_file_name = "abstract/"+str(pmid)+"_abstract.txt"
	abstract = fetch_abstract(pmid)
	if(abstract):
		save_abstract(abstract, abstract_file_name)
		abstract_text = load_text(abstract_file_name)
		
		## Play with tokenization and chunking
		## Get all the commun names in the abstract
		names_found_in_abstract = []
		try:
			tokens = nltk.word_tokenize(abstract.encode('utf8'))
			tagged = nltk.pos_tag(tokens)
			entities = nltk.chunk.ne_chunk(tagged)
		except:
			print "[WARNINGS] => can't perform nlp operation"
			entities = []

		for item in entities:
			try:
				if(item[1] in ["NN", "NNS", "NNP"]):
					if(item[0] not in names_found_in_abstract):
						names_found_in_abstract.append(item[0])
			except:
				## Somethig went wrong
				choucroute = True
		
		## Check validation list				
		for item in names_found_in_abstract:
			for key in validation_keywords.keys():
				keywords_validation_list = validation_keywords[key]
				if(item in keywords_validation_list):
					validation_check[key] = True

		## Check exclusion list
		for item in names_found_in_abstract:
			for key in exclusion_keywords.keys():
				exclusion_validation_list = exclusion_keywords[key]
				if(item in exclusion_validation_list):
					exclusion_check[key] = True


	##--------------##
	## PASS OR FAIL ##
	##--------------##
	## General check phase
	easy_check_passed = False
	smart_check_passed = True

	## Basic check on meta data
	## - check date
	if(int(year) < int(oldest_year_authorized)):
		check_date = False

	## - check language
	if(article_language not in authorized_languages):
		check_language = False

	## Easy Filter
	if(check_date and check_language):
		easy_check_passed = True

	## Complex filter (inclusion)
	if(False in validation_check.values()):
		smart_check_passed = False

	## Complex filter (exclusion)
	if(True in exclusion_check.values()):
		smart_check_passed = False




	## Global check
	if(easy_check_passed and smart_check_passed):
		valid_article = True

	##-------------##
	## SAVING DATA ##
	##-------------##
	## Write and delete files
	if(valid_article):

		## Save meta data in a text file
		## for further use
		title_line = u'>Title;'+unicode(article_title)+u"\n"
		date_line = u'>Date;'+unicode(day)+u"/"+unicode(month)+u"/"+unicode(year)+u"\n"
		#date_line = '>Date;'+str(day.encode('utf8'))+"/"+str(month.encode(utf8))+"/"+str(year.encode("utf8"))+"\n"
		journal_line = u">Journal;"+unicode(journal_name)+u"\n"
		country_line = u">Country;"+unicode(country)+u"\n"
		conflict_of_interest_line = u">Conflict;"+unicode(conflict_of_interest)+u"\n"
		meta_data = open("meta/"+str(pmid)+".csv", "w")
		meta_data.write(title_line.encode('utf8'))
		meta_data.write(date_line.encode('utf8'))
		meta_data.write(journal_line.encode('utf8'))
		meta_data.write(country_line.encode('utf8'))
		meta_data.write(conflict_of_interest_line.encode('utf8'))
		meta_data.close()

	else:
		## Delete the abstract
		try:
			if(abstract):
				os.remove(abstract_file_name)
		except:
			print "[WARNING] => can't delete "+str(abstract_file_name)

	##------------------##
	## RETURN SOMETHING ##
	##------------------##
	## return True if the article pass the 
	## evaluation, else False.
	return (valid_article, easy_check_passed, smart_check_passed)



def get_huge_list_of_artciles(keywords):
	##
	## Create all possible comination of at least two elements
	## from the keywords list. then use these combination
	## to create request (using only the AND operator for now)
	## and screen the pubmed database.
	##
	## return the list of all articles found
	##

	## init variables
	huge_list_of_PMID = []

	## make all combination of at least 2 item in keywords
	combination_list = []
	for x in xrange(2, len(keywords)):
		machin = itertools.combinations(keywords, x)
		for truc in machin:
			combination_list.append(list(truc))

	## create request
	for items_set in combination_list:
		request = ""
		for item in items_set:
			request += item +" AND "
		request = request[:-5]
		
		## screening pubmed
		screened = False
		while(not screened):
			try:
				results_PMID = get_ListOfArticles(request, 9999999)
				screened = True
			except:
				time.sleep(1)
		## increment huge list of pmid
		for pmid in results_PMID:
			if(pmid not in huge_list_of_PMID):
				huge_list_of_PMID.append(pmid)

	return huge_list_of_PMID



def run(request_term):
	##
	## main function, run the bibot programm
	##

	## Dispaly Run information
	print "[INFO] PREPARE FOR RUN"

	## Clean absract and meta folder
	print "[INFO] Cleaning directories"
	for abstract_file in glob.glob("abstract/*.txt"):
		os.remove(abstract_file)
	for meta_data in glob.glob("meta/*.csv"):
		os.remove(meta_data)

	## variables and file initialisation
	print "[INFO] Initialize log file"
	log_file = open("bibot.log", "w")
	
	## Save request term in log file
	print "[INFO] TERMS : "
	request_line = ""
	for item in request_term:
		print "[-]\t"+str(item)
		request_line += str(item)+";"
	request_line = request_line[:-1]
	log_file.write(request_line+"\n")

	## First interrogation of medline, get a
	## huge list of articles possibly relevant and
	## write the numbers of article in lig file
	print "[INFO] Large screening"
	big_list = get_huge_list_of_artciles(request_term)
	Total_number_of_articles = len(big_list)
	log_file.write("Total_number_of_articles;"+str(Total_number_of_articles)+"\n")
	print "[INFO] "+str(Total_number_of_articles) +" articles found"

	## Test each articles retrieved from their pmid
	fetched = 0
	first_fiter_passed = 0
	last_filter_passed = 0
	cmpt = 0

	for article in big_list:

		## try to evaluate the article
		## require a connection to the
		## NCBI Server, if succed go on, 
		## else wait 5 seconds and try again
		article_is_evaluated = False
		while(not article_is_evaluated):
			try:
				#print "|| TRY TO PROCESS ARTICLE "+str(article)+ " ||"
				valid = evaluate_article(article)
				article_is_evaluated = True
			except:
				#print "|| CAN'T REACH NCBI, WAIT FOR 5 SECONDS ||"
				print "[INFO] => CAN'T REACH NCBI, WAIT FOR 5 SECONDS "
				now = datetime.datetime.now()
				time_tag = str(now.hour)+"h:"+str(now.minute)+"m:"+str(now.day)+":"+str(now.month)
				log_file.write("["+str(time_tag)+"];can't reach NCBI, wait for 5 seconds\n")
				time.sleep(5)


		filter_1_status = "FAILED"
		filter_2_status = "FAILED"
		if(valid[0]):
			fetched += 1
		if(valid[1]):
			first_fiter_passed += 1
			filter_1_status = "PASSED"
		if(valid[2]):
			last_filter_passed += 1
			filter_2_status = "PASSED"
		cmpt += 1

		## Display progress on screen
		## and status for each pmid in log file
		print "[RUN] => "+str(cmpt) +" [PROCESSED] || "+ str(fetched) + " [SELECTED] || FIRST FILTERS ["+filter_1_status+ "] || LAST FILTER ["+filter_2_status+ "] || "+str(float((float(cmpt)/float(Total_number_of_articles))*100)) + "% [COMPLETE]"
		log_file.write(">"+str(article)+";First_Filter="+str(filter_1_status)+";Last_filter="+str(filter_2_status)+"\n")

	## close log file
	log_file.close()

	## Save the results files
	## and folders
	now = datetime.datetime.now()
	time_tag = str(now.hour)+"h:"+str(now.minute)+"m:"+str(now.day)+":"+str(now.month)
	abstract_destination  = "SAVE/run_"+str(time_tag)+"/abstract"
	meta_destination = "SAVE/run_"+str(time_tag)+"/meta"
	log_destination = "SAVE/run_"+str(time_tag)+"/bibot.log"
	shutil.copytree("abstract", abstract_destination)
	shutil.copytree("meta", meta_destination)
	shutil.copy("bibot.log", log_destination)



def check_request_terms(request_terms):
	##
	## check request terms provide to the
	## script as an args.
	##
	## -> can be a file containing a semi-col delimited
	## list of terms
	##
	## -> can be a string containing a semi-col delimited
	## list of terms
	##
	## return the list of terms or NA if nothing found.
	##

	request_term_parsed = False
	terms = []

	request_terms = request_terms.decode("utf-8")
	print str(request_terms)

	## check if file exist
	if(os.path.isfile(request_terms)):
		data = open(request_terms, "r")
		for line in data:
			line = line.replace("\n", "")
			line_in_array = line.split(";")

			for elt in line_in_array:
				if(elt not in terms):
					terms.append(str(elt))
		data.close()

		if(len(terms) > 0):
			request_term_parsed = True

	## check if it's semi column seprated list
	elif(";" in str(request_terms)):
		line_in_array = request_terms.split(";")
		for elt in line_in_array:
			if(elt not in terms):
				terms.append(str(elt))

		request_term_parsed = True

	## return a list of terms
	if(request_term_parsed):
		return terms
	else:
		return "NA"



def check_config_file(config_file):
	##
	## check if config_file exist and countain
	## all the required parameters.
	##
	## copy file to config.conf if all is good
	##
	## return a string:
	##			- good
	##			- not a file
	##			- incomplete configuration file
	##

	status = "NA"
	parameter_check = {}
	parameter_check["min year"] = False
	parameter_check["authorized languages"] = False
	parameter_check["validation keywords"] = False

	## check if file exist
	if(os.path.isfile(config_file)):

		## check if all parameters are present
		config_data = open(config_file, "r")
		for line in config_data:
			line = line.replace("\n", "")
			line_in_array = line.split(";")

			if(line_in_array[0] in parameter_check.keys()):
				parameter_check[line_in_array[0]] = True

		config_data.close()

		if(False in parameter_check):
			status = "incomplete configuration file"
		else:
			# rename file if everything is ok
			shutil.copy(config_file, "config.conf")
			status = "good"
	
	else:
		status = "not a file"

	return status


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
		try:
			line = line.decode('utf8')
		except:
			print "Something went wrong"
		line = line.replace("\n", "")
		if(line[0] == ">"):

			line_in_array = line.split(";")
			if(line_in_array[0] == ">Date"):
				date_in_array = line_in_array[1].split("/")
				year = date_in_array[2]

	meta_data.close()

	## return only the year of publication
	return year

def plot_publications_years(meta_data_folder):
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


def plot_country_stats(meta_data_folder):
	"""
	Get the country stat from meta data
	and create a pie chart with these number
	"""

	## init structure
	country_to_count = {}

	## Get data
	meta_file_list = glob.glob(meta_data_folder+"/*.csv")
	for meta_file in meta_file_list:

		meta_data = open(meta_file, "r")
		for line in meta_data:
			line = line.rstrip()
			try:
				line_in_array = line.split(";")
				if(line_in_array[0] == ">Country"):
					country = line_in_array[1]
					if(country not in country_to_count.keys()):
						country_to_count[country] = 1
					else:
						country_to_count[country] += 1
			except:
				print "[WARNING] => Unable to recover country for "+str(meta_file)
			
		meta_data.close()

	## Generate Pie chart 
	plt.figure(1, figsize=(6,6))
	labels = country_to_count.keys()
	fracs = country_to_count.values()
	plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=False, startangle=90)
	plt.savefig("images/country_repartition.png")
	plt.close()


def plot_articles_stats(log_file):
	"""
	IN PROGRESS
	
	Todo : debug shape of pie
	"""

	## Get data
	first_filter_pass_cmpt = 0
	second_filter_pass_cmpt = 0
	pass_both_filter_cmpt = 0
	article_cmpt = 0

	log_data = open(log_file, "r")
	for line in log_data:
		line = line.replace("\n", "")
		line_in_array = line.split(";")
		if(line[0] == ">"):

			article_cmpt += 1
			
			filter_1_info = line_in_array[1].split("=")
			if(filter_1_info[1] == "PASSED"):
				first_filter_pass_cmpt += 1

			filter_2_info = line_in_array[2].split("=")
			if(filter_2_info[1] == "PASSED"):
				second_filter_pass_cmpt += 1


			if(filter_1_info[1] == "PASSED" and filter_2_info[1] == "PASSED"):
				pass_both_filter_cmpt += 1

	log_data.close()

	## Generate pie for data
	fig, (ax1) = plt.subplots(1, ncols=3)

	# plot each pie chart in a separate subplot
	first_filter_failed_cmpt = article_cmpt - first_filter_pass_cmpt
	second_filter_failed_cmpt = article_cmpt - second_filter_pass_cmpt
	failed_one_filter_cmpt = article_cmpt - pass_both_filter_cmpt 
	
	explode = (0, 0.05)

	ax1[0].pie([first_filter_pass_cmpt, first_filter_failed_cmpt], radius=0.5, labels=["Pass", "Failed"], autopct='%1.1f%%', startangle=90)
	ax1[1].pie([second_filter_pass_cmpt, second_filter_failed_cmpt], labels=["Pass", "Failed"], autopct='%1.1f%%', startangle=90)
	ax1[2].pie([pass_both_filter_cmpt, failed_one_filter_cmpt], labels=["Pass", "Failed"], autopct='%1.1f%%', startangle=90)
	
	#plt.show()
	plt.savefig("images/plot_stats.png")
	plt.close()



def write_tex_report():
	"""
	-> Parse information from the log file
	-> Generate figures from the meta data
	-> Write a .tex file

	TODO : - Improve template
		   - Fixe date
	""" 

	## Generate figures
	plot_publications_years("meta")
	plot_country_stats("meta")
	plot_articles_stats("bibot.log")

	## get informations
	keywords_list = []
	total_articles = -1
	selected_articles = 0

	log_file = open("bibot.log", "r")
	cmpt = 0
	for line in log_file:
		line = line.replace("\n", "")
		if(cmpt == 0):
			keywords_list = line.split(";")
		
		elif(cmpt == 1):
			line_in_array = line.split(";")
			total_articles = line_in_array[1]

		elif(line[0] == ">"):

			line_in_array = line.split(";")
			filter_1_info = line_in_array[1].split("=")
			filter_2_info = line_in_array[2].split("=")
			
			if(filter_1_info[1] == "PASSED" and filter_2_info[1] == "PASSED"):
				selected_articles += 1

		cmpt += 1

	log_file.close()

	## create report file
	report_file = open("report.tex", "w")

	## write header
	report_file.write("\\documentclass[a4paper,9pt]{extarticle}\n")
	report_file.write("\\usepackage[utf8]{inputenc}\n")
	report_file.write("\\usepackage[T1]{fontenc}\n")
	report_file.write("\\usepackage{graphicx}\n")
	report_file.write("\\usepackage{xcolor}\n")
	report_file.write("\\usepackage{amsmath,amssymb,textcomp}\n")
	report_file.write("\\everymath{\displaystyle}\n")
	report_file.write("\\usepackage{times}\n")
	report_file.write("\\renewcommand\\familydefault{\sfdefault}\n")
	report_file.write("\\usepackage{tgheros}\n")
	report_file.write("\\usepackage[defaultmono,scale=0.85]{droidmono}\n")
	report_file.write("\\usepackage{multicol}\n")
	report_file.write("\\setlength{\columnseprule}{0pt}\n")
	report_file.write("\\setlength{\columnsep}{20.0pt}\n")
	report_file.write("\\usepackage{geometry}\n")
	report_file.write("\\geometry{\n")
	report_file.write("    a4paper,\n")
	report_file.write("    total={210mm,297mm},\n")
	report_file.write("    left=10mm,right=10mm,top=10mm,bottom=15mm}\n")
	report_file.write("\\linespread{1.3}\n")
	report_file.write("\\makeatletter\n")
	report_file.write("\\renewcommand*{\\maketitle}{%\n")
	report_file.write("\\noindent\n")
	report_file.write("\\begin{minipage}{0.4\\textwidth}\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=6pt,inner sep=10pt,fill=blue!50!black,text width= 0.75\\textwidth] {\\color{white}\\Huge \\@title};\n")
	report_file.write("    \\end{tikzpicture}\n")
	report_file.write("\\end{minipage}\n")
	report_file.write("\\hfill\n")
	report_file.write("\\begin{minipage}{0.55\\textwidth}\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=3pt,inner sep=10pt,draw=blue!50!black,text width= 0.95\\textwidth] {\LARGE \\@author};\n")
	report_file.write("    \\end{tikzpicture}\n")
	report_file.write("\\end{minipage}\n")
	report_file.write("\\bigskip\\bigskip\n")
	report_file.write("}%\n")
	report_file.write("\\makeatother\n")
	report_file.write("\\usepackage[explicit]{titlesec}\n")
	report_file.write("\\newcommand*\\sectionlabel{}\n")
	report_file.write("\\titleformat{\\section}\n")
	report_file.write("{\\gdef\\sectionlabel{}\n")
	report_file.write("\\normalfont\\sffamily\\Large\\bfseries\\scshape}\n")
	report_file.write("{\\gdef\\sectionlabel{\\thesection\\ }}{0pt}\n")
	report_file.write("{\n")
	report_file.write("    \\noindent\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=3pt,inner sep=4pt,fill=blue!50!black,text width= 0.95\\columnwidth] {\\color{white}\\sectionlabel#1};\n")
	report_file.write("    \\end{tikzpicture}\n")		
	report_file.write("}\n")
	report_file.write("\\titlespacing*{\section}{0pt}{15pt}{10pt}\n")	
	report_file.write("\\usepackage{fancyhdr}\n")
	report_file.write("\\makeatletter\n")
	report_file.write("\\pagestyle{fancy}\n")
	report_file.write("\\fancyhead{}\n")
	report_file.write("\\fancyfoot[C]{\\footnotesize \\@date\\ \\ \\@author}\n")	
	report_file.write("\\renewcommand{\\headrulewidth}{0pt}\n")
	report_file.write("\\renewcommand{\\footrulewidth}{0pt}\n")
	report_file.write("\\makeatother\n")
	report_file.write("\\usepackage{tikz}\n")
	report_file.write("\\usetikzlibrary{shapes,arrows}\n")
	report_file.write("\\usepackage{xcolor}\n")
	report_file.write("\\usepackage{chronosys}\n")
	report_file.write("\\usepackage{pgfplots}\n")
	report_file.write("\\usepackage{amsmath,amssymb,textcomp}\n")
	report_file.write("\\everymath{\\displaystyle}\n")

	## write title
	report_file.write("\\title{BIBOT Report}\n")
	report_file.write("\\author{\n")
	report_file.write("\\begin{center}\n")
	report_file.write("   \\vspace{-0.5cm}\n")
	report_file.write("   FIXE DATE\n")
	report_file.write("\\end{center}\n")
	report_file.write("}\n")
	report_file.write("\\date{}\n")

	## write document
	report_file.write("\\begin{document}\n")
	report_file.write("\\maketitle\n")
	report_file.write("\\section{Overview}\n")
	report_file.write("\\subsection{Keywords}\n")

	disposition = "{" + "c" * len(keywords_list) +"}"
	tabular_line = ""
	for keyword in keywords_list:
		tabular_line += keyword + "&"
	tabular_line = tabular_line[:-1]

	report_file.write("\\begin{center}\n")
	report_file.write("\\begin{tabular}"+str(disposition)+"\n")
	report_file.write(str(tabular_line)+"\\\\\n")
	report_file.write("\\end{tabular}\n")
	report_file.write("\\end{center}\n")
	report_file.write("\\subsection{Distributions}\n")
	report_file.write("\\begin{center}\n")
	report_file.write("\\begin{tabular}{cc}\n")
	report_file.write("\\includegraphics[scale=0.5]{images/country_repartition} & \\includegraphics[scale=0.5]{images/years_publications_evolution}\n")
	report_file.write("\\end{tabular}\n")
	report_file.write("\\end{center}\n")
	
	report_file.write("\\section{Filters}\n")
	report_file.write("\\begin{itemize}\n")
	report_file.write("\\item "+str(total_articles)+" articles found\n")
	report_file.write("\\item "+str(selected_articles)+" articles selected\n")
	report_file.write("\\end{itemize}\n")

	report_file.write("\\begin{center}\n")
	report_file.write("\\includegraphics[scale=0.5]{images/plot_stats}\n")
	report_file.write("\\end{center}\n")
	report_file.write("\\end{document}\n")

	## close report file
	report_file.close()



def main(argv):
	##
	## The main function, called when the script
	## is executed.
	##
	## -> parse command line arguments and run
	## 	  BIBOT
	##


	## BIBOT Header
	print "===== BIBOT LIGHT VERSION 1.0 ====="
	print "|       ___  _ ___  ____ ___       |"
	print "|       |__] | |__] |  |  |        |"
	print "|       |__] | |__] |__|  |        |"
	print "|                                  |"
	print "==================================="
                     



	
	## Arguments list
	## init a few values and the help_content
	## variable.
	action = ""
	request_terms = "NA"
	config_file = "NA"
	help_content = "Alpha version\nUsage : bibot.py -a <action> -r <request> -c <configuration>\n"
	readme = open("README.md", "r")
	for line in readme:
		help_content += line
	readme.close

	## check if subfolder exists
	if(not os.path.isdir("abstract")):
		os.mkdir("abstract")
	if(not os.path.isdir("meta")):
		os.mkdir("meta")


	## Retrieve arguments from command line
	## and display help
	## help -h
	## action -a <action>
	## request term -r <rterms>
	## configuration -c <conf>
	try:
		opts, args = getopt.getopt(argv,"ha:r:c:",["action=", "rterms=", "conf="])

	except getopt.GetoptError:
		print help_content
		sys.exit(2)

	for opt, arg in opts:
		
		## Display Help
		if opt == '-h':
			print help_content
			sys.exit()

		## Get action
		elif opt in ("-a", "--action"):
			action = arg
		
		## Get request term
		elif opt in ("-r", "--rterms"):
			request_terms = check_request_terms(arg)
		
		## Get config file (not mandatory)
		elif opt in ("-c", "--conf"):
			config_file = arg
			check = check_config_file(config_file)
			if(check == "not a file"):
				print "[ERROR] => Can't find file "+str(config_file)
				sys.exit(2)
			elif(check == "incomplete configuration file"):
				print "[ERROR] => Missing parameters in configuration file "+str(config_file)	
				sys.exit(2)

	## Run BIBOT
	if(action == "run"):
		if(request_terms != "NA"):
			print "####################################"
			print "============== RUN ================#"
			print "####################################"
			run(request_terms)
		else:
			print "[ERROR] => can't find request terms. Use -r <rterms> option, see help (-h) for details."
			sys.exit(2)

	elif(action == "debug"):
		print "####################################"
		print "==== DEBUG TRACE ==================#"
		print "####################################"
		print "=> Variables :"
		print "<action> => " +str(action)
		print "<rterms> => " +str(request_terms)
		print "<conf> => "+str(config_file)
		print "<help> => " +str(help_content)




##------##
## MAIN ###################################################################
##------##

if __name__ == '__main__':
	
	main(sys.argv[1:])

