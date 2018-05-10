import bibot
import nltk
import glob
import gensim
import numpy

#bibot.plot_country_stats("meta")
#bibot.plot_articles_stats("bibot.log")
#truc = bibot.get_date_from_meta_save("meta/28867810.csv")
#print truc
#bibot.fetch_abstract(28867810)
#bibot.plot_publications_years("meta")
#bibot.write_tex_report()






def sementic_analysis(abstract_file):
	"""
	Add some smart move in this stupid
	programm

	at least, trying
	"""

	sentences_to_investigate = []

	## store abstract in a string
	abstract = open(abstract_file, "r")
	abstract_text = ""
	for line in abstract:
		abstract_text += str(line)
	abstract.close()

	## split into sentences
	sentences = abstract_text.split(". ")
	for sentence in sentences:
		words_in_sentence = sentence.split(" ")		

		## Find Subjects for each sentences
		names_found_in_abstract = []
		try:
			tokens = nltk.word_tokenize(sentence.encode('utf8'))
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
				## Something went wrong
				choucroute = True
		
		print names_found_in_abstract

		## TODO: get global field for each term retrieve,
		## i.e what are we talking about ?

#sementic_analysis("abstract.test")



"""
from textblob import TextBlob
blob = TextBlob("ITP is a two-year graduate program located in the Tisch School of the Arts. Perhaps the best way to describe us is as a Center for the Recently Possible.")
print blob.tags

from nltk.corpus import treebank
t = treebank.parsed_sents('wsj_0001.mrg')[0]
t.draw()
"""


def load_raw_documents(abstract_folder):
	"""
	Load raw documents from abstract folder
	return a list of abstract in strings
	"""

	raw_documents = {}
	abstract_files = glob.glob(str(abstract_folder)+"/*")
	for abstract in abstract_files:

		pmid = abstract.split("/")
		pmid = pmid[-1].split("_")
		pmid = pmid[0]

		abstract_text = ""
		abstract_data = open(abstract, "r")
		for line in abstract_data:
			line = line.replace("\n", "")
			abstract_text += str(line)

		raw_documents[pmid] = abstract_text

	return raw_documents



def create_similarity_map(pmid_to_abstract):
	"""
	create a similarity map between abstract
	"""

	## Extract documents
	## Get the index associated to each pmid
	raw_documents = pmid_to_abstract.values()
	pmid_to_index = {}
	for pmid in pmid_to_abstract:
		abstract = pmid_to_abstract[pmid]
		index = 0
		for document in raw_documents:
			if(abstract == document):
				pmid_to_index[pmid] = index
			index += 1
	
	## Tokenization
	from nltk.tokenize import word_tokenize
	gen_docs = [[w.lower() for w in word_tokenize(text)] for text in raw_documents]

	## Create a Dictionnary
	dictionary = gensim.corpora.Dictionary(gen_docs)

	## Create a Corpus
	corpus = [dictionary.doc2bow(gen_doc) for gen_doc in gen_docs]

	## create a tf - idf model
	tf_idf = gensim.models.TfidfModel(corpus)

	## create a similarity measure object in tf-idf space
	sims = gensim.similarities.Similarity('/home/glorfindel/Spellcraft/BIBOT-light-version/',tf_idf[corpus],num_features=len(dictionary))

	## Create a distance map
	similarity_map = {}
	for pmid in pmid_to_abstract:
		similarity_map[pmid] = {}
		abstract = pmid_to_abstract[pmid]

		query_doc = [w.lower() for w in word_tokenize(abstract)]
		query_doc_bow = dictionary.doc2bow(query_doc)
		query_doc_tf_idf = tf_idf[query_doc_bow]

		index = 0
		for sim in sims[query_doc_tf_idf]:

			pmid_to_test = -1
			for key in pmid_to_index:
				if(pmid_to_index[key] == index):
					pmid_to_test = key

			similarity_map[pmid][pmid_to_test] = sim
			index += 1

	return similarity_map


def find_smallest_distance(similarity_map):
	"""
	-> find smallest distance in similarity map, use to
	define a step for the custerring algorithm.
	"""

	smallest_distance = 1
	for pmid in similarity_map.keys():
		distances = similarity_map[pmid]
		for pmid_to_test in distances.keys():
			if(pmid_to_test != pmid and 1 - distances[pmid_to_test] < smallest_distance):
				smallest_distance = distances[pmid_to_test]
				
	return smallest_distance



def merge_abstract(pmid_list,pmid_to_abstract):
	"""
	-> Merge abstracts associated to the pmid in the pmid
	list and return a new pmid_to_abstract structure with
	the merged abstracts (and new keys for pmid : pmida_pmidb)
	"""

	new_abstract = ""
	new_pmid = ""
	new_pmid_to_abstract = {}

	for pmid in pmid_to_abstract.keys():
		if(pmid in pmid_list):
			new_abstract += str(pmid_to_abstract[pmid])+" "
		else:
			new_pmid_to_abstract[pmid] = pmid_to_abstract[pmid]

	if(new_abstract != ""):
		for item in pmid_list:
			new_pmid += str(item)+"_"
		new_pmid = new_pmid[:-1]
		new_pmid_to_abstract[new_pmid] = new_abstract

	return new_pmid_to_abstract

	

def check_intra_cluster_distance(cluster_id, abstract_folder, intra_dist_treshold):
	"""
	-> Control the variance within the cluster, make sure
	none of the article are more distant than intra_dist_treshold
	"""

	## compute similarity map
	pmid_to_abstract = load_raw_documents(abstract_folder)
	similarity_map = create_similarity_map(pmid_to_abstract)

	## look for distant articles in the same cluster
	valid_cluster = True
	pmid_list = cluster_id.split("_")
	for pmid in pmid_list:
		for pmid_to_test in pmid_list:
			if(pmid_to_test != pmid):

				if(similarity_map[pmid][pmid_to_test] <= intra_dist_treshold):
					valid_cluster = False

	return valid_cluster




def get_distance_distribution(similarity_map):
	"""
	IN PROGRESS

	-> Compute a few stats on the distance distribution,
	The idea is to get the mean distance to other articles
	for each articles, then get the std of this distribution
	and return it to set the max intra cluster distance treshold
	for the find_cluster function


	TODO:
		-> Generate graphical output
	"""

	distances_mean = []

	for pmid in similarity_map.keys():
		
		## get distance
		distances = []
		for pmid_to_test in similarity_map[pmid].keys():
			if(pmid_to_test != pmid):
				distances.append(similarity_map[pmid][pmid_to_test])

		## [TODO] : compute distribution
		## optionnal, just to show a pretty graph


		## compute mean
		mean = numpy.mean(distances)
		distances_mean.append(mean)

	## compute median of means
	## and other different stats
	distribution_median = numpy.median(distances_mean)
	distribution_std = numpy.std(distances_mean)

	return distribution_std




def find_cluster(abstract_folder):
	"""
	IN PROGRESS
	-> Find clusters, "original" algorithm

	TODO:
		-> Work on parameters initialisation, especially
		intra_dist_treshold, seems to be the only one that
		matter actually ...
		-> Work on coefficient for the std ?

	"""

	## Parameters
	number_of_cluster_treshold = 2
	min_similarity_treshold = 0.15
	intra_dist_treshold = 0.1
	something_left_to_do = True

	## Create similarity map
	pmid_to_abstract = load_raw_documents(abstract_folder)
	similarity_map = create_similarity_map(pmid_to_abstract)

	## DEBUG - TEST
	intra_dist_treshold = get_distance_distribution(similarity_map)
	print intra_dist_treshold

	iteration = 0
	while(something_left_to_do):

		## Display a message per iteration
		iteration += 1
		print "["+str(iteration)+"] => "+str(len(similarity_map.keys())) +" clusters"

		## Step 1
		## Find the smallest distance between 2 vectors
		
		print "[STEP1]"
		smallest_distance = find_smallest_distance(similarity_map)
		step = 1 - smallest_distance

		## Step 2
		## Create and extand clusters
		print "[STEP2]"
		possible_extention = False
		clusters = []
		for pmid in similarity_map.keys():
			cluster = []
			distances = similarity_map[pmid]
			for pmid_to_test in distances.keys():
				if(pmid_to_test != pmid and 1 - distances[pmid_to_test] <= step):
					#print str(pmid) +" <=> " +str(pmid_to_test) +" => " +str(distances[pmid_to_test])
					cluster = [pmid, pmid_to_test]	
			clusters.append(cluster)

		valid_clusters = []
		used_entities = []
		valid = True
		for cluster in clusters:

			cluster_id = ""
			for item in cluster:
				cluster_id += str(item)+"_"
				if(item in used_entities):
					valid = False

			cluster_id = cluster_id[:-1]
			if(valid and len(cluster) > 0):
				if(check_intra_cluster_distance(cluster_id, abstract_folder, intra_dist_treshold)):
					valid_clusters.append(cluster)
					possible_extention = True
					for item in cluster:
						used_entities.append(item)


		print "[STEP3]"
		## Merge regrouped abstract
		for cluster in valid_clusters:
			pmid_to_abstract = merge_abstract(cluster, pmid_to_abstract)

		## Compute new similarity map
		similarity_map = create_similarity_map(pmid_to_abstract)


		## Break conditions
		## break if :
		##    -> no more than x entities left (i.e reach the minimum number of clusters treshold)
		##    -> no similarity score above similarity treshold
		##    -> can't expand cluster without breaking the intra dist treshold

		if(len(similarity_map.keys()) <= number_of_cluster_treshold):
			something_left_to_do = False

		no_similarity_above_treshold = True
		for pmid in similarity_map.keys():
			similarities = similarity_map[pmid]
			for pmid_to_test in similarities:
				if(pmid != pmid_to_test and similarities[pmid] >= min_similarity_treshold):
					no_similarity_above_treshold = False

		if(no_similarity_above_treshold):
			something_left_to_do = False

			## DEBUG
			print "[OUT OF CYCLE] => No similarity above treshold"

		if(not possible_extention):
			something_left_to_do = False

			## DEBUG
			print "[OUT OF CYCLE] => No Extention Possible" 

		

	return similarity_map



def analyse_cluster(cluster_to_distance, pmid_to_abstract):
	"""
	IN PROGRESS

	## BUG THERE

	"""

	cluster_to_entities = {}
	for cluster_id in cluster_to_distance.keys():
		pmid_in_cluster = cluster_id.split("_")

		cluster_to_entities[cluster_id] = []

		for pmid in pmid_in_cluster:

			text = pmid_to_abstract[pmid]

			## get entities for this abstract

			## Tokenization
			from nltk.tokenize import word_tokenize
			
			tokenized_text = word_tokenize(text)
			tagged = nltk.pos_tag(tokenized_text)
			entities = nltk.chunk.ne_chunk(tagged)

			for elt in entities:
				if(elt[1] in ["NN", "NNP"]):
					if(elt[0] not in cluster_to_entities[cluster_id]):
						cluster_to_entities[cluster_id].append(elt[0])

		print cluster_to_entities[cluster_id]







## Playing with TF - IDF vector
abstract_folder = "SAVE/run_0h:45m:3:5/abstract"
#abstract_folder = "test/abstract"
pmid_to_abstract = load_raw_documents(abstract_folder)

simi = create_similarity_map(pmid_to_abstract)
get_distance_distribution(simi)

cluster_to_distance = find_cluster(abstract_folder)
analyse_cluster(cluster_to_distance, pmid_to_abstract)





#truc = check_intra_cluster_distance("1_2_5", "test/abstract", 0.2)
#print truc

#print pmid_to_abstract
#new_pmid_to_abstract = merge_abstract(['1','2'],pmid_to_abstract)
#print new_pmid_to_abstract













## Deal with doc creation

## Test Exemple
"""
raw_documents = ["I'm taking the show on the road.",
                 "My socks are a force multiplier.",
             "I am the barber who cuts everyone's hair who doesn't cut their own.",
             "Legend has it that the mind is a mad monkey.",
            "I make my own fun."]

"""



"""

raw_documents = pmid_to_abstract.values()
pmid_to_index = {}

for pmid in pmid_to_abstract:
	abstract = pmid_to_abstract[pmid]
	index = 0
	for document in raw_documents:
		if(abstract == document):
			pmid_to_index[pmid] = index
		index += 1



print("Number of documents:",len(raw_documents))

## Tokenization
from nltk.tokenize import word_tokenize
gen_docs = [[w.lower() for w in word_tokenize(text)] 
            for text in raw_documents]

#print(gen_docs)


## Create a Dictionnary
dictionary = gensim.corpora.Dictionary(gen_docs)
#print(dictionary[5])
#print(dictionary.token2id['road'])
print("Number of words in dictionary:",len(dictionary))
for i in range(len(dictionary)):
    print(i, dictionary[i])


## Create a Corpus
corpus = [dictionary.doc2bow(gen_doc) for gen_doc in gen_docs]
print(corpus)



## create a tf - idf model
tf_idf = gensim.models.TfidfModel(corpus)
print(tf_idf)
s = 0
for i in corpus:
    s += len(i)
print(s)


## create a similarity measure object in tf-idf space
sims = gensim.similarities.Similarity('/home/glorfindel/Spellcraft/BIBOT-light-version/',tf_idf[corpus],num_features=len(dictionary))
print(sims)
print(type(sims))


## Create a distance map
similarity_map = {}
for pmid in pmid_to_abstract:
	similarity_map[pmid] = {}
	abstract = pmid_to_abstract[pmid]

	query_doc = [w.lower() for w in word_tokenize(abstract)]
	query_doc_bow = dictionary.doc2bow(query_doc)
	query_doc_tf_idf = tf_idf[query_doc_bow]

	index = 0
	for sim in sims[query_doc_tf_idf]:

		pmid_to_test = -1
		for key in pmid_to_index:
			if(pmid_to_index[key] == index):
				pmid_to_test = key

		similarity_map[pmid][pmid_to_test] = sim

		index += 1




print similarity_map


query_doc = [w.lower() for w in word_tokenize("Socks are a force for good.")]
print(query_doc)
query_doc_bow = dictionary.doc2bow(query_doc)
print(query_doc_bow)
query_doc_tf_idf = tf_idf[query_doc_bow]
print(query_doc_tf_idf)

print sims[query_doc_tf_idf]
"""