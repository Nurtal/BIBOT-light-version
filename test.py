import bibot
import nltk


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


		print sentence
		

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

sementic_analysis("abstract.test")



"""
from textblob import TextBlob
blob = TextBlob("ITP is a two-year graduate program located in the Tisch School of the Arts. Perhaps the best way to describe us is as a Center for the Recently Possible.")
print blob.tags

from nltk.corpus import treebank
t = treebank.parsed_sents('wsj_0001.mrg')[0]
t.draw()
"""



