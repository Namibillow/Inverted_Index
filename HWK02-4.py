# Nami AKazawa HWK2 Inverted Index

import unicodedata
import functools
import zipfile
import re
from bs4 import BeautifulSoup
from collections import defaultdict
import pprint
from math import log
import math
from collections import OrderedDict
import time
import itertools

start_time = time.time()

_WORD_MIN_LENGTH = 3
_STOP_WORDS = frozenset([
'a', 'about', 'above', 'above', 'across', 'after', 'afterwards', 'again',
'against', 'all', 'almost', 'alone', 'along', 'already', 'also','although',
'always','am','among', 'amongst', 'amoungst', 'amount',  'an', 'and', 'another',
'any','anyhow','anyone','anything','anyway', 'anywhere', 'are', 'around', 'as',
'at', 'back','be','became', 'because','become','becomes', 'becoming', 'been',
'before', 'beforehand', 'behind', 'being', 'below', 'beside', 'besides',
'between', 'beyond', 'bill', 'both', 'bottom','but', 'by', 'call', 'can',
'cannot', 'cant', 'co', 'con', 'could', 'couldnt', 'cry', 'de', 'describe',
'detail', 'do', 'done', 'down', 'due', 'during', 'each', 'eg', 'eight',
'either', 'eleven','else', 'elsewhere', 'empty', 'enough', 'etc', 'even',
'ever', 'every', 'everyone', 'everything', 'everywhere', 'except', 'few',
'fifteen', 'fify', 'fill', 'find', 'fire', 'first', 'five', 'for', 'former',
'formerly', 'forty', 'found', 'four', 'from', 'front', 'full', 'further', 'get',
'give', 'go', 'had', 'has', 'hasnt', 'have', 'he', 'hence', 'her', 'here',
'hereafter', 'hereby', 'herein', 'hereupon', 'hers', 'herself', 'him',
'himself', 'his', 'how', 'however', 'hundred', 'ie', 'if', 'in', 'inc',
'indeed', 'interest', 'into', 'is', 'it', 'its', 'itself', 'keep', 'last',
'latter', 'latterly', 'least', 'less', 'ltd', 'made', 'many', 'may', 'me',
'meanwhile', 'might', 'mill', 'mine', 'more', 'moreover', 'most', 'mostly',
'move', 'much', 'must', 'my', 'myself', 'name', 'namely', 'neither', 'never',
'nevertheless', 'next', 'nine', 'no', 'nobody', 'none', 'noone', 'nor', 'not',
'nothing', 'now', 'nowhere', 'of', 'off', 'often', 'on', 'once', 'one', 'only',
'onto', 'or', 'other', 'others', 'otherwise', 'our', 'ours', 'ourselves', 'out',
'over', 'own','part', 'per', 'perhaps', 'please', 'put', 'rather', 're', 'same',
'see', 'seem', 'seemed', 'seeming', 'seems', 'serious', 'several', 'she',
'should', 'show', 'side', 'since', 'sincere', 'six', 'sixty', 'so', 'some',
'somehow', 'someone', 'something', 'sometime', 'sometimes', 'somewhere',
'still', 'such', 'system', 'take', 'ten', 'than', 'that', 'the', 'their',
'them', 'themselves', 'then', 'thence', 'there', 'thereafter', 'thereby',
'therefore', 'therein', 'thereupon', 'these', 'they', 'thickv', 'thin', 'third',
'this', 'those', 'though', 'three', 'through', 'throughout', 'thru', 'thus',
'to', 'together', 'too', 'top', 'toward', 'towards', 'twelve', 'twenty', 'two',
'un', 'under', 'until', 'up', 'upon', 'us', 'very', 'via', 'was', 'we', 'well',
'were', 'what', 'whatever', 'when', 'whence', 'whenever', 'where', 'whereafter',
'whereas', 'whereby', 'wherein', 'whereupon', 'wherever', 'whether', 'which',
'while', 'whither', 'who', 'whoever', 'whole', 'whom', 'whose', 'why', 'will',
'with', 'within', 'without', 'would', 'yet', 'you', 'your', 'yours', 'yourself',
'yourselves', 'the'])

def parsing_texts(text):
	soup = BeautifulSoup(text,'html.parser')

	#body = soup.find_all('body')

	#links = []
	#for b in body:
		#links.append(b.find('href'))
	#print(links)

	text = [s.extract() for s in soup(['style', 'script','meta','[document]', 'head','title'])]
	text = soup.getText()
	text = text.replace('\n',' ')
	text = " ".join(text.split())
	text = text.lower()
	#print(text)
	return text

def word_split(text):
	"""
    Split a text in words. Returns a list of tuple that contains
    (word, location) location is the starting byte position of the word.
    """
	word_current = []
	word_list = []
	word_index = None
	for i, curr in enumerate(text):
		if curr.isalnum():
			word_current.append(curr)
			word_index = i

		elif word_current:
			word = u''.join(word_current)
			word_list.append((word_index - len(word) +1, word))
			word_current = []

	if word_current:
		word = u''.join(word_current)
		word_list.append((word_index - len(word) +1, word))

	return word_list

def words_eliminate(words):
	"""
	Remove words with length less then a minimum and stopwords.
	"""
	cleaned_words = []
	for index,word in words:
		if len(word) < _WORD_MIN_LENGTH or word in _STOP_WORDS:
			continue
		cleaned_words.append((index, word))
	return cleaned_words

def words_lower(words):
	'''
	lowering the words
	'''
	normalized_words = []
	for index, word in words:
		w_lowered = word.lower()
		normalized_words.append((index, w_lowered))

	return normalized_words

def word_index(text):
	"""
	Just a helper method to process a text.
	It calls word split, normalize and cleanup.
	"""
	words = word_split(text)
	words = words_lower(words)
	words = words_eliminate(words)

	return words

def inverted_index(text):
	"""
	Create an Inverted-Index of the specified text document.
	    {word:[locations]}
	"""
	inverted = {}

	for index, word in word_index(text):
		locations = inverted.setdefault(word,[])
		locations.append(index)

	return inverted


def inverted_index_add(inverted, doc_id, doc_index, f_f):
	"""
	Add Invertd-Index doc_index of the document doc_id to the
	Multi-Document Inverted-Index (inverted),
	using doc_id as document identifier.
	    {word:{doc_id:[locations]}}
	"""
	for word, locations in doc_index.items():
		indices = inverted.setdefault(word,{})
		indices[doc_id] = locations


		f_f[word][doc_id] = len(locations)

	return inverted, f_f

def search_for_and(inverted, query):
	"""
    Returns a set of documents id that contains all the words in your query.
    """
	words = [word for _, word in word_index(query) if word in inverted]
	results = [set(inverted[word].keys()) for word in words]

	###  is used for taking the intersection of two
	### Intersection
	return functools.reduce(lambda x, y: x & y, results) if results else []

def search_for_or(inverted, query):
	"""
	Returns a set of documents id that contains all the words in your query.
	"""
	words = [word for _, word in word_index(query) if word in inverted]
	results = [set(inverted[word].keys()) for word in words]

	return results
    #set.union(*results)

def check_word_validity(inverted,query, ans):
	count = 0
	words = query.split()

	for word in words:
		if word in inverted:
			count +=1

	if ans == 'or':
		if count >= 1:
			return True
	elif ans == 'and':
		if count == len(words):
			return True
	return False

def search_for_phrase(inverted, query, documents):
	### this will get us set of intersectioned files
	words_doc = search_for_and(inverted, query)
	total_freq = 0 ### For the Dictinary Frequency
	#print("words_doc = ", words_doc)
	for doc in words_doc:
		freq = 0 ### For the file frequency
		first_search = 0 ### For counting the total_freq
		text = documents[doc].lower()
		words = query.split()

		for location in inverted[words[0]][doc]:
			l = text.find(query, location,(len(query)+location+1))
			first_search += 1
			#print(l)
			if l == -1:
				continue
			else:
				if(first_search == 1):
					total_freq+= 1
				freq +=1
				print("The file is: ", doc)
				print("Index is: [%d] - [%d]" %(l, l+len(query)-1))
				print('		- %s...' % documents[doc][l:l+40].replace('\n',' '))
		if freq >= 1:
			print("Total frequency of '%s' in %s is: %d \n" % (query, doc, freq ))
	if total_freq >= 1:
		print("This phrase appears %d time(s) in files.\n"% total_freq)


### Step 1: Term Frequency { file: {word: tf, word2: tf}, file2: {word: tf, word2: tf}}
def term_frequency(doc_length, f_f):
	for word in f_f.keys():
		for file in doc_length.keys():
			if file in inverted[word].keys():
				tf[file][word] = float(f_f[word][file]) / float(doc_length[file])

	return tf

### Step2 : IDF  { word : idf, word2 : idf}
def inverse_document_frequency(files_num, word, d_f):
	if word in d_f.keys():
		return 1.0 + log( float(files_num) / d_f[word] )
	else:
		return 1.0

### tf_idf = {word: { file1: tf_idf, file 2: tf_idf}}
def tf_times_idf(tf, idf, query, documents):
	for word in query:
		for file in documents.keys():
			if (word in tf[file].keys()):
				#print(word + " has tf: ", tf[file][word])
				#print(word + " has idf: ", idf[word])
				tf_idf[word][file] = (tf[file][word]) * (idf[word])
				#print(word + " has tf * idf = ", tf_idf[word][file] )
			else:
				tf_idf[word][file] = 0.0
	return tf_idf

def dp(words, file, q_tf_idf, tf_idf):
	total = 0
	for word in words:
		total += (q_tf_idf[word] * tf_idf[word][file])
	return total

def cs(words, file, q_tf_idf, tf_idf):
	query_t = 0
	file_t = 0
	for word in words:
		query_t +=  q_tf_idf[word]**2
		file_t += tf_idf[word][file]**2

	query_t =  math.sqrt(query_t)
	file_t = math.sqrt(file_t)
	return (query_t * file_t)


count = 0
documents = {}
inverted = {}
#frequency of the word in total files
d_f = {}
#frequency of the word in each file: {'apple': {'doc1': 1, 'doc2': 1}, 'banana' : {'doc1' : 2, 'doc3' : 3}}
f_f = defaultdict(dict)
doc_length = {}

### tf = { word: tf, word2: tf2 }
tf = defaultdict(dict)
idf = {}
tf_idf = defaultdict(dict)


with zipfile.ZipFile('rhf.zip') as z:
	for file_path in z.namelist():
		if file_path.endswith(".html"):
			count = count + 1
			with z.open(file_path, 'r') as f:

				text = f.read().decode('ISO-8859-1')
				text = parsing_texts(text)
				#print(text)
				#print('*'* 100)
				documents[file_path] = text
				#print(documents.keys())
				doc_length[file_path] = len(word_index(text))
				#print("length is " , doc_length[file_path])


		if count == 9361:
			break
	### documents {file_path: file_text}
	#pprint.pprint(documents)

### Build Inverted-Index for documents
for doc_id, text in documents.items():
	doc_index = inverted_index(text)
	###doc_index is now {word:[locations]}

	### {word:{doc_id:[locations]}}
	inverted, f_f = inverted_index_add(inverted, doc_id,doc_index, f_f)

for word in f_f.keys():
	d_f[word] = len(f_f[word].keys())

########################################################
total_file_num = 9361 #9361 ## substitute to Count

### Step 1 TF
tf = term_frequency(doc_length, f_f)

###Step 2 IDF
for word in inverted.keys():
		#print("word is ", word)
	idf[word] = inverse_document_frequency(total_file_num, word, d_f)
##############################################################

end_time = time.time()
# print("Length of d_f is ", len(d_f.keys()))
# print("Length of f_f is ", len(f_f.keys()))
# print("Length of inverted is ", len(inverted.keys()))
print("Elapsed time was %g seconds" % (end_time - start_time))
##Elapsed time was 209.074 seconds

# print(f_f)
# print("*"* 100)
# print(d_f)
# print(inverted)
# from itertools import islice
# def take(n, iterable):
# 	return list(islice(iterable,n))

# # n_items = take(10, inverted.items())
# # print(n_items)


type_search = 'B'

ans = 'y'

while type_search != 'Q':
	type_search = input("Which search would you like to do? B for boolean. G for general. Q for quit: ")
	if type_search == 'B':
		while ans == 'y':
			query = input("Search: ")
			print("*" * 100)
			query = query.lower()

			#print("You typed: %s" % query)
			### Either one appears in the file
			if '"' in query:
				query = query.replace('"', '')
				#if check_word_validity(inverted, query,"and"):
				if check_word_validity(inverted, query,"or"):
					search_for_phrase(inverted, query, documents)

				else:
					print("Sorry we couldn't find that phrase!")

			elif " or " in query:
				query = query.replace(" or ", ' ')
				if check_word_validity(inverted, query, 'or'):
					###result_docs contain list of sets [(doc_id, doc_id),(doc_id,doc_id,doc_id)(doc_id)...]
					#result_docs = search_for_or(inverted, query)
					#Underscore will ignore the first value since word_index returns set of [(index, word),(index,word)]

					for word in query.split():
						if word in inverted:
					#for _, word in word_index(query):
							print("The word '%s' appears in total of %d files and files are: " % (word, len(inverted[word].keys())))

							for k, v in inverted[word].items():

								print("'%s' with frequency of %d times and it's locations are: %r" %(k, len(v), v))
							print('\n')
						else:
							print("The word '%s' does not exist"% word)
				else:
					print('Sorry you put non-existing word')

			### BOTH words must be appears in the file
			elif " and " in query:
				searchKey = query
				query = query.replace(" and ", ' ')
				###If NOT is included
				if "(not " in query:
					###FIX ME  ## word (not word2) = word word2
					query = query.replace("not",'')

					#print("(not) is replced ", query.replace("(",'').replace(')',''))
					if check_word_validity(inverted, query.replace("(",'').replace(')',''), 'and'):
						#print(query)
						not_word = re.search(r'\((\s?\w+)\)', query)
						#print(not_word)

						not_word = not_word.group().lower()
						not_word = re.sub(r'[\(\)]',"", not_word) ###left with word
						not_word = not_word.strip()
						#print("not_word" ,not_word)
						###intersection
						#result_docs = search_for_and(inverted, query)

						right_word = re.sub(r"\s?\(.*?\)", "", query)
						right_word = right_word.strip()
						#print("right word", right_word)
						print("Result search for '%s' is following:" %searchKey)

						right_word_docs = set(inverted[right_word].keys())
						not_word_docs = set(inverted[not_word].keys())
						result = right_word_docs.difference(not_word_docs)

						print("The word %s appears %d files:" %(right_word, len(result)))
						for doc in result:
							print("%s with frequency of %d times and it's locations are: %r" %(doc, len(inverted[right_word][doc]), inverted[right_word][doc]))

					else:
						print("Sorry couldn't find it! ")

				else:
					if check_word_validity(inverted, query, 'and'):
						# query = query.split()
						result_docs = search_for_and(inverted, query)
						if len(result_docs) != 0:
							print("The word '%s' appears in total of %d files and files are: %r" % (searchKey, len(result_docs), result_docs))
							for doc in result_docs:
								print("File name : %s " % doc)
								for _, word in word_index(query):
									 print("%s's location(s) in file: %r" %(word, inverted[word][doc]))
								print('\n')
						else:
							print("Sorry couldn't find it! ")

					else:
						print('Sorry you put non-existing word')

			###If there is no BOOLEAN OPERATION and not a phrase but words
			else:
				if check_word_validity(inverted, query,'or'):
					#result_docs = search_for_or(inverted,query)
					#words = word_index(query)
					for word in query.split():
						if word in inverted:
					#for _, word in word_index(query):
							print("The word %s appears in total of %d and files are: " % (word, len(inverted[word].keys())))
							for k, v in inverted[word].items():
								print("%s with frequency of %d times and it's locations are: %r" %(k, len(v), v))
							print('\n')
						else:
							print("The word '%s' does not exist\n"% word)
				else:
					print("Oups! We couldn't find that word(s)!")

			ans = input("\n Search again? y/n: ")

	elif type_search == 'G':

		search_valid = True

		while search_valid == True:
			search = input("type some query: ")
			search = search.split()
			for word in search:
				if word not in inverted.keys():
					# print(word)
					print("Sorry that word doesn't exist. type some valid query. ")
					#search = search.remove(word)
					search = ""
					break
				else:
					search_valid = False


		# ### Step 1 TF
		# tf = term_frequency(doc_length, f_f)

		# ###Step 2 IDF
		# for word in inverted.keys():
		# 		#print("word is ", word)
		# 	idf[word] = inverse_document_frequency(total_file_num, word, d_f)

		### Step 3 TF * IDF
		tf_idf = tf_times_idf(tf, idf, search, documents)

		### Step 4 Calculate the TF*IDF for the query
		d = defaultdict(int)
		query_tf = {}
		query_idf = {}
		query_tf_idf = {}

		for word in search:
			d[word] += 1

		for word in search:
			query_tf[word] = float(d[word]) / len(search)
			#query_idf[word] = 1.0 + (log(float(total_file_num) / len(inverted[word].keys())))
			query_idf[word] = idf[word]
			query_tf_idf[word] = query_tf[word] * query_idf[word]

		### Step 5 Vector Space Model – Cosine Similarity
		cosine = {}

		for doc in documents.keys():

			dot_product = dp(search, doc, query_tf_idf, tf_idf)
			cosine_sim = cs(search, doc, query_tf_idf, tf_idf)

			if cosine_sim == 0:
				continue

			cosine[doc] = dot_product / cosine_sim

		cosine = OrderedDict(sorted(cosine.items(), key=lambda kv:kv[1],reverse=True))
		#print(cosine)

		# print(cosine)
		###ONly prints top 10
		cosine = itertools.islice(cosine.items(),0,10)

		for file, similarity in cosine:
			print(file, similarity)

	else:
		break
	type_search = input("Q to quit. C to continue: ")

		# #print("dot product ", dot_product)
		# #print(cosine_sim)

		# #print("query tf   query idf    query_tf_idf")
		# #print(query_tf, query_idf, query_tf_idf)

		# # print(query_tf['browser'])
		# # print("browser query_idf is ", query_idf['browser'])
		# # print("browser query_tf_idf is ", query_tf_idf['browser'])

		# # print("TF is: ")
		# # pprint.pprint(tf)

		# # print("IDF is :")
		# # pprint.pprint(idf)

		# # print("Tf * idf is:")
		# # pprint.pprint(tf_idf)
























