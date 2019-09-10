# imports usually go at the top of the file
import os  # to run command-line programs from Python
# import pysam
import argparse  # for command-line arguments

parser = argparse.ArgumentParser()  # create ArgumentParser object by calling the ArgumentParser constructor
									# Remember, to use something from a module, prefix that thing with the module name & a dot
                                    # An object is a variable whose data type is a class, which is a type of user-defined data type. (ParseFastQ below is a class.)
									# Create an object by calling the constructor of the class. This is __init__().
									# The constructor may or may not take arguments. Apparently this one does not.
# add command-line arguments
# The arguments to add_argument are the name of the argument and the message you get when you run the program w/ the -h option
parser.add_argument("--fastq", help="insert pooled fastq file here")
parser.add_argument("--clinical", help="insert clinical data file here")
parser.add_argument("--ref", help="insert reference file here")
args = parser.parse_args()  # calling a method (parse_args()) on an object (parser)

# I didn't even really look at this code. That's the nice thing about code reuse & modularity. You don't have to know the inner workings of someone else's code. All you need to know is what it does at a high level, & how to use it. See the example usage in the comment.
class ParseFastQ(object):
	"""Returns a read-by-read fastQ parser analogous to file.readline()"""
	def __init__(self,filePath,headerSymbols=['@','+']):
		"""Returns a read-by-read fastQ parser analogous to file.readline().
		Exmpl: parser.next()
		-OR-
		Its an iterator so you can do:
		for rec in parser:
		... do something with rec ...

		rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
		"""
		if filePath.endswith('.gz'):
			self._file = gzip.open(filePath)
		else:
			self._file = open(filePath, 'rb')
			self._currentLineNumber = 0
			self._hdSyms = headerSymbols

	def __iter__(self):
        while True:
            yield self.__next__()

	def __next__(self):
		"""Reads in next element, parses, and does minimal verification.
		Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
		# ++++ Get Next Four Lines ++++
		elemList = []
		for i in range(4):
			line = self._file.readline()
			self._currentLineNumber += 1 ## increment file position
			if line:
				elemList.append(line.strip('\n'))
			else:
				elemList.append(None)

		# ++++ Check Lines For Expected Form ++++
		trues = [bool(x) for x in elemList].count(True)
		nones = elemList.count(None)
		# -- Check for acceptable end of file --
		if nones == 4:
			raise StopIteration
		# -- Make sure we got 4 full lines of data --
		assert trues == 4,\
		"** ERROR: It looks like I encountered a premature EOF or empty line.\n\
		Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
		# -- Make sure we are in the correct "register" --
		assert elemList[0].startswith(self._hdSyms[0]),\
		"** ERROR: The 1st line in fastq element does not start with '%s'.\n\
		Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber)
		assert elemList[2].startswith(self._hdSyms[1]),\
		"** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
		Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber)
		# -- Make sure the seq line and qual line have equal lengths --
		assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
		Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)

		# ++++ Return fatsQ data as tuple ++++
		return tuple(elemList)

# This code uses a lot of tuples. Example tuple definition: tup=(1,2,3). (I wanted to call this variable tuple, but tuple is a reserved word (keyword).
# function definitions go above the "main" code.
def trim_beg(fastq):      # example tuple passed into trim_beg(): (@arbitrary-seq-111, ATCGATGCACCCC, +, IIIDIIIIIDDFF)
                          # access elements of a tuple by index, just like elements of a list: e.g., fastq[0] is first element of fastq
						  # fastq[1] is a string, whose characters we can access by index - or range of indices.
						  # Select characters 5 thru the end, for sequence & confidence scores: fastq[1][5:], fastq[3][5:]
	fastq=(fastq[0], fastq[1][5:] ,fastq[2], fastq[3][5:])  # Tuple are immutable (can't be modified), so instead of changing elements 1 & 3, we create a new tuple that's identical
	                                                        # except for chopping off the first 5 characters of elements 1 & 3
	return fastq  # return = exit function
	              # This is a value-returning function, so when it exits, it takes something w/ it - in this case, the value of fastq

#@seq13534-419
GCAGTAGCGGTCATAAGTGGTACATTACGAGATTCGGAGTACCATAGATTCGCATGAATCCCTGTGGATACGAGAGTGTGAGATATATGTACGCCAATCCAGTGTGATACCCATGAGATTTAGGACCGATGATGGTTGAGGACCAAGGATTGACCCGATGGATGCAGATTTGACCCCAGATAGAATAAATGCGATGAGATGATTTGGCCGATAGATAGATAGTGTCGTGAGGTGACGTCCGTCACTGGACGAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFDFFDFFDDFDFDFFFFDDFFDDFDDFF##


def trim_end(fastq):
	#while we haven't reached the end of the string, check element i and next element. If both are D/F, exit the loop.
	i=0
	while i < len(fastq[3])-1:
		if fastq[3][i] in ["D", "F"] and fastq[3][i+1] in ["D", "F"]:
			break
		i += 1
	#remove the rest of the string
	fastq=(fastq[0], fastq[1][:i], fastq[2], fastq[3][:i])
	return fastq

def pileup(bam):
	samfile = pysam.AlignmentFile(bam, "rb")
	#Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
	for pileupcolumn in samfile.pileup():
		print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
		#use a dictionary to count up the bases at each position
		ntdict = {}
		ntdict['A']=ntdict['G']=ntdict['C']=ntdict['T']=0
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
				# You can uncomment the below line to see what is happening in the pileup.
				print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
				base = pileupread.alignment.query_sequence[pileupread.query_position]
				########## ADD ADDITIONAL CODE HERE #############
				# Populate the ntdict with the counts of each base
				ntdict[base] += 1
				# This dictionary will hold all of the base read counts per nucletoide per position.
				# Use the dictionary to calculate the frequency of each site, and report it if if the frequency is NOT  100% / 0%.
				#############################################
		print(ntdict)
	samfile.close()

# open clinical data file for reading
clinical_data=open(args.clinical, "r").readlines()[1:]

# remove fastqs directory if it exists
if os.path.isdir("fastqs"):
	#os.rmdir("fastqs")  # "directory is not empty" error, so changed to a Linux command I know will remove a non-empty directory
	os.system("rm -r fastqs")
os.mkdir("fastqs") # create fastqs directory
#os.system("cd fastqs")

# create list of names for use in loop later
names=[]
for line in data:
	names.append(line.split("\t")[0])

fastqfile = ParseFastQ(args.fastq)# replced args with file name
name=None
for fastq in fastqfile:
	for line in clinical_data:
		#print(fastq[1][:5])
		#print(line.split("\t")[2])
		#print("----------------------------------------------------")
		if fastq[1][:5].strip() == line.split("\t")[2].strip():   # These were apparently never equal b/c one of them had a newline. strip() takes care of that.
			name = line.split("\t")[0]
			break
	trim=trim_beg(fastq)
	trim= trim_end(trim)

	file_=open("fastqs/"+name+"_trimmed.fastq", "w+")
	file_.write("".join(trim))  # trim is a tuple, and write() arg must be string. join() to the empty string takes care of that.
	file_.close()

for name in names:

	#os.system("bwa index ${args.ref}") # string interpolation did not work, so used Python string concatenation
	os.system("bwa index "+args.ref)
	#os.system("bwa mem ${args.ref} ${name}_trimmed.fastq > ${name}.sam")
	os.system("bwa mem " + args.ref + " " + name + "_trimmed.fastq > " + name + ".sam")
	os.system("samtools view -bS ${name}.sam > ${name}.bam")
    os.system("samtools sort " + name + ".bam > " + name + ".sorted.bam")
	os.system("samtools index ${name}.sorted.bam")
	os.system("rm ${name}.sam ${name}.bam")
	pileup(name+".sorted.bam")

data.close()
