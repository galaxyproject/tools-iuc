from pprint import pprint
from Bio.Seq import Seq
from BCBio import GFF
from subprocess import call
import os

def dprint(obj):
	pprint(obj)
	x = raw_input("<RETURN>")

class DataInterpreter():
	def __init__(self,files_list=[]):
		self.files_dict = {}
		if files_list == []:
			files_list = self._get_files_from_xml()
		for f in files_list:
			name,ext = os.path.splitext(f)
			if ext == ".bed":
				self.parse_bed(f)
			elif ext in (".bw",".bigWig"):
				self.parse_bigWig(f)
			elif ext == ".gff":
				self.parse_gff3(f)
			else:
				raise ValueError('Unsupported File Format')

	def _get_files_from_xml(self):
		pass

	def parse_gff3(self,f):
		with open(f) as input_handle:
			for rec in GFF.parse(input_handle):
				features_dict = {}
				tmpdict = {'id':rec.id,'seq':rec.seq,'description':rec.description,'dbxrefs':rec.dbxrefs,'features':rec.features,'annotations':rec.annotations,'name':rec.name,'letter_annotations':rec.letter_annotations}
				features_dict[rec.id] = tmpdict
		self.files_dict[f] = features_dict

	def parse_bigWig(self,f):
		features_dict = {}
		locidict = {}
		path = os.path.dirname(f)
		name,ext = os.path.splitext(f)
		newfile = name + '.wig'
		call('./bigWigToWig '+ f +' '+newfile,shell=True) #obviously needs changing
		with open(newfile) as handle:
			i = 0
			currentchrom = ''
			currentspan = 0
			start = 1
			step = 1
			m = 1
			mode = ''
			locidict = {}
			for line in handle.readlines():
				l = line.split()
				if l[0] == 'variableStep':
					mode = 'variable'
					currentchrom = l[1].split('=')[1]
					try:
						currentspan = int(l[2].split('=')[1])
					except IndexError:
						currentspan = 0
				elif l[0] == 'fixedStep':
					mode = 'fixed'
					m = 1
					currentchrom = l[1].split('=')[1]
					start = int(l[2].split()[1])
					step = int(l[3].split()[1])
					try:
						currentspan = int(l[4].split('=')[1]
					except IndexError:
						currentspan = 0
					pass
				elif mode == 'variable':
					locidict[currentchrom+'_'+l[0]+'-'+str(int(l[0])+currentspan)] = {'chromosome':currentchrom,'start':int(l[0]),'end':int(l[0])+currentspan,'val':l[1]}
				elif mode == 'fixed':
					locidict[currentchrom+'_'+str(start+m*step)+'-'+str(start+m*step+currentspan)] = {'chromosome':currentchrom,'start':start+m*step,'end':start+m*step+currentspan,'val':l[0]}
					m+=1
				i+=1
		self.files_dict[f] = locidict

	def parse_bed(self, f):
		bed_standard_fields = ['chromosome','start','end',
				       'name','score','strand',
				       'thickstart','thickend','rgb',
				       'blockcount','blocksizes','blockstarts']
		features_dict = {}
		with open(f) as tmp:
			for l in tmp.readlines():
				data = l.strip().split()	
				tmpdict = dict(zip(bed_standard_fields,data))
				identifier = tmpdict['chromosome']+'_'+tmpdict['start']+'-'+tmpdict['end']
				features_dict[identifier] = tmpdict
		self.files_dict[f] = features_dict

	def integrate_mauve_data(self,data):
		#This method may be required to integrate whatever that parser outputs, or it may be deprecated depending on what that data looks like.
		pass

class CircosPlot():
	def __init__(self):
		self.data_dict = {}

	def append_data(self, interpreter):
		#Add interpreter data to ConfWriter dict. Depends on what kind of data structure we decide on... 
		pass

	def write_base_conf(self):
		pass

	def add_histogram(self):
		pass

	def add_links(self):
		pass

	def add_heatmap(self):
		pass	

	def call_circos(self):
		subprocess.call('circos -conf ' + self.config ,shell=True)


if __name__ == "__main__":
	D = DataInterpreter(['./test-data/ColorByStrandDemo.bed','./test-data/GCF_000146045.2_R64_genomic.gff','./test-data/GSM646567_hg19_wgEncodeUwDgfK562Raw.bigWig'])
	pprint(D.files_dict)
