from pprint import pprint
from Bio.Seq import Seq
from BCBio import GFF

def dprint(obj):
	pprint(obj)
	x = raw_input("<RETURN>")

class DataInterpreter():
	def __init__(self):
		self.files_dict = {}

	def add_file(self, f):
		self.files_dict[f] = ''
		
	def make_files_dict(self):
		pass

	def parse_gff3(self,files_list):
		for f in files_list:
			input_handle = open(f)
			for rec in GFF.parse(input_handle):
				features_dict = {}
				tmpdict = {'id':rec.id,'seq':rec.seq,'description':rec.description,'dbxrefs':rec.dbxrefs,'features':rec.features,'annotations':rec.annotations,'name':rec.name,'letter_annotations':rec.letter_annotations}
				features_dict[rec.id] = tmpdict
			input_handle.close()
			self.files_dict[f] = features_dict

	def parse_bigwig(self):
		pass

	def parse_bed(self, files_list):
		bed_standard_fields = ['chromosome','start','end',
				       'name','score','strand',
				       'thickstart','thickend','rgb',
				       'blockcount','blocksizes','blockstarts']
		for f in files_list:
			features_dict = {}
			tmp = open(f,'r')
			for l in tmp.readlines():
				data = l.strip().split()	
				tmpdict = dict(zip(bed_standard_fields,data))
				identifier = tmpdict['chromosome']+'_'+tmpdict['start']+'-'+tmpdict['end']
				features_dict[identifier] = tmpdict
			tmp.close()
			self.files_dict[f] = features_dict
		pass

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
	D = DataInterpreter()
	D.parse_bed(['./test-data/ColorByStrandDemo.bed'])
	D.parse_gff3(['./test-data/GCF_000146045.2_R64_genomic.gff'])
	dprint(D.files_dict)
	
