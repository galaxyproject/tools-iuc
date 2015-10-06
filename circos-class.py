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
			elif ext in (".bw",".bigWig",'.wig'):
				self.parse_bigWig(f)
			elif ext in (".gff",".gff3"):
				self.parse_gff3(f)
			elif ext == ".fa":
				self.parse_fastA(f)
			else:
				raise ValueError('Unsupported File Format')

	def _get_files_from_xml(self):
		pass

	def parse_gff3(self,f):
		with open(f) as input_handle:
			for rec in GFF.parse(input_handle):
				flat_list = []
				features_dict = {}
				tmpdict = {'id':rec.id,'seq':rec.seq,'description':rec.description,'dbxrefs':rec.dbxrefs,'annotations':rec.annotations,'name':rec.name,'letter_annotations':rec.letter_annotations}
				for obj in rec.features:
					if obj.sub_features == []:
						flat_list.append(obj)
					else:
						flat_list = self._process_obj_subfeatures(obj,flat_list)
				tmpdict['features'] = flat_list
				features_dict[rec.id] = tmpdict
		self.files_dict[f] = features_dict

	def _process_obj_subfeatures(self,obj,flist):
		flist.append(obj)
		l = obj.sub_features
		for element in l:
			if l != []:
				self._process_obj_subfeatures(element,flist)
			else:
				flist.append(element)
		return flist
		
	
	def parse_fastA(self,f):
		self.seq_dict = {}
		seqname = ''
		seqres = ''
		with open(f) as input_handle:
			for line in input_handle.readlines():
				if line[0] == '>':
					if seqname != '':
						self.seq_dict[seqname] = seqres
						seqres = ''
					seqname = line[1:].strip()
				else:
					seqres += line.strip()
		self.seq_dict[seqname] = seqres
		

	def parse_bigWig(self,f):
		features_dict = {}
		locidict = {}
		path = os.path.dirname(f)
		name,ext = os.path.splitext(f)
		newfile = name + '.wig'
		if ext != '.wig':
			call('./bigWigToWig '+ f +' '+newfile,shell=True) #obviously needs changing
		else:
			call('cp '+ f + ' .',shell=True)
		with open(newfile) as handle:
			i = 0
			currentchrom = ''
			currentspan = 0
			start = 1
			step = 1
			m = 1
			mode = ''
			locidict = {}
			tmpdict = {}
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
						currentspan = int(l[4].split('=')[1])
					except IndexError:
						currentspan = 0
				elif mode == 'variable':
					tmpdict[currentchrom+'_'+l[0]+'-'+str(int(l[0])+currentspan)] = {'chromosome':currentchrom,'start':int(l[0]),'end':int(l[0])+currentspan,'val':l[1]}
				elif mode == 'fixed':
					tmpdict[currentchrom+'_'+str(start+m*step)+'-'+str(start+m*step+currentspan)] = {'chromosome':currentchrom,'start':start+m*step,'end':start+m*step+currentspan,'val':l[0]}
					m+=1
				i+=1
			for key in tmpdict:
				if tmpdict[key]['chromosome'] not in locidict.keys():
					locidict[tmpdict[key]['chromosome']] = tmpdict[key]
				else:
					locidict[tmpdict[key]['chromosome']].update({key:tmpdict[key]})
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
	D = DataInterpreter(['./test-data/miro.fa','./test-data/miro.gff3','./test-data/miro.wig'])
	pprint(D.files_dict)
