from pprint import pprint
from jinja2 import Template
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
					locidict[tmpdict[key]['chromosome']] = {}
					locidict[tmpdict[key]['chromosome']].update({key:tmpdict[key]})
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
	def __init__(self,basefilename):
		self.data_dict = {}
		self.basename = basefilename
		self.track_dict = {} #Temporary -- will be obviated by XML at later stage.
		self.master_struct = {'ideograms':{}}
					} #Will be loop-able for jinja2 
		self.last_filled_radius = 1.1
		self.plots = ['<plots>']

	def append_data(self, interpreter):
		#Add interpreter data to ConfWriter dict. Depends on what kind of data structure we decide on... 
		self.raw_dict = interpreter.files_dict
		self.seq_dict = interpreter.seq_dict

	def create_base_conf(self,spacing_opt='0.005r',radius_opt='0.9r',thickness_opt='20p',fill_opt='yes',color_opt=[]):
		self._make_karyotype()
		boilerplate = ['karyotype = '+self.karyotype_filename,
			       '<ideogram>','<spacing>','default = '+spacing_opt,
			       '</spacing>','radius = '+radius_opt,'thickness = '+thickness_opt,'fill = '+fill_opt,'</ideogram>','<image>','<<include etc/image.conf>>','file*= '+self.basename,'</image>','<<include etc/colors_fonts_patterns.conf>>','<<include etc/housekeeping.conf>>']
		self.mainconf = self.basename + '.conf'
		with open(self.mainconf,'w') as conf:
			for line in boilerplate:
				conf.write(line+'\n')

	def create_base_conf_template(self):
		T = Template("""karyotype = {{karyotype_filename}}
				<ideogram>
				<spacing>
				</spacing>
				{% if zooms in self.master_struct %}
					{% for zoom in self.master_struct[zooms] %}
					{% endfor %}
				{% endif %}
				{% if highlights in self.master_struct %}
					{% for highlight in self.master_struct[highlights] %}
						{% if rules in self.master_struct[highlights][highlight]%}
							{% for rule in self.master_struct[highlights][highlight][rules] %}
							{% endfor %}
						{% endif %}
					{% endfor %}
				{% endif %}
				{% if plots in self.master_struct %}
					{% for plot in self.master_struct[plots] %}
					{% endfor %}
				{% endif %}
				{% if links in self.master_struct %}
					{% for link in self.master_struct[links] %}
					{% endfor %}
				{% endif %}
				{% if ticks in self.master_struct %}
					{% for tick in self.master_struct[ticks] %}
					{% endfor %}
				{% endif %}
				{% for ideogram in self.master_struct[ideograms] %}
				</ideogram>
				<image>
				<<include etc/image.conf>>
				file* = {{filename}}
				</image>
				<<include etc/colors_fonts_patterns.conf>>
				<<include etc/houskeeping.conf>>
				""")
		T.render(self.master_struct)

	def _make_karyotype(self,colors_list=[]):
		if colors_list == []:
			colors_list = ['chr1','chr2','chr3','chr4','chr5','chr6',
				       'chr7','chr8','chr9','chr10','chr11','chr12',
				       'chr13','chr14','chr15','chr16','chr17','chr18'
				       'chr19','chr20','chr21','chr22','chrx','chry'
					] #Temporary -- will update to include more diverse list of colors, Brewer...
		self.karyotype_filename = self.basename + '-karyotype.txt'
		i = 0
		with open(self.karyotype_filename,'w') as karyotype:
			for element in self.seq_dict:
				karyotype.write('chr - '+element+' '+element+' 0 '+str(len(self.seq_dict[element]))+' '+colors_list[i]+'\n')
				if i < len(colors_list)-1:
					i+=1
				else:
					raise ValueError('Colors list exhausted')

	def add_four_elem_track(self,type_='four-element',filename=''):
		#Generic for heatmap, scatterplot, and histogram
		if type_ not in self.track_dict:
			key = type_
		else:
			i = 1
			while type_ +'-'+str(i) in self.track_dict:
				i+=1
			key = type_+'-'+str(i)
		if filename == '':
			f = self.basename+'-'+key+'.txt'
		else:
			f = filename
		self.track_dict[key] = []
		self._four_element_processing(f,key)
		self._add_plot(type_,f)

	def _four_element_processing(self,file_,key):
		#Heatmaps, scatterplots, and histograms share the same 4-element format...
		with open(file_,'w') as handle:
			for f in self.raw_dict:
				if os.path.splitext(f)[1] in ('.wig','.bw','.bigWig'):
					for chromosome in self.raw_dict[f]:
						for feature in self.raw_dict[f][chromosome]:
							datastring = chromosome+' '+str(self.raw_dict[f][chromosome][feature]['start'])+' '+str(self.raw_dict[f][chromosome][feature]['end'])+' '+str(self.raw_dict[f][chromosome][feature]['val'])+'\n'
							handle.write(datastring)
							self.track_dict[key].append(datastring.split())

	def update_master_conf(self):
		self.plots.append('</plots>')
		with open(self.mainconf,'a') as conf:
			for element in self.plots:
				conf.write(element+'\n')
	
	def _add_plot(self,block_type,filename,extend_bin='no',
		      max_ = '1.0', min_ = '0.0', glyph='rectangle', glyph_size='8', color_list='spectral-5-div', 
		      stroke_color='dred', thickness='1',color='red'):
		r1 = str(self.last_filled_radius + 0.01)
		r2 = str(self.last_filled_radius + 0.02)
		self.last_filled_radius += 0.02
		if block_type == 'hist':
			block = ['<plot>',
				 'type = histogram',
				 'file = '+filename,
				 'r1 = '+r1,
				 'r2 = '+r2,
				 'extend_bin = '+extend_bin,
				 '</plot>']
		elif block_type == 'heat':
			block = ['<plot>',
				 'show = yes',
				 'type = heatmap',
				 'file = '+filename,
				 'color = '+color_list,
				 'r0 = '+r1,
				 'r1 = '+r2,
				 '</plot>']
		elif block_type == 'scatter':
			block = ['<plot>',
				 'show = yes',
				 'type = scatter',
				 'file = '+filename,
				 'r0 = '+r1,
				 'r1 = '+r2,
				 'max = '+max_,
				 'min = '+min_,
				 'glyph = '+glyph,
				 'glyph_size = '+glyph_size,
				 'color = '+color,
				 'stroke_color = '+stroke_color,
				 'stroke_thickness = '+thickness,
				 '</plot>']
		else:
			raise ValueError('Unsupported plot type.')
		for line in block:
			self.plots.append(line)
					
	def add_links(self):
		pass

	def call_circos(self):
		call('circos -conf ' + self.mainconf ,shell=True)

class XML_parser():
	#Reads the XML file and manages user-specified options
	def __init__(self):
		pass	

if __name__ == "__main__":
	D = DataInterpreter(['./test-data/miro.fa','./test-data/miro.gff3','./test-data/miro.wig'])
	C = CircosPlot('miro')
	C.append_data(D)
	C.create_base_conf()
	C.add_four_elem_track('scatter')
	C.update_master_conf()
