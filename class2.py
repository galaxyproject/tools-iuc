from pprint import pprint
from jinja2 import Template, Environment, PackageLoader
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
		#TODO Associate each file in the dict with a user-given name in the XML allowing relevant data to be aggregated together.
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

class CircosObj():
	def __init__(self):
		pass
	
	def set_param(self,param,val):
		vars(self)[param] = val

class TopLevelObj(CircosObj):
	def __init__(self):
		self.llo = []
		pass

	def add_rule(self):
		pass

class Highlight(TopLevelObj):
	def __init__(self):
		self.file_ = None
		self.fill_color = None
		self.ideogram = None
		self.r0 = None
		self.r1 = None
		self.stroke_color = None
		self.stroke_thickness = None
		self.z = None
	
	def make_data_file(self):
		pass

class Ideogram(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.band_stroke_color = None
		self.band_stroke_thickness = None
		self.band_transparency = None
		self.fill = None
		self.fill_bands = None
		self.fill_color = None
		self.label_case = None
		self.label_font = None
		self.label_format = None
		self.label_parallel = None
		self.label_radius = None
		self.label_size = None
		self.label_with_tag = None
		self.radius = None
		self.show = None
		self.show_bands = None
		self.show_label = None
		self.stroke_color = None
		self.stroke_thickness = None
		self.thickness = None
		self.add_spacing()

	def add_spacing(self):
		self.llo += [Spacing()]

class Plot(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)

class Tick(TopLevelObj):
	def __init__(self):
		self.color = None
		self.format_ = None
		self.grid = None
		self.grid_color = None
		self.grid_end = None
		self.grid_start = None
		self.grid_thickness = None
		self.label_offset = None
		self.label_separation = None
		self.label_size = None
		self.multiplier = None
		self.radius = None
		self.show_label = None
		self.size = None
		self.skip_first_label = None
		self.skip_last_label = None
		self.spacing = None
		self.thickness = None
		self.tick_separation = None

class Image(TopLevelObj):
	def __init__(self):
		self.24bit* = None
		self.auto_alpha_colors* = None
		self.auto_alpha_steps* = None
		self.angle_offset* = None
		self.angle_orientation* = None
		self.background* = None
		self.dir* = None
		self.file* = None
		self.radius* = None

class LowLevelObj(CircosObj):
	def __init__(self):
		pass

class Pairwise(LowLevelObj):
	def __init__(self):
		self.chromosome_pair = ['','']
		self.spacing = None

class Spacing(LowLevelObj):
	def __init__(self):
		self.axis_break = None
		self.axis_break_at_edge = None
		self.axis_break_style = None
		self.break_ = None
		self.default = None
		self.llo = [Break_Style() for num in range(int(axis_break_style))]
		

class Break_Style(LowLevelObj):
	def __init__(self):
		self.fill_color = None
		self.stroke_color = None
		self.stroke_thickness = None
		self.thickness = None

class Rule(LowLevelObj):
	def __init__(self):
		self.condition = None
		self.show = None
		self.show_bands = None
		self.show_ticks = None

class CircosPlot():
	def __init__(self,basefilename):
		self.basename = basefilename
		self.data_dict = {}
		self.last_filled_radius = 1.1
		self.master_struct = {'ideograms':{}} #Will be loop-able for jinja2 
		self.plots = ['<plots>']
		self.track_dict = {} #Temporary -- will be obviated by XML at later stage.
		self.add_ideogram()
		self.add_top_lvl_params()

	def add_top_lvl_params():
		self.anglestep = None
		self.beziersamples = None
		self.chromosomes_display_default = None
		self.chromosomes = None
		self.chromosomes_breaks = None
		self.chromosomes_order = None
		self.chromosomes_radius = None
		self.chromosomes_reverse = None
		self.chromosomes_units = None
		self.debug = None
		self.imagemap = None
		self.minslicestep = None
		self.show_ticks = None
		self.show_tick_labels = None
		self.units_nounit = None
		self.units_ok = None
		self.warnings = None

	def append_data(self, interpreter):
		#Add interpreter data to ConfWriter dict. Depends on what kind of data structure we decide on... 
		self.raw_dict = interpreter.files_dict
		self.seq_dict = interpreter.seq_dict
		self._make_karyotype()

	def add_ideogram(self,spacing_opt='0.005r',radius_opt='0.5r',thickness_opt='20p',fill_opt='yes',color_opt=[]):
		l = len(self.master_struct['ideograms']) + 1
		ideogramid = 'ideogram-'+str(l)
		self.master_struct['ideograms'][ideogramid] = {}
		self.add_spacing(spacing_opt,ideogramid)
		block = {'radius':radius_opt,
			 'thickness':thickness_opt,
			 'fill':fill_opt}
		self.master_struct['ideograms'][ideogramid].update(block)

	def add_spacing(self,spacing_opt,ideogramid):
		self.master_struct['ideograms'][ideogramid].update({'spacing':{'default':spacing_opt}})

	def create_base_conf_template(self):
		E = Environment(loader=PackageLoader('dummy','Templates'))
		T = E.get_template('template.conf')
		with open(self.basename + '.conf','w') as f:
			rendering = T.render(master_struct = self.master_struct, filename = self.basename + '.png', karyotype_filename = self.karyotype_filename)
			f.write(rendering)

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
		if 'plots' not in self.master_struct:
			self.master_struct['plots'] = {}
		l = len(self.master_struct['plots']) + 1
		plotname = 'plot-'+type_+'-'+str(l)
		self.master_struct['plots'][plotname] = {}
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
		self.track_dict[key] = '' 
		self._four_element_processing(f,key)
		self._add_plot(type_,f,plotname)

	def _four_element_processing(self,file_,key):
		#Heatmaps, scatterplots, and histograms share the same 4-element format...
		with open(file_,'w') as handle:
			for f in self.raw_dict:
				if os.path.splitext(f)[1] in ('.wig','.bw','.bigWig'):
					for chromosome in self.raw_dict[f]:
						for feature in self.raw_dict[f][chromosome]:
							datastring = chromosome+' '+str(self.raw_dict[f][chromosome][feature]['start'])+' '+str(self.raw_dict[f][chromosome][feature]['end'])+' '+str(self.raw_dict[f][chromosome][feature]['val'])+'\n'
							handle.write(datastring)

	def update_master_conf(self):
		self.plots.append('</plots>')
		with open(self.mainconf,'a') as conf:
			for element in self.plots:
				conf.write(element+'\n')
	
	def _add_plot(self,block_type,filename,plotname,extend_bin='no',
		      max_ = '1.0', min_ = '0.0', glyph='rectangle',
		      glyph_size='8', color_list='spectral-9-div', 
		      stroke_color='dred', thickness='1',color='red',log_base_opt='1'):
		r1 = str(self.last_filled_radius + 0.100) + 'r'
		r2 = str(self.last_filled_radius + 0.200) + 'r'
		self.last_filled_radius += 0.20
		if block_type == 'hist':
			block = {'type':'histogram',
				 'file':filename,
				 'r0':r1,
				 'r1':r2,
				 'extend_bin':extend_bin}
				 
		elif block_type == 'heat':
			block = {'show':'yes',
				 'type':'heatmap',
				 'file':filename,
				 'color':color_list,
				 'stroke_thickness':thickness,
				 'r0':r1,
				 'r1':r2,
				 'scale_log_base':log_base_opt
				 }
		elif block_type == 'scatter':
			block = {'show':'yes',
				 'type':'scatter',
				 'file':filename,
				 'r0':r1,
				 'r1':r2,
				 'max':max_,
				 'min':min_,
				 'glyph':glyph,
				 'glyph_size':glyph_size,
				 'color':color,
				 'stroke_color':stroke_color,
				 'stroke_thickness':thickness,
				 }
		else:
			raise ValueError('Unsupported plot type.')
		self.master_struct['plots'][plotname] = block
					
	def add_links(self):
		pass

class XML_parser():
	#Reads the XML file and manages user-specified options
	def __init__(self):
		pass	

if __name__ == "__main__":
	I = Ideogram()
	pprint(vars(I))
