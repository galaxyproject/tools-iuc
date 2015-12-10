from pprint import pprint
from jinja2 import Template, Environment, PackageLoader
from BCBio import GFF
from subprocess import call
from Bio.Seq import Seq
import sys
import copy
import os
import re
import xml
import brewer2mpl
import argparse

COLORS_LIST = ['red','green','blue','purple','orange',
	       'yellow','black','grey','white','pred','pgreen','pblue',
               'ppurple','pyellow','porange']

def dprint(obj):
	pprint(obj)
	x = raw_input("<RETURN>")

def test_routine():
	A = argparse.ArgumentParser()
	A.add_argument('xml',nargs=1)
	A.add_argument('-k',nargs='*')
	A.add_argument('-H',nargs='*')
	A.add_argument('-I',nargs='*')
	A.add_argument('-s',nargs='*')
	A.add_argument('-b',nargs='*')
	args = A.parse_args()
	X = args.xml[0] 
	T = xml.etree.ElementTree.parse(X)
	R = T.getroot()
	D = DataInterpreter(R,[])
	D.files_dict = {'karyotype':args.k,
			'scatter':args.s,
			'heatmap':args.H,
			'histogram':args.b,
			'highlight':args.I}
	D._parse_files()
	P = CircosPlot('plot',D.object_list,D.karyotype,D.files_dict)
	#P.files_list = D.txtlist
	P.write_conf()

class DataInterpreter():
	def __init__(self,xmlroot,files_dict={}):
		self.xmlroot = xmlroot
		self._create_object_list()
		self.files_dict = {} 
		if self.files_dict == {}:
		#FIXME
			files_dict = self._get_files_from_xml()
		self._parse_files()

	def _create_object_list(self):
		self.object_list = []
		current_obj = None 
		current_sub = None
		subobjind= 0
		recognized_objects = {'ideogram':Ideogram(),'image':Image(),'link':Link(),
                                      'tick':Tick(),'zoom':Zoom(),'highlight':Highlight(),
				      'plot':Plot()}
		recognized_subobjects = {'pairwise':Pairwise(),'rule':Rule(),'break_style':Break_Style(),'spacing':Spacing()}
		for element in self.xmlroot.iter():
			if element.tag in recognized_objects:
				if current_sub is not None:
					current_obj.llo.append(copy.copy(current_sub))
				if current_obj is not None:
					self.object_list.append(copy.copy(current_obj))
				current_obj = recognized_objects[element.tag]
				subobjind = 0
			elif current_obj is not None:
				if element.tag in recognized_subobjects:
					subobjind = 1
					if current_sub is not None:
						current_obj.llo.append(copy.copy(current_sub))
					current_sub = recognized_subobjects[element.tag]
				if subobjind == 0:
					setattr(current_obj,element.tag,element.text)
				elif subobjind == 1:
					setattr(current_sub,element.tag,element.text)
		self.object_list.append(copy.copy(current_obj))
	
	def _get_files_from_xml(self):
		#TODO When the time comes to create the final Galaxy tool...
		pass

	def integrate_mauve_data(self,data):
		#This method may be required to integrate whatever that parser outputs, or it may be deprecated depending on what that data looks like.
		pass

	def _make_karyotype(self,f,filename="karyotype"):
		self.karyotype = filename+'.txt'
		i = 0
		seqlen = 0
		space = ' '
		with open(f) as input_handle:
			g = open(self.karyotype,'w')
			for line in input_handle:
				if line[0] == '>':
					if i != 0:
						g.write(space.join(['chr','-',cname,str(i),str(0),str(seqlen),COLORS_LIST[i]+'\n']))
						seqlen = 0
					i+=1
					cname = line[1:].strip()
				else:
					seqlen += len(line.strip())
			g.write(space.join(['chr','-',cname,str(i),str(0),str(seqlen),COLORS_LIST[i]+'\n']))
			g.close()
		return self.karyotype

	def parse_bed(self, key,f,index):
		bed_standard_fields = ['chromosome','start','end',
				       'name','score','strand',
				       'thickstart','thickend','rgb',
				       'blockcount','blocksizes','blockstarts']
		features_dict = {}
		with open(self.files_dict[key]) as tmp:
			for l in tmp:
				data = l.strip().split()	
				tmpdict = dict(zip(bed_standard_fields,data))
				identifier = tmpdict['chromosome']+'_'+tmpdict['start']+'-'+tmpdict['end']
				features_dict[identifier] = tmpdict
		self.files_dict[f] = features_dict

	def parse_bigWig(self,sourcekey,f,index):
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
			for line in handle:
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
					tmpdict[currentchrom+'_'+str(start+m*step)+'-'+str(start+m*step+currentspan)] = {'chromosome':currentchrom,'start':start+m*step,
															 'end':start+m*step+currentspan,'val':l[0]}
					m+=1
				i+=1
			for key in tmpdict:
				if tmpdict[key]['chromosome'] not in locidict.keys():
					locidict[tmpdict[key]['chromosome']] = {}
					locidict[tmpdict[key]['chromosome']].update({key:tmpdict[key]})
				else:
					locidict[tmpdict[key]['chromosome']].update({key:tmpdict[key]})
		if sourcekey == 'scatter':
			filename = 'scatter'+'-'+str(index)+'.txt'
		elif sourcekey == 'heatmap':
			filename = 'heatmap'+'-'+str(index)+'.txt'
		elif sourcekey == 'histogram':
			filename = 'histogram'+'-'+str(index)+'.txt'
		g = open(filename,'w')
		for genome in locidict:
			for feature in locidict[genome]:
				g.write(locidict[genome][feature]['chromosome']+ ' ' +str(locidict[genome][feature]['start']) +\
					' ' +str(locidict[genome][feature]['end']) + ' ' +locidict[genome][feature]['val']+'\n')
		g.close()
		return filename

	def _parse_files(self):
		tmp = {'karyotype':[],'scatter':[],'heatmap':[],'histogram':[],'highlight':[]}
		i = 1
		for key in self.files_dict:
			if self.files_dict[key] is not None:
				for f in self.files_dict[key]:
				#TODO Each parse method needs to update self.files_dict to a text file compatible with Circos.
					ext = os.path.splitext(f)[1]
					if key == "karyotype":
						tmp[key] = self._make_karyotype(f)
					elif key in ("scatter","heatmap","histogram","highlight"):
						if ext in ('.gff','.gff3'):
							tmp[key].append(self.parse_gff3(key,f,i))
						elif ext in ('.wig','.bigWig','.bw'):
							tmp[key].append(self.parse_bigWig(key,f,i))
						elif ext == ".bed":
							tmp[key].append(self.parse_bed(key,f,i))
						else:
							raise Exception('Unsupported File Format')
						i+=1
					else:
						raise Exception('Unsupported Data Type')
		self.files_dict = tmp

	def parse_gff3(self,key,f,index):
		#FIXME Seems there's no nontrivial way to represent gff3 file info in Circos... 
		#For now...
		filename = str(key)+'-'+str(index)+'.txt'
		with open(f,'r') as input_handle:
			for rec in GFF.parse(input_handle):
				flat_list = []
				features_dict = {}
				tmpdict = {'id':rec.id,'seq':rec.seq,'description':rec.description,
				           'dbxrefs':rec.dbxrefs,'annotations':rec.annotations,
                                           'name':rec.name,'letter_annotations':rec.letter_annotations}
				for obj in rec.features:
					if obj.sub_features == []:
						flat_list.append(obj)
					else:
						flat_list = self._process_obj_subfeatures(obj,flat_list)
				tmpdict['features'] = flat_list
				features_dict[rec.id] = tmpdict
			g = open(filename,'w')
			for element in tmpdict['features']:
				g.write(' '.join([rec.id,str(int(element.location.start)),str(int(element.location.end))+'\n']))
			g.close()
		return filename

		
	def _process_obj_subfeatures(self,obj,flist):
		flist.append(obj)
		l = obj.sub_features
		for element in l:
			if l != []:
				self._process_obj_subfeatures(element,flist)
			else:
				flist.append(element)
		return flist
	
class CircosObj():
	def __init__(self):
		pass
	

class TopLevelObj(CircosObj):
	def __init__(self):
		self.llo = []
		pass

class Highlight(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'highlight'
		self.type = 'highlight'

class Karyotype(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'karyotype'
	
class Ideogram(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'ideogram'

class Plot(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'plot'

class Zoom(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'zoom'

class Link(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'link'

class Tick(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'tick'

class Image(TopLevelObj):
	def __init__(self):
		TopLevelObj.__init__(self)
		self.__name__ = 'image'
		self.boilerplate = '<<include etc/image.conf>>'

class LowLevelObj(CircosObj):
	def __init__(self):
		pass

class Pairwise(LowLevelObj):
	def __init__(self):
		self.chromosome_pair = ['','']
		self.__name__ = 'pairwise'

class Spacing(LowLevelObj):
	def __init__(self):
		self.__name__ = 'spacing'
		

class Break_Style(LowLevelObj):
	def __init__(self):
		self.__name__ = 'break_style'

class Rule(LowLevelObj):
	def __init__(self):
		self.__name__ = 'rule'
	
	def associate_parent(self,parent):
		self.Parent = parent #Not to be added to conf file.

class CircosPlot():
	def __init__(self,basefilename,obj_list,karyotype,files_dict):
		self.files_dict = files_dict
		self.karyotype = karyotype
		self.basename = basefilename
		self.data_dict = {}
		self.last_filled_radius = 1.1
		self.master_struct = {'ideograms':{}} #Will be loop-able for jinja2 
		self.plots = ['<plots>']
		self.track_dict = {} #Temporary -- will be obviated by XML at later stage.
		self.add_ideogram()
		self.files_list = []
		self.add_top_lvl_params()
		self.obj_list = obj_list
		self.boilerplate = ['<<include colors_fonts_patterns.conf>>\n',
				    '<<include housekeeping.conf>>\n'
					]

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

	def add_ideogram(self,spacing_opt='0.005r',radius_opt='0.5r',thickness_opt='20p',fill_opt='yes',color_opt=[]):
		l = len(self.master_struct['ideograms']) + 1
		ideogramid = 'ideogram-'+str(l)
		self.master_struct['ideograms'][ideogramid] = {}
		self.add_spacing(spacing_opt,ideogramid)
		block = {'radius':radius_opt,
			 'thickness':thickness_opt,
			 'fill':fill_opt}
		self.master_struct['ideograms'][ideogramid].update(block)

	def add_links(self):
		pass

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

	def add_spacing(self,spacing_opt,ideogramid):
		self.master_struct['ideograms'][ideogramid].update({'spacing':{'default':spacing_opt}})

	def add_top_lvl_params(self):
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

	def create_base_conf_template(self):
		E = Environment(loader=PackageLoader('dummy','Templates'))
		T = E.get_template('template.conf')
		with open(self.basename + '.conf','w') as f:
			rendering = T.render(master_struct = self.master_struct, filename = self.basename + '.png', karyotype_filename = self.karyotype_filename)
			f.write(rendering)

	def _four_element_processing(self,file_,key):
		#Heatmaps, scatterplots, and histograms share the same 4-element format...
		with open(file_,'w') as handle:
			for f in self.raw_dict:
				if os.path.splitext(f)[1] in ('.wig','.bw','.bigWig'):
					for chromosome in self.raw_dict[f]:
						for feature in self.raw_dict[f][chromosome]:
							datastring = chromosome+' '+str(self.raw_dict[f][chromosome][feature]['start'])+' '+str(self.raw_dict[f][chromosome][feature]['end'])+' '+str(self.raw_dict[f][chromosome][feature]['val'])+'\n'
							handle.write(datastring)

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

	def update_master_conf(self):
		self.plots.append('</plots>')
		with open(self.mainconf,'a') as conf:
			for element in self.plots:
				conf.write(element+'\n')

	def write_conf(self):
		current_obj = None
		prev_obj = None
		f = open(self.basename+'.conf','w') 
		for elem in self.boilerplate:
			f.write(elem)
		f.write('karyotype = '+self.karyotype)
		for obj in self.obj_list:
			if obj.__name__ in ('ideogram','rule','image','plot','zoom','highlight','tick','link','karyotype'):
				if current_obj is not None:
					f.write('</'+current_obj+'>'+'\n')
				prev_obj = current_obj
				current_obj = obj.__name__
				if current_obj != prev_obj and prev_obj in ('plot','highlight'):
					f.write('</'+prev_obj+'s>\n')
				if current_obj in ('plot','highlight') and prev_obj != current_obj:
					f.write('<'+current_obj+'s>\n')
				f.write('\n<'+current_obj+'>'+'\n')
				if current_obj == 'ideogram':
					f.write('<spacing>\ndefault = 0.25r\n</spacing>\n')
				if current_obj in ('plot','highlight') and 'file' not in dir(obj):
					obj.file = self.files_dict[obj.type][0]
					self.files_dict[obj.type].pop(0)
			for att in dir(obj):
				if att[0:2] != '__' and att not in ['llo','boilerplate'] and getattr(obj,att) is not None and hasattr(att,'__call__') == False:
					if getattr(obj,att)[:2] == '__':
						val = str(getattr(obj,att)).strip()[len(str(getattr(obj,att)).strip())-6:]
						h = [val[0:2],val[2:4],val[4:6]]
						rgb = []
						for num in h:
							rgb.append(str(int(str(num),16)))
						s = ','
						val = s.join(rgb)
					else:
						val = str(getattr(obj,att)).strip()	
					if att[0:4] == obj.__name__[0:4]:
						if att[5:] in ['thickness','radius','r0','r1']:
							f.write(att[5:]+' = '+val.strip()+'r\n')
						else:	
							f.write(att[5:]+' = '+val.strip()+'\n')
					else:
						f.write(att+' = '+val.strip()+'\n')
				elif att == 'boilerplate':
					f.write(str(getattr(obj,att)).strip()+'\n')
		f.write('</'+current_obj+'>'+'\n')
		if current_obj in ('plot','highlight'):
			f.write('</'+current_obj+'s>\n')
		f.close()
	
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

if __name__ == "__main__":
	test_routine()

#TODO
"""
Add support for sub-objects.
LAST: Write paper as Application Note...
"""
