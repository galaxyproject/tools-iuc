from pprint import pprint
from jinja2 import Template, Environment, PackageLoader
from BCBio import GFF
from subprocess import call
from Bio.Seq import Seq
import itertools
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
	for argument in ['-k','-H','-I','-s','-b','-l']:
		A.add_argument(argument,nargs='*')
	args = A.parse_args()
	X = args.xml[0] 
	T = xml.etree.ElementTree.parse(X)
	R = T.getroot()
	D = DataInterpreter(R,[])
	D.files_dict = {'karyotype':args.k,
			'scatter':args.s,
			'heatmap':args.H,
			'histogram':args.b,
			'highlight':args.I,
			'link':[args.l[0]]}
	D.assoc_chromosomes = args.l[1:]
	D._parse_files()
	P = CircosPlot('christmas',D.object_list,D.karyotype,D.files_dict)
	P.write_conf()

class DataInterpreter():
	def __init__(self,xmlroot,files_dict={}):
		self.xmlroot = xmlroot
		self._create_object_list()
		self.files_dict = {} 
		if self.files_dict == {}:
		#FIXME
			files_dict = self._get_files_from_xml()
		#self._parse_files()

	def _create_object_list(self):
		self.object_list = []
		current_obj = None 
		current_sub = None
		subobjind= 0
		recognized_objects = ['ideogram','image','link','tick','zoom','highlight','plot']
		recognized_subobjects = ['pairwise','rule','break_style','spacing']
		for element in self.xmlroot.iter():
			if element.tag in recognized_objects:
				if current_sub is not None:
					current_obj.llo.append(copy.copy(current_sub))
				if current_obj is not None:
					self.object_list.append(copy.copy(current_obj))
				current_obj = CircosObj(element.tag) 
				subobjind = 0
			elif current_obj is not None:
				if element.tag in recognized_subobjects:
					subobjind = 1
					if current_sub is not None:
						current_obj.llo.append(copy.copy(current_sub))
					current_sub = CircosObj(element.tag)
				if subobjind == 0:
					setattr(current_obj,element.tag,element.text)
				elif subobjind == 1:
					setattr(current_sub,element.tag,element.text)
		self.object_list.append(copy.copy(current_obj))
	
	def _get_files_from_xml(self):
		#TODO When the time comes to create the final Galaxy tool...
		pass

	def _make_karyotype(self,f,filename="karyotype"):
		self.karyotype = filename+'.txt'
		self.chromosomes = []
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
					self.chromosomes.append(cname)
				else:
					seqlen += len(line.strip())
			g.write(space.join(['chr','-',cname,str(i),str(0),str(seqlen),COLORS_LIST[i]+'\n']))
			g.close()
		return self.karyotype

	def parse_backbone(self,key,f,index):
		filename = str(key)+'-'+str(index)+'.txt'
		chrom_list = [self.chromosomes[int(c)] for c in self.assoc_chromosomes]
		with open(f,'r') as input_handle:
			g = open(filename,'w')
			for line in input_handle:
				links = []
				i = 0
				data = line.split()
				for element in range(0,len(data),2):
					if data[element] != '0' and data[element][0]!='s':
						links.append((chrom_list[element/2],data[element],data[element+1]))
				for n in itertools.combinations([i for i in range(len(links))],2):
					for j in n:
						g.write(' '.join([word for word in links[j]]))
						g.write(' ')
					g.write('\n')
			g.close()
		return filename 

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
		tmp = {'karyotype':[],'scatter':[],'heatmap':[],'histogram':[],'highlight':[],'link':[]}
		i = 1
		tmp['karyotype'] = self._make_karyotype(self.files_dict['karyotype'][0])
		for key in self.files_dict:
			if self.files_dict[key] is not None:
				for f in self.files_dict[key]:
				#TODO Each parse method needs to update self.files_dict to a text file compatible with Circos.
					ext = os.path.splitext(f)[1]
					if key == "karyotype":
						continue
					elif key in ("link","scatter","heatmap","histogram","highlight"):
						if ext in ('.gff','.gff3'):
							tmp[key].append(self.parse_gff3(key,f,i))
						elif ext in ('.wig','.bigWig','.bw'):
							tmp[key].append(self.parse_bigWig(key,f,i))
						elif ext == ".bed":
							tmp[key].append(self.parse_bed(key,f,i))
						elif ext == ".backbone":
							tmp[key].append(self.parse_backbone(key,f,i))
						else:
							raise Exception('Unsupported File Format')
						i+=1
					else:
						raise Exception('Unsupported Data Type')
		self.files_dict = tmp

	def parse_gff3(self,key,f,index):
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
	def __init__(self,name):
		self.llo = []
		self.__name__ = name
		self.type = name
		if self.__name__ == 'image':
			self.boilerplate = '<<include etc/image.conf>>'
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
				    '<<include housekeeping.conf>>\n']

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
				if current_obj != prev_obj and prev_obj in ('link','plot','highlight'):
					f.write('</'+prev_obj+'s>\n')
				if current_obj in ('plot','highlight','link') and prev_obj != current_obj:
					f.write('<'+current_obj+'s>\n')
				f.write('\n<'+current_obj+'>'+'\n')
				if current_obj == 'ideogram':
					f.write('\n'.join(['<spacing>','default = 0.25r','</spacing>\n']))
				if current_obj in ('plot','highlight','link') and 'file' not in dir(obj):
					obj.file = self.files_dict[obj.type][0]
					self.files_dict[obj.type].pop(0)
			for att in dir(obj):
				if att[0:2] != '__' and att not in ['llo','boilerplate','rules'] and getattr(obj,att) is not None and hasattr(att,'__call__') == False:
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
		if current_obj in ('plot','highlight','link'):
			f.write('</'+current_obj+'s>\n')
		f.close()

if __name__ == "__main__":
	test_routine()

#TODO
"""
Add support for sub-objects.
"""
