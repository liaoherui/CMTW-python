#####	author:liaoherui	############
####	mail:liaoherui@mail.dlut.edu.cn	######
import re
import os
import numpy as np
import pandas as pd
import math
import getopt
import sys
import time

########### Get Option ##############
opts,args=getopt.getopt(sys.argv[1:],"hg:s:k:")
gl=''	 #gene_abundance_list
sl=''	#species_abundance_list
kl=''	#ko_abundance_list

for opt,arg in opts:
	if opt=='-h':
		print 'If you want to convert the data,you will need to give:\n\tgene/species/ko profile dir.\nAnd use this script as:\n\t python '+os.path.basename(sys.argv[0])+' -g(s/k) [gene/species/ko  dir ] '
	elif opt=='-g':
		gl=arg
	elif opt=='-s':
		sl=arg
	elif opt=='-k':
		kl=arg
########## Initialize anno hash ######
if not gl=='':
	anno={}
	f=open('/mnt/osf1/data/MSV_config/profile/Other_profile/annotation.summary','r')
	while True:
		line=f.readline()
       	 	if not line:break
	        if re.search('^ID',line):continue
	        eles=line.split('\t')
	        anno[eles[0]]=eles[1]
	f.close()
####### Gene related function ########
def stat_gene(gl,anno):
	pwd=os.getcwd()
	for filename in os.listdir(gl):
		if not re.search('\.xls',filename):continue
		pro=re.split('_',filename)
		pro=pro[0]
		if not os.path.exists(pwd+'/Result/Gene/'+pro):
			os.makedirs(pwd+'/Result/Gene/'+pro,0755)
		od=pwd+'/Result/Gene/'+pro
		fg1=pd.DataFrame(pd.read_table(gl+'/'+filename,sep='\t'))
		fg2=open(gl+'/'+filename,'r')
		while True:
			line=fg2.readline()
			line=line.strip()
			ele=line.split('\t')
			ele=ele[0:]  #first line of gene profile
			#print ele
			if True:break
		if re.search('Head',line):
			fg1.rename(columns={'Head':'Gene'},inplace=True)
		else:
			fg1.rename(columns={'Unnamed: 0':'Gene'},inplace=True)
		#print fg1.columns.tolist()
		fg1=fg1.set_index('Gene')
		for e in ele:
			fg1=fg1.sort_values(e)
			total=0
			for i in fg1[e]:
				total+=i
			og1=open(od+'/'+e+'.xls','w+')
			og1.write('Name\tPercent\n')
			gn=list(fg1[e][-100:].index)
			count=0
			for n in fg1[e][-100:]:
				value=n/total
				value=str('%.2f%%' % (value*100))
				if re.search('0.00%',value):continue
				num=gn[count]
				if str(num) in anno:
					og1.write(anno[str(num)]+'\t'+value+'\n')
					count+=1
				else:
					og1.write(num+'\t'+value+'\n')
					count+=1

###### Species related function #######
def stat_taxo(sl):
	pwd=os.getcwd()
	for filename1 in os.listdir(sl):
		if re.search('genus',filename1):
			d=pwd+'/Result/Taxo/genus'
			if not os.path.exists(d):
				os.makedirs(d,0755)
		if re.search('phylum',filename1):
                        d=pwd+'/Result/Taxo/phylum'
                        if not os.path.exists(d):
                                os.makedirs(d,0755)
		if re.search('species',filename1):
                        d=pwd+'/Result/Taxo/species'
                        if not os.path.exists(d):
                                os.makedirs(d,0755)
		le=sl+'/'+filename1
		for filename2 in os.listdir(le):
			pre=re.split('_',filename2)	
			pre=pre[0]
			if not os.path.exists(d+'/'+pre):
				os.makedirs(d+'/'+pre,0755)
			f1=pd.DataFrame(pd.read_table(le+'/'+filename2,sep='\t'))
			f2=open(le+'/'+filename2,'r')
			while True:
				line=f2.readline()
				line=line.strip()
				ele=line.split('\t')	
				ele=ele[1:]
				#print ele
				#exit()
				if True:break
			f1.rename(columns={'#Header':'Taxo'},inplace=True)
			f1=f1.set_index('Taxo')
			#print f1['SRR3131692'][0:10]
			#exit()
			for e in ele:
				f1=f1.sort_values(e)
				total=0
				for i in f1[e]:
					total+=i
				ot=open(d+'/'+pre+'/'+e+'.xls','w+')
				tn=list(f1[e][-100:].index)
				count=0
				for n in f1[e][-100:]:
					value=n/total
					value=str('%.2f%%' % (value*100))
					if re.search('0.00%',value):continue
					name=tn[count]
					na=name.split()
					name='_'.join(na)
					ot.write(name+'\t'+value+'\n')
					count+=1
##### KO related function ########
def stat_ko(kl):
	pwd=os.getcwd()
	anno={}
	t=0
	fa=open('profile/anno_path.xls','r')
	while True:
		line=fa.readline()
		if True:break
	while True:
		line=fa.readline()
		line=line.strip()
		if not line:break
		ele=line.split('\t')
		if ele[-1] not in anno:
			anno[ele[-1]]=ele[-2]
		else:
			print ele[-1]+'has more than one  level3'
			anno[ele[-1]]+=ele[-2]
	for filename in os.listdir(kl):
		km={}
		mm={}  #key -> map, value -> ko plus matrix
		###  convert ko to pathway profile ######
		f=open(kl+'/'+filename,'r')
		if not os.path.exists('KO_to_Path'):
			os.makedirs('KO_to_Path',0755)
		name=re.sub('KO','Path',filename)
		o=open('KO_to_Path/'+name,'w+')
		while True:
			line=f.readline()
			o.write(line)
			if True:break
		while True:
			line=f.readline()
			line=line.strip()
			if not line:break
			ele=line.split('\t')	
			km[ele[0]]=[]
			for e in ele[1:]:
				km[ele[0]].append(float(e))
			km[ele[0]]=np.array(km[ele[0]])
		fp=open('profile/prokaryote_ko_map.tab','r')	
		while True:
			line=fp.readline()
			line=line.strip()
			if not line:break
			ele=line.split('\t')
			if ele[0] not in km:
				continue
			em=ele[1].split()
			for e in em:
				ma='map'+str(e)
				if ma not in mm:
					mm[ma]=km[ele[0]]
				else:
					mm[ma]+=km[ele[0]]
		for key in mm:		
			o.write(key)
			for e in list(mm[key]):
				o.write('\t'+str(e))
			o.write('\n')
		o.close()
		##### stat the pathway profile #####
		for filename in os.listdir(pwd+'/KO_to_Path'):
			#new=open(pwd+'/KO_to_Path/'+filename,'r')	
			pro=re.split('_',filename)	
			pro=pro[0]
			if not os.path.exists('Result/KO/'+pro):
				os.makedirs('Result/KO/'+pro,0755)
			ok='Result/KO/'+pro
			fk1=pd.DataFrame(pd.read_table('KO_to_Path/'+filename,sep='\t'))
			fk2=open('KO_to_Path/'+filename,'r')
			while True:
				line=fk2.readline()
				line=line.strip()
				ele=line.split('\t')
				ele=ele[1:]
				if True:break
			fk1.rename(columns={'#Header':'Map'},inplace=True)
			fk1=fk1.set_index('Map')
			for e in ele:
				fk1=fk1.sort_values(e)
				total=0
				for i in fk1[e]:
					total+=i
				ok1=open(ok+'/'+e+'.xls','w+')
				kn=list(fk1[e][-100:].index)
				count=0
				for n in fk1[e][-100:]:
					value=n/total
					value=str('%.2f%%' % (value*100))
					if re.search('0.00%',value):continue
					num=kn[count]
					if str(num) in anno:
						ok1.write(anno[str(num)]+'\t'+value+'\n')
						count+=1
					else:
						ok1.write('NA'+'\t'+value+'\n')
						count+=1
if not gl=='':	
	stat_gene(gl,anno)
if not sl=='':
	stat_taxo(sl)
if not kl=='':
	stat_ko(kl)
