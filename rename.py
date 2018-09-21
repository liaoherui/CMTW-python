import re
import os

path='Result/KO_TWO/ERP004605'
#h={}
os.makedirs('ERP004605',0755)
for filename in os.listdir(path):
	f=open(path+'/'+filename,'r')
	o=open('ERP004605/'+filename,'w+')
	while True:
		line=f.readline()
		line=line.strip()
		if not line:break
		ele=line.split('\t')
		ele=ele[1:]
		n='\t'.join(ele)
		o.write(n+'\n')
		

exit()
for filename in os.listdir(path):
	#for filename2 in os.listdir(path+'/'+filename):
	h={}
	name=re.split('\.',filename)
	name=name[0]
	f=open(path+'/'+filename,'r')
	'''
	if not os.path.exists('Output/Gene/'+filename):
		os.makedirs('Output/Gene/'+filename,0755)
	'''
	f1=open('/mnt/osf1/data/MSV_config/profile/Phenotype_profile/All.xls','r')	
	while True:
		line=f1.readline()
		line=line.strip()
		if not line:break
		if not re.search('ERP002469',line):continue		
		ele=line.split('\t')
		#print ele[31]
		disease=re.sub(' ','',ele[31])
		oldname=disease+'_'+ele[0]
		newname=ele[1]
		h[oldname]=newname
	if name in h:
		o=open('Output/species/ERP002469/'+h[name]+'.xls','w+')
	else:
		print name+' is wrong ! '
		continue
	while True:
		line=f.readline()
		line=line.strip()
			#if True:break
		if not line:break
			#ele=line.split('\t')	
			#ele=ele[1:]
			#n='\t'.join(ele)
		o.write(line+'\n')
			#print ele 
			#exit()
		'''
		while True:
			line=f.readline()
			if not line:break
			o.write(line)
		'''

			
			
