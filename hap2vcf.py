import argparse
from numpy import *

def hap2vcf(input, chr, datatype, output):
	f1 = open(input+".hap", 'r')
	f2 = open(input+".ind", 'r')
	f3 = open(input+".pos", 'r')
	g = open(output+".vcf", 'w')

	hapdata = []
	for a in f1:
		if a[0] == '9':
			hapdata.append(hapdata[-1])	
		else:
			temphap = []
			for i in range(0, len(a)-1):
				temphap.append(a[i])
			hapdata.append(temphap)

	hapdataM = mat(hapdata)
	hapdataM = hapdataM.T

	indlist = []
	for a in f2:
		a = a.split()
		indlist.append(a[0])

	posinfo = {}
	poslist = []
	for a in f3:
		a = a.split()
		poslist.append(a[0])
		posinfo[a[0]] = (a[1], a[2])

	if (hapdataM.shape[0] != len(poslist)):
		print("Error: position number not equal to that in hap file!")
		exit()
	if (int(hapdataM.shape[1]/2) != len(indlist)):
		print("Error: individual number not equal to that in hap file!")
		exit()

	g.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">\n")
	g.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
	for i in range(0, len(indlist)-1):
		g.write(indlist[i]+"\t")
	g.write(indlist[len(indlist)-1]+"\n")

	for i in range(0, len(poslist)):
		pos = poslist[i]
		g.write(chr+"\t"+pos+"\t.\t"+posinfo[pos][0]+"\t"+posinfo[pos][1]+"\t.\t.\t.\tGT\t")
		if datatype == "A":
			for j in range(0, len(indlist)):
				for k in [0,1]:
					if hapdataM[i,2*j+k] == posinfo[pos][0]:
						g.write('0')
					elif hapdataM[i,2*j+k] == posinfo[pos][1]:
						g.write('1')
					else:
						g.write(hapdataM[i,2*j+k])
					if k == 0:
						g.write('|')
				if j != len(indlist)-1:
					g.write("\t")
			g.write("\n")
		else:
			for j in range(0, len(indlist)):
				g.write(hapdataM[i,2*j]+'|'+hapdataM[i,2*j+1])
				if j != len(indlist)-1:
			 		g.write("\t")
			g.write("\n")	

	print("Done")

	f1.close()
	f2.close()
	f3.close()
	g.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", type = str, required = True, \
						help = "Prefix of hap, ind, pos file")
	parser.add_argument("--chr", type = str, required = True, \
						help = "chromosome label, 1,2,..,22,X")
	parser.add_argument("--datatype", type = str, required = True, \
						help = "hap datatype, A for AGCT, 0 for 01")
	parser.add_argument("--output", type = str, required = True, \
						help = "Prefix of output vcf file")

	args = parser.parse_args()
	if args.chr not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']:
		print("Error: not correct chr!")
		exit()

	if args.datatype not in ['A', '0']:
		print("Error: not correct datatype!")
		exit()

	
	hap2vcf(args.input, args.chr, args.datatype, args.output)
	
