import argparse
from numpy import *

def vcf2hap(input, chr, output):
	f = open(input+".vcf", 'r')
	g1 = open(output+".hap", 'w')	
	g2 = open(output+".ind", 'w')
	g3 = open(output+".pos", 'w')

	A = []
	if chr != 'X':
		for a in f:
			if a[0] == "#" and a[1] == "C":
				a = a.split()
				for i in range(9, len(a)):
					g2.write(a[i]+"\n")
			elif a[0] != "#":
				a = a.split()
				g3.write(a[1]+"\t"+a[3]+"\t"+a[4]+"\n")
				A1 = []
				for i in range(9, len(a)):
					A1.append(a[i][0])
					A1.append(a[i][2])
				A.append(A1)
	else:
		for a in f:
			if a[0] == "#" and a[1] == "C":
				a = a.split()
				for i in range(9, len(a)):
					g2.write(a[i]+"\n")
			elif a[0] != "#":
				a = a.split()
				g3.write(a[1]+"\t"+a[3]+"\t"+a[4]+"\n")
				A1 = []
				for i in range(9, len(a)):
					A1.append(a[i][0])
					A1.append('9')
				A.append(A1)

	AM = mat(A)
	AM = AM.T


	for i in range(0, AM.shape[0]):
		for j in range(0, AM.shape[1]):
			g1.write(str(AM[i,j]))
		g1.write("\n")

	print("Done")

	f.close()
	g1.close()
	g2.close()
	g3.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", type = str, required = True, \
						help = "Prefix of input vcf")
	parser.add_argument("--chr", type = str, required = True, \
						help = "chromosome label, 1,2,...,22,X")
	parser.add_argument("--output", type = str, required = True, \
						help = "Prefix of output file(.hap, .ind, .pos)")

	args = parser.parse_args()
	if args.chr not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']:
		print("Error: not correct chr!")
		exit()

	vcf2hap(args.input, args.chr, args.output)
	
