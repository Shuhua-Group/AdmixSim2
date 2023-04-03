from functools import reduce
import argparse
import pandas
import gc
import math
from math import sqrt
from scipy import integrate
from scipy.special import gamma
import numpy as np

## 1989-Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism
def calculate_Dtajima(s, n, pi):
	## s: number of segregating sites
	## n: number of DNA sequences

	## Tajima's D
	a1 = reduce(lambda x,y: x+1.0/y, range(1,n))
	a2 = reduce(lambda x,y: x+1.0/pow(y,2), range(1,n))
	b1 = (n+1.0)/(3.0*n-3.0)
	b2 = 2.0*(pow(n,2)+n+3.0)/(9.0*n*(n-1.0))
	c1 = b1-1.0/a1
	c2 = b2-(n+2.0)/(a1*n)+a2/pow(a1,2)
	e1 = c1/a1
	e2 = c2/(pow(a1,2)+a2)
	D = (pi-s/a1)/sqrt(e1*s+e2*s*(s-1.0))
	## P value for Tajima's D, assuming that D follows the beta distrbution
	if n%2 == 0:
		Dmax = (n/(2.0*(n-1))-1.0/a1)/sqrt(e2)
	else:
		Dmax = ((n+1.0)/(2.0*n)-1.0/a1)/sqrt(e2)
	Dmin = (2.0/n-1.0/a1)/sqrt(e2)
	a = Dmin
	b = Dmax
	alpha = -(1.0+a*b)*b/(b-a)
	beta = (1.0+a*b)*a/(b-a)
	func = lambda d: gamma(alpha+beta)*pow(b-d,alpha-1.0)*pow(d-a,beta-1.0)/(gamma(alpha)*gamma(beta)*pow(b-a,alpha+beta-1.0))
	pvalue = 2*min(integrate.quad(func, Dmin, D)[0], integrate.quad(func, D, Dmax)[0])

	return D, pvalue	


def theta_pi_k(hapcount, s, n):
	pi = (hapcount['0'].values * hapcount['1'].values).sum()*1.0/(n*(n-1.0)/2.0)
	k = s*1.0/reduce(lambda x,y: x+1.0/y, range(1,n))

	return pi, k


### Nei, M., and Tajima, F. (1981). DNA POLYMORPHISM DETECTABLE BY RESTRICTION ENDONUCLEASES. Genetics 97, 145-163
def haplotype_diversity(haps):
	## Haplotype Diversity (H), H = N/(N-1)*(1-sigma(x^2))
	## x is the haplotype frequency of each haplotype
	## N is the sample size (haplotypes)
	## This measure of gene diversity is analogous to the heterozygosity at a single locus
	haplist = haps.apply(lambda x: "".join(list(x)), axis=0)	#assemble each hap to string
	nsample = haps.shape[1]
	sigmax2 = reduce(lambda x,y: x+pow(y,2), [0]+[z*1.0/nsample for z in list(haplist.value_counts())])
	nhap = len(set(haplist))
	H = nsample*1.0/(nsample-1.0)*(1.0-sigmax2)

	return nhap, H	 


def calculate_one_region_stat(hap_df, hapcount_df, nseq, chromid, start, end):
	## info in the given region
	haps = hap_df[(hap_df['#CHROM']==str(chromid)) & (hap_df['POS']>=int(start)) & (hap_df['POS'] <= int(end))].copy()
	hapcount = hapcount_df[(hapcount_df['#CHROM']==str(chromid)) & (hapcount_df['POS']>=int(start)) & (hapcount_df['POS']<=int(end))].copy()
	
	if haps.empty:
		nmarker = 0
		thetaPI = 0
		thetaK = 0
		seg = 0
		nhap = 0
		H = 0
		Dtajima = np.nan
		DtajimaP = np.nan
		return nmarker, thetaPI, thetaK, seg, nhap, H, Dtajima, DtajimaP	
	else:
		nmarker = hapcount.shape[0]
		haps.drop(['#CHROM','POS'], axis = 1, inplace = True)
		hapcount.drop(['#CHROM','POS'], axis = 1, inplace = True)
		
		## check non-biallelic (or missing) sites + homozygotes
		site2rm = list(hapcount[((hapcount['0']+hapcount['1']) < nseq) | (hapcount['1']==0) | (hapcount['0']==0)].index)
		if len(site2rm) == nmarker:
			thetaPI = 0
			thetaK = 0
			seg = 0
			nhap = 0
			H = 0
			Dtajima = np.nan
			DtajimaP = np.nan
			return nmarker, thetaPI, thetaK, seg, nhap, H, Dtajima, DtajimaP
		else:
			if len(site2rm) > 0:
				haps.drop(site2rm, inplace = True)
				hapcount.drop(site2rm, inplace = True)
			else:
				pass
			
		seg = hapcount.shape[0]		## num of segregating site
		thetaPI, thetaK = theta_pi_k(hapcount, seg, nseq)
		nhap, H = haplotype_diversity(haps)
		Dtajima, DtajimaP = calculate_Dtajima(seg, nseq, thetaPI) 
		return nmarker, thetaPI, thetaK, seg, nhap, H, Dtajima, DtajimaP


## remain required geno data
## convert to one snp a row one haplotype a column format
## count alleles
def convert_vcf(vcf, regionfile, window_shift, info):
	if window_shift == 'target_region':
		windowsize = 5000
	else:
		windowsize = int(window_shift.split('@')[0])

	region = pandas.read_csv(regionfile,sep='\s+',header=None,usecols=[0,1,2,3], names=['regionID','chr', 'start','end'])
	region['chr'] = region['chr'].astype(str)
	region['start'] = region['start']-windowsize
	region['end'] = region['end']+windowsize
	region.sort_values(by=['chr','start','end'], ascending=True, inplace=True)
	region.reset_index(inplace=True, drop=True)

	##	merge regions	
	if region.shape[0] == 1:
		pass
	else:
		for index in list(region.index)[:-1]:
			chrom1,start1,end1 = list(region.loci[index])[1:]
			chrom2,start2,end2 = list(region.loci[index+1])[1:]	
		
			if ((chrom2 == chrom1) & (start2 <= end1+1)):
				new_end = max(end1, end2)
				region.loci[index+1,'start'] = start1
				region.loci[index+1,'end'] = new_end
				region.drop(index,inplace=True)
			else:
				continue

	
	## extract vcf, and convert format
	mlist = [samplename for samplename in list(info.keys()) if info[samplename]==1]
	flist = [samplename for samplename in list(info.keys()) if info[samplename]==2]

	geno = pandas.concat(list(region.apply(lambda x: vcf[(vcf['#CHROM']==x['chr']) & (vcf['POS']>=x['start']) & (vcf['POS']<=x['end'])].copy(), axis=1)),ignore_index=True)

	if geno.empty:
		haps = pandas.DataFrame()
		hapcount = pandas.DataFrame()
		nseq = 0
		return haps, hapcount, nseq
	else:
		## convert format, for male individuals
		if ((len(mlist) > 0) & ('X' in list(geno['#CHROM']))):	
			nseq = len(mlist)+2*len(flist)
			geno[mlist] = geno[mlist].applymap(lambda x: x[0]).astype('category')
		else:
			nseq = 2*len(mlist)+2*len(flist)

		if nseq < 4:
			haps = pandas.DataFrame()
			hapcount = pandas.DataFrame()
			return haps, hapcount, nseq
		else:
			haps = pandas.DataFrame(geno.apply(lambda x: '|'.join(x[2:]).split('|'),axis=1).tolist(),columns=['hap'+str(x) for x in range(1,nseq+1)]).astype('category')
			count0 = (haps=='0').apply(sum,axis=1)
			count1 = (haps=='1').apply(sum,axis=1)
			hapcount = pandas.concat([count0, count1],axis = 1).astype('int64')
			hapcount.rename(columns=lambda x:str(x),inplace=True)
		
			haps['#CHROM'] = hapcount['#CHROM'] = geno['#CHROM'].values
			haps['POS'] = hapcount['POS'] = geno['POS'].values
		
			## compress
			haps = haps.astype({'#CHROM':'category','POS':'int32'})
			hapcount = hapcount.astype({'#CHROM':'category','POS':'int32'})
		
			return haps,hapcount, nseq


def split_window(regionID,chromID,start,end,window_shift):
	windowsize = int(window_shift.split('@')[0])
	stepsize = int(window_shift.split('@')[1])
	overlapsize = windowsize - stepsize

	length = end-start+1
	bin_num = max(int(math.ceil((length-overlapsize)*1.0/stepsize)),1)
	ex_len = bin_num*stepsize+overlapsize
	ex_start = int(max(start-(ex_len-length)/2.0, 1.0))
	ex_end = int(end+(ex_len-length)/2.0)

	region = pandas.DataFrame(columns = ['regionID', 'chr', 'start', 'end'])
	region['regionID'] = [regionID]*bin_num
	region['chr'] = chromID
	region['start'] = [ex_start+num*stepsize for num in range(bin_num)]
	region['end'] = region['start']+windowsize-1

	return region


def make_regions(regionfile, window_shift):
	region = pandas.read_csv(regionfile, sep='\s+', header=None, usecols=[0,1,2,3], names=['regionID','chr','start','end'])
	region['chr'] = region['chr'].astype(str)
	
	if window_shift == 'target_region':
		pass
	else:
		region['tmp'] = region.apply(lambda x: split_window(x['regionID'],x['chr'],x['start'],x['end'],window_shift), axis=1)
		region = pandas.concat(list(region['tmp']), ignore_index=True)
	region.sort_values(by=['chr','start','end'], ascending=True, inplace=True)
	
	return region


def fdr(pvaluelist):
	## numpy.array format
	## should be sorted, decreasing order (ascending Pvalues)
	n = len(pvaluelist)
	pvalues = pvaluelist[~np.isnan(pvaluelist)]
	if len(pvalues) <= 1:
		return list(pvaluelist)
	else:
		num = len(pvalues)
		adj_pvalues= pvalues*num/range(1,num+1)
		if adj_pvalues[-1] > 1.0:
			adj_pvalues[-1] = 1.0
		for i in range(num-2,-1,-1):
			adj_pvalues[i] = min(adj_pvalues[i+1], adj_pvalues[i])
		adj_pvalues = list(adj_pvalues)+[np.nan]*(n-num)
		return adj_pvalues


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", type = str, required = True, \
						help = "vcf file")
	parser.add_argument("--region", type = str, required = True, \
						help = "region file, 4 columns: <region ID> <chrom ID> <start pos> <end pos>, no header line, tab or space seperated")
	parser.add_argument("--samples", type = str, required = False, default = 'all', \
						help = "/included/sample/ID/list, 1 or 2 column: <sample ID> <gender 1/2, optional>, no header")
	parser.add_argument("--window_shift", type = str, required = False, default = 'target_region', \
						help = "windowsize@increment, for example, 50000@10000.")
	parser.add_argument("--out", type = str, required = False, default = 'out', \
						help = "Prefix of output file")
	args = parser.parse_args()

	## sample info	
	with open(args.vcf) as f:
		headerline = 0
		line = f.readline()
		while line[:2] == "##":
			headerline += 1
			line = f.readline()
		
	allsamplelist = line.strip().split('\t')[9:]
	if args.samples == 'all':
		samplelist = allsamplelist
		sampleinfo = {samplename:2 for samplename in samplelist}
	else:
		sampleinfo = pandas.read_csv(args.samples, header = None, sep ='\s+')
		sampleinfo = sampleinfo[smapleinfo[0].isin(allsamplelist)]
		if sampleinfo.shape[0] == 0:
			print("Error: No sample included.")
			exit()
		else:
			pass
	
		samplelist = list(sampleinfo[0])
		if sampleinfo.shape[1] == 1:
			sampleinfo = dict(zip(samplelist, [2]*len(samplelist)))	
		else:
			sampleinfo[1] = sampleinfo[1].astype(int)
			sampleinfo = dict(zip(samplelist, list(sampleinfo[1])))

	## input	
	## read str using Categorical dtypes, to save memory
	datatype = dict(zip(['#CHROM','POS']+samplelist, ['category','int32']+['category']*len(samplelist)))
	vcfdata = pandas.read_csv(args.vcf, sep = '\t', skiprows = range(headerline), usecols=['#CHROM','POS']+samplelist, dtype=datatype)

	## convert		
	hapdata, hapcountdata, nseq = convert_vcf(vcfdata, args.region, args.window_shift, sampleinfo)
	if nseq < 4:
		print('No enough sequences, at least 4 sequences are required.')
		exit()
	else:
		del(vcfdata)
		gc.collect()


	## get output
	result = make_regions(args.region, args.window_shift)
	result['#sequence'] = nseq
	result['tmp'] = result.apply(lambda x: calculate_one_region_stat(hapdata, hapcountdata, nseq, x['chr'],x['start'], x['end']),axis = 1)	

	result = pandas.concat([result[['regionID', 'chr', 'start', 'end', '#sequence']], pandas.DataFrame(result['tmp'].tolist(), columns=['#marker','ThetaPI', 'ThetaK', '#segregating', '#haplotype', 'Hap_diversity', 'Dtajima', 'Dtajima_P'], index=list(result.index))],axis=1)

	result.sort_values(by='Dtajima_P', ascending=True,inplace=True)
	result['Dtajima_adj.P'] = fdr(result['Dtajima_P'].values)
	result.sort_values(by=['chr','start','end'],ascending=True,inplace=True)

	result.to_csv(args.out+".stat", sep = "\t", index=None, na_rep='NA')

	print('Done.')

if __name__ == '__main__':
	main() 
