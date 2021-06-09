'''
SEGMENT
Initialization: left end, right end, label of original individual, mutations' position, alleles of markers under selection.
Function:inserting mutations, adding labels, updating makers under selection and their status, separating segments and extract parts.
'''
class SEGMENT:
	def __init__(self, left, right, lab, mut, sel):
		self.left = left 	#int,left end(inside)
		self.right = right 	#int, right end(inside)
		self.mut = mut		#list, mutations carried by the segment
		self.lab = tuple(lab) 		#list, a list of individuals' labels (each label is a tuple)
		self.sel = sel		#dict, keys are physical positions and values are alleles.
	#insert mutation 
	def update_mut(self,x):
		bisect.insort(self.mut, x)
	#add label
	def update_lab(self,x):
		self.lab = tuple(x)
	#update makers under selection and their status
	def add_sel(self,x,y):
		self.sel[x] = y
	#separate segments and extract parts
	def separate(self,x,y):
		return SEGMENT(x, y, self.lab, [i for i in self.mut if i>=x and i<=y], {key:value for key,value in self.sel.items() if key>=x and key<=y})

'''
HAPLOTYPE
Initialization: status, segments
Function: generating mutation and inserting them to the corresponding segments
'''
class HAPLOTYPE:
	def __init__(self, status, seg):
		self.status = status 	#int, haplotype status，0 normal，1 abnormal（exist in male individual while simulating X chromosome, no recombination,s no mutations, no selctions)
		self.seg = seg 			#list, consist of classes of all segment
	# generating novel mutations and insert them to corresponding segments
	def add_mut(self):
		global gl_mutation
		if self.status == 1: 	#return without any processes for abnormal haplotype
			return
		k = PoissonRandom(gl_mut_dis[-1])
		n = 0
		start = [i.left for i in self.seg]
		while n<k:
			x = PositionRandom(gl_mut_dis[-1], 1)[0]
			y = Map_pos(gl_mut_dis, x)
			if gl_pos[y[1]] == y[0]:
				continue
			z = bisect.bisect_left(gl_mutation, y[0])
			if z != len(gl_mutation):
				if gl_mutation[z] == y[0]:
					continue	#Poisson distribution used for counts of novel mutations.
			bisect.insort(gl_mutation, y[0])
			m = bisect.bisect_right(start, y[0])
			self.seg[m-1].update_mut(y[0])
			n += 1
		return

	#assemble alleles or ancestral origins of all makers under selections
	def assemle_sel(self,x,y):
		if self.status == 1:
			return []
		tempsel = {key:value for i in self.seg for key, value in i.sel.items()}			#combine all dicts
		if y == 0:																		
			return [tempsel[i] for i in x]																		
		else:
			return [i.lab[0] for i in self.seg for j in i.sel.items() if j[0] in x]	
'''
INDIVIDUAL
initialization, haplotypes, fitness, label, parents, sex
function: calculating fitness, adding mutation, adding label, generating a haplotype of a gamete, cutting segment
'''
class INDIVIDUAL:
	def __init__(self, hap, coe, index, parents, sex):
		self.hap = hap 			#list, two haplotypes
		self.coe = coe 			#float, fitness
		#pop,generation,number
		self.index = index		#tuple，individual's label(population name，generation name，individual name）
		self.parents = parents 	#list ,labels of parents
		self.sex = sex			#character，individual's sex，0 for unknown, 1 for male, 2 for female
	def update_coe(self, x):
		for i in x:
			for hapsel in [k.assemle_sel(i[0], i[2]) for k in self.hap]:	#assembled segments on postions under selection
				if hapsel in i[1]:										#for all conditions
					if self.sex == '2':										#add the second fitness for female
						self.coe += i[3][1]
					else:													#other scenarios
						self.coe += i[3][0]
				if i[4] == '2':											#if dominant model
					break
	#adding novel muatations for two haplotypes
	def add_mut(self):
		for i in self.hap:
			i.add_mut()
		
	#adding labels for all segments of two haplotypes
	def update_lab(self):
		for i,j in enumerate(self.hap):
			for k in j.seg:
				k.update_lab(self.index + tuple([i+1]))

	#extract segments between two adjacent recombination postions from a list segments of a haplotype
	def Cutsegment(self, seglist, x, y):
		a = bisect.bisect_right([i.left for i in seglist], x) - 1 	
		b = bisect.bisect_right([i.left for i in seglist], y) - 1 	
		if a == b:
			parttempseg = copy.deepcopy([seglist[a].separate(x+1,y)])
		else:
			if x == seglist[a].right or x == gl_pos[0]-1:
				k = []
			else:
				k = [seglist[a].separate(x+1,seglist[a].right)]
			parttempseg = copy.deepcopy(k + seglist[a+1:b] + [seglist[b].separate(seglist[b].left,y)]) 	
		return parttempseg

	def Recombination(self):
		if gl_rec == '0' or (gl_chr == 'X' and self.sex == '1'):
			return copy.deepcopy(self.hap[np.random.randint(0,2)])
		m = PoissonRandom(gl_gen_dis[-1])
		n = 0
		nbk = [gl_pos[0]-1, gl_pos[-1]]
		while n<m :
			x = PositionRandom(gl_gen_dis[-1], 1)
			y = Map_pos(gl_gen_dis, x)[0]
			bisect.insort(nbk, y)
			n += 1
										
		k = np.random.randint(0,2)												#determine the index of haplotype which the segment before the first segement comes from
		tempseg = []
		for i  in range(len(nbk)-1):											
			tempseg += self.Cutsegment(self.hap[(i+k)%2].seg, nbk[i], nbk[i+1])				
		H = HAPLOTYPE(0,tempseg)										
		return H

'''
GENERATION
variable: population name, generation name, positions,conditions,type,coefficient, mode of selection, Ne, ancestral proportion
function: setting Ne, setting ancestral proportion, adding individuals, adding selection contents, summarizing frequencies of haplotypes satisfying each selection conditions
 		calculating fitness, adding labels, adding novel mutations, summarizing all present nutations in this generation, sampling ancestral populations, sampling individuals in a populations.
'''
class GENERATION:
	def __init__(self, popname, generation):
		self.popname = popname 			#characters，population name
		self.generation = generation 	#int, generation
		self.Sel = [] 					#list，a list of selection content(tuple, (postions, conditions, type, coefficient, mode))
		self.Ind = {} 					#dict，individual name: INDIVIDUAL
		self.Ne = -9 					#list, one or two value
		self.Proportion = {} 			#dict, population name: proportion
	
	def set_Ne(self,x):
		if self.Ne != -9 and x !=self.Ne: 	#检查是否已经设置，若已设置是否相同，不同则报错
			print('Error: Contradiction in setting Ne of population '+self.popname+' in generation '+str(self.generation)+': '+','.join(self.Ne)+' '+','.join(x))
			exit()
		self.Ne = x
	
	def set_Proportion(self,x):
		if self.Proportion != {} and x != self.Proportion:	#检查是否已经设置，若已设置是否相同，不同则报错
			print('Error: Contradiction in setting proportion of population '+self.popname+' in generation '+str(self.generation)+': ')
			print(self.Proportion)
			print(x)
			exit()
		self.Proportion = x
	
	def add_Ind(self,x,y):
		self.Ind[x] = y
	
	def add_Sel(self,x):
		if x in self.Sel: 	
			return
		for i in [j for j in self.Sel if j[0] == x[0]]: 
			if i != x:
				print('Warning: Selections happend on indentical postion with different arguments: ')
				print(x, i)
		self.Sel.append(x)

	def cal_freq(self, x, y):
		if len(x) == 0:
			return 0
		else:
			return round(x.count(y)/len(x),4)
	def Sel_Freq(self,x):
		tempfreq = []
		for i in x:
			if gl_sex == 0:
				tempseg = [k.assemle_sel(i[0], i[2]) for j in self.Ind.values() for k in j.hap]
				tempfreq.append([(self.cal_freq(tempseg, y),'/','/') for y in i[1]])
			else:
				tempseg1 = [k.assemle_sel(i[0], i[2]) for j in self.Ind.values() if j.sex == '1' for k in j.hap]
				tempseg2 = [k.assemle_sel(i[0], i[2]) for j in self.Ind.values() if j.sex == '2' for k in j.hap]
				tempseg = tempseg1 + tempseg2
				tempfreq.append([(self.cal_freq(tempseg, y), self.cal_freq(tempseg1, y), self.cal_freq(tempseg2, y)) for y in i[1]])
		return tempfreq

	def update_coe(self):
		for i in self.Ind.values():
			i.update_coe(self.Sel)
		self.Selfreq = self.Sel_Freq(self.Sel)

	def update_lab(self):
		for i in self.Ind.values():
			i.update_lab()
		
	def add_mut(self): 
		for i in self.Ind.values():
			i.add_mut()

	def summary_mut(self):
		tempmut = []
		for i in self.Ind.values():
			for j in i.hap:
				tempmut += [x for k in j.seg for x in k.mut]
		return tempmut

	def IndSampling(self, x, y, w):		
		while True:
			tempind = [i for i in self.Ind.values() if i.sex in x]
			if tempind != []:
				break
			if self.Ne == -9:
				print('Error: missing some sex in the input population '+str(self.popname))
				exit()
			else:
				print('Missing some sex in the generation '+str(self.generation)+' of the population '+str(self.popname))
				print('Re-simulating the generation '+str(self.generation)+' of the population '+str(self.popname))
				self.Reproduce()
				if gl_sel == '1':
					self.update_coe()
				if (self.popname,self.generation) in list(zip(gl_outpop,gl_outgen)):
					print('record',(self.popname,self.generation))
					gl_out[(self.popname,self.generation)] = copy.deepcopy(gl_pop[self.popname].generation[self.generation])
				if gl_mut == '1':
					self.add_mut()

		tempcoe = np.array([1+i.coe for i in tempind])
		sumcoe = tempcoe.sum()
		probcoe = np.array(tempcoe/sumcoe)
		z = np.random.choice(tempind, size = y, replace = w, p = probcoe)

		return z

	def PopSampling(self):
		temppop = list(self.Proportion.keys())
		sumprob1 = sum([i[0] for i in self.Proportion.values()])
		sumprob2 = sum([i[1] for i in self.Proportion.values()])
		x1 = np.random.choice(temppop, sum(self.Ne), True, [self.Proportion[i][0]/sumprob1 for i in temppop])
		x2 = np.random.choice(temppop, sum(self.Ne), True, [self.Proportion[i][1]/sumprob2 for i in temppop])

		return x1, x2

	def Reproduce(self):
		if sum(self.Ne) == 0:
			return 
		fatherpop, motherpop = self.PopSampling()
		fathertemppop = [i for i,j in self.Proportion.items() if j[0] != 0]
		mothertemppop = [i for i,j in self.Proportion.items() if j[1] != 0]
		fatherind = {x: gl_pop[x].generation[self.generation - 1].IndSampling(['0','1'], sum(fatherpop == x), True) for x in fathertemppop if self.generation - 1 in gl_pop[x].generation.keys()}
		motherind = {x: gl_pop[x].generation[self.generation - 1].IndSampling(['0','2'], sum(motherpop == x), True) for x in mothertemppop if self.generation - 1 in gl_pop[x].generation.keys()}
		fatheridx = {x:0 for x in fathertemppop}
		motheridx = {x:0 for x in mothertemppop}
		if len(self.Ne) == 1:
			k = 0
			fathermut = []
			mothermut = []
			offspringmut = []
			for i in range(self.Ne[0]):				
				father = fatherind[fatherpop[k]][fatheridx[fatherpop[k]]]
				fatheridx[fatherpop[k]] += 1
				mother = motherind[motherpop[k]][motheridx[motherpop[k]]]
				motheridx[motherpop[k]] += 1
				hap2 = father.Recombination()
				hap1 = mother.Recombination()
				if gl_chr == 'X':
					if hap2.status == 1:
						tempsex = '1'
					else:
						tempsex = '2'
				elif gl_sex == 1:
					tempsex = str(np.random.randint(1,3))
				else:
					tempsex = '0'
				k += 1
				self.add_Ind(k, INDIVIDUAL([hap1,hap2], 0, (self.popname, self.generation, k), [father.index, mother.index],tempsex))
		else:
			k = 0
			for i in range(self.Ne[0]):
				father = fatherind[fatherpop[k]][fatheridx[fatherpop[k]]]
				fatheridx[fatherpop[k]] += 1
				mother = motherind[motherpop[k]][motheridx[motherpop[k]]]
				motheridx[motherpop[k]] += 1
				if gl_chr == 'X': 
					hap2 = father.hap[1]
				else:
					hap2 = father.Recombination()
				hap1 = mother.Recombination()
				k += 1
				self.add_Ind(k, INDIVIDUAL([hap1,hap2], 0, (self.popname, self.generation, k), [father.index, mother.index],'1'))

			for i in range(self.Ne[1]):
				father = fatherind[fatherpop[k]][fatheridx[fatherpop[k]]]
				fatheridx[fatherpop[k]] += 1
				mother = motherind[motherpop[k]][motheridx[motherpop[k]]]
				motheridx[motherpop[k]] += 1
				if gl_chr == 'X': 
					hap2 = father.hap[0]
				else:
					hap2 = father.Recombination()
				k += 1
				self.add_Ind(k, INDIVIDUAL([hap1,hap2], 0, (self.popname, self.generation, k), [father.index, mother.index],'2'))
		return

	def Delete(self):
		del self.Ind

class POPULATION:
	def __init__(self, popname):
		self.popname = popname
		self.generation = {}
	def set_generation(self,x,y):
		self.generation[x] = y

	def Founder(self):
		k = min(list(self.generation.keys()))
		self.set_generation(k-1, GENERATION(self.popname, k-1))
		for j in gl_indinfo[self.popname]:
			if gl_chr == 'X' and j[2] == '1':
				self.generation[k-1].add_Ind(j[0], INDIVIDUAL([HAPLOTYPE(0, [SEGMENT(gl_pos[0], gl_pos[-1], [], [], {})]),HAPLOTYPE(1, [SEGMENT(0, 0, [], [], {})])], 0, (self.popname, k-1, j[0]), (0,0), j[2]))
			else:
				self.generation[k-1].add_Ind(j[0], INDIVIDUAL([HAPLOTYPE(0, [SEGMENT(gl_pos[0], gl_pos[-1], [], [], {})]),HAPLOTYPE(0, [SEGMENT(gl_pos[0], gl_pos[-1], [], [], {})])], 0, (self.popname, k-1, j[0]), (0,0), j[2]))
		self.generation[k-1].update_lab()

def Check_value_int(x,y):
	try: 
		int(x)
	except:
		print('Error: Invalid '+y+': '+x)
		exit()

def Check_value_float(x,y):
	try: 
		float(x)
	except:
		print('Error: Invalid '+y+': '+x)
		exit()

def Map_pos(x, y):
        k = bisect.bisect_left(x, y)
        if k == 0:
                z = gl_pos[0]
        else:
                z = gl_pos[k-1] + int((gl_pos[k] - gl_pos[k-1])*(y - x[k-1])/(x[k] - x[k-1]))   
        return z, k

def assemble_pos(x):
	x = x.split(',')
	temppos = []
	for i in x:
		if '-' not in i:
			temppos.append(int(i))
		else:
			temppos += list(gl_pos[bisect.bisect_left(gl_pos, int(i.split('-')[0])):bisect.bisect_right(gl_pos, int(i.split('-')[1]))])
	return temppos

def PoissonRandom(x):
        return np.random.poisson(x)

def PositionRandom(x, y):
        return x*np.random.rand(y)
	
def Prepare_Get_Check_Arg():
        opts, args = getopt.getopt(sys.argv[1:],'hi:o:p:g:n:',['help','no-rec','no-mut','no-sel','hap=','snv=','mod=','ind=','out=','rec-rate=','mut-rate=','chr='])
        x = [i[0] for i in opts]
        if '-h' in x or '--help' in x or opts == []:
                Prepare_Print_Help()
                exit()
        return opts, args

def Print_Header():
        print('===================================================')
        print('AdmixSim v2.0 [Alpha/2020-05-30]')
        print('===================================================')
        print()
        return

def Prepare_Print_Help():
        print('--hap [filename]\t(required) the haplotype file containing 2*N haplotypes for N ancestral individuals')
        print('--mod [filename]\t(required) the model file describing admixture events and selection events')
        print('--snv [filename]\t(required) the snp file providing a genetic distance and a mutation rate of simulated each site')
        print('--ind [filename]\t(required) the information of population and gender of n ancestral individuals')
        print('-o/--out [prefix]\t(required) the prefix of output files (default = current directory)')
        print('-p [name1,name2,...]\t(optional) names of populations for output (default=the final generated population)')
        print('-g [generation1,generation2,...] generations of populations corresponding to \'-p\'for output (default = the last generation in the modfile)')
        print('-n [number1,number2,...]\t(optional) numbers of populations corresponding to \'-p\'for output (default = 10)')
        print('--no-rec\tno recombination in the simulation process')
        print('--no-mut\tno mutation in the simulation process')
        print('--no-sel\tno selection in the simulation process')
        print('--rec-rate\ta general recombination rate across the simulated genome while ignoring genetic distances in the snp file')
        print('--mut-rate\ta general mutation rate across the simulated genome while ignoring mutation rates in the snp file')
        print('--chr [character] \t\tset the type of simluated genomes. \'A\' for autosomes, \'X\' for X chromosomes')
        print('-h/--help\tprint this help')
        return

def Prepare_Set_Value(opts):
	global gl_hapfile, gl_snpfile, gl_modfile, gl_indfile, gl_outpop, gl_outgen, gl_outnum, gl_outprefix, gl_rec, gl_mut, gl_sel, gl_recrate, gl_mutrate, gl_rateunit, gl_chr

	gl_rateunit = 10**-3 *3
	gl_snpfile = ''
	gl_indfile = ''
	gl_outpop = []
	gl_outgen = []
	gl_outnum = []
	gl_recrate = -9
	gl_mutrate = -9
	gl_outprefix = './out'
	gl_rec = '1'            
	gl_mut = '1'            
	gl_sel = '1'            
	gl_chr = 'A'            

	for opt_name,opt_value in opts:
		if opt_name == '-i':
			gl_hapfile = opt_value + '.hap'
			gl_modfile = opt_value + '.mod'
			gl_snpfile = opt_value + '.snv'
			gl_indfile = opt_value + '.ind'
		if opt_name == '--hap':
		    gl_hapfile = opt_value
		if opt_name == '--snv':
		    gl_snpfile = opt_value
		if opt_name == '--mod':
		    gl_modfile = opt_value
		if opt_name == '--ind':
		    gl_indfile = opt_value
		if opt_name in ('-o', '--out'):
		    gl_outprefix = opt_value
		if opt_name == '-p':
		    gl_outpop = opt_value.split(',')
		if opt_name == '-g':
			tempopt = opt_value.split(',')
			for i in tempopt: 
				Check_value_int(i, 'output generation')
			gl_outgen = list(map(int, tempopt))
		if opt_name == '-n':
			tempopt = opt_value.split(',')
			for i in tempopt: 
				Check_value_int(i, 'output individual number')
			gl_outnum = list(map(int, tempopt))
		if opt_name == '--no-rec':
		    gl_rec = '0'
		if opt_name == '--no-mut':
		    gl_mut = '0'
		if opt_name == '--no-sel':
		    gl_sel = '0'
		if opt_name == '--rec-rate':
		    Check_value_float(opt_value, 'global recombination rate')
		    gl_recrate = float(opt_value)
		if opt_name == '--mut-rate':
		    Check_value_float(opt_value, 'global mutation rate')
		    gl_mutrate = float(opt_value)
		if opt_name == '--chr':
			if opt_value not in ['A','X']:
				print('Invalid chromosome: '+opt_value)
				exit()
			gl_chr = opt_value

	if gl_outpop != []:
		tempout = list(zip(gl_outpop, gl_outgen))
		if len(tempout) != len(list(set(tempout))):
			print('Warning: Repeat output populaitons and generations fuond: '+' '.join([','.join(list(i)) for i in tempout if tempout.count(i) > 1]))

def Prepare_Check_Read_indfile():

		global gl_indinfo, gl_sex

		gl_indinfo = {}
		print('Reading indfile...')
		f = open(gl_indfile, 'r')
		A = [i.split() for i in f.readlines()]
		f.close()
		if not set([i[2] for i in A]).issubset(set(['0','1','2'])):
			print('Error: Invalid sex information: '+' '.join(list(set([i[2] for i in A]).difference(set(['0','1','2'])))))
			exit()
		if '0' in [i[2] for i in A]:
			gl_sex = 0
		else:
			gl_sex = 1
		for j in set([i[1] for i in A]):
			tempid = [k[0] for k in A if k[1] == j]
			if len(set(tempid)) < len(tempid):
				print('Error: Different individuals from population '+j+' take same id: '+' '.join([k for k in tempid if tempid.count(k) > 1]))
				exit()
		for j in set([i[1] for i in A]):
			gl_indinfo[j] = [i + [A.index(i)] for i in A if i[1] == j]
		print(str(len(A))+' individuals found in the indfile')
		return

def Check_snpfile_value(x):
	Check_value_int(x[0], 'physical position')
	if int(x[0]) <= 0:
		print('Error: Invalid physical position: '+x[0])
		exit()
	if gl_recrate == -9 and gl_rec == '1':
		Check_value_float(x[1], 'genetic distance')
		if float(x[1]) < 0:
			print('Error: Invalid genetic distance: '+x[1])
			exit()
	if gl_mutrate == -9 and gl_mut != '0':
		Check_value_float(x[2], 'mutation rate')
		if float(x[2]) < 0:
			print('Error: Invalid mutation rate: '+x[2])
			exit()

def Prepare_Check_Read_snpfile():
 
		global gl_pos, gl_gen_dis, gl_mut_dis

		gl_pos = []
		gl_gen_dis = []
		gl_mut_dis = []

		print('Reading snvfile...')
		f = open(gl_snpfile, 'r')
		a = f.readline().split()

		Check_snpfile_value(a)
		gl_pos.append(int(a[0]))

		gl_mut_dis.append(0)    
		if gl_recrate == -9:
		        gl_gen_dis.append(float(a[1]))
		else:
		        gl_gen_dis.append(0)    
		while True:
			a = f.readline().split()        
			if not a:
				break
			Check_snpfile_value(a)
			if int(a[0]) <= gl_pos[-1]:
				print('Error: Invalid physical positions: ' + str(gl_pos[-1])+' '+str(a[0]))
				exit()
			gl_pos.append(int(a[0]))

			if gl_recrate == -9:
				if float(a[1]) < gl_gen_dis[-1]:
					print('Error: Invalid genetic distances: ' + str(gl_gen_dis[-1])+' '+str(a[1]))
					exit()
				gl_gen_dis.append(float(a[1]))
			else:
				gl_gen_dis.append(gl_gen_dis[-1]+gl_recrate*(gl_pos[-1]-gl_pos[-2]))    
			if gl_mutrate == -9:
			    gl_mut_dis.append(gl_mut_dis[-1]+float(a[2])*(gl_pos[-1]-gl_pos[-2]))  
			else:
			    gl_mut_dis.append(gl_mut_dis[-1]+gl_mutrate*(gl_pos[-1]-gl_pos[-2])) 
		f.close()
		gl_pos = np.array(gl_pos, dtype = np.uint32)
		gl_gen_dis = np.array(gl_gen_dis, dtype = np.float64)
		gl_mut_dis = np.array(gl_mut_dis, dtype = np.float64)
		print(str(len(gl_pos)) + ' positions found in the snvfile')
		print('gl_mut_dis: ' +str(gl_mut_dis[-1]))
		return

def Prepare_Check_Read_modfile():

		global gl_pop, gl_pos, gl_indinfo, gl_outpop, gl_outgen, gl_outnum, gl_mutation, gl_popname

		gl_mutation = []

		print('Reading modfile...')
		f = open(gl_modfile, 'r')
		A = [i.strip() for i in f.readlines() if i.strip() != '']
		A = [i.split() for i in A if i[0] != '#']
		f.close()

		modelindex = [A.index(i) for i in A if i[0][0]=='*']+[len(A)]
		model = [A[modelindex[i]:modelindex[i+1]] for i in range(len(modelindex)-1)]
		print(str(len(model)) + ' modules found in the modfile')
		for i in model:
			if len(list(set([len(j) for j in i]))) >1 :
				print('Error: Exist lines with different columns in model '+i[0][0])
				exit()
		gl_popname = list(set([i for j in model for i in j[0][1].split(',')]))
		gl_pop = {x:POPULATION(x) for x in gl_popname}
		if not set(gl_indinfo.keys()).issubset(gl_pop.keys()):
			print('Warning: Input populations are not in model: '+' '.join([i for i in gl_indinfo.keys() if i not in gl_pop.keys()]))

		for i in [x for i in model for x in i[0][0][1:].split('-')]:
			Check_value_int(i, 'generation in model')
		for i in model:
			x = i[0][0][1:].split('-')
			if len(list(range(int(x[0]), int(x[1])+1))) != len(i) - 1:
				print('Error: different length given between generations and rows in model: '+i[0][0][1:])
				exit()
			for j in range(int(x[0]), int(x[1])+1):
				for k in i[0][1].split(','):
					if j not in gl_pop[k].generation.keys():
						gl_pop[k].set_generation(j, GENERATION(k, j))
		for i in gl_pop.keys():
			if i in gl_indinfo.keys():
				gl_pop[i].Founder()
		
		for i in model:
			temppop = i[0][1].split(',')
			tempgen = list(range(int(i[0][0][1:].split('-')[0]), int(i[0][0][1:].split('-')[1])+1))
			newpop = [j for j in temppop if min(list(gl_pop[j].generation.keys())) in tempgen]
			if len(newpop)>=2:
				print('Error: Too many new populations given in model: '+' '.join(newpop))
				exit()
			for j in range(1, len(i)):
				tempNe = [x.split(',') for x in i[j][0].split(':')]
				if len(tempNe) > 2:
					print('Error: Too many : in Ne in models: '+i[0][1])
					exit()
				elif gl_sex == 0 and len(tempNe) == 2:
					print('Error: Input sex-specific Ne in the modfile while no-sex individuals exsit in the indfile')
					exit()
				for m in [y for x in tempNe for y in x]:
					Check_value_int(m, 'Ne in models')
					if int(m) < 0:
						print('Error: Invalid Ne in models: '+m)
						exit()
				for k,m in enumerate(tempNe):
					if len(m) == 1:
							tempNe[k] = [m[0] for x in temppop]
					elif len(m) != len(temppop):
							print('Error: Different length given between populations and Ne: ' + i[0][1]+' '+i[j][0])
							exit()
				if len(tempNe) == 1:
					tempNe = list(zip(tempNe[0]))
				else:
					tempNe = list(zip(tempNe[0], tempNe[1]))
				for k,m in enumerate(temppop):
						gl_pop[m].generation[tempgen[j-1]].set_Ne(list(map(int,list(tempNe[k]))))

				if len(newpop) == 0:
					for k in temppop:
						gl_pop[k].generation[tempgen[j-1]].set_Proportion({k:[1,1]})
				else:
					for k in temppop:
						if k in newpop:
							continue
						gl_pop[k].generation[tempgen[j-1]].set_Proportion({k:[1,1]})
					tempprop = [k.split(',') for k in i[j][1].split(':')]
					if len(tempprop) > 2:
						print('Error: Too many : in proportion in models: '+i[0][1])
						exit()
					elif gl_sex == 0 and len(tempprop) == 2:
						print('Error: Input sex-specific proportion in the modfile while no-sex individuals exsit in the indfile') 
						exit()
					if len(tempprop[0]) != len(temppop) or len(tempprop[-1]) != len(temppop):
						print('Error: Different length given between populations and proportions: '+ i[0][1]+' '+i[j][1])
						exit()
					for k in tempprop[0] + tempprop [-1]:
						Check_value_float(k, 'proportion in model')
						if float(k) < 0:
							print('Error: Invalid proportion in model: '+k)
							exit()
					tempprop = list(zip(tempprop[0],tempprop[-1]))
					if j==1 and (float(tempprop[temppop.index(newpop[0])][0])!=0 or float(tempprop[temppop.index(newpop[0])][-1]) != 0):
						print('Error: Invalid proportion in founder generation of a new population: '+ ' '.join(tempprop[temppop.index(newpop[0])]))
					gl_pop[newpop[0]].generation[tempgen[j-1]].set_Proportion({m:list(map(float, tempprop[k])) for k,m in enumerate(temppop)})
			if gl_sel == '1':
				tempsel = []
				for j in i[0][2:]:
					tempj = j.split(':')
					jpos = assemble_pos(tempj[0])
					jcon = tempj[1].split(',')
					if set(jcon).issubset(set(gl_pop.keys())):
						jcon = [[x for k in jpos] for x in jcon]
						jlab = 1
					elif set(list(''.join(jcon))).issubset(set(['A','T','G','C'])) or set(list(''.join(jcon))).issubset(set(['0','1'])):
						jcon = [list(x) for x in jcon]
						jlab = 0
					else:
						print('Error: Invalid selection condition: '+tempj[1])
						exit()
					if len(tempj) == 2:
						jmod = '1'
					elif len(tempj) == 3:
						if tempj[2] not in ['1','2']:
							print('Error: Invalid selection mode: '+tempj[2])
							exit()
						jmod = tempj[2]
					tempsel.append([jpos, jcon, jlab, jmod])

				for j in range(1, len(i)):
					for k, m in enumerate(i[j][2:]):
						tempcoe = m.split(':')
						if len(tempcoe) > 2:
							print('Error: Too many : in selection coefficent in models: '+i[0][1])
							exit()
						elif gl_sex == 0 and len(tempcoe) == 2:
							print('Error: Input sex-specific selection coefficient in the modfile while no-sex individuals exsit in the indfile')
							exit()
						tempcoe = list(zip(tempcoe[0].split(','), tempcoe[-1].split(',')))
						if len(tempcoe) == 1:
							for n in temppop:
								Check_value_float(tempcoe[0][0], 'selection coefficient in model')
								Check_value_float(tempcoe[0][1], 'selection coefficient in model')
								jcoe = tempsel[k][:3] + [list(map(float,tempcoe[0]))] + tempsel[k][3:]
								gl_pop[n].generation[tempgen[j-1]].add_Sel(jcoe)
						elif len(tempcoe) == len(temppop): 
							for n,p in enumerate(tempcoe):
								Check_value_float(p[0], 'selection coefficient in model')
								Check_value_float(p[1], 'selection coefficient in model')
								jcoe = tempsel[k][:3] + [list(map(float,p))] + tempsel[k][3:]
								gl_pop[temppop[n]].generation[tempgen[j-1]].add_Sel(jcoe)
						else:
							print('Error: Invalid number of selection coeffients in model')
							exit()
		for i in gl_pop.keys():
			existgen = list(gl_pop[i].generation.keys())
			if i in gl_outpop:
				x = [k for j,k in enumerate(gl_outgen) if gl_outpop[j] == i]
				if min(x) < min(existgen):
					print('Error: The output generation does not exist in models: '+str(min(x)))
					exit()
				maxgen = max(x + existgen)
			else:
				maxgen = max(existgen)
			for j in range(min(existgen), maxgen+1):
				if j not in existgen:
					gl_pop[i].set_generation(j, GENERATION(i, j))
					gl_pop[i].generation[j].set_Ne(gl_pop[i].generation[j-1].Ne)
					gl_pop[i].generation[j].set_Proportion(gl_pop[i].generation[j-1].Proportion)
					for k in gl_pop[i].generation[j-1].Sel:
						gl_pop[i].generation[j].add_Sel(k)

		start = min([min(list(i.generation.keys())) for i in gl_pop.values()])
		end = max([max(list(i.generation.keys())) for i in gl_pop.values()])

		for i in range(start, end+1):
			for j in gl_pop.keys():
				if i in gl_pop[j].generation.keys():
					for k in gl_pop[j].generation[i].Proportion.keys():
						if float(sum(gl_pop[j].generation[i].Proportion[k]))!=0:
							if i not in gl_pop[k].generation.keys():
								gl_pop[k].set_generation(i, GENERATION(k, i))
								gl_pop[k].generation[i].set_Ne(gl_pop[k].generation[i-1].Ne)
								gl_pop[k].generation[i].set_Proportion(gl_pop[k].generation[i-1].Proportion)
								for m in gl_pop[k].generation[i-1].Sel:
									gl_pop[k].generation[i].add_Sel(m)

		if gl_outpop == []:
			tempsort = [(k, min(list(gl_pop[k].generation.keys()))) for k in gl_pop.keys()]
			tempsort.sort(key = lambda x: x[1])
			gl_outpop.append(tempsort[-1][0])
			gl_outgen.append(max(list(gl_pop[gl_outpop[0]].generation.keys())))
			gl_outnum = [10]
		return

def Check_file_row(x):
	k = 0 
	f = open(x,'r')
	while True:
		a = f.readline().strip()
		if not a:
			break
		k += 1
	f.close()
	return k

def Check_character(x):
	tempchar = []
	k = 0
	for i in x:
		k += 1
		if i not in tempchar:
			tempchar.append(i)
	return k, set(tempchar)

def Prepare_Check_Read_hapfile():
		
		global gl_pos, gl_gen_dis, gl_mut_dis, gl_mutation, gl_type

		if Check_file_row(gl_hapfile) < 2*Check_file_row(gl_indfile):
			print('Error: Less rows in the hapfile than expectation according to the indfile')
			exit()
		elif Check_file_row(gl_hapfile) > 2*Check_file_row(gl_indfile):
			print('Error: More rows in the hapfile than expectation according to the indfile')
			exit()

		tempind = [x for i in gl_indinfo.values() for x in i]
		tempind.sort(key = lambda x: x[-1])

		tempselpos = list(set([z for x in gl_pop.values() for y in x.generation.values() for z in [u for w in y.Sel for u in w[0]]]))
		print('Reading hapfile...')
		gl_type = 0
		k = 0
		f = open(gl_hapfile, 'r')
		while True:
			a = f.readline().strip()
			if not a:
			        break
			length, charset = Check_character(a)
			if length != len(gl_pos):
				print('Error: Inconsistence length between postions in the mapfile and line '+str(k+1)+' in the hapfile')
				exit()
			x = tempind[int(k/2)]
			y = min(list(gl_pop[x[1]].generation.keys()))
			if gl_type == 0:
				if charset.issubset(set(['A','G','T','C'])):
					gl_type = ['A','G','T','C']
				elif charset.issubset(set(['0','1'])):
					gl_type = ['0','1']
				else:
					print('Error: Invalid characters in line '+str(k+1)+' in the hapfile')
					exit()
				for i in tempselpos:
					gl_pop[x[1]].generation[y].Ind[x[0]].hap[k%2].seg[0].add_sel(i, a[bisect.bisect_left(gl_pos, i)])
			elif x[2] == '1' and gl_chr == 'X' and k%2 == 1:
				if not charset.issubset(set(['9'])):
					print('Error: Invalid characters in line '+str(k+1)+' in the hapfile')
					exit()
			else:
				if not charset.issubset(gl_type):
					print('Error: Invalid characters in line '+str(k+1)+' in the hapfile')
					exit()
				for i in tempselpos:
					gl_pop[x[1]].generation[y].Ind[x[0]].hap[k%2].seg[0].add_sel(i, a[bisect.bisect_left(gl_pos, i)])
			k += 1
		f.close()
		print('-'.join(gl_type)+' type found in the hapfile')
		return

def update_gl_mutation(x):
	tempmut = []
	for i in gl_pop.values():
		if x in i.generation.keys():
			tempmut += i.generation[x].summary_mut()
	return list(set(tempmut))


def TimeFlying():
	global gl_mutation, gl_mut, gl_sel, gl_out
	start = min([min(list(i.generation.keys())) for i in gl_pop.values()])
	end = max([max(list(i.generation.keys())) for i in gl_pop.values()])
	gl_out = {}
	print('Time is flying!')
	for i in range(start, end+1):
		for j in gl_pop.values():
			if i not in j.generation.keys():
				continue
			if j.generation[i].Ne == -9:
				if (j.popname,i) in list(zip(gl_outpop,gl_outgen)):
					print('record',(j.popname,i))
					gl_out[(j.popname,i)] = copy.deepcopy(j.generation[i])
				if gl_mut == '1':
					j.generation[i].add_mut()
				continue
			j.generation[i].Reproduce()
			if gl_sel == '1':
				j.generation[i].update_coe()
			if (j.popname,i) in list(zip(gl_outpop,gl_outgen)):
				print('record',(j.popname,i))
				gl_out[(j.popname,i)] = copy.deepcopy(j.generation[i])
			if gl_mut == '1':
				j.generation[i].add_mut()
			print('generation', i, j.popname, 'finish')
		if gl_mut == '1':
			gl_mutation = update_gl_mutation(i)
		for j in gl_pop.keys():
			if i-2 in gl_pop[j].generation.keys():
				gl_pop[j].generation[i-2].Delete()
		gc.collect()

def Trace_mut(x,y):
	x = gl_pop[x[0]].generation[x[1]].Ind[x[2]]
	if y not in [k for i in x.hap for j in i.seg for k in j.mut]:
		return []
	else:
		if x.parents == (0,0):
			return [x.index]
		a = Trace_mut(x.parents[0], y)
		b = Trace_mut(x.parents[1], y)
		if a == [] and b == []:
			return [x.index]
		elif a == []:
			return b
		elif b == []:
			return a
		else:
			return a+b

def calculate_gd(x):
	y = bisect.bisect_left(gl_pos,x)
	if x == gl_pos[y]:
		return gl_gen_dis[y]
	else:
		return ((x-gl_pos[y-1])*gl_gen_dis[y] + (gl_pos[y]-x)*gl_gen_dis[y-1])/(gl_pos[y]-gl_pos[y-1])

def Output():
	tempoutind = {}
	temppos = []
	tempind = [x for i in gl_indinfo.values() for x in i]
	tempind.sort(key = lambda x: x[-1])
	tempmatchind = [tuple(x[:2]) for x in tempind]
	for i in list(zip(gl_outpop,gl_outgen,gl_outnum)):
		x = list(gl_out[i[:2]].IndSampling(['0','1','2'],i[2],False))
		x.sort(key = lambda y:y.index[2])
		tempoutind[i] = x
		for j in tempoutind[i]:
			for k in j.hap:
				for m in k.seg:
					temppos += list(m.mut)
	temppos = list(set(temppos))
	temppos.sort(key = lambda x:int(x))
	tempposindex = [bisect.bisect(gl_pos, i) for i in temppos]
	haplen = len(gl_pos)+len(temppos)
	
	k = 0
	fsnp = open(gl_outprefix+'.snv','w')
	for i,j in enumerate(gl_pos):
		if k >= len(temppos):
			tempgd = str(round(gl_gen_dis[i],8))
			if i == 0:
				tempmt = '0'
			else:
				tempmt = str(round((gl_mut_dis[i] - gl_mut_dis[i-1])/(gl_pos[i]-gl_pos[i-1]),8))
			fsnp.write(str(j)+'\t'+tempgd+'\t'+tempmt+'\tF\n')
		elif i != tempposindex[k]:
			tempgd = str(round(gl_gen_dis[i],8))
			if i == 0:
				tempmt = '0'
			else:
				tempmt = str(round((gl_mut_dis[i] - gl_mut_dis[i-1])/(gl_pos[i]-gl_pos[i-1]),8))
			fsnp.write(str(j)+'\t'+tempgd+'\t'+tempmt+'\tF\n')
		else:
			while i == tempposindex[k]:
				tempgd = str(round((gl_gen_dis[i]-gl_gen_dis[i-1])*(temppos[k]-gl_pos[i-1])/(gl_pos[i]-gl_pos[i-1])+gl_gen_dis[i-1], 8))
				tempmt = str(round((gl_mut_dis[i] - gl_mut_dis[i-1])/(gl_pos[i]-gl_pos[i-1]),8))
				fsnp.write(str(temppos[k])+'\t'+tempgd+'\t'+tempmt+'\tT,'+gl_type[0]+'/'+gl_type[1]+'\n')
				k += 1
				if k == len(temppos):
					break
			tempgd = str(round(gl_gen_dis[i],8))
			tempmt = str(round((gl_mut_dis[i] - gl_mut_dis[i-1])/(gl_pos[i]-gl_pos[i-1]),8))
			fsnp.write(str(j)+'\t'+tempgd+'\t'+tempmt+'\tF\n')				
	fsnp.close()

	fsel = open(gl_outprefix+'.sel','w')
	start = min([min(list(i.generation.keys())) for i in gl_pop.values()])
	end = max([max(list(i.generation.keys())) for i in gl_pop.values()])
	for i in range(start, end+1):
		for j in list(gl_pop.keys()):
			if i not in gl_pop[j].generation.keys():
				continue
			for k,m in enumerate(gl_pop[j].generation[i].Sel):
				for p,q in enumerate(m[1]):
					fsel.write(j+'\t'+str(i)+'\t'+','.join(list(map(str, m[0])))+'\t'+','.join(q)+'\t'+'\t'.join(list(map(str, list(gl_pop[j].generation[i].Selfreq[k][p]))))+'\n')
	fsel.close()

	finhap = open(gl_hapfile, 'r')
	hapdata = [i.strip() for i in finhap.readlines()]
	finhap.close()
	fhap = open(gl_outprefix + '.hap', 'w')
	find = open(gl_outprefix + '.ind', 'w')
	fseg1 = open(gl_outprefix + '.recseg', 'w')
	fseg2 = open(gl_outprefix + '.indseg', 'w')
	fseg3 = open(gl_outprefix + '.popseg', 'w')
	t = 0
	for i in list(zip(gl_outpop,gl_outgen,gl_outnum)):
		for j in tempoutind[i]:
			find.write(str(j.index[2])+'\t'+str(j.index[0])+'_'+str(j.index[1])+'\t'+str(j.sex)+'\n')
			for k in j.hap:
				fseg1.write('\t'.join(list(map(str,list(j.index)+[j.hap.index(k)+1])))+'\n')
				fseg2.write('\t'.join(list(map(str,list(j.index)+[j.hap.index(k)+1])))+'\n')
				fseg3.write('\t'.join(list(map(str,list(j.index)+[j.hap.index(k)+1])))+'\n')
				if k.status == 1:
					temphap = ''.join(['9' for i in range(haplen)])
				else:
					fseg1.write(str(k.seg[0].left)+'\t')
					fseg2.write(str(k.seg[0].left)+'\t')
					fseg3.write(str(k.seg[0].left)+'\t'+str(round(calculate_gd(k.seg[0].left),8))+'\t')
					for m,n in enumerate(k.seg):
						if m == 0:
							continue
						fseg1.write(str(k.seg[m-1].right)+'\t'+'\t'.join(list(map(str,list(k.seg[m-1].lab))))+'\n')
						fseg1.write(str(n.left)+'\t')
						if n.lab[0] != k.seg[m-1].lab[0]:
							fseg2.write(str(k.seg[m-1].right)+'\t'+'\t'.join(list(map(str,list(k.seg[m-1].lab))))+'\n')
							fseg2.write(str(n.left)+'\t')
							fseg3.write(str(k.seg[m-1].right)+'\t'+str(round(calculate_gd(k.seg[m-1].right),8))+'\t'+str(k.seg[m-1].lab[0])+'\n')
							fseg3.write(str(n.left)+'\t'+str(round(calculate_gd(n.left),8))+'\t')
						elif n.lab[2]!= k.seg[m-1].lab[2]:
							fseg2.write(str(k.seg[m-1].right)+'\t'+'\t'.join(list(map(str,list(k.seg[m-1].lab))))+'\n')
							fseg2.write(str(n.left)+'\t')
					fseg1.write(str(k.seg[-1].right)+'\t'+'\t'.join(list(map(str,list(k.seg[-1].lab))))+'\n')
					fseg2.write(str(k.seg[-1].right)+'\t'+'\t'.join(list(map(str,list(k.seg[-1].lab))))+'\n')
					fseg3.write(str(k.seg[-1].right)+'\t'+str(round(calculate_gd(k.seg[-1].right),8))+'\t'+str(k.seg[-1].lab[0])+'\n')

					tempseg = ''
					temphapmut = []
					for m in k.seg:
						tempseg += hapdata[tempmatchind.index((m.lab[2],m.lab[0]))*2+m.lab[3]-1][bisect.bisect_left(gl_pos, m.left): bisect.bisect_right(gl_pos, m.right)]
						temphapmut += list(m.mut)
					if tempposindex == []:
						temphap = tempseg
					else:
						temphap = tempseg[:tempposindex[0]]
						k = 0
						for u,v in enumerate(temppos):
							if temphapmut == []:
								mutallele = gl_type[0]
							elif v == temphapmut[k]:
								mutallele = gl_type[1]
								k = (k+1)%len(temphapmut)
							else:
								mutallele = gl_type[0]
							if v != temppos[-1]:
								temphap += mutallele + tempseg[tempposindex[u]:tempposindex[u+1]]
							else:
								temphap += mutallele + tempseg[tempposindex[u]:]
				fhap.write(temphap+'\n')
				t += 1
	fhap.close()
	find.close()
	fseg1.close()
	fseg2.close()
	fseg3.close()

def Prepare_Print():
	global flog
	print('haplotype file: ' + gl_hapfile)
	print('genetic map file: ' + gl_snpfile)
	print('individual file: ' + gl_indfile)
	print('model file: ' + gl_modfile)
	print('length(bp): ' + str(gl_pos[-1]))
	print('input individuals (population:number): ' + ' '.join([i+':'+str(len(j)) for i,j in gl_indinfo.items()]))
	print('output individuals (population:generation:number): ' + ' '.join([i+':'+str(j)+':'+str(k) for i,j,k in zip(gl_outpop,gl_outgen,gl_outnum)]))
	print('output prefix: ' + gl_outprefix)
	print('recombination state: ' + gl_rec)
	print('mutation state ' + gl_mut)
	print('selection state ' + gl_sel)
	print('gl_recrate', gl_recrate)
	print('gl_mutrate', gl_mutrate)
	print('gl_sex',gl_sex)
	print('\n===================================================\n')
	flog = open(gl_outprefix + '.log', 'w')
	flog.write('haplotype file: ' + gl_hapfile + '\n')
	flog.write('genetic map file: ' + gl_snpfile + '\n')
	flog.write('individual file: ' + gl_indfile + '\n')
	flog.write('model file: ' + gl_modfile + '\n')
	flog.write('length(bp): ' + str(gl_pos[-1]) + '\n')
	flog.write('input individuals (population:number): ' + ' '.join([i+':'+str(len(j)) for i,j in gl_indinfo.items()]) + '\n')
	flog.write('output individuals (population:generation:number): ' + ' '.join([i+':'+str(j)+':'+str(k) for i,j,k in zip(gl_outpop,gl_outgen,gl_outnum)]) + '\n')
	flog.write('output prefix: ' + gl_outprefix + '\n')
	flog.write('recombination state: ' + gl_rec + '\n')
	flog.write('mutation state ' + gl_mut + '\n')
	flog.write('selection state ' + gl_mut + '\n')
	flog.write('gl_recrate ' + str(gl_recrate) + '\n')
	flog.write('gl_mutrate ' + str(gl_mutrate) + '\n')
#	flog.close()

def Prepare():
        Print_Header()
        opts, args = Prepare_Get_Check_Arg()
        Prepare_Set_Value(opts)

        Prepare_Check_Read_indfile()
        Prepare_Check_Read_snpfile()
        Prepare_Check_Read_modfile()
        Prepare_Check_Read_hapfile()
        Prepare_Print()

if __name__ == '__main__':
		import time
		starttime = time.time()
		import sys, getopt, bisect, copy, time, gc
		import numpy as np
		Prepare()
		TimeFlying()
		Output()
		print('Simulation Completed!')
		endtime = time.time()
		flog.write(str(starttime)+'\n')
		flog.write(str(endtime)+'\n')
		flog.close()
