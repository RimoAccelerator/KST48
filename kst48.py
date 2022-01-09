import os, sys, re, numpy, math

dE_thresh = 0.000050
rms_thresh = 0.0025
max_dis_thresh = 0.004
max_dis_thresh = 0.004
max_g_thresh = 0.0007
rms_g_thresh = 0.0005
max_steps = 100
max_step_size = 0.1
n_atom = 0
reduced_factor = 0.5

list_element = []
header_a = ''
header_b = ''
tail1 = ''
tail2 = ''
constraints = []
scans = []
geom = []

#constraints = [['R',1,2,1.0]]
#constraints = [['A',2,1,3,90],['R',1,2,1.2],['R',1,3,1.2]]
prog = 'gaussian'
prog_comm = 'g16'
method = ''
td1 = ''
td2 = ''
mp2 = False
delete_gbw = False

def buildHeader(nprocs, mem, charge, mult1, mult2, method, prog, td1='', td2=''):
	header_1, header_2 = ['', '']
	if prog == 'gaussian':
		header_1 = f'%chk=a.chk\n%nprocs={nprocs} \n%mem={mem} \n# {method} {td1} nosymm\n\n Title Card \n\n{charge} {mult1}'
		header_2 = f'%chk=b.chk\n%nprocs={nprocs} \n%mem={mem} \n# {method} {td2} nosymm\n\n Title Card \n\n{charge} {mult2}'
	elif prog == 'orca':
		header_1 = f'%pal nprocs {nprocs} end\n%maxcore {mem} \n! {method} \n\n *xyz {charge} {mult1}'
		header_2 = f'%pal nprocs {nprocs} end\n%maxcore {mem} \n! {method} \n\n *xyz {charge} {mult2}'
		if '***' in header_1:
			header_1 = header_1.replace('***','JOBS/a.gbw')
		if '***' in header_2:
			header_2 = header_2.replace('***','JOBS/b.gbw')
	return [header_1, header_2]

def buildInitJob(nprocs, mem, charge, mult1, mult2, method, prog, run_mode, td1='', td2=''):
	header_1, header_2 = ['','']
	if run_mode == 'normal' or run_mode == 'read' or run_mode == 'noread':
		header_1, header_2 = buildHeader(nprocs, mem, charge, mult1, mult2, method, prog, td1=td1, td2=td2)
	elif run_mode == 'stable'  or run_mode == 'inter_read':
		if prog == 'gaussian':
			header_1, header_2 = buildHeader(nprocs, mem, charge, mult1, mult2, method + ' stable=opt', prog)
		elif prog == 'orca':
			header_1, header_2 = buildHeader(nprocs, mem, charge, mult1, mult2, method + '\n %scf stabperform true StabRestartUHFifUnstable true end \n', prog)
		else:
			raise Exception('This program not supported for STABLE mode!')
	else:
		raise Exception('Unknown mode!')
	if '"***"' in header_1:
		header_1 = header_1.replace('"***"', f'"JOBS/a.gbw"')
	if '"***"' in header_2:
		header_1 = header_2.replace('"***"', f'"JOBS/b.gbw"')
	return [header_1, header_2]

def runPrePoint(header_1, header_2, prog, run_mode='normal'):
	#os.system('cp *.chk JOBS; cp *.gbw JOBS')
	if prog == 'gaussian':
		writeGjf(geom, header_1, tail1, f'JOBS/pre_A.gjf')
		writeGjf(geom, header_2, tail2, f'JOBS/pre_B.gjf')
		os.system(f'{prog_comm} JOBS/pre_B.gjf')
		if run_mode == 'inter_read':
			os.system(f'cp b.chk a.chk')
			os.system(f'sed -i "s/#/# guess=(read,mix)/g" JOBS/pre_A.gjf')
		os.system(f'{prog_comm} JOBS/pre_A.gjf')
		
	elif prog == 'orca':
		writeORCA(geom, header_1, tail1, f'JOBS/pre_A.inp')
		writeORCA(geom, header_2, tail2, f'JOBS/pre_B.inp')
		os.system(f'{prog_comm} JOBS/pre_B.inp > JOBS/pre_B.out')
		os.system(f'cp JOBS/pre_B.gbw JOBS/b.gbw')
		if run_mode == 'inter_read':
			os.system(f'cp JOBS/pre_B.gbw JOBS/a.gbw')
			print('Note: The inter_read mode is set for ORCA. In Gaussian, the program automatically add guess=mix for state A, \
but this will not be done for ORCA. If you want to converge to correct OSS wavefunction from a triplet wavefunction, guess=(read,mix) \
is always beneficial. So please do not forget to add relevant convergen controlling in your ORCA tail part.')
		os.system(f'{prog_comm} JOBS/pre_A.inp > JOBS/pre_A.out')
		os.system(f'cp JOBS/pre_A.gbw JOBS/a.gbw')

	#elif prog == 'xtb':
	#	writeXYZ(geom, f'JOBS/pre_A.xyz')
	#	writeXYZ(geom, f'JOBS/pre_B.xyz')
	#	os.system(f'{prog_comm} JOBS/pre.xyz --engrad')
	#	os.system(f'{prog_comm} JOBS/pre.xyz --engrad')

def inputParser(path):
	global geom
	global list_element
	global n_atom
	global prog
	global td1
	global td2
	global mp2
	global prog_comm
	global tail1
	global tail2
	global reduced_factor
	charge, mult1, mult2, method, nprocs, mem = ['', '', '', '', '', '']
	command = {'gau':'', 'orca':'', 'xtb':''}
	run_mode = 'normal'
	with open(path) as f:
		isGeom = False
		isTail1 = False
		isTail2 = False
		isConst = False
		for l in f:
			l = l.lower().strip()
			if '#' in l:
				l = l.split('#')[0]
			if '*geom' in l:
				isGeom = True
			elif l == '*' and isGeom:
				isGeom = False
			elif '*tail1' in l:
				isTail1 = True
			elif '*tail2' in l:
				isTail2 = True
			elif l == '*' and isTail2:
				isTail2 = False
			elif l == '*' and isTail1:
				isTail1 = False
			elif '*constr' in l:
				isConst = True
			elif l == '*' and isConst:
				isConst = False
			elif isGeom and (re.match("\s*\S+\s*\-*[0-9]+",l)):
				geom.extend(l.split()[1:])
				list_element.append(l.split()[0])
				n_atom += 1
			elif isTail1:
				tail1 += l
				tail1 += '\n'
			elif isTail2:
				tail2 += l
				tail2 += '\n'
			elif isConst and l != '':
				l_splitted = l.split()
				if(l_splitted[0] != 's'):
					constraints.append(l.split())
				else:
					if l_splitted[1] == 'r':
						scans.append([l_splitted[1:4], l_splitted[4:]])
					elif l_splitted[1] == 'a':
						scans.append([l_splitted[1:5], l_splitted[5:]])
			if '=' in l and not isTail1 and not isTail2:
				parameter = l.split('=')[1].strip()
				if 'mem' in l:
					mem = parameter
				elif 'nprocs' in l:
					nprocs = parameter
				elif 'charge' in l:
					charge = parameter
				elif 'method' in l:
					method = parameter
				elif 'mult1' in l:
					mult1 = parameter
				elif 'mult2' in l:
					mult2 = parameter
				elif 'dE_thresh' in l:
					dE_thresh = float(parameter)
				elif 'rms_thresh' in l:
					rms_thresh = float(parameter)
				elif 'max_dis_thresh' in l:
					max_dis_thresh = float(parameter)
				elif 'max_g_thresh' in l:
					max_g_thresh = float(parameter)
				elif 'rms_g_thresh' in l:
					rms_g_thresh = float(parameter)
				elif 'max_steps' in l:
					max_steps = float(parameter)
				elif 'max_step_size' in l:
					max_step_size = float(parameter)
				elif 'program' in l:
					prog = parameter
				elif 'gau_comm' in l:
					command['gaussian'] = parameter
				elif 'orca_comm' in l:
					command['orca'] = parameter
				elif 'xtb_comm' in l:
					command['xtb'] = parameter
				elif 'mode' in l:
					run_mode = parameter
				elif 'td1' in l:
					td1 = parameter
				elif 'td2' in l:
					td2 = parameter
				elif 'mp2' in l and parameter == 'true':
					mp2 = True
				elif 'reduced_factor' in l:
					reduced_factor = float(parameter)
	geom = numpy.mat(geom)
	prog_comm = command[prog]
	return [nprocs, mem, charge, mult1, mult2, method, run_mode]

def modifyMethod(prog,method,run_mode):
	if prog == 'gaussian':
		method += ' force '
	elif prog == 'orca':
		method += ' engrad '
	# if noread mode is set, do not read the initial wavefunctions
	if run_mode != 'noread':
		if prog == 'gaussian':
			method += ' guess=read '
		elif prog == 'orca':
			method += '\n!moread \n %moinp "***"\n'
	return method

def readForceAndGeomForGaussian(path):
	forceArr = []
	geomArr = []
	with open(path) as f:
		isForce = False
		isGeom = False
		E = 0
		archive_part = ''
		isArchive = False
		for l in f.readlines():
			if 'Input orientation' in l:
				isGeom = True
			elif 'Distance matrix' in l:
				isGeom = False
			elif 'Forces (Hartrees/Bohr)' in l:
				isForce = True
			elif 'Cartesian Forces:  Max' in l:
				isForce = False
			elif 'SCF DONE' in l.upper():
				E = float(l.split('=')[1].upper().split('A.U.')[0])
			elif 'E(TD-HF/TD-DFT)' in l.upper():
				E = float(l.split('=')[1])
			elif mp2 and 'MP2=' in l.upper(): #read the MP2 energy in the archive part of a log file
				archive_part += l.upper().strip()
				isArchive = True
			elif isArchive:
				archive_part += l.upper().strip() # sometimes the MP2 energy may be separated by a \n. Two lines have to be combined
				isArchive = False
				E = float(rchive_part.split('MP2=')[1].split('=')[0])
			elif isForce and (re.match("\s*[0-9]+\s*[0-9]+\s*\-*[0-9]+",l) != None):
				forceArr.extend(l.split()[2:])
			elif isGeom and (re.match("\s*[0-9]+\s*[0-9]+\s*[0-9]+\s*\-*[0-9]+",l) != None):
				geomArr.extend(l.split()[3:])
	geomArr = [float(i) for i in geomArr]
	#forceArr = [float(i)/0.529 for i in forceArr]
	forceArr = [float(i) for i in forceArr]
	forceArr = addConst(geomArr, forceArr, constraints)
	return [geomArr, forceArr, E]

def readForceAndGeomForORCA(path):
	forceArr = []
	geomArr = []
	path_splitted = path.split('.')
	path = ''.join(path_splitted[:-1])+'.engrad' #replace A.xxx into A.engrad
	with open(path) as f:
		isForce = False
		isGeom = False
		isEnergy = False
		E = 0
		for l in f.readlines():
			if 'The atomic numbers and current coordinates in Bohr' in l:
				isGeom = True
			elif '#' in l and len(geomArr)>3:
				isGeom = False
			elif 'The current gradient' in l:
				isForce = True
			elif '#' in l and len(forceArr)>3:
				isForce = False
			elif 'current total energy in Eh' in l:
				isEnergy = True
			elif '#' in l and E != 0:
				isEnergy = False
			elif isEnergy and (re.match("\s*\-*[0-9]+",l) != None):
				E = float(l.strip())
			elif isForce and (re.match("\s*\-*[0-9]+",l) != None):
				forceArr.append(l.strip())
			elif isGeom and (re.match("\s*[0-9]+\s*\-*[0-9]+",l) != None):
				geomArr.extend(l.split()[1:])
	geomArr = [float(i)*0.52918 for i in geomArr] #change Bohr to angstrom
	forceArr = [-float(i) for i in forceArr] #ORCA outputs gradients. Here it is adapted to gaussian, which is the force
	forceArr = addConst(geomArr, forceArr, constraints)
	return [geomArr, forceArr, E]

def readForceAndGeom(path):
	if prog == 'gaussian':
		return readForceAndGeomForGaussian(path)
	elif prog == 'orca':
		return readForceAndGeomForORCA(path)
	else:
		raise Exception('Unsupported program!')

# read constraints and add it to the gradients according to a harmonic term
def addConst(geomArr, forceArr, constrList): #constr list:[['r', 1 , 2, 1.5], ['a', 1, 2, 3]], starting from 1
	constrR = []
	constrA = []
	constrD = []
	newForce = forceArr[:]
	applyForce = 10
	for i in constrList:
		if i[-1] == -1: # skip invalid constrain. This feature is used in the PES-scan part
			continue
		if i[0] == 'r':
			constrR.append([ int(i[1]) - 1, int(i[2]) - 1, float(i[3]) ])
		if i[0] == 'a':
			constrA.append([ int(i[1]) - 1, int(i[2]) - 1, int(i[3]) - 1, math.cos(float(i[4])/180.0*math.pi)])
	for i in constrR: # [atom_a, atom_b, dist]
		x1, y1, z1= geomArr[i[0] * 3 : i[0] * 3 + 3]
		x2, y2, z2= geomArr[i[1] * 3 : i[1] * 3 + 3]
		r12 = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
		newForce[i[0] * 3 + 0]  += applyForce * (x1-x2) * (1 - i[2] / r12)
		newForce[i[0] * 3 + 1]  += applyForce * (y1-y2) * (1 - i[2] / r12)
		newForce[i[0] * 3 + 2]  += applyForce * (z1-z2) * (1 - i[2] / r12)
		newForce[i[1] * 3 + 0]  -= applyForce * (x1-x2) * (1 - i[2] / r12)
		newForce[i[1] * 3 + 1]  -= applyForce * (y1-y2) * (1 - i[2] / r12)
		newForce[i[1] * 3 + 2]  -= applyForce * (z1-z2) * (1 - i[2] / r12)
	for i in constrA:
		# constraint for angle ABC: lambda * ((BA dot BC)/|BA||BC| - cosAngle)^2
		# its deriv: lambda * (cos_current - target) * d(cos)
		xA, yA, zA= geomArr[i[0] * 3 : i[0] * 3 + 3]
		xB, yB, zB= geomArr[i[1] * 3 : i[1] * 3 + 3]
		xC, yC, zC= geomArr[i[2] * 3 : i[2] * 3 + 3]
		rBA = math.sqrt((xA - xB)**2 + (yA - yB)**2 + (zA - zB)**2)
		rBC = math.sqrt((xC - xB)**2 + (yC - yB)**2 + (zC - zB)**2)
		BAdotBC = (xA - xB)*(xC - xB) + (yA - yB)*(yC - yB) + (zA - zB)*(zC - zB)
		cosAngle_current = BAdotBC / rBA / rBC
		cosTarget = i[3]
		dCos_dxA = (rBA*(xC-xB) - BAdotBC*(xA-xB)) / (rBC * rBA**2)
		dCos_dxC = (rBC*(xA-xB) - BAdotBC*(xC-xB)) / (rBA * rBC**2)
		dCos_dyA = (rBA*(yC-yB) - BAdotBC*(yA-yB)) / (rBC * rBA**2)
		dCos_dyC = (rBC*(yA-yB) - BAdotBC*(yC-yB)) / (rBA * rBC**2)
		dCos_dzA = (rBA*(zC-zB) - BAdotBC*(zA-zB)) / (rBC * rBA**2)
		dCos_dzC = (rBC*(zA-zB) - BAdotBC*(zC-zB)) / (rBA * rBC**2)
		dCos_dxB = (rBA*rBC*(2*xB-xA-xC) + BAdotBC*xB*(rBA**2 + rBC**2)/rBA/rBC) / (rBC**2 * rBA**2)
		dCos_dyB = (rBA*rBC*(2*yB-yA-yC) + BAdotBC*yB*(rBA**2 + rBC**2)/rBA/rBC) / (rBC**2 * rBA**2)
		dCos_dzB = (rBA*rBC*(2*zB-zA-zC) + BAdotBC*zB*(rBA**2 + rBC**2)/rBA/rBC) / (rBC**2 * rBA**2)
		newForce[i[0] * 3 + 0]  += applyForce * (cosAngle_current - cosTarget) * dCos_dxA
		newForce[i[0] * 3 + 1]  += applyForce * (cosAngle_current - cosTarget) * dCos_dyA
		newForce[i[0] * 3 + 2]  += applyForce * (cosAngle_current - cosTarget) * dCos_dzA
		newForce[i[1] * 3 + 0]  += applyForce * (cosAngle_current - cosTarget) * dCos_dxB
		newForce[i[1] * 3 + 1]  += applyForce * (cosAngle_current - cosTarget) * dCos_dyB
		newForce[i[1] * 3 + 2]  += applyForce * (cosAngle_current - cosTarget) * dCos_dzB
		newForce[i[2] * 3 + 0]  += applyForce * (cosAngle_current - cosTarget) * dCos_dxC
		newForce[i[2] * 3 + 1]  += applyForce * (cosAngle_current - cosTarget) * dCos_dyC
		newForce[i[2] * 3 + 2]  += applyForce * (cosAngle_current - cosTarget) * dCos_dzC
	return newForce

def writeGjf(geom, header, tail, name):
	f = open(name,'w+')
	f.write(header)
	f.write('\n')
	for i in range(n_atom):
		f.write('{ele}  {x}  {y}  {z}'.format(ele = list_element[i], x = geom[0, i*3], y = geom[0, i*3+1], z = geom[0, i*3+2]))
		f.write('\n')
	f.write('\n')
	f.write(tail)
	f.write('\n')
	f.close()

def writeORCA(geom, header, tail, name):
	f = open(name,'w+')
	f.write(header)
	f.write('\n')
	for i in range(n_atom):
		f.write('{ele}  {x}  {y}  {z}'.format(ele = list_element[i], x = geom[0, i*3], y = geom[0, i*3+1], z = geom[0, i*3+2]))
		f.write('\n')
	f.write('*\n')
	f.write(tail)
	f.write('\n')
	f.close()

def writeXYZ(geom, name):
	f = open(name,'w+')
	f.write(n_atom)
	f.write('\n\n')
	for i in range(n_atom):
		f.write('{ele}  {x}  {y}  {z}'.format(ele = list_element[i], x = geom[0, i*3], y = geom[0, i*3+1], z = geom[0, i*3+2]))
		f.write('\n')
	f.write('\n')
	f.close()

def runEachStep(geom, nstep, header_a, header_b, tail1, tail2):
	if prog == 'gaussian':
		writeGjf(geom, header_a, tail1, f'JOBS/{nstep}_A.gjf')
		writeGjf(geom, header_b, tail2, f'JOBS/{nstep}_B.gjf')
		os.system(f'{prog_comm} JOBS/{nstep}_B.gjf')
		os.system(f'{prog_comm} JOBS/{nstep}_A.gjf')
	elif prog == 'orca':
		writeORCA(geom, header_a, tail1, f'JOBS/{nstep}_A.inp')
		writeORCA(geom, header_b, tail2, f'JOBS/{nstep}_B.inp')
		os.system(f'{prog_comm} JOBS/{nstep}_B.inp > JOBS/{nstep}_B.log')
		os.system(f'{prog_comm} JOBS/{nstep}_A.inp > JOBS/{nstep}_A.log')
		if delete_gbw: 
			os.system(f'rm -rf JOBS/{nstep}_B.gbw')
			os.system(f'rm -rf JOBS/{nstep}_A.gbw')
		else:
			os.system(f'mv JOBS/{nstep}_B.gbw  JOBS/b.gbw')
			os.system(f'mv JOBS/{nstep}_A.gbw  JOBS/a.gbw')
	elif prog == 'xtb':
		writeXYZ(geom, header_a, tail1, f'JOBS/{nstep}_A.xyz')
		writeXYZ(geom, header_b, tail2, f'JOBS/{nstep}_B.xyz')
		os.system(f'{prog_comm} JOBS/{nstep}_B.xyz')
		os.system(f'{prog_comm} JOBS/{nstep}_A.xyz')

def getG(nstep):
	_, f1, E1 = readForceAndGeom(f'JOBS/{nstep}_A.log')
	_, f2, E2 = readForceAndGeom(f'JOBS/{nstep}_B.log')
	f1 = -numpy.array(f1)
	f2 = -numpy.array(f2)
	xVec = (f1-f2) # Please check this sign here
	xVecNorm = xVec / numpy.linalg.norm(xVec)
	fVec = (E1-E2) * xVecNorm # In Chem Phys Lett (1994) 223: 269 xVecNorm is used, although in Harvey's work xVec is used herein
	#fVec = (E1-E2)*xVec #if use this, f will be rather small, resulting in inefficient optimization with repsect to dE
	gVec = ( f1 - numpy.dot(xVecNorm, f1)*xVecNorm ) #*0.001 may help BGFS
	return fVec + gVec, E1, E2 #return an python array
	#return f1, E1, E2 # for debugging: to output force for the 1st state

def BFGSpropagation(X0, Gk, Bk):
	#rho = 0.01 #set rho = 0.01 for minimum opt, and 15 for MECP opt
	rho = 15
	dk = -numpy.linalg.solve(Bk, Gk)
	if numpy.linalg.norm(dk) > 0.1:
		dk = dk * 0.1 / numpy.linalg.norm(dk) 
	XNew = X0 + rho * dk  
	return XNew

def GDIISpropagation(Xs, Gs, Hs): #Produce a new geometry based on the GDIIS algorith, see https://manual.q-chem.com/5.3/A1.S7.html
	dimension = len(Xs)
	if len(Hs) != len(Xs):
		raise Exception('Runtime exception: H and X numbers are different.')
	EMat = numpy.mat(numpy.zeros(shape=(dimension, n_atom*3, )))
	H_mean = Hs[0]
	for i in range(1,dimension):
		H_mean += Hs[i]
	H_mean /= dimension
	for i in range(dimension):
		#print(Hs[i])
		#EMat[i] = (Hs[i].I * numpy.mat(Gs[i]).T).flatten() # different updation model of Hessian matrix.
		EMat[i] = (H_mean.I * numpy.mat(Gs[i]).T).flatten()
		#EMat[i] = (Gs[i].T).flatten() # Now use a unity matrix instead of the approximated Hessian
	onesBlock = numpy.mat(numpy.ones(dimension))
	BMat = numpy.block([[EMat * EMat.T, onesBlock.T],\
						[onesBlock, numpy.mat([0])] ]) #construct the Bmat: B 1|1 0, Bij = ei*ej
	y = numpy.append(numpy.zeros(dimension),1)
	c = numpy.linalg.solve(BMat, y)
	c = numpy.delete(c, -1) #delete the last element, lambda, in c vector
	XNew_prime = numpy.mat(numpy.zeros(n_atom*3))
	GNew_prime = numpy.mat(numpy.zeros(n_atom*3))
	HNew_prime = numpy.mat(numpy.zeros(shape=(n_atom*3,n_atom*3)))
	for i in range(dimension):
		XNew_prime += Xs[i]*c[i]
		GNew_prime += Gs[i]*c[i]
		HNew_prime += Hs[i]*c[i]
	XNew = XNew_prime - GNew_prime
	#XNew = XNew_prime - (H_mean.I * numpy.mat(GNew_prime).T).flatten()
	XNew = XNew_prime - (HNew_prime.I * numpy.mat(GNew_prime).T).flatten()
	return XNew

def GEDIISpropagation(Xs, Gs, Hs, Es): #Produce a new geometry based on the GEDIIS algorith (J. Chem. Theory Comput. 2006, 2, 835-839)
# Note that the energy to be minimized should not be E1 or E2, but produced from the pritimive function of G (not implemented yet)
	dimension = len(Xs)
	if len(Hs) != len(Xs):
		raise Exception('Runtime exception: H and X numbers are different.')
	EMat = numpy.mat(numpy.zeros(shape=(dimension, dimension)))
	for i in range(dimension):
		EMat[i, i] = 0
		for j in range(i + 1, dimension):
			EMat[i, j] = -1 * numpy.dot(Gs[i] - Gs[j], numpy.array((Xs[i] - Xs[j])).flatten()) # Gs is force rather than gradient
			EMat[j, i] = EMat[i, j] 
	onesBlock = numpy.mat(numpy.ones(dimension))
	BMat = numpy.block([[EMat, onesBlock.T],\
						[onesBlock, numpy.mat([0])] ])
	y = numpy.append(-1 *numpy.array(Es), 1)
	c = numpy.linalg.solve(BMat, y)
	c = numpy.delete(c, -1) #delete the last element, lambda, in c vector
	XNew_prime = numpy.mat(numpy.zeros(n_atom * 3))
	GNew_prime = numpy.mat(numpy.zeros(n_atom * 3))
	for i in range(dimension):
		XNew_prime += Xs[i] * c[i]
		GNew_prime += Gs[i] * c[i]
	XNew = XNew_prime + GNew_prime
	return XNew

def MaxStep(X0, XNew, factor = 1):
	dX = XNew - X0
	dX *= factor
	stepsize_norm = numpy.linalg.norm(dX)
	if  stepsize_norm > max_step_size:
		print(f'current stepsize: {stepsize_norm} is reduced to max_size')
		dX = dX * max_step_size / numpy.linalg.norm(dX)
	return X0 + dX

def HessianUpdator(Bk, yk, sk):
	#10.1002/wcms.34
	Bk_BFGS = Bk + float(1.0/(sk*yk.T)) * ( float(1+(yk*Bk*yk.T)/(sk*yk.T))*(sk.T*sk)- (sk.T*yk*Bk+Bk*yk.T*sk))
	Bk_SR1 = Bk + (yk.T - Bk*sk.T)*(yk.T-Bk*sk.T).T/float((yk.T-Bk*sk.T).T*sk.T)# slow convergence
	Bk_PSB = Bk + ((yk.T - Bk*sk.T)*sk + sk.T*(yk.T-Bk*sk.T).T)/(sk*sk.T) \
	- (float(sk*(yk.T-Bk*sk.T))*sk.T*sk)/float(sk*sk.T)**2
	
	phi = float((yk.T - Bk*sk.T).T*sk.T)**2 / \
	float((yk.T-Bk*sk.T).T*(yk.T-Bk*sk.T)) / float(sk*sk.T)
	Bk_Bofill = phi*Bk_SR1 + (1-phi)*Bk_PSB # slow convergence
	Bk = Bk_PSB
	eigval, eigvec = numpy.linalg.eig(Bk)
	eigval = numpy.real(eigval)
	print(eigval)
	numimag = 0
	#Bk_BFGS is always positive finite, while Bk_PBS is not.
	for i in eigval:
		if i<0:
			print(i)
			numimag += 1
	if numimag > 0:
		print('Negative eigenvals found for Bk_PSB, Bk_BFGS used instead')
		Bk = Bk_BFGS
	return Bk

def BFGS(X0, G0, B0, nstep):  
	print(f'\nNow Entering BFGS Step {nstep}')
	XNew = BFGSpropagation(X0, G0, B0)
	XNew = MaxStep(X0, XNew)
	runEachStep(XNew, nstep + 1, header_a, header_b, tail1, tail2)
	sk = XNew - X0 
	GNew, E1, E2 = getG(nstep + 1)
	yk = numpy.mat(GNew - G0)
	Bk = HessianUpdator(B0, yk, sk)
	#Bk_BFGS = Bk - (Bk * sk.T * sk * Bk) / float(sk * Bk * sk.T) + (yk.T * yk) / float(yk * sk.T)
	return [XNew, GNew, Bk, E1, E2]

def GDIIS(Xs, Gs, Bs, nstep, Es, flag = 'gdiis'):  
	print(f'\nNow Entering GDIIS Step {nstep}')
	XNew = GDIISpropagation(Xs, Gs, Bs)
	if flag == 'gediis':
		XNew2 = GEDIISpropagation(Xs, Gs, Bs, Es)
		XNew = XNew * 0.5 + XNew2 * 0.5
	factor = 1
	if numpy.linalg.norm(Gs) < rms_g_thresh*10:
		factor = reduced_factor
	XNew = MaxStep(Xs[-1], XNew, factor)
	runEachStep(XNew, nstep + 1, header_a, header_b, tail1, tail2)
	sk = XNew - Xs[-1]
	GNew, E1, E2 = getG(nstep + 1)
	yk = numpy.mat(GNew - Gs[-1])
	Bk = HessianUpdator(Bs[-1], yk, sk)
	return [XNew, GNew, Bk, E1, E2]

def Converged(E1, E2, X0, X1, G1):
	dE = E2 - E1
	rms = numpy.linalg.norm(X0 - X1)/math.sqrt(n_atom*3)
	max_dis = 0
	for i in range(0, n_atom*3, 3):
		dist = math.sqrt((X0[0,i]-X1[0,i])**2 + (X0[0,i+1]-X1[0,i+1])**2 + (X0[0,i+2]-X1[0,i+2])**2)
		if dist > max_dis:
			max_dis = dist
	rms_g = numpy.linalg.norm(G1)/math.sqrt(n_atom*3)
	max_g = 0
	for i in range(0, n_atom*3, 3):
		g = abs(G1[i])
		if g > max_g:
			max_g = g
	dE_isConv, rms_isConv, max_dis_isConv, max_g_isConv, rms_g_isConv = ['NO', 'NO', 'NO', 'NO', 'NO']
	conv_flag = 0
	if abs(rms) < rms_thresh:
		rms_isConv = 'YES'
		conv_flag += 1
	if abs(dE) < dE_thresh:
		dE_isConv = 'YES'
		conv_flag += 1
	if abs(max_dis) < max_dis_thresh:
		max_dis_isConv = 'YES'
		conv_flag += 1
	if abs(max_g) < max_g_thresh:
		max_g_isConv = 'YES'
		conv_flag += 1
	if abs(rms_g) < rms_g_thresh:
		rms_g_isConv = 'YES'
		conv_flag += 1
	print(f'E1 = {E1}\nE2 = {E2}')
	print(f'deltaE                    {dE:5f}     {dE_thresh:5f}     {dE_isConv}')
	print(f'RMS Gradient              {rms_g:5f}     {rms_g_thresh:5f}     {rms_g_isConv}')
	print(f'Maximium Gradient         {max_g:5f}     {max_g_thresh:5f}     {max_g_isConv}')
	print(f'RMS Displacement          {rms:5f}     {rms_thresh:5f}     {rms_isConv}')
	print(f'Maximium Displacement     {max_dis:5f}     {max_dis_thresh:5f}     {max_dis_isConv}')
	if conv_flag == 5:
		return True
	else:
		return False

n_step = 0 
def runOpt(X0, flag='hybrid'):
	B0 = numpy.eye(numpy.shape(X0)[0])
	X0 = numpy.mat(X0)
	G0, E1, E2 = getG(0)
	n_step = 0
	Xs, Bs, Gs, Es = [[], [], [], []]
	X1, G1, B1, E1, E2 = [None, None, None, 0, 0]
	E0 = E1
	while True:
		if flag == 'pure' or n_step < 3:
			X1, G1, B1, E1, E2 = BFGS(X0, G0, B0, n_step)
		elif flag == 'hybrid':
			X1, G1, B1, E1, E2 = GDIIS(Xs, Gs, Bs, n_step, Es) 
		elif flag == 'gediis':
			X1, G1, B1, E1, E2 = GDIIS(Xs, Gs, Bs, n_step, Es, flag='gediis') 
		if len(Xs) > 3:
			Xs.pop(0)
			Bs.pop(0)
			Gs.pop(0)
			Es.pop(0)
		Xs.append(X1)
		Bs.append(B1)
		Gs.append(G1)
		Es.append(E1)
		n_step += 1
		if Converged(E1, E2, X0, X1, G1):
			return [n_step, X1]
		if n_step > max_steps:
			raise Exception('Maximium number of steps exceeded.')
		X0 = X1
		G0 = G1
		B0 = B1

def main():
	global header_a
	global header_b
	global delete_gbw
	global geom
	global constraints
	print('****KST48 PROGRAM: a powerful tool for MECP locating****')
	print('****By Yumiao Ma, BSJ Institute, 2022/01/09****\n')
	_, inp = sys.argv
	nprocs, mem, charge, mult1, mult2, method, run_mode = inputParser(sys.argv[1])
	header_1, header_2 = buildInitJob(nprocs, mem, charge, mult1, mult2, method, prog, run_mode)
	print(f'Note: This program is now running on {run_mode} mode')
	if run_mode == 'noread':
		delete_gbw = True

	# the preparation phase: run single points or stability calcs to obtain the wavefunction
	print('****Initialization: running the first single point calculations according to the mode****')
	runPrePoint(header_1, header_2, prog, run_mode)
	if run_mode == 'stable' or run_mode == 'inter_read':
		if prog == 'orca':
			print('In an RHF calculation in ORCA, it will not restart automatically if an unstability is found. \
Remember to write UKS when you are handling singlet state! ')
			print('RI is unsupported for stability analysis. It is recommended to MANUALLY obtain the correct wavefunction, \
and then use the read model of KST48, rather than the stable mode, in order to use RI.')
			os.system(f'cp JOBS/pre_A.gbw JOBS/a.gbw')
			os.system(f'cp JOBS/pre_B.gbw JOBS/b.gbw')
		run_mode = 'read'

	# the main loop 
	scan_step = 1
	if run_mode == 'read' or run_mode == 'normal' or run_mode == 'noread':
		method = modifyMethod(prog, method, run_mode)
		print('****Initialization OK, now entering main loop****')
		print('****Before that, please check the keywords list****')
		header_a, header_b = buildInitJob(nprocs, mem, charge, mult1, mult2, method, prog, run_mode, td1=td1, td2=td2)
		print(f'Header A:\n {header_a}')
		print(f'Header B:\n {header_b}')
		print('****If everything is OK, then go to the main loop****\n')
		if constraints != []:
			print(f'Note: This is a contrained optimization. The constraints are added \
				by a harmonic potential rather than the Lagrange method, and therefore the final geometry \
				parameter may not strictly equal to what you wrote in the input file.')

		
		if scans == []: #normal opt without scan
			runEachStep(geom, 0, header_a, header_b, tail1, tail2) # run JOBS/0_A and 0_B, to obtain the first gradient
			geom, force, E = readForceAndGeom('JOBS/0_A.log')
			conv_step, geom = runOpt(geom, flag = 'hybrid')
			print('****Congrats! MECP has converged****')
			if prog == 'gaussian':
				os.system(f'cp JOBS/{conv_step}_A.gjf .')
				os.system(f'cp JOBS/{conv_step}_B.gjf .') #.chk file is already at current folder
			elif prog == 'orca':
				os.system(f'cp JOBS/{conv_step}_A.inp .')
				os.system(f'cp JOBS/{conv_step}_B.inp .')
				os.system(f'cp JOBS/a.gbw .; cp JOBS/b.gbw .')
			os.system(f'cp JOBS/{conv_step}_A.log .')
			os.system(f'cp JOBS/{conv_step}_B.log .')
		# Run PES Scan
		else: #scans: [ [[r,A,B], [start, num, size] ], ... ]
			initial_cons_num = len(constraints)
			scanVars = [[], []]
			scan_para1 = [float(i) for i in scans[0][1]]
			if len(scans) == 1:
				scans.append([[], [0, 0, 0]]) #add one invalid scan to the 2nd dimension if only 1-D scan is requested
			scan_para2 = [float(i) for i in scans[1][1]]
			for i  in range(int(scan_para1[1])):
				scanVars[0].append(scan_para1[0] + i*scan_para1[2])
			for i  in range(int(scan_para2[1])):
				scanVars[1].append(scan_para2[0] + i*scan_para2[2])
			if scanVars[1] == []:
				scanVars[1].append(-1) # make sure it contains something  for the following loop to run
			for i in scanVars[0]:
				new_const = scans[0][0][:]
				new_const.append(i)
				constraints.append(new_const)
				for j in scanVars[1]:
					new_const = scans[1][0][:]
					new_const.append(j)
					constraints.append(new_const)
					print(f'constraints after pop and append')
					print(constraints)

					print(f'****Scan Cycle {i:4f}_{j:4f}****')
					runEachStep(geom, 0, header_a, header_b, tail1, tail2) # run JOBS/0_A and 0_B, to obtain the first gradient
					geom, force, E = readForceAndGeom('JOBS/0_A.log')
					conv_step, geom = runOpt(geom, flag = 'hybrid')
					numpy.mat(geom)
					print('****Congrats! MECP has converged****\n')
					if prog == 'gaussian':
						os.system(f'cp JOBS/{conv_step}_A.gjf {i:4f}_{j:4f}.gjf')
						os.system(f'cp JOBS/{conv_step}_B.gjf {i:4f}_{j:4f}.gjf') #.chk file is already at current folder
					elif prog == 'orca':
						os.system(f'cp JOBS/{conv_step}_A.inp {i:4f}_{j:4f}.inp')
						os.system(f'cp JOBS/{conv_step}_B.inp {i:4f}_{j:4f}.inp')
						#os.system(f'cp JOBS/a.gbw .; cp JOBS/b.gbw .')
					os.system(f'cp JOBS/{conv_step}_A.log {i:4f}_{j:4f}.log')
					os.system(f'cp JOBS/{conv_step}_B.log {i:4f}_{j:4f}.log')
					constraints.pop(-1)

				constraints.pop(-1)
		



main()
