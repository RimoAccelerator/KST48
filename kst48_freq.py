import numpy, math, re

def readForceAndGeomForGaussian(path):
    forceArr = []
    geomArr = []
    with open(path) as f:
        isForce = False
        isGeom = False
        E = 0
        archivePart = ''
        #isArchive = False
        for l in f.readlines():
            if 'Input orientation' in l:
                isGeom = True
                geomArr = []
            elif 'Distance matrix' in l or 'Rotational constants' in l:
                isGeom = False
            elif 'Forces (Hartrees/Bohr)' in l:
                isForce = True
            elif 'Cartesian Forces:  Max' in l:
                isForce = False
            elif 'SCF DONE' in l.upper():
                E = float(l.split('=')[1].upper().split('A.U.')[0])
            elif 'extrapolated energy' in l.lower():
                E = float(l.split('=')[1])
            elif 'E(TD-HF/TD-DFT)' in l.upper():
                E = float(l.split('=')[1])
            elif isForce and (re.match("\\s*[0-9]+\\s+[0-9]+\\s*\\-*[0-9]+", l) is not None):
                forceArr.extend(l.split()[2:])
            elif isGeom and (re.match("\\s*[0-9]+\\s+[0-9]+\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                geomArr.extend(l.split()[3:])
    geomArr = [float(i) for i in geomArr]
    # forceArr = [float(i)/0.529 for i in forceArr]
    forceArr = [float(i) for i in forceArr]
    return [geomArr, forceArr, E]

def readHessianForGaussian(path):
    N = 0
    isHessianCartesian = False
    hessian = numpy.zeros(0)
    with open(path) as f:
        jFrom = 0
        jTo = 0
        for l in f:
            if 'NAtoms=' in l and hessian.shape[0] == 0:
                N = int(l.split()[1]) # number of atoms
                hessian = numpy.zeros((3 * N, 3 * N))
            if 'Force constants in Cartesian coordinates' in l:
                isHessianCartesian = True
            elif isHessianCartesian:
                if not l.split()[0].isdigit():
                    break
                if not '.' in l: #header like '  1             2             3             4             5'
                    jTo = int(l.split()[-1]) - 1
                    jFrom = int(l.split()[0]) - 1
                else: # data line like '2  0.000000D+00  0.549637D+00'
                    lSplitted = l.split()
                    thisI = int(lSplitted[0]) - 1
                    thisJ = jFrom
                    for datum in lSplitted[1:]:
                        hessian[thisI, thisJ] = float(datum.replace('D', 'E'))
                        hessian[thisJ, thisI] = float(datum.replace('D', 'E'))
                        thisJ += 1
        if not isHessianCartesian: #No Hessian part found
            raise Exception('No force constant matrix found. Please check your file. Do not forget #p.')
    return hessian

def readGradFromFchk(path):
	N = 0
	isGradient = False
	g = []
	with open(path) as f:
		for l in f:
			if 'Number of atoms' in l and len(g) == 0:
				N = int(l.split('I')[1]) # number of atoms
			if 'Cartesian Gradient' in l:
			 	isGradient = True
			elif isGradient:
				if not l.split()[0][0] == '-' and not l.split()[0][0].isdigit():
					break
				for datum in l.split():
					g.append(float(datum))
		if not isGradient: #No Hessian part found
			raise Exception('No gradient found. Please check your file. ')
	if len(g) != N * 3:
		raise Exception('Something wrong occured when reading gradient.')
	return g

def readGeomFromFchk(path):
	N = 0
	isGeom = False
	g = []
	with open(path) as f:
		for l in f:
			if 'Number of atoms' in l and len(g) == 0:
				N = int(l.split('I')[1]) # number of atoms
			if 'Current cartesian coordinates' in l:
				isGeom = True
			elif isGeom:
				if not l.split()[0][0] == '-' and not l.split()[0][0].isdigit():
					break
				for datum in l.split():
					g.append(float(datum))
		if not isGeom: #No Hessian part found
			raise Exception('No geometry found. Please check your file. ')
	if len(g) != N * 3:
		raise Exception('Something wrong occured when reading geometry.')
	return g

def readHessianForGaussianFromFchk(path):
	N = 0
	isHessianCartesian = False
	hessian = numpy.zeros(0)
	thisI = 0
	thisJ = 0
	with open(path) as f:
		for l in f:
			if 'Number of atoms' in l and hessian.shape[0] == 0:
				N = int(l.split('I')[1]) # number of atoms
				hessian = numpy.zeros((3 * N, 3 * N))
			if 'Cartesian Force Constants' in l:
				isHessianCartesian = True
			elif isHessianCartesian:
				if not l.split()[0][0] == '-' and not l.split()[0][0].isdigit():
					break
				for datum in l.split():
					hessian[thisI, thisJ] = float(datum)
					hessian[thisJ, thisI] = float(datum)
					thisJ += 1
					if thisJ > thisI:
						thisI += 1
						thisJ = 0
		if not isHessianCartesian: #No Hessian part found
		    raise Exception('No force constant matrix found. Please check your file. Do not forget #p.')
	return hessian

def readMassForGaussian(path):
    masses = []
    with open(path) as f:
        for l in f:
            if 'has atomic number' in l:
                masses.append(float(l.split()[-1]))
    return masses

def readMassForGaussianFromFchk(path):
	N = 0
	isMass = False
	m = []
	with open(path) as f:
		for l in f:
			if 'Number of atoms' in l and len(m) == 0:
				N = int(l.split('I')[1]) # number of atoms
			if 'Vib-AtMass' in l:
				isMass = True
			elif isMass:
				if not l.split()[0][0] == '-' and not l.split()[0][0].isdigit():
					break
				for datum in l.split():
					m.append(float(datum))
		if not isMass: #No Hessian part found
			raise Exception('No mass found. Please check your file. ')
	if len(m) != N:
		raise Exception('Something wrong occured when reading mass.')
	return m

def buildInverseSqrtM(masses):
    inverseSqrtM = numpy.zeros((3 * len(masses), 3 * len(masses)))
    for i in range(len(masses)):
        thisInverseSqrtM = 1 / math.sqrt(masses[i])
        inverseSqrtM[i * 3 + 0, i * 3 + 0] = thisInverseSqrtM
        inverseSqrtM[i * 3 + 1, i * 3 + 1] = thisInverseSqrtM
        inverseSqrtM[i * 3 + 2, i * 3 + 2] = thisInverseSqrtM
    return inverseSqrtM

def buildProjectionForMotion(masses, geom): #build the projecting matrix for translation and rotation. Input both python arrays
    DMat = numpy.zeros((len(masses) * 3, 6))
    for i in range(len(masses)):
        DMat[i * 3 + 0, 0] = math.sqrt(masses[i])
        DMat[i * 3 + 1, 1] = math.sqrt(masses[i])
        DMat[i * 3 + 2, 2] = math.sqrt(masses[i])
        DMat[i * 3 + 1, 3] = -math.sqrt(masses[i]) * geom[i * 3 + 2]
        DMat[i * 3 + 2, 3] = +math.sqrt(masses[i]) * geom[i * 3 + 1]
        DMat[i * 3 + 0, 4] = +math.sqrt(masses[i]) * geom[i * 3 + 2]
        DMat[i * 3 + 2, 4] = -math.sqrt(masses[i]) * geom[i * 3 + 0]
        DMat[i * 3 + 0, 5] = -math.sqrt(masses[i]) * geom[i * 3 + 1]
        DMat[i * 3 + 1, 5] = +math.sqrt(masses[i]) * geom[i * 3 + 0]
    DMat, R = numpy.linalg.qr(DMat) 
    return numpy.eye(len(masses) * 3) - DMat @ DMat.T

def buildProjectionForCP(f1, f2): #build the projecting matrix for translation and rotation. Input both numpy arrays
    x1 = numpy.array([f1 - f2]).T
    x1 /= numpy.linalg.norm(x1)
    #return x1 @ x1.T
    return numpy.eye(len(x1)) - x1 @ x1.T

def buildHessianForCP(f1, f2, HMat1, HMat2): #build the projecting matrix for translation and rotation. Input both numpy arrays
    x1 = f1 - f2
    factor = numpy.linalg.norm(f1) / numpy.linalg.norm(f1 - f2)
    H = HMat1 * (1 - factor) + HMat2 * factor
    #H = numpy.linalg.norm(f1) * HMat2 - numpy.linalg.norm(f2) * HMat1
    #H /= numpy.linalg.norm(f1 - f2)
    return H

def checkCosine(f1, f2):
    cosine = numpy.dot(f1, f2) / numpy.linalg.norm(f1) / numpy.linalg.norm(f2)
    maxError = 0.05
    print(f'Cosine for the angle formed by f1 and f2 is {cosine}')
    if abs(abs(cosine) - 1) > maxError:
        print(f'Warning: f1 and f2 not colinear. ')

def writeFreqsAndModesOneBatch(freqs, modes, startIndex): #print Gaussian freq info for less than 3 modes
    output = ''
    headerLines = [' ', \
    ' ', \
    ' Frequencies --   ', \
    ' Red. masses --   ', \
    ' IR Inten    --   ', \
    '  Atom  AN      ']
    for i in range(len(freqs)):
        headerLines[0] += '                      '
        headerLines[0] += str(i + 1 + startIndex)
        headerLines[1] += '                      '
        headerLines[1] += 'A'
        headerLines[2] += f'{freqs[i]:9.4f}              '
        headerLines[3] += '   0.0000              '
        headerLines[4] += '   0.0000              '
        headerLines[5] += 'X      Y      Z        '
    output += '\n'.join(headerLines)
    output += '\n'
    NAtoms = int(len(modes[0]) / 3)
    modeLines = [''] * NAtoms
    for n in range(NAtoms):
        modeLines[n] = f'{n + 1:6d}   1  '
    for i in range(len(freqs)):
        for n in range(NAtoms):
            modeLines[n] += f'{modes[i][n * 3 + 0].real:7.2f}'
            modeLines[n] += f'{modes[i][n * 3 + 1].real:7.2f}'
            modeLines[n] += f'{modes[i][n * 3 + 2].real:7.2f}'
            modeLines[n] += '  '
    output += '\n'.join(modeLines)
    output += '\n'
    return output


def writeFreqsAndModes(freqs, modes):
    output = '''
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:\n'''
    threeFreqs = []
    threeModes = []
    thisFreqIndex = -1
    startIndex = 0
    while True:
        if thisFreqIndex == len(freqs):
            break
        for i in range(3):
            thisFreqIndex += 1
            if thisFreqIndex == len(freqs):
                break
            if abs(freqs[thisFreqIndex]) > 5: #only freqs greater than 5 cm-1 will be written
                startIndex += 1
                threeFreqs.append(freqs[thisFreqIndex])
                threeModes.append(modes[:, thisFreqIndex])
        if len(threeFreqs) > 0:
            output += writeFreqsAndModesOneBatch(threeFreqs, threeModes, startIndex - len(threeFreqs))
        threeFreqs = []
        threeModes = []
    return output

def getReplacedLog(path, data):
    output = ''
    isFreq = False
    with open(path) as f:
        for l in f.readlines():
            if not isFreq and 'Harmonic frequencies (cm**-1), IR intensities (KM/Mole)' in l:
                isFreq = True
                output += data
            if '-------------------' in l and isFreq:
                isFreq = False
            if not isFreq:
                output += l
    return output

def main():
    pathFreqA = input('Freq fchk file for state A?')
    HCartesianA = readHessianForGaussianFromFchk(pathFreqA)
    
    pathFreqB = input('Freq fchk file for state B?')
    HCartesianB = readHessianForGaussianFromFchk(pathFreqB)
    
    geomArr = readGeomFromFchk(pathFreqA)
    f1Arr = readGradFromFchk(pathFreqA)
    f2Arr = readGradFromFchk(pathFreqB)
    #geomArr, f1Arr, E1 = readForceAndGeomForGaussian(pathFreqA)
    #geomArr, f2Arr, E2 = readForceAndGeomForGaussian(pathFreqB)
    f1 = numpy.array(f1Arr)
    f2 = numpy.array(f2Arr)
    
    checkCosine(f1,f2) # f1 and f2 must be colinear for a true MECP
    HEff = buildHessianForCP(f1, f2, HCartesianA, HCartesianB)
    
    PMat1 = buildProjectionForCP(f1,f2)
    HEff = PMat1 @ HEff @ PMat1 # Now the direction along d(E1-E2)/dq is projected out
    
    masses = readMassForGaussianFromFchk(pathFreqA)
    inverseSqrtM = buildInverseSqrtM(masses)
    HEff = inverseSqrtM @ HEff @ inverseSqrtM # Transfomred to mass-weighted coord
    PMat2 = buildProjectionForMotion(masses, geomArr)
    HEff = PMat2 @ HEff @ PMat2
    
    eigs, modes = numpy.linalg.eig(HEff)
    freqs = []
    for eig in eigs:
        freq = math.sqrt(abs(eig)) / 2 / math.pi * 3.24318E4
        if eig < 0:
            freq = -freq
        freqs.append(freq)
    freqs = numpy.array(freqs)
    sortedIndices = numpy.argsort(freqs)
    freqs = freqs[sortedIndices]
    modes = modes[:, sortedIndices]
    print('Frequencies for the HEff of CP:')
    print(freqs)

    print("Now it's time to output something.\
    Give me a freq log file so that I could replace the vibration information\
    and output a new file to kst48_freq.out. Note that only the frequencies and modes are replaced, and all the other things are remained.\
    It can be read by GoodVibes to obtain free energy correction.")
    pathLog = input('Your freq .log file?')
    freqSection = writeFreqsAndModes(freqs, modes)
    #print(freqSection)
    with open('kst48_freq.out', 'w') as f:
        f.write(getReplacedLog(pathLog, freqSection))

main()