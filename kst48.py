import os
import sys
import re
import numpy
import math

THRESH_dE = 0.000050
THRESH_RMS = 0.0025
THRESH_MAX_DIS = 0.004
THRESH_MAX_G = 0.0007
THRESH_RMS_G = 0.0005
MAX_STEPS = 100
MAX_STEP_SIZE = 0.1
NUM_ATOM = 0
REDUCED_FACTOR = 0.5

LIST_ELEMENT = []
HEADER_A = ''
HEADER_B = ''
TAIL1 = ''
TAIL2 = ''
CONSTRAINTS = []
SCANS = []
GEOM = []

# constraints = [['R',1,2,1.0]]
# constraints = [['A',2,1,3,90],['R',1,2,1.2],['R',1,3,1.2]]
PROG = 'gaussian'
PROG_COMM = 'g16'
METHOD = ''
TD1 = ''
TD2 = ''
MP2 = False
DELETE_GBW = False

LAMBDAS = []


def buildHeader(
        NProcs,
        Mem,
        Charge,
        Mult1,
        Mult2,
        Method,
        Prog,
        Td1='',
        Td2=''):
    header_1, header_2 = ['', '']
    if Prog == 'gaussian':
        header_1 = f'%chk=a.chk\n%nprocs={NProcs} \n%mem={Mem} \n# {Method} {Td1} nosymm\n\n Title Card \n\n{Charge} {Mult1}'
        header_2 = f'%chk=b.chk\n%nprocs={NProcs} \n%mem={Mem} \n# {Method} {Td2} nosymm\n\n Title Card \n\n{Charge} {Mult2}'
    elif Prog == 'orca':
        header_1 = f'%pal nprocs {NProcs} end\n%maxcore {Mem} \n! {Method} \n\n *xyz {Charge} {Mult1}'
        header_2 = f'%pal nprocs {NProcs} end\n%maxcore {Mem} \n! {Method} \n\n *xyz {Charge} {Mult2}'
        if '***' in header_1:
            header_1 = header_1.replace('***', 'JOBS/a.gbw')
        if '***' in header_2:
            header_2 = header_2.replace('***', 'JOBS/b.gbw')
    return [header_1, header_2]


def buildInitJob(
        NProcs,
        Mem,
        Charge,
        Mult1,
        Mult2,
        Method,
        Prog,
        RunMode,
        Td1='',
        Td2=''):
    header_1, header_2 = ['', '']
    if RunMode == 'normal' or RunMode == 'read' or RunMode == 'noread':
        header_1, header_2 = buildHeader(
            NProcs, Mem, Charge, Mult1, Mult2, Method, Prog, Td1=Td1, Td2=Td2)
    elif RunMode == 'stable' or RunMode == 'inter_read':
        if Prog == 'gaussian':
            header_1, header_2 = buildHeader(
                NProcs, Mem, Charge, Mult1, Mult2, Method + ' stable=opt', Prog)
        elif Prog == 'orca':
            header_1, header_2 = buildHeader(NProcs, Mem, Charge, Mult1, Mult2, Method +
                                             '\n %scf stabperform true StabRestartUHFifUnstable true end \n', Prog)
        else:
            raise Exception('This PROGram not supported for STABLE mode!')
    else:
        raise Exception('Unknown mode!')
    if '"***"' in header_1:
        header_1 = header_1.replace('"***"', f'"JOBS/a.gbw"')
    if '"***"' in header_2:
        header_1 = header_2.replace('"***"', f'"JOBS/b.gbw"')
    return [header_1, header_2]


def runPrePoint(Header_1, Header_2, Prog, RunMode='normal'):
    # os.system('cp *.chk JOBS; cp *.gbw JOBS')
    if Prog == 'gaussian':
        writeGjf(GEOM, Header_1, TAIL1, f'JOBS/pre_A.gjf')
        writeGjf(GEOM, Header_2, TAIL2, f'JOBS/pre_B.gjf')
        os.system(f'{PROG_COMM} JOBS/pre_B.gjf')
        if RunMode == 'inter_read':
            os.system(f'cp b.chk a.chk')
            os.system(f'sed -i "s/#/# guess=(read,mix)/g" JOBS/pre_A.gjf')
        os.system(f'{PROG_COMM} JOBS/pre_A.gjf')

    elif Prog == 'orca':
        writeORCA(GEOM, Header_1, TAIL1, f'JOBS/pre_A.inp')
        writeORCA(GEOM, Header_2, TAIL2, f'JOBS/pre_B.inp')
        os.system(f'{PROG_COMM} JOBS/pre_B.inp > JOBS/pre_B.out')
        os.system(f'cp JOBS/pre_B.gbw JOBS/b.gbw')
        if RunMode == 'inter_read':
            os.system(f'cp JOBS/pre_B.gbw JOBS/a.gbw')
            print(
                'Note: The inter_read mode is set for ORCA. In Gaussian, the PROGram automatically add guess=mix for state A, \
but this will not be done for ORCA. If you want to converge to correct OSS wavefunction from a triplet wavefunction, guess=(read,mix) \
is always beneficial. So please do not forget to add relevant convergence controlling in your ORCA tail part.')
        os.system(f'{PROG_COMM} JOBS/pre_A.inp > JOBS/pre_A.out')
        os.system(f'cp JOBS/pre_A.gbw JOBS/a.gbw')

    # elif PROG == 'xtb':
    #	writeXYZ(geom, f'JOBS/pre_A.xyz')
    #	writeXYZ(geom, f'JOBS/pre_B.xyz')
    #	os.system(f'{prog_comm} JOBS/pre.xyz --engrad')
    #	os.system(f'{prog_comm} JOBS/pre.xyz --engrad')


def inputParser(path):
    global GEOM
    global LIST_ELEMENT
    global NUM_ATOM
    global PROG
    global TD1
    global TD2
    global MP2
    global PROG_COMM
    global TAIL1
    global TAIL2
    global REDUCED_FACTOR
    global MAX_STEPS
    global MAX_STEP_SIZE
    charge, mult1, mult2, method, nprocs, mem = ['', '', '', '', '', '']
    command = {'gau': '', 'orca': '', 'xtb': ''}
    runMode = 'normal'
    with open(path) as f:
        isGEOM = False
        isTAIL1 = False
        isTAIL2 = False
        isConst = False
        for l in f:
            l_bk = l
            l = l.lower().strip()
            if '#' in l:
                l = l.split('#')[0]
            if '*geom' in l:
                isGEOM = True
            elif l == '*' and isGEOM:
                isGEOM = False
            elif '*tail1' in l:
                isTAIL1 = True
            elif '*tail2' in l:
                isTAIL2 = True
            elif l == '*' and isTAIL2:
                isTAIL2 = False
            elif l == '*' and isTAIL1:
                isTAIL1 = False
            elif '*constr' in l:
                isConst = True
            elif l == '*' and isConst:
                isConst = False
            elif isGEOM and (re.match("\\s*\\S+\\s*\\-*[0-9]+", l)):
                GEOM.extend(l.split()[1:])
                LIST_ELEMENT.append(l.split()[0])
                NUM_ATOM += 1
            elif isTAIL1:
                TAIL1 += l
                TAIL1 += '\n'
            elif isTAIL2:
                TAIL2 += l
                TAIL2 += '\n'
            elif isConst and l != '':
                l_splitted = l.split()
                if (l_splitted[0] != 's'):
                    CONSTRAINTS.append(l.split())
                else:
                    if l_splitted[1] == 'r':
                        SCANS.append([l_splitted[1:4], l_splitted[4:]])
                    elif l_splitted[1] == 'a':
                        SCANS.append([l_splitted[1:5], l_splitted[5:]])
            if '=' in l and not isTAIL1 and not isTAIL2:
                parameter = '='.join(l.split('=')[1:]).strip()
                parameter_bk = '='.join(l_bk.split('=')[1:]).strip()
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
                elif 'thresh_de' in l:
                    THRESH_dE = float(parameter)
                elif 'thresh_rms' in l:
                    THRESH_RMS = float(parameter)
                elif 'thresh_max_dis' in l:
                    THRESH_MAX_DIS = float(parameter)
                elif 'thresh_max_g' in l:
                    THRESH_MAX_G = float(parameter)
                elif 'thresh_rms_g' in l:
                    THRESH_RMS_G = float(parameter)
                elif 'max_steps' in l:
                    MAX_STEPS = float(parameter)
                elif 'max_step_size' in l:
                    MAX_STEP_SIZE = float(parameter)
                elif 'program' in l:
                    PROG = parameter
                elif 'gau_comm' in l:
                    command['gaussian'] = parameter_bk
                elif 'orca_comm' in l:
                    command['orca'] = parameter_bk
                elif 'xtb_comm' in l:
                    command['xtb'] = parameter_bk
                elif 'mode' in l:
                    runMode = parameter
                elif 'td1' in l:
                    TD1 = parameter
                elif 'td2' in l:
                    TD2 = parameter
                elif 'mp2' in l and parameter == 'true':
                    MP2 = True
                elif 'reduced_factor' in l:
                    REDUCED_FACTOR = float(parameter)
    GEOM = numpy.mat(GEOM)
    PROG_COMM = command[PROG]
    return [nprocs, mem, charge, mult1, mult2, method, runMode]


def modifyMETHOD(Prog, Method, RunMode):
    if Prog == 'gaussian':
        Method += ' force '
    elif Prog == 'orca':
        Method += ' engrad '
    # if noread mode is set, do not read the initial wavefunctions
    if RunMode != 'noread':
        if Prog == 'gaussian':
            Method += ' guess=read '
        elif Prog == 'orca':
            Method += '\n!moread \n %moinp "***"\n'
    return Method


def readForceAndGeomForGaussian(path):
    forceArr = []
    geomArr = []
    with open(path) as f:
        isForce = False
        isGeom = False
        E = 0
        archivePart = ''
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
            elif MP2 and 'MP2=' in l.upper():  # read the MP2 energy in the archive part of a log file
                archivePart += l.upper().strip()
                isArchive = True
            elif isArchive:
                # sometimes the MP2 energy may be separated by a \n. Two lines
                # have to be combined
                archivePart += l.upper().strip()
                isArchive = False
                E = float(archivePart.split('MP2=')[1].split('=')[0])
            elif isForce and (re.match("\\s*[0-9]+\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                forceArr.extend(l.split()[2:])
            elif isGeom and (re.match("\\s*[0-9]+\\s*[0-9]+\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                geomArr.extend(l.split()[3:])
    geomArr = [float(i) for i in geomArr]
    # forceArr = [float(i)/0.529 for i in forceArr]
    forceArr = [float(i) for i in forceArr]
    forceArr = addConstLag(geomArr, forceArr, CONSTRAINTS)
    if len(geomArr) == NUM_ATOM * 3:
        geomArr.extend(LAMBDAS)
    return [geomArr, forceArr, E]


def readForceAndGeomForORCA(path):
    forceArr = []
    geomArr = []
    splittedPath = path.split('.')
    # replace A.xxx into A.engrad
    path = ''.join(splittedPath[:-1]) + '.engrad'
    with open(path) as f:
        isForce = False
        isGeom = False
        isEnergy = False
        E = 0
        for l in f.readlines():
            if 'The atomic numbers and current coordinates in Bohr' in l:
                isGeom = True
            elif '#' in l and len(geomArr) > 3:
                isGeom = False
            elif 'The current gradient' in l:
                isForce = True
            elif '#' in l and len(forceArr) > 3:
                isForce = False
            elif 'current total energy in Eh' in l:
                isEnergy = True
            elif '#' in l and E != 0:
                isEnergy = False
            elif isEnergy and (re.match("\\s*\\-*[0-9]+", l) is not None):
                E = float(l.strip())
            elif isForce and (re.match("\\s*\\-*[0-9]+", l) is not None):
                forceArr.append(l.strip())
            elif isGeom and (re.match("\\s*[0-9]+\\s*\\-*[0-9]+", l) is not None):
                geomArr.extend(l.split()[1:])
    geomArr = [float(i) * 0.52918 for i in geomArr]  # change Bohr to angstrom
    # ORCA outputs gradients. Here it is adapted to gaussian, which is the
    # force
    forceArr = [-float(i) for i in forceArr]
    forceArr = addConstLag(geomArr, forceArr, CONSTRAINTS)
    if len(geomArr) == NUM_ATOM * 3:
        geomArr.extend(LAMBDAS)
    return [geomArr, forceArr, E]


def readForceAndGeom(path):
    if PROG == 'gaussian':
        return readForceAndGeomForGaussian(path)
    elif PROG == 'orca':
        return readForceAndGeomForORCA(path)
    else:
        raise Exception('Unsupported program!')


# read constraints and add it to the gradients according to a harmonic term; which is deprecated because of the slow convergence
# constr list:[['r', 1 , 2, 1.5], ['a', 1, 2, 3]], starting from 1
def addConst(geomArr, forceArr, constrList):
    constrR = []
    constrA = []
    constrD = []
    newForce = forceArr[:]
    applyForce = 10
    for i in constrList:
        if i[-1] == -1:  # skip invalid constrain. This feature is used in the PES-scan part
            continue
        if i[0] == 'r':
            constrR.append([int(i[1]) - 1, int(i[2]) - 1, float(i[3])])
        if i[0] == 'a':
            constrA.append([int(i[1]) - 1, int(i[2]) - 1, int(i[3]) - 1, math.cos(float(i[4]) / 180.0 * math.pi)])
    for i in constrR:  # [atom_a, atom_b, dist]
        x1, y1, z1 = geomArr[i[0] * 3: i[0] * 3 + 3]
        x2, y2, z2 = geomArr[i[1] * 3: i[1] * 3 + 3]
        r12 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        newForce[i[0] * 3 + 0] += applyForce * (x1 - x2) * (1 - i[2] / r12)
        newForce[i[0] * 3 + 1] += applyForce * (y1 - y2) * (1 - i[2] / r12)
        newForce[i[0] * 3 + 2] += applyForce * (z1 - z2) * (1 - i[2] / r12)
        newForce[i[1] * 3 + 0] -= applyForce * (x1 - x2) * (1 - i[2] / r12)
        newForce[i[1] * 3 + 1] -= applyForce * (y1 - y2) * (1 - i[2] / r12)
        newForce[i[1] * 3 + 2] -= applyForce * (z1 - z2) * (1 - i[2] / r12)
    for i in constrA:
        # constraint for angle ABC: lambda * ((BA dot BC)/|BA||BC| - cosAngle)^2
        # its deriv: lambda * (cos_current - target) * d(cos)
        xA, yA, zA = geomArr[i[0] * 3: i[0] * 3 + 3]
        xB, yB, zB = geomArr[i[1] * 3: i[1] * 3 + 3]
        xC, yC, zC = geomArr[i[2] * 3: i[2] * 3 + 3]
        rBA = math.sqrt((xA - xB) ** 2 + (yA - yB) ** 2 + (zA - zB) ** 2)
        rBC = math.sqrt((xC - xB) ** 2 + (yC - yB) ** 2 + (zC - zB) ** 2)
        BAdotBC = (xA - xB) * (xC - xB) + (yA - yB) * (yC - yB) + (zA - zB) * (zC - zB)
        cosAngle_current = BAdotBC / rBA / rBC
        cosTarget = i[3]
        dCos_dxA = (rBA * (xC - xB) - BAdotBC * (xA - xB)) / (rBC * rBA ** 2)
        dCos_dxC = (rBC * (xA - xB) - BAdotBC * (xC - xB)) / (rBA * rBC ** 2)
        dCos_dyA = (rBA * (yC - yB) - BAdotBC * (yA - yB)) / (rBC * rBA ** 2)
        dCos_dyC = (rBC * (yA - yB) - BAdotBC * (yC - yB)) / (rBA * rBC ** 2)
        dCos_dzA = (rBA * (zC - zB) - BAdotBC * (zA - zB)) / (rBC * rBA ** 2)
        dCos_dzC = (rBC * (zA - zB) - BAdotBC * (zC - zB)) / (rBA * rBC ** 2)
        dCos_dxB = (rBA * rBC * (2 * xB - xA - xC) + BAdotBC * xB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        dCos_dyB = (rBA * rBC * (2 * yB - yA - yC) + BAdotBC * yB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        dCos_dzB = (rBA * rBC * (2 * zB - zA - zC) + BAdotBC * zB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        newForce[i[0] * 3 + 0] += applyForce * (cosAngle_current - cosTarget) * dCos_dxA
        newForce[i[0] * 3 + 1] += applyForce * (cosAngle_current - cosTarget) * dCos_dyA
        newForce[i[0] * 3 + 2] += applyForce * (cosAngle_current - cosTarget) * dCos_dzA
        newForce[i[1] * 3 + 0] += applyForce * (cosAngle_current - cosTarget) * dCos_dxB
        newForce[i[1] * 3 + 1] += applyForce * (cosAngle_current - cosTarget) * dCos_dyB
        newForce[i[1] * 3 + 2] += applyForce * (cosAngle_current - cosTarget) * dCos_dzB
        newForce[i[2] * 3 + 0] += applyForce * (cosAngle_current - cosTarget) * dCos_dxC
        newForce[i[2] * 3 + 1] += applyForce * (cosAngle_current - cosTarget) * dCos_dyC
        newForce[i[2] * 3 + 2] += applyForce * (cosAngle_current - cosTarget) * dCos_dzC
    return newForce


# add constraints with the Lagrange method
# constr list:[['r', 1 , 2, 1.5], ['a', 1, 2, 3]], starting from 1
def addConstLag(geomArr, forceArr, constrList):
    global LAMBDAS
    constrR = []
    constrA = []
    constrD = []
    newForce = forceArr[:]
    applyForce = 10
    for i in constrList:
        if i[-1] == -1:  # skip invalid constrain. This feature is used in the PES-scan part
            continue
        if i[0] == 'r':
            constrR.append([int(i[1]) - 1, int(i[2]) - 1, float(i[3])])
        if i[0] == 'a':
            constrA.append([int(i[1]) - 1, int(i[2]) - 1, int(i[3]) - 1, math.cos(float(i[4]) / 180.0 * math.pi)])
    totalConst = len(constrR) + len(constrA)
    if totalConst == 0:
        return newForce
    index = 0
    dL_dlambda = []
    dCondition_dR = numpy.zeros(shape=(totalConst, NUM_ATOM * 3))
    for j in range(len(constrR)):  # [atom_a, atom_b, dist]
        i = constrR[index]
        x1, y1, z1 = geomArr[i[0] * 3: i[0] * 3 + 3]
        x2, y2, z2 = geomArr[i[1] * 3: i[1] * 3 + 3]
        r12 = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        dCondition_dR[index, i[0] * 3 + 0] = + (x1 - x2) / r12
        dCondition_dR[index, i[0] * 3 + 1] = + (y1 - y2) / r12
        dCondition_dR[index, i[0] * 3 + 2] = + (z1 - z2) / r12
        dCondition_dR[index, i[1] * 3 + 0] = - (x1 - x2) / r12
        dCondition_dR[index, i[1] * 3 + 1] = - (y1 - y2) / r12
        dCondition_dR[index, i[1] * 3 + 2] = - (z1 - z2) / r12
        dL_dlambda.append(r12 - i[2])
        index += 1
    for i in constrA:
        # constraint for angle ABC: lambda * ((BA dot BC)/|BA||BC| - cosAngle)^2
        # its deriv: lambda * (cos_current - target) * d(cos)
        xA, yA, zA = geomArr[i[0] * 3: i[0] * 3 + 3]
        xB, yB, zB = geomArr[i[1] * 3: i[1] * 3 + 3]
        xC, yC, zC = geomArr[i[2] * 3: i[2] * 3 + 3]
        rBA = math.sqrt((xA - xB) ** 2 + (yA - yB) ** 2 + (zA - zB) ** 2)
        rBC = math.sqrt((xC - xB) ** 2 + (yC - yB) ** 2 + (zC - zB) ** 2)
        BAdotBC = (xA - xB) * (xC - xB) + (yA - yB) * (yC - yB) + (zA - zB) * (zC - zB)
        cosAngle_current = BAdotBC / rBA / rBC
        cosTarget = i[3]
        dCos_dxA = (rBA * (xC - xB) - BAdotBC * (xA - xB)) / (rBC * rBA ** 2)
        dCos_dxC = (rBC * (xA - xB) - BAdotBC * (xC - xB)) / (rBA * rBC ** 2)
        dCos_dyA = (rBA * (yC - yB) - BAdotBC * (yA - yB)) / (rBC * rBA ** 2)
        dCos_dyC = (rBC * (yA - yB) - BAdotBC * (yC - yB)) / (rBA * rBC ** 2)
        dCos_dzA = (rBA * (zC - zB) - BAdotBC * (zA - zB)) / (rBC * rBA ** 2)
        dCos_dzC = (rBC * (zA - zB) - BAdotBC * (zC - zB)) / (rBA * rBC ** 2)
        dCos_dxB = (rBA * rBC * (2 * xB - xA - xC) + BAdotBC * xB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        dCos_dyB = (rBA * rBC * (2 * yB - yA - yC) + BAdotBC * yB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        dCos_dzB = (rBA * rBC * (2 * zB - zA - zC) + BAdotBC * zB * (rBA ** 2 + rBC ** 2) / rBA / rBC) / (
                    rBC ** 2 * rBA ** 2)
        dCondition_dR[index, i[0] * 3 + 0] = dCos_dxA
        dCondition_dR[index, i[0] * 3 + 1] = dCos_dyA
        dCondition_dR[index, i[0] * 3 + 2] = dCos_dzA
        dCondition_dR[index, i[1] * 3 + 0] = dCos_dxB
        dCondition_dR[index, i[1] * 3 + 1] = dCos_dyB
        dCondition_dR[index, i[1] * 3 + 2] = dCos_dzB
        dCondition_dR[index, i[2] * 3 + 0] = dCos_dxC
        dCondition_dR[index, i[2] * 3 + 1] = dCos_dyC
        dCondition_dR[index, i[2] * 3 + 2] = dCos_dzC
        dL_dlambda.append(cosAngle_current - cosTarget)
        index += 1
    if LAMBDAS == []:  # now generate lambda = C_T*G/(C.T*C), G=original gradient, C=dConstraint/dx
        LAMBDAS = [(-float(dCondition_dR[i, :] * numpy.mat(forceArr).T) \
                    / numpy.dot(dCondition_dR[0, :], dCondition_dR[0, :].T)) \
                   * (len(constrR) + len(constrA)) for i in range(totalConst)]
    # LAMBDAS = [0.5] * totalConst
    print(LAMBDAS)
    dCondition_dR = numpy.mat(LAMBDAS) * dCondition_dR
    for i in range(len(newForce)):
        for j in range(totalConst):
            newForce[i] += dCondition_dR[j, i]
    newForce.extend(dL_dlambda)
    return newForce


def writeGjf(Geom, Header, Tail, Name):
    f = open(Name, 'w+')
    f.write(Header)
    f.write('\n')
    for i in range(NUM_ATOM):
        f.write('{ele}  {x}  {y}  {z}'.format(
            ele=LIST_ELEMENT[i], x=Geom[0, i * 3], y=Geom[0, i * 3 + 1], z=Geom[0, i * 3 + 2]))
        f.write('\n')
    f.write('\n')
    f.write(Tail)
    f.write('\n')
    f.close()


def writeORCA(Geom, Header, Tail, Name):
    f = open(Name, 'w+')
    f.write(Header)
    f.write('\n')
    for i in range(NUM_ATOM):
        f.write('{ele}  {x}  {y}  {z}'.format(
            ele=LIST_ELEMENT[i], x=Geom[0, i * 3], y=Geom[0, i * 3 + 1], z=Geom[0, i * 3 + 2]))
        f.write('\n')
    f.write('*\n')
    f.write(Tail)
    f.write('\n')
    f.close()


def writeXYZ(Geom, Name):
    f = open(Name, 'w+')
    f.write(NUM_ATOM)
    f.write('\n\n')
    for i in range(NUM_ATOM):
        f.write('{ele}  {x}  {y}  {z}'.format(
            ele=LIST_ELEMENT[i], x=Geom[0, i * 3], y=Geom[0, i * 3 + 1], z=Geom[0, i * 3 + 2]))
        f.write('\n')
    f.write('\n')
    f.close()


def runEachStep(Geom, NStep, Header_A, Header_B, Tail1, Tail2):
    if PROG == 'gaussian':
        writeGjf(Geom, Header_A, Tail1, f'JOBS/{NStep}_A.gjf')
        writeGjf(Geom, Header_B, Tail2, f'JOBS/{NStep}_B.gjf')
        os.system(f'{PROG_COMM} JOBS/{NStep}_B.gjf')
        os.system(f'{PROG_COMM} JOBS/{NStep}_A.gjf')
    elif PROG == 'orca':
        writeORCA(Geom, Header_A, Tail1, f'JOBS/{NStep}_A.inp')
        writeORCA(Geom, Header_B, Tail2, f'JOBS/{NStep}_B.inp')
        os.system(f'{PROG_COMM} JOBS/{NStep}_B.inp > JOBS/{NStep}_B.log')
        os.system(f'{PROG_COMM} JOBS/{NStep}_A.inp > JOBS/{NStep}_A.log')
        if DELETE_GBW:
            os.system(f'rm -rf JOBS/{NStep}_B.gbw')
            os.system(f'rm -rf JOBS/{NStep}_A.gbw')
        else:
            os.system(f'mv JOBS/{NStep}_B.gbw  JOBS/b.gbw')
            os.system(f'mv JOBS/{NStep}_A.gbw  JOBS/a.gbw')
    elif PROG == 'xtb':
        writeXYZ(Geom, Header_A, Tail1, f'JOBS/{NStep}_A.xyz')
        writeXYZ(Geom, Header_B, Tail2, f'JOBS/{NStep}_B.xyz')
        os.system(f'{PROG_COMM} JOBS/{NStep}_B.xyz')
        os.system(f'{PROG_COMM} JOBS/{NStep}_A.xyz')

def getG(NStep):
    _, f1, E1 = readForceAndGeom(f'JOBS/{NStep}_A.log')
    _, f2, E2 = readForceAndGeom(f'JOBS/{NStep}_B.log')
    f1 = -numpy.array(f1)
    f2 = -numpy.array(f2)
    xVec = (f1 - f2)  # Please check this sign here
    xVecNorm = xVec / numpy.linalg.norm(xVec)
    # In Chem Phys Lett (1994) 223: 269 xVecNorm is used, although in Harvey's
    # work xVec is used herein
    fVec = (E1 - E2) * xVecNorm
    # fVec = (E1-E2)*xVec #if use this, f will be rather small, resulting in
    # inefficient optimization with repsect to dE
    gVec = (f1 - numpy.dot(xVecNorm, f1) * xVecNorm)  # *0.001 may help BGFS
    return fVec + gVec, E1, E2  # return an python array
    # return f1, E1, E2 # for debugging: to output force for the 1st state

# return f1, E1, E2 # for debugging: to output force for the 1st state

def propagationBFGS(X0, Gk, Bk):
    global LAMBDAS
    # rho = 0.01 #set rho = 0.01 for minimum opt, and 15 for MECP opt
    rho = 15
    dk = -numpy.linalg.solve(Bk, Gk)
    if numpy.linalg.norm(dk) > 0.1:
        dk = dk * 0.1 / numpy.linalg.norm(dk)
    XNew = X0 + rho * dk
    if len(LAMBDAS) != 0:
        LAMBDAS = XNew[0, -len(LAMBDAS):].tolist()[0]
    return XNew

def propagationGDIIS(Xs, Gs, Hs):
	# Produce a new geometry based on the GDIIS algorith, see https://manual.q-chem.com/5.3/A1.S7.html
    global LAMBDAS
    dimension = len(Xs)
    if len(Hs) != len(Xs):
        raise Exception('Runtime exception: H and X dimensions are different.')
    EMat = numpy.mat(numpy.zeros(shape=(dimension, NUM_ATOM * 3 + len(LAMBDAS))))
    H_mean = Hs[0]
    for i in range(1, dimension):
        H_mean += Hs[i]
    H_mean /= dimension
    for i in range(dimension):
        # print(Hs[i])
        # EMat[i] = (Hs[i].I * numpy.mat(Gs[i]).T).flatten()
        EMat[i] = (H_mean.I * numpy.mat(Gs[i]).T).flatten()
        # EMat[i] = (Gs[i].T).flatten()
    onesBlock = numpy.mat(numpy.ones(dimension))
    # construct the Bmat: B 1|1 0, Bij = ei*ej
    BMat = numpy.block([[EMat * EMat.T, onesBlock.T],
                        [onesBlock, numpy.mat([0])]])
    y = numpy.append(numpy.zeros(dimension), 1)
    c = numpy.linalg.solve(BMat, y)
    c = numpy.delete(c, -1)  # delete the last element, lambda, in c vector
    XNew_prime = numpy.mat(numpy.zeros(NUM_ATOM * 3 + len(LAMBDAS)))
    GNew_prime = numpy.mat(numpy.zeros(NUM_ATOM * 3 + len(LAMBDAS)))
    HNew_prime = numpy.mat(numpy.zeros(shape=(NUM_ATOM * 3 + len(LAMBDAS), NUM_ATOM * 3 + len(LAMBDAS))))
    for i in range(dimension):
        XNew_prime += Xs[i] * c[i]
        GNew_prime += Gs[i] * c[i]
        HNew_prime += Hs[i] * c[i]
    # XNew = XNew_prime - GNew_prime
    XNew = XNew_prime - (H_mean.I * numpy.mat(GNew_prime).T).flatten()
    # XNew = XNew_prime - (HNew_prime.I * numpy.mat(GNew_prime).T).flatten()
    LAMBDAS = XNew[0, NUM_ATOM * 3 - len(LAMBDAS) + 1:].tolist()[0]
    return XNew


def propagationGEDIIS(Xs, Gs, Hs,
                      Es):  # Produce a new geometry based on the GEDIIS algorith (J. Chem. Theory Comput. 2006, 2, 835-839)
    # Note that the energy to be minimized should not be E1 or E2, but produced from the pritimive function of G (not implemented yet)
    dimension = len(Xs)
    if len(Hs) != len(Xs):
        raise Exception('Runtime exception: H and X numbers are different.')
    EMat = numpy.mat(numpy.zeros(shape=(dimension, dimension)))
    for i in range(dimension):
        EMat[i, i] = 0
        for j in range(i + 1, dimension):
            EMat[i, j] = -1 * numpy.dot(Gs[i] - Gs[j],
                                        numpy.array((Xs[i] - Xs[j])).flatten())  # Gs is force rather than gradient
            EMat[j, i] = EMat[i, j]
    onesBlock = numpy.mat(numpy.ones(dimension))
    BMat = numpy.block([[EMat, onesBlock.T], \
                        [onesBlock, numpy.mat([0])]])
    y = numpy.append(-1 * numpy.array(Es), 1)
    c = numpy.linalg.solve(BMat, y)
    c = numpy.delete(c, -1)  # delete the last element, lambda, in c vector
    XNew_prime = numpy.mat(numpy.zeros(NUM_ATOM * 3))
    GNew_prime = numpy.mat(numpy.zeros(NUM_ATOM * 3))
    for i in range(dimension):
        XNew_prime += Xs[i] * c[i]
        GNew_prime += Gs[i] * c[i]
    XNew = XNew_prime + GNew_prime
    return XNew


def MaxStep(X0, XNew, factor=1):
    dX = XNew - X0
    dX *= factor
    stepsize_norm = numpy.linalg.norm(dX)
    if stepsize_norm > MAX_STEP_SIZE:
        print(f'current stepsize: {stepsize_norm} is reduced to max_size')
        dX = dX * MAX_STEP_SIZE / numpy.linalg.norm(dX)
    return X0 + dX


def HessianUpdator(Bk, yk, sk):
    # 10.1002/wcms.34
    Bk_BFGS = Bk + float(1.0 / (sk * yk.T)) * (
                float(1 + (yk * Bk * yk.T) / (sk * yk.T)) * (sk.T * sk) - (sk.T * yk * Bk + Bk * yk.T * sk))
    Bk_SR1 = Bk + (yk.T - Bk * sk.T) * (yk.T - Bk * sk.T).T / float((yk.T - Bk * sk.T).T * sk.T)  # slow convergence
    Bk_PSB = Bk + ((yk.T - Bk * sk.T) * sk + sk.T * (yk.T - Bk * sk.T).T) / (sk * sk.T) \
             - (float(sk * (yk.T - Bk * sk.T)) * sk.T * sk) / float(sk * sk.T) ** 2

    phi = float((yk.T - Bk * sk.T).T * sk.T) ** 2 / \
          float((yk.T - Bk * sk.T).T * (yk.T - Bk * sk.T)) / float(sk * sk.T)
    Bk_Bofill = phi * Bk_SR1 + (1 - phi) * Bk_PSB  # slow convergence
    Bk = Bk_PSB
    # Bk = Bk_BFGS
    # Bk = Bk_Bofill
    eigval, eigvec = numpy.linalg.eig(Bk)
    eigval = numpy.real(eigval)
    numimag = 0
    # Bk_BFGS is always positive finite, while Bk_PBS is not.
    for i in eigval:
        if i < 0:
            # print(i)
            numimag += 1
    # if numimag > 0:
    #	print('Negative eigenvals found for Bk_PSB, Bk_BFGS used instead')
    #	Bk = Bk_BFGS
    return Bk


def BFGS(X0, G0, B0, nstep):
    print(f'\nNow Entering BFGS Step {nstep}')
    XNew = propagationBFGS(X0, G0, B0)
    XNew = MaxStep(X0, XNew)
    runEachStep(XNew, nstep + 1, HEADER_A, HEADER_B, TAIL1, TAIL2)
    sk = XNew - X0
    GNew, E1, E2 = getG(nstep + 1)
    yk = numpy.mat(GNew - G0)
    Bk = HessianUpdator(B0, yk, sk)
    # Bk_BFGS = Bk - (Bk * sk.T * sk * Bk) / float(sk * Bk * sk.T) + (yk.T * yk) / float(yk * sk.T)
    return [XNew, GNew, Bk, E1, E2]


def GDIIS(Xs, Gs, Bs, nstep, Es, flag='gdiis'):
    print(f'\nNow Entering GDIIS Step {nstep}')
    XNew = propagationGDIIS(Xs, Gs, Bs)
    if flag == 'gediis':
        XNew2 = GEDIISpropagation(Xs, Gs, Bs, Es)
        XNew = XNew * 0.5 + XNew2 * 0.5
    factor = 1
    if numpy.linalg.norm(Gs) < THRESH_RMS_G * 10:
        factor = REDUCED_FACTOR
    XNew = MaxStep(Xs[-1], XNew, factor)
    runEachStep(XNew, nstep + 1, HEADER_A, HEADER_B, TAIL1, TAIL2)
    sk = XNew - Xs[-1]
    GNew, E1, E2 = getG(nstep + 1)
    yk = numpy.mat(GNew - Gs[-1])
    Bk = HessianUpdator(Bs[-1], yk, sk)
    return [XNew, GNew, Bk, E1, E2]


def isConverged(E1, E2, X0, X1, G1):
    dE = E2 - E1
    rms = numpy.linalg.norm(X0 - X1) / math.sqrt(NUM_ATOM * 3)
    maxDis = 0
    for i in range(0, NUM_ATOM * 3, 3):
        dist = math.sqrt((X0[0,
                             i] - X1[0,
                                     i])**2 + (X0[0,
                                                  i + 1] - X1[0,
                                                              i + 1])**2 + (X0[0,
                                                                               i + 2] - X1[0,
                                                                                           i + 2])**2)
        if dist > maxDis:
            maxDis = dist
    rmsG = numpy.linalg.norm(G1) / math.sqrt(NUM_ATOM * 3)
    maxG = 0
    for i in range(0, NUM_ATOM * 3, 3):
        g = abs(G1[i])
        if g > maxG:
            maxG = g
    dE_isConv, rms_isConv, max_dis_isConv, max_g_isConv, rms_g_isConv = [
        'NO', 'NO', 'NO', 'NO', 'NO']
    conv_flag = 0
    if abs(rms) < THRESH_RMS:
        rms_isConv = 'YES'
        conv_flag += 1
    if abs(dE) < THRESH_dE:
        dE_isConv = 'YES'
        conv_flag += 1
    if abs(maxDis) < THRESH_MAX_DIS:
        max_dis_isConv = 'YES'
        conv_flag += 1
    if abs(maxG) < THRESH_MAX_G:
        max_g_isConv = 'YES'
        conv_flag += 1
    if abs(rmsG) < THRESH_RMS_G:
        rms_g_isConv = 'YES'
        conv_flag += 1
    print(f'E1 = {E1}\nE2 = {E2}')
    print(
        f'deltaE                    {dE:5f}     {THRESH_dE:5f}     {dE_isConv}')
    print(
        f'RMS Gradient              {rmsG:5f}     {THRESH_RMS_G:5f}     {rms_g_isConv}')
    print(
        f'Maximium Gradient         {maxG:5f}     {THRESH_MAX_G:5f}     {max_g_isConv}')
    print(
        f'RMS Displacement          {rms:5f}     {THRESH_RMS:5f}     {rms_isConv}')
    print(
        f'Maximium Displacement     {maxDis:5f}     {THRESH_MAX_DIS:5f}     {max_dis_isConv}')
    if conv_flag == 5:
        return True
    else:
        return False

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
        if isConverged(E1, E2, X0, X1, G1):
            return [n_step, X1]
        if n_step > MAX_STEPS:
            raise Exception('Maximium number of steps exceeded.')
        X0 = X1
        G0 = G1
        B0 = B1


def main():
    global HEADER_A
    global HEADER_B
    global DELETE_GBW
    global GEOM
    global CONSTRAINTS
    print('****KST48 PROGRAM: a powerful tool for MECP locating****')
    print('****By Yumiao Ma, BSJ Institute, 2022/01/09****\n')
    _, inp = sys.argv
    os.system('mkdir JOBS')
    nprocs, mem, charge, mult1, mult2, method, runMode = inputParser(
        sys.argv[1])
    header_1, header_2 = buildInitJob(
        nprocs, mem, charge, mult1, mult2, method, PROG, runMode)
    print(f'Note: This program is now running on {runMode} mode')
    if runMode == 'noread':
        DELETE_GBW = True

    # the preparation phase: run single points or stability calcs to obtain the wavefunction
    print('****Initialization: running the first single point calculations according to the mode****')
    if runMode != 'read':
        runPrePoint(header_1, header_2, PROG, runMode)
        if runMode == 'stable' or runMode == 'inter_read':
            if PROG == 'orca':
                print('In an RHF calculation in ORCA, it will not restart automatically if an unstability is found. \
	Remember to write UKS when you are handling singlet state! ')
                print('RI is unsupported for stability analysis. It is recommended to MANUALLY obtain the correct wavefunction, \
	and then use the read model of KST48, rather than the stable mode, in order to use RI.')
                os.system(f'cp JOBS/pre_A.gbw JOBS/a.gbw')
                os.system(f'cp JOBS/pre_B.gbw JOBS/b.gbw')
            runMode= 'read'

    # the main loop
    scan_step = 1
    if runMode in ['read', 'normal', 'noread']:
        method = modifyMETHOD(PROG, method, runMode)
        print('****Initialization OK, now entering main loop****')
        print('****Before that, please check the keywords list****')
        HEADER_A, HEADER_B = buildInitJob(
            nprocs, mem, charge, mult1, mult2, method, PROG, runMode, Td1=TD1, Td2=TD2)
        print(f'Header A:\n {HEADER_A}')
        print(f'Header B:\n {HEADER_B}')
        print('****If everything is OK, then go to the main loop****\n')
        if SCANS == []:  # normal opt without scan
            runEachStep(GEOM, 0, HEADER_A, HEADER_B, TAIL1, TAIL2)
            # run JOBS/0_A and 0_B, to obtain the first gradient
            GEOM, force, E = readForceAndGeom('JOBS/0_A.log')
            conv_step, GEOM = runOpt(GEOM, flag='hybrid')
            print('****Congrats! MECP has converged****')
            if PROG == 'gaussian':
                os.system(f'cp JOBS/{conv_step}_A.gjf .')
                os.system(f'cp JOBS/{conv_step}_B.gjf .')  # .chk file is already at current folder
            elif PROG == 'orca':
                os.system(f'cp JOBS/{conv_step}_A.inp .')
                os.system(f'cp JOBS/{conv_step}_B.inp .')
                os.system(f'cp JOBS/a.gbw .; cp JOBS/b.gbw .')
            os.system(f'cp JOBS/{conv_step}_A.log .')
            os.system(f'cp JOBS/{conv_step}_B.log .')
        # Run PES Scan
        else:  # SCANS: [ [[r,A,B], [start, num, size] ], ... ]
            initial_cons_num = len(CONSTRAINTS)
            scanVars = [[], []]
            scan_para1 = [float(i) for i in SCANS[0][1]]
            if len(SCANS) == 1:
                SCANS.append([[], [0, 0, 0]])  # add one invalid scan to the 2nd dimension if only 1-D scan is requested
            scan_para2 = [float(i) for i in SCANS[1][1]]
            for i in range(int(scan_para1[1])):
                scanVars[0].append(scan_para1[0] + i * scan_para1[2])
            for i in range(int(scan_para2[1])):
                scanVars[1].append(scan_para2[0] + i * scan_para2[2])
            if scanVars[1] == []:
                scanVars[1].append(-1)  # make sure it contains something  for the following loop to run
            for i in scanVars[0]:
                new_const = SCANS[0][0][:]
                new_const.append(i)
                CONSTRAINTS.append(new_const)
                for j in scanVars[1]:
                    new_const = SCANS[1][0][:]
                    new_const.append(j)
                    CONSTRAINTS.append(new_const)
                    print(f'constraints after pop and append')
                    print(CONSTRAINTS)
                    print(f'****Scan Cycle {i:4f}_{j:4f}****')
                    runEachStep(GEOM, 0, HEADER_A, HEADER_B, TAIL1, TAIL2)
                    GEOM, force, E = readForceAndGeom('JOBS/0_A.log')
                    conv_step, GEOM = runOpt(GEOM, flag='hybrid')
                    print('****Congrats! MECP has converged****\n')
                    if PROG == 'gaussian':
                        os.system(f'cp JOBS/{conv_step}_A.gjf {i:4f}_{j:4f}.gjf')
                        os.system(
                            f'cp JOBS/{conv_step}_B.gjf {i:4f}_{j:4f}.gjf')  # .chk file is already at current folder
                    elif PROG == 'orca':
                        os.system(f'cp JOBS/{conv_step}_A.inp {i:4f}_{j:4f}.inp')
                        os.system(f'cp JOBS/{conv_step}_B.inp {i:4f}_{j:4f}.inp')
                    # os.system(f'cp JOBS/a.gbw .; cp JOBS/b.gbw .')
                    os.system(f'cp JOBS/{conv_step}_A.log {i:4f}_{j:4f}.log')
                    os.system(f'cp JOBS/{conv_step}_B.log {i:4f}_{j:4f}.log')
                    CONSTRAINTS.pop(-1)
                CONSTRAINTS.pop(-1)


main()
