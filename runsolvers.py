#!/usr/bin/env python

import sys,shutil,os,re
import subprocess
from glob import glob
from time import time
from ipdb import set_trace
from phd import soliton
from phd import postprocess as pproc

local_computer = os.uname()[1].split('.')[0]  #gives back 'dagda', 'medea' and so on

compNameShort={'omeyocan':'om','dagda':'da','medea':'me','shamash':'sh'}
solverShort = {'OpenFOAM': 'OF', 'Truchas':'Tr', 'TruchasEnSight':'Tr', 'Thetis':'Th','Gerris':'Ge'}
inputExt = {'Gerris': '.gfs','OpenFOAM': '.foam','Thetis': '.don','Truchas': '.inp', 'TruchasEnSight': '.inp'}

def logInScreen(job,computer='',ending='',solver=None):
    """
    Changes the name of the currently used screen, for information purposes.
    """
    solver = solver or pproc.determineSolver()
    solver = solverShort[solver]
    screenCommand = 'screen -X title \"{}\"'.format(' '.join([solver,job,compNameShort[computer],ending]))
    if os.getenv('STY'): #env var giving the name of the screen session
        subprocess.call(screenCommand,shell=True)
        
   
def updateCaseDir(caseDir):
    """
    Removes the -inProgress ending of the case directory name, after the simulation is finished.
    """
    os.chdir('../')
    newCaseDir = caseDir.replace('-inProgress','')
    if os.path.exists(newCaseDir):
        shutil.rmtree(newCaseDir)
        print 'Removed old directory',newCaseDir
    shutil.move(caseDir,newCaseDir)
    os.chdir(newCaseDir)
    return newCaseDir


def copyRemote(path,target, remove=True,force_remove=False,verbose=False):
    """
    Copy the given path to the same path on the target computer
    """
    if local_computer == target:
        sys.exit('Error: Target computer is the same as local computer')
        
    if type(path) == list:
        for p in path:
            copyRemote(p,target,remove=True,force_remove=force_remove,verbose=verbose)
    else:
        path = os.path.abspath(path)
        dest = ':'.join([target,path])

        subprocess.call('ssh {} mkdir -p {}'.format(target,os.path.split(path)[0] ),shell=True)
        if not subprocess.call('ssh {} test -e {}'.format(target,path), shell=True) and remove:
            if force_remove or raw_input('Remove {} ?\n'.format(dest)).lower() in ['y','yes']:
                print 'REMOVING',dest
                subprocess.call('ssh {} rm -r {}'.format(target,path), shell=True)
        print 'Copying {} from {} to {}'.format(path,local_computer,target)
        quiet = '' if verbose else 'q'
        subprocess.call('{} {} {}'.format('scp -r'+quiet,path,dest), shell=True)
        #subprocess.call('rsync -a'+quiet+'e ssh '+path+' '+dest+' --progress',shell=True)  #alternative, not really better here
        return 0

inputPatterns = {  'Gerris': ['*.gfs'],
                'OpenFOAM': ['*.foam','0','constant','system'],
                'Thetis': ['*.don','caract.par','defaut.don','licence.dat','thetis.dat'],
                'Truchas': ['*.inp',] }

def copyInputRemote(source,target):
    """
    Copies the input found below the source directory to the target maching.
    Retains the directory structure.
    """
    source = os.path.abspath(source)
    solver = pproc.determineSolver(source)
    for i in range(10):
        if glob(os.path.join(source,i*'*/','*'+inputExt[solver])):
            for pattern in inputPatterns[solver]:
                for path in glob(os.path.join(source,i*'*/',pattern)):
                    copyRemote(path, target,force_remove=True)
            break

def removeInput(path):
    """
    Removes input files in the given directory
    """
    path = os.path.abspath(path)
    solver = pproc.determineSolver(path)
    for pattern in inputPatterns[solver]:
        for obj in glob(os.path.join(path,pattern)):
            subprocess.call('rm -r '+obj,shell=True)

def backupInput(source,dest):
    """
    Copies all input files found under the source directory, preserving the directory structure (files won't get overwritten).
    """
    cases = ['-1','-2','-3','-4','-5','e1','e2','e3','e4','e5']  #e1 from case1 etc.
    for root, dirs, files in os.walk(source):
        #print root, dirs, files
        for file in files:
            if (os.path.splitext(file)[1] in ['.inp','.gfs','.foam','.don'] or file in ['caract.par','licence.dat','thetis.dat']) and root[-2:] not in cases:
                print 'Copying',os.path.join(root,file),os.path.join(dest,root,file)
                subprocess.call('mkdir -p '+os.path.join(dest,root),shell=True)
                shutil.copy(os.path.join(root,file),os.path.join(dest,root,file))
        for dir in dirs:
            if dir in ['0','constant','system'] and root[-2:] not in cases:
                print 'Copying',os.path.join(root,dir),os.path.join(dest,root,dir)
                subprocess.call('mkdir -p '+os.path.join(dest,root),shell=True)
                shutil.copytree(os.path.join(root,dir),os.path.join(dest,root,dir))


def runGerrisCase(testFile, testDir, testName, case):
    os.chdir(testDir)
    caseDir = os.path.join(testDir,'case'+str(case)+'-inProgress')
    caseName = testName+'-'+str(case)
    
    if os.path.exists(caseDir):
        shutil.rmtree(caseDir)
    os.mkdir(caseDir)
    
    mesh = open(testFile).readlines()
    for i,l in enumerate(mesh):
        l = l.split()
        if len(l)==2 and l[0] in ['GfsRefine', 'Refine']  and l[1].isdigit():
            print 'Changed Refine from',l[1],'to',int(l[1])+case-1
            l[1] = str(int(l[1])+case-1)+'\n'
            mesh[i] = ' '.join(l)
        if len(l)==3 and l[0] == 'Define' and l[1]== 'REFLEVEL':
            print 'Changed REFLEVEL from',l[2],'to',int(l[2])+case-1
            l[2] = str(int(l[2])+case-1)+'\n'
            mesh[i] = ' '.join(l)
        if len(l)==3 and l[0] == 'Define' and l[1] == 'NAME':
            l[2] = caseName+'\n'
            mesh[i] = ' '.join(l)
    
    open(os.path.join(caseDir,caseName+'.gfs'),'w').write(''.join(mesh))
    for f in glob('*.gts') + glob('*.cgd'): #copy additional input files
        shutil.copy(f,caseDir)
    
    os.chdir(caseDir)
    
    if '.cgd' in open(os.path.join(caseDir,caseName+'.gfs')).read() and ('propagation' in testDir or 'runup' in testDir):
        soliton.initialize()
    
    t0=time()
    subprocess.call('gerris2D -m '+caseName+'.gfs 2>&1 | tee '+caseName+'.out',shell=True)  #2>&1 streams stderr to stdout, so that tee can catch it
    open(caseName+'_timing.dat','w').write(str(time()-t0)+'\n')
    caseDir = updateCaseDir(caseDir)
    return caseDir
    
def chooseFoamSolver(caseDir='./'):
    """
    Chooses a proper OpenFOAM solver to be used, according to the path.
    Returns a name of the chosen solver.
    """
    caseDir = os.path.abspath(caseDir)
    if glob('constant/dynamicMeshDict'):
        if glob('./0/alpha1'):
            solver = 'interDyMFoam'
        else:
            solver = 'pimpleDyMFoam'
    else:
        if glob('./0/alpha1'):
            solver = 'interFoam'
        elif 'pimple' in caseDir:
            #solver = 'pimpleInviscidFoam'
            solver = 'pimpleFoam'
        elif 'simple' in caseDir:
            solver = 'simpleFoam'
        else:
            solver = 'icoFoam'
    return solver

def runFoamCase(testFile, testDir, testName, case,fenton=0):
    inpext = '.foam'
    os.chdir(testDir)
    caseName = testName+'-'+str(case)
    caseDir = os.path.join(testDir,caseName)
    
    if os.path.exists(caseDir):
        shutil.rmtree(caseDir)
        print 'Removed old directory',caseDir
    os.mkdir(caseDir)
    for f in ['0','system','constant','initialDisplacement','generateSolid.py']:
        if os.path.exists(f): subprocess.call('cp -r {0} {1}'.format(f,caseName),shell=True)
    shutil.copy(testFile,os.path.join(caseName,os.path.splitext(testName)[0]+'-'+str(case)+inpext))
    mesh = open('constant/polyMesh/blockMeshDict').readlines()
    for i,l in enumerate(mesh):
        if l.split() and l.split()[0] == 'hex':
            params = re.split('[()]+', l)[3]
            paramsInt = [str(int(p)*2**(case-1)) if int(p)!=1 else p for p in params.split()] #the structure 'a if b else c' is called a 'ternary operator'
            mesh[i] = mesh[i].replace(params,' '.join(paramsInt))
            print 'In mesh file, replaced line:'
            print l,
            print 'with:'
            print mesh[i],
    open(os.path.join(caseName,'constant/polyMesh/blockMeshDict'),'w').write(''.join(mesh))
    
    foamDotFile = subprocess.check_output('. $HOME/.bashrc; echo $foamDotFile',shell=True).splitlines()[-1]

    if not foamDotFile: sys.exit('foamDotFile not defined in .bashrc')
    if not os.path.exists(foamDotFile): sys.exit('foamDotFile does not exist: {}'.format(foamDotFile))
    os.chdir(caseName)
    for f in glob('constant/triSurface/*.stl'):
        subprocess.call('. '+foamDotFile+'; surfaceFeatureExtract -includedAngle 150 '+f+' features',shell=True)
    subprocess.call('cp ./0/alpha1.org ./0/alpha1; rm constant/polyMesh/[c-f,h-z]* -f; rm 0.* -rf; rm [1-9]* -rf',shell=True)
    subprocess.call('. '+foamDotFile+'; blockMesh',shell=True)
    
    if os.path.exists('./system/snappyHexMeshDict'): subprocess.call('. '+foamDotFile+'; snappyHexMesh -overwrite',shell=True)
    if os.path.exists('./generateSolid.py'): subprocess.call('. '+foamDotFile+'; ./generateSolid.py',shell=True)
    if os.path.exists('./initialDisplacement'): subprocess.call('. '+foamDotFile+'; cd initialDisplacement; ./run.sh; cd ..',shell=True)
    
    if 'runup' in testDir or 'propagation' in testDir:
        if fenton:
            subprocess.call('. '+foamDotFile+'; initializeFoamFenton',shell=True)
        else:
            soliton.initialize()
    elif glob('system/setFieldsDict'):
        subprocess.call('. '+foamDotFile+'; setFields',shell=True)
    
    t0=time()
    if subprocess.call('. '+foamDotFile+'; ' + chooseFoamSolver(),shell=True): sys.exit('runsolvers.py: Exiting because of the error of the solver')
    open(caseName+'_timing.dat','w').write(str(time()-t0)+'\n')

    if 'cylinder' in caseDir:
        subprocess.call('. '+foamDotFile+'; foamToVTK -excludePatches \'( \".*\" )\' -fields \'( alpha1 U p )\' ',shell=True)
    else:
        subprocess.call('. '+foamDotFile+'; foamToVTK -excludePatches \'( \".*\" )\' -fields \'( alpha1 U )\' ',shell=True) #-ascii #-noPointValues
    subprocess.call('mkdir timesteps; mv [1-9]* timesteps; mv 0.* timesteps; mv ./VTK/* ./; rm VTK -r',shell=True)
    return caseDir

def runThetisCase(testFile, testDir, testName, case,fenton=0):
    inpext = '.don'
    os.chdir(testDir)
    caseDir = os.path.join(testDir,'case'+str(case)+'-inProgress')
    caseName = testName+'-'+str(case)
    
    if os.path.exists(caseDir):
        shutil.rmtree(caseDir)
        print 'Removed old directory',caseDir
    os.mkdir(caseDir)
    
    for file in ['caract.par','defaut.don','licence.dat']:
        shutil.copy(file,caseDir)
    
    mesh = open(testFile).readlines()
    for i,l in enumerate(mesh):
        l = re.split('[ ]+',l.strip())
        if len(l)==3 and l[0] == 'MAILLAGE':
            Nx = str(int(l[1])*2**(case-1))
            Ny = str(int(l[2])*2**(case-1))
            mesh[i] = 'MAILLAGE   '+Nx+' '+Ny+'\n'
        if len(l)==3 and l[0] == 'DIM_MAX':
            L = l[1].replace('D0','')
            H = str(float(l[2].replace('D0',''))+1)
        
    open(os.path.join(caseDir,caseName+inpext),'w').write(''.join(mesh))
    open(os.path.join(caseDir,'thetis.dat'),'w').write('REPERTOIRE_DON ./\nFICHIER_DONNEE {file}.don\nFICHIER_DEFAUT defaut.don\nFICHIER_CARACT caract.par\n'.format(file=caseName))
    
    os.chdir(caseDir)
    if not 'stillwater-flat' in testDir:
        if (not 'rotated' in testFile) and fenton:
            amp = pproc.determineAmplitude(testDir)
            subprocess.call('initialConditionThetis '+ str(amp) +' 15 '+ L +' '+ H +' '+ Nx +' '+ Ny,shell=True)
        else:
            soliton.initialize()
    t0=time()
    subprocess.call('thetis',shell=True) # | tee '+caseName+'.out'
    open(caseName+'_timing.dat','w').write(str(time()-t0)+'\n')
    caseDir = updateCaseDir(caseDir)
    return caseDir
    
def runTruchas200Case(testFile, testDir, testName, case):
    #runs the old version of Truchas, with the moving solid module
    os.chdir(testDir)
    caseDir = os.path.join(testDir,'case'+str(case)+'-inProgress')
    caseName = testName+'-'+str(case)
    
    if os.path.exists(caseDir):
        shutil.rmtree(caseDir)
        print 'Removed old directory',caseDir
    os.mkdir(caseDir)
    
    mesh = open(testFile).readlines()
    for i,l in enumerate(mesh):
        l = re.split('[ =,]+',l.strip())
        if len(l)==4 and l[0] == 'Ncell':
            print 'In mesh file, replaced line:'
            print mesh[i],
            for j in range(1,len(l)):
                if int(l[j]) != 1:
                    l[j] = str(int(l[j])*2**(case-1))
            mesh[i] = '    Ncell             = '+l[1]+', '+l[2]+', '+l[3]+'\n'
            print 'with:'
            print mesh[i]
    
    os.chdir(caseDir)
    open(caseName+'.inp','w').write(''.join(mesh)) #first located and computed in testDir, then moved to caseDir
    t0=time()
    subprocess.call('truchas-moving-2.0 '+caseName+'.inp',shell=True) # | tee '+caseName+'.out'
    
    open(caseName+'_timing.dat','w').write(str(time()-t0)+'\n')
    caseDir = updateCaseDir(caseDir)
    return caseDir

def parseTruchasOutput(dir='./'):
    print "Parsing results..."
    subprocess.call('cd '+dir+';TBrookParse > /dev/null;rm *.bin -f;cd '+os.getcwd(),shell=True)
 
def runTruchasCase(testFile, testDir, testName, case,fenton=0):

    os.chdir(testDir)
    caseDir = os.path.join(testDir,'case'+str(case)+'-inProgress')
    if os.path.exists(caseDir):
        shutil.rmtree(caseDir)
        print 'Removed old directory',caseDir
    os.mkdir(caseDir)
    
    caseName = testName+'-'+str(case)
    
    mesh = open(testFile).readlines()
    for i,l in enumerate(mesh):
        l = re.split('[ =,]+',l.strip())
        if len(l)==4 and l[0] == 'Ncell':
            print 'In mesh file, replaced line:'
            print mesh[i],
            for j in range(1,len(l)):
                if int(l[j]) != 1:
                    l[j] = str(int(l[j])*2**(case-1))
            mesh[i] = '    Ncell             = '+l[1]+', '+l[2]+', '+l[3]+'\n'
            print 'with:'
            print mesh[i]
    
    os.chdir(caseDir)
    open(caseName+'.inp','w').write(''.join(mesh)) #first located and computed in testDir, then moved to caseDir
    soliton.initialize()
    t0=time()
    subprocess.call('$SOLVERS/truchas-2.6.0/bin/t-linux.x86_64.g95.serial.opt-2.6.0 '+caseName+'.inp',shell=True) # | tee '+caseName+'.out'
    
    for f in glob(caseName+'_output/*'):
        shutil.move(f,caseDir)
    shutil.rmtree(caseName+'_output')
    
    open(caseName+'_timing.dat','w').write(str(time()-t0)+'\n')
    parseTruchasOutput()
    caseDir = updateCaseDir(caseDir)
    return caseDir

run={'Gerris':runGerrisCase, 'OpenFOAM':runFoamCase, 'Thetis':runThetisCase,'Truchas':runTruchasCase,'TruchasEnSight':runTruchas200Case}

def runCase(job,case,results_collector=''):
    """
    Runs a single simulation on a local machine.
    job is a string of the format 'path_to_the_case_file[!name_of_the_remote_compouter]'
    If the optional remote_computer is provided, the results are sent to the remote computer and removed from the local machine.
    """
    testFile = os.path.abspath(job)
    solver = pproc.determineSolver(testFile)
    if not os.path.exists(testFile) or os.path.splitext(testFile)[1]!=inputExt[solver]:
        sys.exit("Please provide a {} {} case file as an argument".format(solver,inputExt[solver]))
        
    testDir,testName = os.path.split( os.path.splitext(testFile)[0] )
    
    print 'Running {}, case {}, on {}'.format(testFile,case,local_computer)
    sys.stdout.flush()
    caseDir = run[solver](testFile, testDir, testName, case)

    if 'propagation' in testDir or 'runup' in testDir or 'stillwater' in testDir:
        pproc.postprocess(caseDir)
    elif 'moving' in testDir or 'cylinder' in testDir:
        pproc.postprocess1phase(caseDir)
    
    os.chdir(testDir)

    if results_collector:
        r = copyRemote(caseDir,results_collector,force_remove=True) #copying results from a single case dir

        #cleanup
        if not r:
            shutil.rmtree(caseDir)
        removeInput(testDir)
        if not os.listdir(testDir): #empty dir
            os.chdir(os.path.split(testDir)[0])
            shutil.rmtree(testDir)
    
def runSeries(job,cases=[1,2,3],results_collector=''):
    """
    The main public function of the module.
    'job' parameter is a string formatted as [remoteComputerName:]pathToInputFile
    The function runs a series of simulations for automatically refined resolutions.
    The jobs are sent in forms of single cases to remote servers, if the name of the computer is provided in the job name.
    The function should only be runned on the main computer collecting data.
    """
    remote_computer = ''
    if ':' in job:
        remote_computer,testFile = job.split(':')
    else:
        testFile = job
    testFile = os.path.abspath( testFile )
    testDir = os.path.split(testFile)[0]
    solver = pproc.determineSolver(testFile)
    testName = os.path.split( os.path.splitext(testFile)[0] )[1]
    for case in cases:
        logInScreen('{} {}/{}'.format(testName,case,len(cases)),computer=remote_computer or local_computer,solver=solver)
        if remote_computer:
            copyInputRemote(testDir,remote_computer)
            subprocess.call('ssh -t {} runsolvers.py {} {} --results_collector {}'.format(remote_computer,testFile,case,local_computer),shell=True)  #thanks to -t, output appears immediately
        else:
            runCase(testFile,case,results_collector=results_collector)

    if not results_collector: #do the plotting on the local machine, to take into account result of the previous simulations
        if 'runup' in testFile or 'propagation' in testFile:
            pproc.plotCollections(testDir)
        elif 'moving' in testFile or 'cylinder' in testFile:
            xsections.cylinder()
 

#remote computers at UiO: 'omeyocan','shamash','keto','kingu','nefele'

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Runs simulation on a local or remote machine')
    parser.add_argument('job', type=str, help='Path to the input file, proceeded with an optional [remote_computer:]')
    parser.add_argument('cases', type=int, default=[1,2,3], nargs='*',help='List of integers giving the refinement. Default is 1 2 3')
    parser.add_argument('--results_collector', type=str, help='Remote computer collecting the results')

    args = parser.parse_args()

    runSeries(args.job,args.cases,results_collector=args.results_collector)
        
