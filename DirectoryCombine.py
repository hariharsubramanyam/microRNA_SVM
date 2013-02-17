#!/usr/bin/env python

# Assumptions Made:
    # All input file names will end with .txt
    # Every input file name (excluding the .txt) will be of the form <number><A or B> (ex. 60A.txt, 12B.txt, 52A.txt)
    # There are N assays in the A-group whose name includes "U6" and N assays in the B-group whose name includes "U6"  
    # Minimum variant assays are chosen by finding assays which
        # Belong to group A
        # Have all numerical CT values
        # Minimize the chosen metric (--compMetric or -c)
    # We normalize between patients using the bottom 10% of least variant assays
    # Patients groups have been specified in the --groupfile (or -g) argument. Or, the information is available in the ./PatientGroups/patgroups.txt
    # We will run SVM only on "A" mirs that have all their CT values

import argparse     # for parsing command line arguments
import sys          # for running command line programs
import re           # for pattern matching in strings
import os           # for changing directories
import subprocess   # for running command line programs
import shlex        # for string splitting
from subprocess import *    # for using subprocess functions by writing myFunction() instead of subprocess.myFunction()
from os import listdir      # for listing directories
from os.path import isfile, join    # for checking if files exist and combining multiple strings into one 
from mx.Misc.Cache import DENOM     # for command line programs

def main():
    # Parse command line arguments (see README for description of arguments)
    aParser = argparse.ArgumentParser(description="Aggregate PCR files and optionally run SVM")
    aParser.add_argument("inDir",help="Directory containing PCR files")
    aParser.add_argument("-l","--libsvmdir", help="Directory containing libsvm (ex. user/stuff/libsvm-3.12)")
    aParser.add_argument("-g","--groupfile", help="File explaining how are patients broken into groups")
    aParser.add_argument("-t","--testingDir",help="Directory containing patients to be classified")
    aParser.add_argument("-c","--compMetric",help="Metric used to normalize between patients (0 = chi-squared (default), 1 = variance*mean)")
    
    argNamespace = aParser.parse_args()

    # map command line arguments to local variables
    inDir = argNamespace.inDir
    libsvmdir = argNamespace.libsvmdir
    groupfile = argNamespace.groupfile
    scriptDir = os.path.dirname(os.path.realpath(__file__))
    compMetric = argNamespace.compMetric
    testingDir = argNamespace.testingDir
    # Find the text files
    print "Identifying text files..."
    sys.stdout.flush() 
    txtfiles = [f for f in listdir(inDir) if isfile(join(inDir,f)) and (".txt" == f[-4:])]

    # Get array of A and B assays    
    print "Combining text file data..."
    sys.stdout.flush()
    (AAssays,BAssays) = filesToAssayArray(txtfiles,inDir)

    # Print out table
    print "Producing raw output..."
    sys.stdout.flush()
    fileRawOutput = open(scriptDir + "/Output/A_RawOutput.txt","w")
    fileRawOutput.write(assayOutputString((AAssays,BAssays),ABMerged=True))
    fileRawOutput.close()
    
    # Combine A and B assays
    print "Merging A and B assays..."
    sys.stdout.flush()
    (assays, adjuster) = combineABAssays(AAssays,BAssays)
    
    # Print out another table
    print "Producing combined output..."
    sys.stdout.flush()
    fileCombineOutput = open(scriptDir + "/Output/B_CombinedOutput.txt","w")
    fileCombineOutput.write(assayOutputString(assays,ABMerged=False))
    fileCombineOutput.close()
        
    
    # Find least variant assays according to the metric that the user inputted
    print "Finding least variant assays..."
    sys.stdout.flush()
    minChiSqAssays(assays)
    if int(compMetric) == 1:
        (minAssays,metricName) = minVarianceTimesMeanAssays(assays)
    else:
        (minAssays,metricName) = minChiSqAssays(assays)
    
    # Print least variant assays into table
    print "Producing least variant assay output with metric = " + metricName + "..."
    sys.stdout.flush()
    fileMinAssayOutput = open(scriptDir + "/Output/C_MinAssays.txt","w")
    fileMinAssayOutput.write(minAssayOutputString(minAssays,metricName))
    fileMinAssayOutput.close()
        
    # Reorganize the data from a list of assays to a list of patients
    print "Reorganizing data structures..."
    sys.stdout.flush()
    patients = assaysToPatients(assays,txtfiles)
    
    # Normalize patients
    print "Normalizing patients with respect to least variant assays..."
    sys.stdout.flush()
    normalizeMirNames = [minAssays[x].name for x in range(0,int(len(minAssays)*0.1))]
    for patient in patients:
        patient.normalizeWRTMirNames(normalizeMirNames)
    
    # Reorganize patients into assays so we can make output table
    print "Reorganizing data structures again..."
    sys.stdout.flush()
    assays = patientsToAssays(patients,assays)
    
    # Make output table
    print "Producing normalized output..."
    sys.stdout.flush()
    fileNormalizedOutput = open(scriptDir + "/Output/D_NormalizedOutput.txt","w")
    fileNormalizedOutput.write(assayOutputString(assays,ABMerged=False))
    fileNormalizedOutput.close()
    
    # If the user does not want to run an SVM, stop the program here
    if not libsvmdir:
        print "\nProgram finished running!\n"
        print "Output files have been written to " + scriptDir + "/" + "Output\n"
        print "Next time, if you want to run an SVM, run this script as\npython DirectoryCombine.py <inputdir> <libsvmdir>\nlibsvmdir is the path to the libsvm directory (ex. /user/home/libsvm-3.12). This is an optional parameter\n\n"
        return
    
    # Read group file and label patients based on their group
    print "\nLabeling patients for SVM..."
    sys.stdout.flush()
    groups = labelPatients(patients, groupfile)
    
    # Get all candidate mirs
    print "\nChoosing mirs from A group with all CTs..."
    sys.stdout.flush()
    lvMirNames = [x.name for x in assays if x.group == "A"]
    
    
    for x in xrange(0,len(groups)-1):
        for y in xrange(x+1,len(groups)):
            # Generate testing data string
            outString = ""
            for patient in patients:
                if not(patient.label is None):
                    if patient.label == groups[x] or patient.label == groups[y]:
                        outString += patient.toSVMFormatStringForMirs(lvMirNames) + "\n" 
            outString = outString.strip()
            # Write data string to file
            currDir = os.getcwd()
            os.chdir(libsvmdir)
            dataFile = open(libsvmdir+"/PatientData","w")
            dataFile.write(outString)
            dataFile.close()
            if (not (testingDir is None)):
                testingFile = open(libsvmdir+"/PatientDataTesting","w")
                testingFile.write(outString)
                testingFile.close()
                
            os.chdir(libsvmdir+"/tools")
            
            if (not (testingDir is None)):
                os.system("python fselect.py ../PatientData ../PatientDataTesting")
                copyFileToPath(libsvmdir + "/tools/PatientData.model", scriptDir + "/Output/Model" + str(groups[x]) + "_" + str(groups[y]) + ".model")
                copyFileToPath(libsvmdir + "/PatientDataTesting",scriptDir + "/Output/Testing" + str(groups[x]) + "_" + str(groups[y]))
                copyFileToPath(libsvmdir + "/PatientDataTesting", scriptDir + "/Output/Prediction" + str(groups[x]) + "_" + str(groups[y]) + ".pred")
            else:
                os.system("python fselect.py ../PatientData")
                
            copyFileToPath(libsvmdir + "/PatientData",scriptDir + "/Output/Training" + str(groups[x]) + "_" + str(groups[y]))
            copyFileToPath(libsvmdir + "/tools/PatientData.select", scriptDir + "/Output/Select" + str(groups[x]) + "_" + str(groups[y]) + ".select")
            os.chdir(currDir)
            
                    
    print "Done"

def assayOutputString(assays,ABMerged=False):
    # Produces a tabular representation of an array of assays
    # ABMerged = True means the assays have been normalized
    # assays is a tuple of the form (AAssays, BAssays)
    outString = "DirectoryCombine.py\n\n"
    
    # If the assays are merged, treat all assays as A assays
    if ABMerged:
        (AAssays,BAssays) = assays
    else:
        (AAssays,BAssays) = (assays,[])

    # Get experiment names for A and B assays
    AExps = [AAssays[0].cts[x].experiment for x in xrange(0,len(AAssays[0].cts))]
    if ABMerged:
        BExps = [BAssays[0].cts[x].experiment for x in xrange(0,len(BAssays[0].cts))]
    else:
        BExps = []
    
    # Which has more experiments, A or B?
    exps = AExps
    if len(BExps) > len(AExps):
        exps = BExps
    
    # Sort experiments in ascending order
    exps.sort()
    
    # Create column headers
    outString += ("%-10s%-30s") % ("Well","Target_Name")
    for exp in exps:
        outString += ("%-30s")%(exp)
        
     
    outString += "\n\n"
    
    # For all the AAssays
    for x in xrange(0,len(AAssays)):
        # Print assay well and name
        outString += ("%-10s%-30s") % (x+1,AAssays[x].name + "A")
        # Print all CT values
        for exp in exps:
            ct = AAssays[x].getCTForExperiment(exp)
            if ct:  # If the CT does not exist, print an X. If it is NOT "Undetermined", round it to 3 decimal places
                outString += ("%-30s") % ("Undetermined" if ct.value == "Undetermined" else round(float(ct.value),3))
            else:
                outString += ("%-30s") % ("X")
        outString += "\n"
    
    # Repeat above steps for BAssays
    # NOTE: If ABMerged = True, BAssays = [], so this code would not execute at all
    for x in xrange(0,len(BAssays)):
        outString += ("%-10s%-30s") % (x+1,BAssays[x].name + "B")
        for exp in exps:
            ct = BAssays[x].getCTForExperiment(exp)
            if ct:
                outString += ("%-30s") % ("Undetermined" if ct.value == "Undetermined" else round(float(ct.value),3))
            else:
                outString += ("%-30s") % ("X")
        outString += "\n"
    
    # Return the table
    return outString

def assaysToPatients(assays,txtfiles):
    # Converts the assays from txtfiles into Patient objects
    
    # Get patient numbers
    expnames = [x.replace(".txt","").replace("A","").replace("B","") for x in txtfiles]
    expnames = list(set(expnames))
    
    # Create patients from number 
    patients = [Patient(x) for x in expnames]
    
    # Add Mir CTs to patients
    for assay in assays:
        for ct in assay.cts:
            for patient in patients:
                if patient.name == ct.experiment:
                    patient.mirs.append(Mir(assay.name,ct.value,assay.group))
    return patients

def combineABAssays(AAssays,BAssays,refExp = "01"):
    # Merge A and B assays using the U6 reference assay
    # Get U6 assay CTs from A and B assays
    AU6 = [float(x.getCTForExperiment("01").value) for x in AAssays if "U6" in x.name and x.getCTForExperiment("01").value <> "Undetermined"]
    BU6 = [float(x.getCTForExperiment("01").value) for x in BAssays if "U6" in x.name and x.getCTForExperiment("01").value <> "Undetermined"]
    # Find difference in the averages
    adjuster = sum(AU6)/len(AU6) - sum(BU6)/len(BU6)
    # Adjust the CT by this difference
    for assay in BAssays:
        for ct in assay.cts:
            if ct.value <> "Undetermined":
                ct.value = float(ct.value) + adjuster
    # Return merged assays
    mergedRes = AAssays + BAssays
    return (mergedRes,adjuster)

def copyFileToPath(filePath,targetPath):
    sourceFile = open(filePath,"r")
    lines = sourceFile.readlines()
    sourceFile.close()
    targetFile = open(targetPath,"w")
    targetFile.writelines(lines)
    targetFile.close()

def fileNamesForLetter(txtfiles,inDir, letter):
    # Returns all file names for group A or B
    fileNames = [inDir + "/" + f for f in txtfiles if f[-5:-4] == letter]
    return fileNames

def filesToAssayArray(txtfiles,inDir):
    # Turns files into array of assays
    
    (AFileNames, BFileNames) =  (fileNamesForLetter(txtfiles,inDir, "A"),fileNamesForLetter(txtfiles,inDir, "B"))
    AAssays = []
    BAssays = []
    AFileNames.sort()
    BFileNames.sort()
    
    # Get the A and B assay names
    (targetnames,cts) = getTargetNamesAndCTs(open(AFileNames[0],"r").readlines(),AFileNames[0])
    for x in xrange(0,len(targetnames)):
        AAssays.append(Assay(targetnames[x],"A"))
    

    (targetnames,cts) = getTargetNamesAndCTs(open(BFileNames[0],"r").readlines(),BFileNames[0])
    for x in xrange(0,len(targetnames)):
        BAssays.append(Assay(targetnames[x],"B"))
    
    # Append the cts to each assay
    for AFileName in AFileNames:
        (targetnames,cts) = getTargetNamesAndCTs(open(AFileName,"r").readlines(),AFileName)
        for x in xrange(0,len(targetnames)):
            for AAssay in AAssays:
                if AAssay.name == targetnames[x]:
                    AAssay.cts.append(CT(cts[x],AFileName.replace(".txt","").replace("A","").replace(inDir+"/","")))
    
    print "Finished the As"

    for BFileName in BFileNames:
        (targetnames,cts) = getTargetNamesAndCTs(open(BFileName,"r").readlines(),BFileName)
        for x in xrange(0,len(targetnames)):
            for BAssay in BAssays:
                if BAssay.name == targetnames[x]:
                    BAssay.cts.append(CT(cts[x],BFileName.replace(".txt","").replace("B","").replace(inDir+"/","")))
    return (AAssays,BAssays)

def getTargetNamesAndCTs(lines, fname):
    # Returns all the assay names and ctvalues in a file
    resLine  = -1   # 1 means this is the line that contains the [Results] marker
    tInt = -1
    cInt = -1
    targetNames = []
    cts = []
    rowLen = -1
    u6 = 1
    for x in xrange(0,len(lines)):
        lines[x]=lines[x].replace('\n','').replace('\r','')
        if resLine <> -1:   # If this is the line after resLine, get column information
            if x == resLine + 1:
                spl = lines[x].split("\t")
                rowLen = len(spl)
                for y in xrange(0,len(spl)):
                    if spl[y] == "Target Name":
                        tInt = y
                    elif spl[y] == "CT":
                        cInt = y
            elif x >= resLine + 2 and tInt <> -1 and cInt <> -1:    # If this is a later line, get CT and TargetName
                spl = lines[x].split("\t")
                if len(spl) < rowLen and x>resLine+2:  # If this is not a row, ignore it
                    continue
                if "U6" in spl[tInt]:
                    spl[tInt] += str(u6)    # Give U6 assays distinct names
                    u6 += 1
                targetNames.append(spl[tInt])
                cts.append(spl[cInt])
        elif "[Results]" in lines[x]:
            resLine = x
    if tInt == -1 or cInt == -1:
        raise Exception(fname + " does not have the CT and Target Name columns")
    return (targetNames, cts)
   
def histogramDictionary(L, binSz, low, high):
    # Create a histogram in the form of a dictionary
    d = {}
    for x in L:
        myPos = int((x - low)/binSz)
        if myPos in d:
            d[myPos] += 1
        else:
            d[myPos] = 1
    return d

def drawHistogramString(D, binSz, low, high):
    # Draw histogram
    outString = ""
    for key in sorted(D.keys()):
        outString += str(int(key)*binSz + low)
        outString += "|"
        for x in xrange(0,int(D[key])):
            outString += "*"
        outString += "\n"
    return outString

def labelPatients(patients, groupfile):
    if groupfile is None:
        groupfile = os.path.dirname(os.path.realpath(__file__)) + "/PatientGroups/patgroups.txt"
    try:
        gf = open(groupfile,"r")
    except:
        raise Exception("Could not open patient group data")
    groupData = []
    for line in gf.readlines():
        spl = line.split()
        groupData.append(tuple((spl[0],tuple([str(spl[x]) for x in range(1,len(spl))]))))
    for p in patients:
        for gDat in groupData:
            if str(p.name) in gDat[1]:
                p.label = int(gDat[0])
    return tuple([int(gd[0]) for gd in groupData]) # group numbers

def minVarianceTimesMeanAssays(assays):
    minAssays = [assay for assay in assays if assay.group == "A" and assay.hasAllCTs()]
    for assay in minAssays:
        assay.comparisonMetric = 0.0
        mn = 0.0
        for ct in assay.cts:
            mn += float(ct.value)
        mn /= len(assay.cts)
        for ct in assay.cts:
            assay.comparisonMetric += (float(ct.value)-mn)**2
        assay.comparisonMetric /= (len(assay.cts) -1.0)
        assay.comparisonMetric *= mn
    minAssays.sort(key=lambda assay: assay.comparisonMetric, reverse=False)
    return (minAssays,"Variance*Mean")

def minChiSqAssays(assays):
    consideredAssays = [assay for assay in assays if assay.group == "A" and assay.hasAllCTs()]
    exps = tuple(set([ct.experiment for ct in consideredAssays[0].cts]))
    for assay in consideredAssays:
        if len(assay.cts) <> len(exps):
            raise Exception("Error: CT dimension mismatch for least variant assays")
    expected = []
    for exp in exps:
        currSum = 0.0
        for assay in consideredAssays:
            currSum += float(assay.getCTForExperiment(exp).value)
        currSum /= (1.0*len(consideredAssays))
        expected.append(tuple((exp,currSum)))
    for assay in consideredAssays:
        assay.comparisonMetric = 0
        for ele in expected:
            assay.comparisonMetric += (float(assay.getCTForExperiment(ele[0]).value) - ele[1])**2/ele[1]
    consideredAssays.sort(key=lambda assay: assay.comparisonMetric, reverse=False)
    return (consideredAssays,"Chi-Squared")

def minAssayOutputString(minAssays, metricName):
    outString = "DirectoryCombine.py\n\n"
    outString += str(len(minAssays)) + " assays were A assays that had a numerical CT for every patient\n\n"
    outString += ("%5s     %-30s%-20s\n\n")%("Well","Target Name",metricName)
    for x in xrange(0,len(minAssays)):
        outString += ("%5s     %-30s%-20s\n")%(x+1,minAssays[x].name,round(float(minAssays[x].comparisonMetric),3))
    
    for assay in minAssays:
        L = [float(ct.value) for ct in assay.cts]
        outString += "\n\n" + assay.name + "\n\n"
        outString += drawHistogramString(histogramDictionary(L, 1, 0, 40), 1, 0, 40)
    return outString

def patientsToAssays(patients,oldAssays):
    assays = [Assay(x.name,x.group) for x in oldAssays]
    for patient in patients:
        for mir in patient.mirs:
            for assay in assays:
                if assay.name == mir.name:
                    assay.cts.append(CT(mir.ct, patient.name))
    return assays

def populateAssaysForNames(AAssayNames,BAssayNames,AFileNames,BFileNames,inDir):
    AAssays = [Assay(AAssayNames[x],"A") for x in range(0,len(AAssayNames))]
    BAssays = [Assay(BAssayNames[x],"B") for x in range(0,len(BAssayNames))]
    for AFileName in AFileNames:
        AFile = open(AFileName)
        assayCounter = 0
        for line in AFile.readlines():
            if line[0].isdigit():
                spl = line.split("\t")
                if spl[1] <> AAssays[assayCounter].name:
                    raise Exception("Error")        
                AAssays[assayCounter].cts.append(CT(spl[2],AFileName.replace(inDir+"/", "").replace(".txt","").replace("A", "")))             
                assayCounter += 1
        AFile.close()
    for BFileName in BFileNames:
        BFile = open(BFileName)
        assayCounter = 0
        for line in BFile.readlines():
            if line[0].isdigit():
                spl = line.split("\t")
                if spl[1] <> BAssays[assayCounter].name:
                    raise Exception("Error")        
                BAssays[assayCounter].cts.append(CT(spl[2],BFileName.replace(inDir+"/", "").replace(".txt","").replace("B", "")))             
                assayCounter += 1
        BFile.close()
    return (AAssays,BAssays)
    
class CT:
    def __init__(self, val, exp):
        self.value = val
        self.experiment = exp
        self.group = ""

class Assay:
    def __init__(self,nm,grp):
        self.name = nm
        self.group = grp
        self.cts = []
        self.comparisonMetric = -1
    
    def getCTForExperiment(self,exp):
        for x in self.cts:
            if x.experiment == exp:
                return x
        return None
    def hasAllCTs(self):
        for ct in self.cts:
            if ct.value == "Undetermined":
                return False
        return True

class Patient:
    def __init__(self,nm):
        self.name = nm
        self.mirs = []
        self.label = None
    def getMirForName(self,nm):
        for mir in self.mirs:
            if mir.name == nm:
                return mir
        return None
    def normalizeWRTMirNames(self,mirnames):
        mirCTs = [float(x.ct) for x in self.mirs if x.name in mirnames and x.ct <> "Undetermined"]
        avgCT = sum(mirCTs)/(1.0*len(mirCTs))
        for mir in self.mirs:
            if mir.ct <> "Undetermined":
                mir.ct = float(mir.ct) - avgCT
    def toSVMFormatString(self):
        outString = str(self.label) + " "
        for x in xrange(0,len(self.mirs)):
            if self.mirs[x].letter == "B":  # Skip B-assays
                continue;
            outString += str(x+1) + ":"
            if self.mirs[x].ct == "Undetermined":
                outString += "40 "
            else:
                outString += str(self.mirs[x].ct) + " "
        outString = outString.strip()
        return outString
    def toSVMFormatStringForMirs(self,mirNames):
        if (self.label is None):
            self.label = -1
        outString = str(self.label) + " "
        x = 1
        for mir in self.mirs:
            if mir.name in mirNames:
                outString += str(x) + ":"
                x += 1
                if mir.ct == "Undetermined":
                    outString += "40 "
                else:
                    outString += str(mir.ct) + " "
        return outString
            
class Mir:
    def __init__(self,nm,c,lt):
        self.name = nm
        self.ct = c
        self.letter = lt

class SVMHelper:
    def __init__(self,libsvmdir, accMode):
        self.accMode = accMode
        self.libsvmdir = libsvmdir
        is_win32 = (sys.platform == 'win32')
        if not is_win32:
            self.svmscale_exe = libsvmdir + "/svm-scale"
            self.svmtrain_exe = libsvmdir + "/svm-train"
            self.svmpredict_exe = libsvmdir + "/svm-predict"
            self.grid_py = libsvmdir + "/tools/grid.py"
            self.gnuplot_exe = "/usr/bin/gnuplot"
            self.libsvmtoolsdir = libsvmdir + "/tools"
            self.sep = "/"
        else:
        # example for windows
            self.svmscale_exe = libsvmdir + r"\windows\svm-scale.exe"
            self.svmtrain_exe = libsvmdir + r"\windows\svm-train.exe"
            self.svmpredict_exe = libsvmdir + r"\windows\svm-predict.exe"
            self.gnuplot_exe = r"c:\tmp\gnuplot\binary\pgnuplot.exe"
            self.grid_py = libsvmdir + r"\tools\grid.py"
            self.libsvmtoolsdir = libsvmdir + r"\tools"
            self.sep = "\\"
        assert os.path.exists(self.svmscale_exe),"svm-scale executable not found"
        assert os.path.exists(self.svmtrain_exe),"svm-train executable not found"
        assert os.path.exists(self.svmpredict_exe),"svm-predict executable not found"
        assert os.path.exists(self.gnuplot_exe),"gnuplot executable not found"
        assert os.path.exists(self.grid_py),"grid.py not found"
        
    def scaleFile(self,origFileName,scaledFileName):
        currDir = os.getcwd()
        os.chdir(self.libsvmdir)
        os.system( self.svmscale_exe + " -l -1 -u 1 " + origFileName + " > " + scaledFileName)
        os.chdir(currDir)

    def getModelLines(self,fileName):
        currDir = os.getcwd()
        os.chdir(self.libsvmdir)
        fullFileName = self.libsvmdir + self.sep + fileName + ".scale"
        os.system(self.svmtrain_exe + " " + fullFileName + " " + self.libsvmdir + self.sep + fileName + ".model")
        modFile = open(self.libsvmdir + self.sep + fileName + ".model","r")
        lines = modFile.readlines()
        os.chdir(currDir)
        return lines
        
    def pickParamsAndCV(self,fileName):
        currDir = os.getcwd()
        # cmd = '{0} -svmtrain "{1}" -gnuplot "{2}" "{3}"'.format(self.grid_py, self.svmtrain_exe, self.gnuplot_exe, self.scaled_file)
        cmd = '{0} -svmtrain "{1}" "{2}"'.format(self.grid_py, self.svmtrain_exe, fileName)
        f = Popen(cmd, shell = True, stdout = PIPE).stdout
        line = ''
        while True:
            last_line = line
            line = f.readline()
            if not line: break
        c,g,rate = map(float,last_line.split())
        os.chdir(currDir)
        return (c,g,rate)
    
    def featureselect(self,l1,l2,patients,mirsAndFScores,maxMirs):
        mirNames = [m[0] for m in mirsAndFScores]
        datFileName = self.libsvmdir + self.sep + "PatientData"
        
        maxCV = 0.0
        bestMirs = [x for x in mirNames]
        
        if maxMirs >= len(mirNames):
            maxMirs = len(mirNames) - 1
        
        mirNames = [x for x in mirNames[0:maxMirs]]
        while(len(mirNames) > 0):
            outString = self.patientsToDataString(patients, l1, l2, mirNames)
            fileMySVMData = open(datFileName,"w")
            fileMySVMData.write(outString)
            fileMySVMData.close()
            
            currDir = os.getcwd()    
            
            os.chdir(self.libsvmdir)
            self.scaleFile(datFileName, datFileName + ".scale")
            
            os.chdir(self.libsvmtoolsdir)
            (c,g,rate) = self.pickParamsAndCV(datFileName + ".scale")
            if float(rate) > maxCV:
                maxCV = float(rate)
                bestMirs = [x for x in mirNames]
            # check for accMode = 0
            elif float(rate) < maxCV and self.accMode == 0:
                break
            mirNames = [x for x in mirNames[0:len(mirNames)-1]]
            os.chdir(currDir)
        print "Best CV = " + str(maxCV) + " with num mirs = " + str(len(bestMirs))
        return bestMirs
    
    def featureselectV2(self,l1,l2,patients,mirsAndFScores,maxMirs):
        mirNames = [m[0] for m in mirsAndFScores]
        datFileName = self.libsvmdir + self.sep + "PatientData"
        
        maxCV = 0.0
        bestMirs = []
        
        if maxMirs >= len(mirNames):
            maxMirs = len(mirNames) - 1
        
    
        consideredMirs = [x for x in mirNames[0:maxMirs]]
        mirNames = [consideredMirs[0]]
        tried = 1
        while(tried < len(consideredMirs)):
            outString = self.patientsToDataString(patients, l1, l2, mirNames)
            fileMySVMData = open(datFileName,"w")
            fileMySVMData.write(outString)
            fileMySVMData.close()
            
            currDir = os.getcwd()    
            
            os.chdir(self.libsvmdir)
            self.scaleFile(datFileName, datFileName + ".scale")
            os.chdir(self.libsvmtoolsdir)
            (c,g,rate) = self.pickParamsAndCV(datFileName + ".scale")
            if float(rate) > maxCV:
                maxCV = float(rate)
                bestMirs = [x for x in mirNames]
            elif float(rate) < maxCV:
                mirNames.pop()
            tried += 1
            mirNames.append(consideredMirs[tried-1])
            os.chdir(currDir)
        mirNames = mirNames[1:]
        outString = self.patientsToDataString(patients, l1, l2, mirNames)
        fileMySVMData = open(datFileName,"w")
        fileMySVMData.write(outString)
        fileMySVMData.close()
            
        currDir = os.getcwd()    
            
        os.chdir(self.libsvmdir)
        self.scaleFile(datFileName, datFileName + ".scale")
        os.chdir(self.libsvmtoolsdir)
        (c,g,rate) = self.pickParamsAndCV(datFileName + ".scale")
            
        if float(rate) > maxCV:
            maxCV = float(rate)
            bestMirs = [x for x in mirNames]
        os.chdir(currDir)
        print "Best CV = " + str(maxCV) + " with num mirs = " + str(len(bestMirs))
        return bestMirs
    
    def patientsToDataString(self,patients,l1,l2,mirNames):
        outString = ""
        for patient in patients:
            if not(patient.label is None):
                if patient.label == l1 or patient.label == l2:
                    outString += patient.toSVMFormatStringForMirs(mirNames) + "\n" 
        outString = outString.strip()
        return outString
        
    def mirsAndFScores(self,patients,possibleMirNames,l1,l2):
        # Returns list of mirs and fscores [(mirname,fscore),(mirname,fscore)...]
        mirsAndFScores = {} # key: mir name, value: fscore
        for mirName in possibleMirNames:
            (totSum, posSum, negSum) = (0.0, 0.0, 0.0)
            (totCount,posCount,negCount) = (0.0, 0.0, 0.0)
            for pat in patients:
                if pat.label <> l1 and pat.label <> l2:
                    continue
                val = pat.getMirForName(mirName).ct
                if val == "Undetermined":
                    val = 40
                totSum += float(val)
                totCount += 1
                if pat.label == l1:
                    posSum += float(val)
                    posCount += 1
                elif pat.label == l2:
                    negSum += float(val)
                    negCount += 1
            totMean = totSum/totCount
            posMean = posSum/posCount
            negMean = negSum/negCount
            numerator = (posMean-totMean)**2 + (negMean-totMean)**2
            (posSqSum,negSqSum) = (0.0, 0.0)
            for pat in patients:
                if pat.label <> l1 and pat.label <> l2:
                    continue
                val = pat.getMirForName(mirName).ct
                if val == "Undetermined":
                    val = 40
                if pat.label == l1:
                    posSqSum += (val - posMean)**2
                elif pat.label == l2:
                    negSqSum += (val - negMean)**2
            denominator = 1.0/(posCount-1)*posSqSum + 1.0/(negCount-1)*negSqSum
            if denominator == 0:
                numerator = 0
                denominator = 1
            theFScore = numerator/denominator
            mirsAndFScores[mirName] =  theFScore
        tupMirsAndFScores = [(mirsAndFScores[k],k) for k in mirsAndFScores.keys()]
        tupMirsAndFScores.sort(reverse=True)
        tupMirsAndFScores = [(x[1],x[0]) for x in tupMirsAndFScores]
        return tupMirsAndFScores
    
    def testPatientsLines(self,patients,possibleMirNames,l1,l2,modelFilePath):
        currDir = os.getcwd()
        os.chdir(self.libsvmdir)
        dataString = self.patientsToDataString(patients, l1, l2, possibleMirNames)
        dataFilePath = self.libsvmdir + self.sep + "PatientData"
        # Write testing data to file
        fileMySVMData = open(dataFilePath,"w")
        fileMySVMData.write(dataString)
        fileMySVMData.close()
        # Scale data
        self.scaleFile(dataFilePath, dataFilePath + ".scale")
        # Test Data
        os.system(self.svmpredict_exe + " " + dataFilePath + ".scale " +  modelFilePath + " " + dataFilePath + ".pred")
        # Read Prediction
        predFile = open(dataFilePath + ".pred","r")
        predLines = predFile.readlines()
        predFile.close()
        os.chdir(currDir)
        return predLines
    
main()
