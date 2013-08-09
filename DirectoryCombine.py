#!/usr/bin/env python

# Assumptions Made:
    # Script is in directory containing libsvm directory, patgroups.txt, training directory, and (optional) testing directory
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
#from mx.Misc.Cache import DENOM     # for command line programs

def main():
    
    '''
    ###########################
    CMD LINE ARGS & INPUT FILES
    ###########################
    '''

    aParser = argparse.ArgumentParser(description="Aggregate PCR files and optionally run SVM")
    aParser.add_argument("-c","--compMetric",help="Metric used to normalize between patients (0 = chi-squared (default), 1 = variance*mean)")
    aParser.add_argument("-f","--feature_set_size",help="Number of micro-RNAs to look for when determining most distinguishing micro-RNAs")
    argNamespace = aParser.parse_args()

    scriptDir = os.path.dirname(os.path.realpath(__file__))
    libsvmdir = checkForLibSVM()
    groupfile = checkForGroupFile()
    trainingDir = checkForTrainingData()
    testingPatientsFile = checkForTestingPatientsFile()
    feature_set_size = argNamespace.feature_set_size
    compMetric = argNamespace.compMetric
    if compMetric is None:
        compMetric = 0
    if feature_set_size is None:
        feature_set_size = -1

    if not(libsvmdir or groupfile):
        print "\n"
        print "Please ensure that the libsvm directory and the patient group file (patgroups.txt) are both present, then rerun the script"
    if not(trainingDir):
        print "\n"
        print "Please ensure that the script is placed in a directory which CONTAINS the directory of training data (call it training_data)"
    
    '''
    ###########################
    NORMALIZATION
    ###########################
    '''

    # Find the text files
    print "Identifying text files..."
    sys.stdout.flush() 
    txtfiles = [f for f in listdir(trainingDir) if isfile(join(trainingDir,f)) and (".txt" == f[-4:])]
    
    # Get array of A and B assays    
    print "Combining text file data..."
    sys.stdout.flush()
    (AAssays,BAssays) = filesToAssayArray(txtfiles,trainingDir)
    
    # Print out table
    print "Producing raw output..."
    sys.stdout.flush()
    fileRawOutput = open(scriptDir + "/Output/A_RawOutput.txt","w")
    fileRawOutput.write(assayOutputString((AAssays,BAssays),ABMerged=True))
    fileRawOutput.close()
    
    # Combine A and B assays
    if len(BAssays) > 0:
        print "Merging A and B assays..."
        sys.stdout.flush()
        (assays, adjuster) = combineABAssays(AAssays,BAssays)

        # Print out another table
        print "Producing combined output..."
        sys.stdout.flush()
        fileCombineOutput = open(scriptDir + "/Output/B_CombinedOutput.txt","w")
        fileCombineOutput.write(assayOutputString(assays,ABMerged=False))
        fileCombineOutput.close()
    else:
        assays = AAssays
        adjuster = 0
        
    
    # Find least variant assays according to the metric that the user inputted
    print "Finding least variant assays..."
    sys.stdout.flush()
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
    
    # Determine which files need to be tested
    print "Determining files to be tested..."
    sys.stdout.flush()
    testPatients = [x.replace("\n","") for x in open(testingPatientsFile).readlines()]
    
    # Get all candidate mirs
    print "Choosing mirs from A group with all CTs..."
    sys.stdout.flush()
    lvMirNames = [x.name for x in assays if x.group == "A"]
    
    
    for x in xrange(0,len(groups)-1):
        for y in xrange(x+1,len(groups)):
            # Generate training data string
            outString = ""
            for patient in patients:
                if not(patient.label is None):
                    if patient.label == groups[x] or patient.label == groups[y]:
                        outString += patient.toSVMFormatStringForMirs(lvMirNames) + "\n" 
            outString = outString.strip()
            # Write data string to file
            dataFile = open(libsvmdir + "/PatientData","w")
            dataFile.write(outString)
            dataFile.close()

            # Generate testing data string
            outString = ""
            for patient in patients:
                if patient.name in testPatients:
                    outString += patient.toSVMFormatStringForMirs(lvMirNames) + "\n"
            outString = outString.strip()
            testingFile = open(libsvmdir+"/PatientDataTesting","w")
            testingFile.write(outString)
            testingFile.close()
            
            currDir = os.getcwd()
            os.chdir(scriptDir + "/" + libsvmdir+"/tools")
            
            if (not (testingPatientsFile is None)):
                if feature_set_size is None:
                    os.system("python fselect.py ../PatientData ../PatientDataTesting")
                else:
                    os.system("python fselect.py ../PatientData " + str(feature_set_size) + " ../PatientDataTesting")
                os.chdir(currDir)
                modelFile = [f for f in listdir(scriptDir+"/"+libsvmdir+"/tools") if isfile(join(scriptDir+"/"+libsvmdir+"/tools",f)) and (".model" in f)][0]
                predFile = [f for f in listdir(scriptDir+"/"+libsvmdir+"/tools") if isfile(join(scriptDir+"/"+libsvmdir+"/tools",f)) and (".pred" in f)][0]
                copyFileToPath(scriptDir+"/"+libsvmdir+"/tools/"+modelFile, scriptDir + "/Output/Model" + str(groups[x]) + "_" + str(groups[y]) + ".model")
                copyFileToPath(scriptDir + "/" + libsvmdir + "/PatientDataTesting",scriptDir + "/Output/Testing" + str(groups[x]) + "_" + str(groups[y]))
                copyFileToPath(scriptDir+"/"+libsvmdir+"/tools/"+predFile, scriptDir + "/Output/Prediction" + str(groups[x]) + "_" + str(groups[y]) + ".pred")
            else:
                if feature_set_size is None:
                    os.system("python fselect.py ../PatientData")
                else:
                    os.system("python fselect.py ../PatientData " + str(feature_set_size))
            
            os.chdir(currDir)
            copyFileToPath(libsvmdir + "/PatientData",scriptDir + "/Output/Training" + str(groups[x]) + "_" + str(groups[y]))
            copyFileToPath(libsvmdir + "/tools/PatientData.select", scriptDir + "/Output/Select" + str(groups[x]) + "_" + str(groups[y]) + ".select")
            os.chdir(currDir)
            
                    
    print "Done"

def assayOutputString(assays,ABMerged=False,tabDelimited=True):
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
    if ABMerged and len(BAssays) > 0:
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
    if tabDelimited:
        outString += "Well\tTarget_Name"
    else:
        outString += ("%-10s%-30s") % ("Well","Target_Name")
    
    for exp in exps:
        if tabDelimited:
            outString += "\t" + exp
        else:
            outString += ("%-30s")%(exp)
        
     
    outString += "\n\n"
    
    # For all the AAssays
    for x in xrange(0,len(AAssays)):
        # Print assay well and name
        if tabDelimited:
            outString += str(x+1) + "\t" + AAssays[x].name + "A"
        else:
            outString += ("%-10s%-30s") % (x+1,AAssays[x].name + "A")
        
        # Print all CT values
        for exp in exps:
            ct = AAssays[x].getCTForExperiment(exp)
            if tabDelimited:
                if ct:
                    outString += "\t" + ("Undetermined" if ct.value == "Undetermined" else str(round(float(ct.value),3)))
                else:
                    outString += "\t" + "X"
            else:
                if ct:  # If the CT does not exist, print an X. If it is NOT "Undetermined", round it to 3 decimal places
                    outString += ("%-30s") % ("Undetermined" if ct.value == "Undetermined" else round(float(ct.value),3))
                else:
                    outString += ("%-30s") % ("X")
        outString += "\n"
    
    # Repeat above steps for BAssays
    # NOTE: If ABMerged = True, BAssays = [], so this code would not execute at all
    for x in xrange(0,len(BAssays)):
        if tabDelimited:
            outString += str(x+1) + "\t" + AAssays[x].name + "B"
        else:
            outString += ("%-10s%-30s") % (x+1,BAssays[x].name + "B")
        for exp in exps:
            ct = BAssays[x].getCTForExperiment(exp)
            if tabDelimited:
                if ct:
                    outString += "\t" + ("Undetermined" if ct.value == "Undetermined" else str(round(float(ct.value),3)))
                else:
                    outString += "\t" + "X"
            else:
                if ct:  # If the CT does not exist, print an X. If it is NOT "Undetermined", round it to 3 decimal places
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

def checkForLibSVM():
    libsvmdir = None
    try:
        libsvmdir = [f for f in listdir('.') if 'libsvm' in f][0]
        print "Found", libsvmdir
        sys.stdout.flush()
        return libsvmdir
    except:
        print "Could not find libsvm..."
        print "Please ensure that this script (DirectoryCombine.py) is in a directory which CONTAINS the libsvm directory (ex. libsvm-3.17)"
        sys.stdout.flush()
        return None

def checkForGroupFile():
    patgroups = None
    try:
        patgroups = [f for f in listdir('.') if 'patgroups.txt' in f][0]
        print "Found patgroups.txt..."
        sys.stdout.flush()
        return patgroups
    except:
        print "Could not find patgroups.txt"
        print "Please ensure that this script (DirectoryCombine.py) is in a directory which CONTAINS the patient groups file (patgroups.txt)"
        sys.stdout.flush()
        return None

def checkForTestingPatientsFile():
    testingPatientsFile = None
    try:
        testingPatientsFile = [f for f in listdir('.') if 'testing_patients.txt' in f][0]
        print "Found testing_patients.txt..."
        sys.stdout.flush()
        return testingPatientsFile
    except:
        print "Could not find testing_patients.txt"
        print "Please ensure that this script (DirectoryCombine.py) is in a directory which CONTAINS the testing patients file (testing_patients.txt)"
        sys.stdout.flush()
        return None

def checkForTrainingData():
    trainingDir = None
    try:
        trainingDir = [f for f in listdir('.') if 'patient_data' in f][0]
        print "Found training_data..."
        sys.stdout.flush()
        return trainingDir
    except:
        print "Could not find training_data"
        print "Please ensure that this script (DirectoryCombine.py) is in a directory which CONTAINS the directory of training data files (called training_data)"
        sys.stdout.flush()
        return None


def combineABAssays(AAssays,BAssays):
    # Pick a reference patient to compute adjuster. That is, which experiment do A and B have in common (ex. if 03 was a experiment, there would be a 03A.txt and 03B.txt)
    aCts = AAssays[0].cts
    bCts = BAssays[0].cts
    refExp = None
    for aCt in aCts:
        if not(refExp is None):
            break
        for bCt in bCts:
            if bCt.experiment == aCt.experiment:
                refExp = aCt.experiment
                break

    # Merge A and B assays using the U6 reference assay
    # Get U6 assay CTs from A and B assays
    AU6 = [float(x.getCTForExperiment(refExp).value) for x in AAssays if "U6" in x.name and x.getCTForExperiment(refExp).value != "Undetermined"]
    BU6 = [float(x.getCTForExperiment(refExp).value) for x in BAssays if "U6" in x.name and x.getCTForExperiment(refExp).value != "Undetermined"]
    # Find difference in the averages
    adjuster = sum(AU6)/len(AU6) - sum(BU6)/len(BU6)
    # Adjust the CT by this difference
    for assay in BAssays:
        for ct in assay.cts:
            if ct.value != "Undetermined":
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

def filesToAssayArray(txtFiles,trainingDir):
    ''' Turns files into list of assays '''
    
    # Get files that end with *A.txt and *B.txt
    AFileNames = [trainingDir + "/" + f for f in txtFiles if f[-5:-4] == "A"]
    BFileNames = [trainingDir + "/" + f for f in txtFiles if f[-5:-4] == "B"]
    AAssays = []
    BAssays = []
    
    # IF there are A files
    if len(AFileNames) > 0:
        # Get the target names (the mir names) and CT values from the first file
        (targetnames,cts) = getTargetNamesAndCTs(open(AFileNames[0],"r").readlines(),AFileNames[0])
        # Add this to list of A assays
        AAssays += [Assay(x,"A") for x in targetnames]

        # For all the A files,
        for AFileName in AFileNames:
            (targetnames,cts) = getTargetNamesAndCTs(open(AFileName,"r").readlines(),AFileName)
            # for each target name within the file
            for x in xrange(0,len(targetnames)):
                # find that target name in the list of A assays
                for AAssay in AAssays:
                    if AAssay.name == targetnames[x]:
                        # and add the ct value to the A assay's list of CTs (NOTE, the second argument to the CT constructor is experiment number)
                        AAssay.cts.append(CT(cts[x],AFileName.replace(".txt","").replace("A","").replace(trainingDir+"/","")))
    # Repeat the above procedure for B files
    if len(BFileNames) > 0:
        (targetnames,cts) = getTargetNamesAndCTs(open(BFileNames[0],"r").readlines(),BFileNames[0])
        BAssays = [Assay(x,"B") for x in targetnames]
    
        for BFileName in BFileNames:
            (targetnames,cts) = getTargetNamesAndCTs(open(BFileName,"r").readlines(),BFileName)
            for x in xrange(0,len(targetnames)):
                for BAssay in BAssays:
                    if BAssay.name == targetnames[x]:
                        BAssay.cts.append(CT(cts[x],BFileName.replace(".txt","").replace("B","").replace(trainingDir+"/","")))
    return (AAssays,BAssays)

def getTargetNamesAndCTs(lines,fname):
    ''' Use state machine to extract target names and CTs from a file'''
    # params: lines = list of the lines in the file, fname = name of file

    # states
    LOOKING_FOR_DATA = 0
    READING_HEADER = 1
    READING_DATA = 2
    state = LOOKING_FOR_DATA

    target_name_col = -1    # column of the data which contains the target name
    ct_col = -1             # column of the data which contians the CT
    num_cols = -1           # number of columns of data

    u6 = 1

    # list of target names and cts from the file
    targetNames = []
    cts = []

    # for every line
    for x in xrange(0,len(lines)):
        # if we see [Results] then the following line is the header for the table (i.e. the listing of column names)
        if state == LOOKING_FOR_DATA:
            if "[Results]" in lines[x]:
                state = READING_HEADER
        # read the header to find out target_name_col, ct_col, and num_cols
        elif state == READING_HEADER:
            split_line = lines[x].split("\t")
            num_cols = len(split_line)
            for y in xrange(0,len(split_line)):
                if split_line[y] == "Target Name":
                    target_name_col = y
                elif split_line[y] == "CT":
                    ct_col = y
            if target_name_col == -1 or ct_col == -1:   # panic if we can't find the ct_col or target_name_col
                raise Exception(fname + " does not have the CT and Target Name columns")
            state = READING_DATA    # otherwise, change state because the following line is a data line
        # split the line and get the target name (from the target_name_col) and ct (from the ct_col)
        elif state == READING_DATA:
            split_line = lines[x].split("\t")
            if len(split_line) != num_cols:
                continue
            if "U6" in split_line[target_name_col]: # Give each U6 assay a different name
                split_line[target_name_col] += str(u6)
                u6 += 1
            targetNames.append(split_line[target_name_col])
            cts.append(split_line[ct_col])
    return (targetNames, cts)

'''
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
        if resLine != -1:   # If this is the line after resLine, get column information
            if x == resLine + 1:
                spl = lines[x].split("\t")
                rowLen = len(spl)
                for y in xrange(0,len(spl)):
                    if spl[y] == "Target Name":
                        tInt = y
                    elif spl[y] == "CT":
                        cInt = y
            elif x >= resLine + 2 and tInt != -1 and cInt != -1:    # If this is a later line, get CT and TargetName
                spl = lines[x].split("\t")
                if len(spl) != rowLen:  # If this is not a row, ignore it
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
'''   


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
        groupfile = os.path.dirname(os.path.realpath(__file__)) + "/patgroups.txt"
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
        if len(assay.cts) != len(exps):
            raise Exception(("%s has %s cts, but there are %s patients")%(str(assay.name),str(len(assay.cts)),str(len(exps))))
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

def populateAssaysForNames(AAssayNames,BAssayNames,AFileNames,BFileNames,trainingDir):
    AAssays = [Assay(AAssayNames[x],"A") for x in range(0,len(AAssayNames))]
    BAssays = [Assay(BAssayNames[x],"B") for x in range(0,len(BAssayNames))]
    for AFileName in AFileNames:
        AFile = open(AFileName)
        assayCounter = 0
        for line in AFile.readlines():
            if line[0].isdigit():
                spl = line.split("\t")
                if spl[1] != AAssays[assayCounter].name:
                    raise Exception("Error")        
                AAssays[assayCounter].cts.append(CT(spl[2],AFileName.replace(trainingDir+"/", "").replace(".txt","").replace("A", "")))             
                assayCounter += 1
        AFile.close()
    for BFileName in BFileNames:
        BFile = open(BFileName)
        assayCounter = 0
        for line in BFile.readlines():
            if line[0].isdigit():
                spl = line.split("\t")
                if spl[1] != BAssays[assayCounter].name:
                    raise Exception("Error")        
                BAssays[assayCounter].cts.append(CT(spl[2],BFileName.replace(trainingDir+"/", "").replace(".txt","").replace("B", "")))             
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
        mirCTs = [float(x.ct) for x in self.mirs if x.name in mirnames and x.ct != "Undetermined"]
        avgCT = sum(mirCTs)/(1.0*len(mirCTs))
        for mir in self.mirs:
            if mir.ct != "Undetermined":
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
        #assert os.path.exists(self.gnuplot_exe),"gnuplot executable not found"
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
                if pat.label != l1 and pat.label != l2:
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
                if pat.label != l1 and pat.label != l2:
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
