#!/usr/bin/env python -i

## ProminentEnds version ag00 12-24-22
## ProminentEnds finds k-mers and surrounding segments which are candidates for DNA termini
## The program is reference independent, so the input is one or more short-read sequencing runs (generally Illumina)
## Output is a list of K-mers which are candidates for DNA termini that appear frequently at the beginnings of reads but more rarely than expected internally
## Because standard Tagmentation ("Nextera") often captures ends near but not at the end (with some variability), the program
## focuses on K-mers that are often the beginning of a read but where the presence downstream by a specified amount is not common.
## So for every observed K-mer that is observed at the start of one or more reads we count reads that start with that given K-mer (surrogate for possible end capture)
## This is compared to a count of internal reads (starting this "internal" count several nucleotides in since the end-biased tagmentation is somewhat sloppy around the terminus
## Arbitrary set points are present in the initial block of the program (FirstInternal is an important parameter, this is where the program starts counting "internal' instances of a K-mer
## All of the values in the initial block can be reset in the command line
## The program requires VSG_ModuleEN for the command line interface
## Note that the program yields candidates for termini in a mixed sequence pool, not definitive ends (which would require additional validation)
## For sequences that exist both as internal and terminal, sensitivity will depend on both the parameters below and on the details of the fragment capture and sequencing.
## This is a work in progress-- check back frequently to see if there are updated versions

## Example Syntax
##    python ProminentEnds##.py data=<mydata.fastq/mydata.fastq.gz>
## Output is a table with the following fields

## Columns in any output 
##   0: DataSet                  ## Sequence Read file where the K-mer was found
##   1: EndRatioRank             ## Rank among the candidate K-mer termini
##   2: Reads_Starting_With_KMer ## Number of reads that start with that K-mer
##   3: Reads_With_Later_Kmer    ## Number of reads with that K-mer internal to the read (after position X)
##   4: Positions_in_read        ## K-mer occurence by position in read.  (negative=antisense instances), positioned relative to end or read (so -1 is an antisense copy of the k-mer at the very end of a read).  Numbers in parentheses are deduplicated counts
##   5: Read-Consensus-upstream  ## Consensus sequence upstream of K-mer based on any internal reads
##   6: K-mer                    ## Sequence of K-mer 
##   7: Read-Consensus-downstream## Consensus sequence downstream of K-mer from read data

## Optional columns (numbers depend on which are included based on command line choices)
##   Instances                   ## List of all reads containing a given K-mer (only if recordallinstances is set to True)
##   Reference Position(s)       ## Identity and Positions of K-mer start in the reference (e.g. genome) specified by user
##   Ref-Consensus-upstream      ## Consensus sequence upstream of K-mer in the reference file
##   K-mer                       ## Sequence of K-mer (repeat of column 6 to allow facile copy paste of entire region around k-mer)
##   Ref-Consensus-downstream    ## Consensus sequence downstream of K-mer in the reference file
##   Sample_ID                   ## ID of original sample (if a RunTable is available)

## Basic parameters in identifying terminus-enriched k-mers
klen1 = 25   ## (kmer=) ProminentEnds search for k-mers of this length that are preferentially at the beginnings of reads
MinStartCount1 = 10  ## (minstart=) ProminentEnds will only look at k-mers that are present at the beginning of MinStartCount1 reads (default 10) 
FirstInternal1 = 8 ## (firstinternal=) To obtain a measure of the overall incidence of a given k-mer ProminentEnds will look for that k-mer at positions >firstinternal nucleotides into the read.  So k-mers between 2 and firstinternal-1 bases into the read are effectively ignored in deciding whether a k-mer is likely toi be terminal since they could be either a truly internal k-mer or a due to variable capture relative to the true end  

MinStartToInternalRatio1 = 1.0 ## (minstarttointernalratio=) This value sets the minimum ratio of starts to internal k-mer occurences.  1.0 is a very arbitrary ratio and this can be adjusted if needed.  Note that the terminal k-mers are in one position while the internal k-mers can be in 100 or more positions, so a null hypothesis for this value might be very low (e.g. 0.01).  But preferences in Tagmentation, capture and sequencing will tend to produce many false positives as this value lowers below about 1.0
IgnorePairedEnd1 = False ## (ignorepairedend=) Setting this to True will only use the indicated R1 file (not the paired file)
## Input data (files can be fasta or fastq and can be gzipped)
DataFile0 = '' ## (data=)  Test With 'ERR1938563' This can be a single fasta/fastq file, a group of files, or a wildcard ('*') designation.  Each file is handled separately and candidate k-mers are listed in a single output file by source file
RunsToSamples1 = '' ## (runstosamples=) An optional file that maps run numbers into sample names. General format is "SampleName<tab>Run,Run,Run" with one sample on each line
StopAfterRead1 = 0 ## (stopafterread=)  Instructs ProminentEnds to stop after a specified number of reads.  0=Analyze all reads
## Output formatting and details
OutputFileBase0 = 'default' ## (out=) This will indicate the name for the output file  (default will assign a name based on time and date)
delimiter1 = '\r' ## (delimiter=) The line endings used for output
ExampleCount1 = 4 ## (examplecount=)  How many exemplary reads to return for each k-mer (not implemented)
MinReported1 = 10  ## (minreported=) ProminentEnds will return at least this many candidate k-mers
MaxReported1 = 100 ## (maxreported=) ProminentEnds will return at most this many candidate k-mers
RecordAllInstances1 = False ## (recordallinstances=)  Setting this to true reports all relevant reads for each K-mer (massive additions to the output table but useful for some downstream assembly tasks)

ReferenceFile0 = '' ## (ref=) A file or list of files with reference sequences for restriction and/or mapping of k-mers-- test with 'BdAdV1.fasta' 
RestrictToReference1 = False ## (restricttoreference=) Setting this to true will return only k-mers mapping to the reference file (default "False" maps k-mers to reference file but doesn't resttrict to only mapped k-mers)
ReferenceCircular1 = False ## (referencecircular=) Setting this to true will allow permuted k-mers from the reference
MaxAllRefLenForPrefilter1 = 10000000 ## (maxallreflengforprefilter=) This allows a faster filtering for a short ref.  Memory and speed tradeoffs here but no difference in results (10M is a reasonable limit on a modern system)
## Assembly of local consensus sequences
## ProminentEnds assembles a very rough sequence around each identified k-mer based on the consensus from reads containing the k-mer.
## This is a simple "plurailty" algorithm with each position relative to the k-mer assigned based on how many reads would have that base in a simple context of a single sequence containing the k-mer
## Where there is only one read or any ambiguity, a small letter is used, where there is no majority or strong plurality an "N" is added
## The threshholds for calling a definitive base or a tentative but likely base are as set below
DefinitiveCall1 = 0.95 ## (definitivecall=) Threshold percentage to definitively call as base in the consensus local assembly
TentativeCall1 = 0.5 ## (tentatitvecall=) Threshold percentage for a tentative (small letter) call of a base in the reported local consensus

## Removal of linker-containing (artefactual) fake-ends
TetritisFile1 = 'illuminatetritis.fa' ## (tetritisfile=).  A file with various linkers used in Illumina or other sequencing
tetritisK1 = 12 ## (tetritisK=) This is the K-mer length used in removing linker-derived k-mers that are commonly different between datasets

## Other unlikely-to-be-needed settings
regularization1 = 1 ## (regularization=).  This is a relatively minor regularization value for determining the ratio between starting and internal k-mer occurences
PairedIDStart1 = 10 ## (pairedidstart=) This is the position in the paired read used to deduplicate individual paired end combinations.
PairedIDEnd1 = 18  ## (pairedidend=)
FastQDumpProgram1 = 'fastq-dump' ## (fastqdump=) or 'fasterq-dump' Full path of a program that will download files from SRA if needed.
Threads1 = 1# 8 ## (threads=) Number of threads to execute in parallel.  Only applies for multiple files being analzed
ReportGran1 = 1000000 ## (reportgranularity=) How often to report progress during analysis
KsInSummary1 = 10 ## (ksinsummary=) How many K-mers to include in the summary blurb for each dataset

import gzip
from collections import Counter
from VSG_ModuleEN import *
from random import randint, choice
import subprocess
from itertools import chain, zip_longest
import glob
vcommand()

class StartK():
    def __init__(self, k):
        self.k = k
        self.a = antisense(k)
        self.instances = Counter()
        self.readpos = {}
        self.readcontext = {}
        self.refpos = []
        self.refcontext = {}
        
def antisense(s):
    return s.replace('G','c').replace('C','g').replace('A','t').replace('T','a')[::-1].upper()
antibase = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}

def find1(n):
    ''' find a file on this machine'''
    o = []
    for r,d,f in os.walk('/'): 
        if n in f:
            o.append(os.path.join(r, n))
    return o
def versioner1(v88):
    s88 = str(v88).split('version')[-1].split('.')
    v1 =[]
    for v88 in s88:
        v1.append(int(''.join(list(filter(str.isdigit, v88)))))
    return v1
def findfastqdump(candidate):
    ver = 0
    try:
        ver = subprocess.check_output([candidate,'-V'])
        return candidate
    except:
        pass
    if os.path.isfile('FastQDumpLocation.txt'):
        candidate = open('FastQDumpLocation.txt', mode='rt').read()
        try:
            ver = subprocess.check_output([candidate,'-V'])
            return candidate
        except:
            pass
    vLog('Looking for fastq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    vLog('or reinstall and allow execution (chmod +X <path>)')
    vLog('Note this can take some real time, get a cup of coffee or tea')
    NewCands = find1('fastq-dump')
    NewItems = []
    for candidate in NewCands:
        try:
            ver = subprocess.check_output([candidate,'-V'])
            NewItems.append([versioner1(ver),candidate])
        except:
            pass
    if NewItems:
        candidate = sorted(NewItems, reverse=True)[0][1]
        open('FastQDumpLocation.txt', mode='w').write(candidate)
        try:
            open('~/FastQDumpLocation.txt', mode='w').write(candidate)
        except:
            pass
        return candidate
    vLog('Unable to find fast-q dump.  Recommend that you reinstall this from ncbi or download fastq/fasta files directly')
    return ''
    
def ContextToConsensus1(cx,sx):
    ups1 = ''
    j1 = -1
    while j1 in cx:
        myD1 = cx[j1]
        t1 = sum(myD1.values())
        b1,x1 = myD1.most_common()[0]
        if x1 > DefinitiveCall1*t1 and x1>1:
            ups1 += b1.upper()
        elif x1 > TentativeCall1*t1:
            ups1 += b1.lower()
        else:
            ups1 += 'N'
        j1-=1
    ups1 = ups1[::-1]
    j1 = 0
    dwn1 = ''
    while j1 in cx:
        myD1 = cx[j1]
        t1 = sum(myD1.values())
        b1,x1 = myD1.most_common()[0]
        if x1 > DefinitiveCall1*t1 and x1>1:
            dwn1+=b1.upper()
        elif x1 > TentativeCall1*t1:
            dwn1+=b1.lower()
        else:
            dwn1+='N'
        j1 += 1
    if not(ups1): ups1=' '
    if not(dwn1): dwn1=' '
    return (ups1,sx,dwn1)
def DictionaryDedupStats1(D):
    DeDupD = {}
    CountD = {}
    for d in D:
        DeDupD[d] = len(D[d])
        CountD[d] = sum(D[d].values())
    SortedKeys = sorted(list(D.keys()),key=lambda x:CountD[x], reverse=True)
    Report = '{'
    for d in SortedKeys:
        Report += str(d)+':'+str(CountD[d])+'('+str(DeDupD[d])+'), '
    if len(Report)>3:
        Report = Report[:-2]+'}'
    return Report

## Set Up RunsToSamples if a file is provided, the code below allows a number of different file structures
RunsToSamplesD1 = {}
if os.path.isfile(RunsToSamples1):
    IDFile1 = open(RunsToSamples1, mode='rt')
    for L0 in IDFile1:
        L1 = L0.split('\t')
        if len(L1)>=2:
            RunsToSamplesD1[L1[0].strip()] = L1[1].strip()
            RunsToSamplesD1[L1[1].strip()] = L1[0].strip()
            for lx11 in L1[1].split(','):
                RunsToSamplesD1[lx11.strip()] = L1[0]
            for lx11 in L1[0].split(','):
                RunsToSamplesD1[lx11.strip()] = L1[1]
            for lx11 in L1[1:]:
                RunsToSamplesD1[lx11] = L1[0]
    IDFile1.close()

if type(DataFile0)==str:
    if '*' in DataFile0:
        DataFile1 = list(glob.glob(DataFile0))
    else:
        DataFile1 = DataFile0.strip().strip('[').strip(']').strip('(').strip(')').split(',')
else:
    DataFile1 = DataFile0

if type(ReferenceFile0)==str:
    if '*' in ReferenceFile0:
        ReferenceFile1 = list(glob.glob(ReferenceFile0))
    else:
        ReferenceFile1 = ReferenceFile0.strip().strip('[').strip(']').strip('(').strip(')').split(',')
else:
    ReferenceFile1 = ReferenceFile0
    
rD1 = {}
arD1 = {}
rL1 = {}
for rf1 in ReferenceFile1:
    if os.path.isfile(rf1):
        rD1.update(vFastAToDict(rf1))
for n1 in rD1:
    rD1[n1] = rD1[n1].upper()
    arD1[n1] = antisense(rD1[n1]).upper()
    rL1[n1] = len(rD1[n1])
    if ReferenceCircular1:
        rD1[n1] = rD1[n1]+rD1[n1][:klen1-1]  
        arD1[n1] = arD1[n1]+arD1[n1][:klen1-1]

AllRefLen1 = sum(rL1.values())
rfD1 = Counter()
if RestrictToReference1 and AllRefLen1<MaxAllRefLenForPrefilter1:
    for rn1 in rD1:
        rS1 = rD1[rn1]
        arS1 = arD1[rn1]
        for ri1 in range(len(rS1)-klen1+1):
            rk1 = rS1[ri1:ri1+klen1]
            rfD1[rk1]+=1
            ark1 = arS1[ri1:ri1+klen1]
            rfD1[ark1]+=1
    




def myStr1(x):
    mst = str(x).split('Counter(')[-1].strip(')').strip(']').strip('[')
    if mst:
        return mst
    else:
        return ' '
## make a dictionary of kmers (Length tetritisK1) to get rid of anything with a k-mer from an illumina linker
tetritis=open(TetritisFile1,mode='rt').readlines()
T1=Counter()
for l1 in tetritis:
    l2=l1.strip().upper()
    a2=antisense(l2)
    for i in range(len(l2)-tetritisK1+1):
        T1[l2[i:i+tetritisK1]]=1
        T1[a2[i:i+tetritisK1]]=1

TempFileUID1 = ''
for i in range(4):
    TempFileUID1 += choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm0123456789')

if OutputFileBase0 == 'default':
    OutputFileBase0 = 'PromEndOut_'+vnow ##+'_'+TempFileUID1
OutFileName1 = OutputFileBase0+'_Results.tsv'
OutFileName2 = OutputFileBase0+'_Summary.tsv'
OutFileName1 = OutputFileBase0+'_Results.tsv'
OutFileName2 = OutputFileBase0+'_Summary.tsv'
OutFile1 = open(OutFileName1,mode='w')
OutFile2 = open(OutFileName2,mode='w')
Headers1 =['DataSet','EndRatioRank_For_The_Dataset','KMer-Starts_Read','KMer-Later_In_Read','Positions_In_Read_(Deduplicated_Count_In_Parentheses)','Read_Upstream_Consensus','KMer','Read_Downstream_Consensus']
Headers2 = ['Run','R1File','R2File','Sample','Reads','LengthFilteredReads','UniqueStarts','AbundanceFilteredStarts','RatioFilteredStarts','BestRatioHits']
if RecordAllInstances1:
    Headers1.append('Instances')
if rD1:
    Headers1.append('Ref_Position')
    Headers1.append('Ref_Upstream_Consensus')
    Headers1.append('KMer') ## repeat of column 6 included to allow facile copying of entire reference sequence by selecting over three adjacent columns
    Headers1.append('Read_Downstream_Consensus')
if RunsToSamplesD1:
    Headers1.append('Sample_ID')
def HeaderTranspose(hT2):
    hT0 = '<!--\tOutput_Key\t\t-->'+delimiter1
    hT0 += '<!--\tColumnNumber\tColumnHeader\t-->'+delimiter1
    for iT2,nT2 in enumerate(hT2):
        nT2 = nT2.strip()
        if not(nT2.startswith('<!')):
            hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'   -->'+delimiter1
    hT3 = sys.argv
    hT4 = '<!--\tResults From Running ProminentEnds End-Finder    -->'+delimiter1
    hT4 += '<!--\t    Date/Time of Run = '+vnow+'    -->'+delimiter1
    for hT5 in hT3:
        if hT5:
            hT4 += '<!--\t'+hT5+'      -->'+delimiter1
    hT4 += '<!--\t\t\t-->'+delimiter1
    return hT4+hT0+'<!--\t\t\t-->'+delimiter1
RunDetails1 = delimiter1.join(['<!--'+x+'-->' for x in vSysLogInfo1.splitlines()])
AbbrevHeader1 = '|'.join(vSysLogInfo1.replace('\t',' ').splitlines())+'\t'+vnow
##OutFile1.write('<!-- Results from ProminentEnds Run: '+vnow+' -->'+delimiter1)
OutFile1.write(HeaderTranspose(Headers1)+'<!--\t\t\t-->'+delimiter1)
## OutFile1.write(RunDetails1+delimiter1)
OutFile1.write('<!--\t\t-->'+delimiter1)
OutFile1.write('\t'.join(Headers1)+'\t'+AbbrevHeader1+delimiter1)
OutFile2.write(HeaderTranspose(Headers2)+'<!--\t\t\t-->'+delimiter1)
## OutFile1.write(RunDetails1+delimiter1)
OutFile2.write('<!--\t\t-->'+delimiter1)
OutFile2.write('\t'.join(Headers2)+'\t'+AbbrevHeader1+delimiter1)
OutFile1.close()
OutFile2.close()

for fn1 in DataFile1:
    if not(os.path.isfile(fn1)):
        FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
        break
DataFileCount1 = len(DataFile1)
def ProminentEndsProcess1(Param1,MyQueue1):
    fnum1,Fn1 = Param1[0]+1,Param1[1]
    if not(os.path.isfile(Fn1)) and not('.' in Fn1):
        if not(os.path.isfile(Fn1+'_1.fasta.gz')):
            if Threads1 == 1:
                vLog(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
                vLog("Preparing to download sequence read set "+Fn1+" from NCBI")
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                '--split-files',
                                                '--fasta',
                                                '0',
                                                '--origfmt',
                                                '--gzip',
                                                '--outdir',
                                                './',
                                                Fn1])
            vLog("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        Fn1 = Fn1+'_1.fasta.gz'
        if not(os.path.isfile(Fn1)):
            vLog('Looks Like fast(er)q-dump failed for '+Fn1)
            return ''
    Fn11 = os.path.basename(Fn1).split('.')[0]
    F2 = iter(())
    F2Open1 = False
    if '_1.' in Fn1:
        Fn2 = Fn1[::-1].replace('.1_','.2_',1)[::-1]
        Fn11a = Fn11[::-1].split('1_',1)[-1][::-1]
    elif '_R1' in Fn1:
        Fn2 = Fn1[::-1].replace('1R_','2R_',1)[::-1]
        Fn11a = Fn11[::-1].split('1R_',1)[-1][::-1]
    elif 'R1' in Fn1:
        Fn2 = Fn1[::-1].replace('1R','2R',1)[::-1]
        Fn11a = Fn11[::-1].split('1R',1)[-1][::-1]
    if Fn1.lower().endswith('fasta') or Fn1.lower().endswith('fasta.gz') or Fn1.lower().endswith('fa') or Fn1.lower().endswith('fa.gz'):
        LineGran1 = 2
    else:
        LineGran1 = 4
    if Fn1.lower().endswith('gz'):       
        F1 = gzip.open(Fn1, mode='rt')
        if not(IgnorePairedEnd1) and os.path.isfile(Fn2):
            F2 = gzip.open(Fn2, mode='rt')
            F2Open1 = True
            Fn11 = Fn11a
    else:
        F1 = open(Fn1, mode='rt')    
        if not(IgnorePairedEnd1) and os.path.isfile(Fn2):
            F2 = open(Fn2, mode='rt')
            F2Open1 = True
            Fn11 = Fn11a
    Fn11 = Fn11.split('_')[0]
    if Fn11 in RunsToSamplesD1:
        Sample1 = RunsToSamplesD1[Fn11]
    else:
        Sample1 = Fn11
    FileMnemonic1 = Fn11 + ' ('+str(fnum1)+'/'+str(DataFileCount1)+')'
    Results1 = ''
    ##[Run,Fn1,Fn2,Sample,Reads,LengthFiltered,UniqueStarts,StartsMeetingAbundance,StartMeetingRatio]
    if F2Open1:
        Summary1 = [Fn11,os.path.basename(Fn1),os.path.basename(Fn2),Sample1,0,0,0,0,0,[]]
    else:
        Summary1 = [Fn11,os.path.basename(Fn1),' ',Sample1,0,0,0,0,0,[]]        
    S1 = Counter()
    for i0,L0 in enumerate(chain(F1,F2)):
        if i0%LineGran1!=1: continue
        if StopAfterRead1 and (i0>StopAfterRead1*LineGran1): break
        if (Threads1==1) and (i0%(LineGran1*ReportGran1)==1):
            vLog(FileMnemonic1+" Collecting Starts: "+str(i0//LineGran1))
        if len(L0)>klen1:
            k1 = L0[:klen1]
            if RestrictToReference1 and rfD1 and not(k1 in rfD1): continue
            S1[k1] += 1
    Summary1[4] = i0
    Summary1[5] = sum(S1.values())
    Summary1[6] = len(S1)
    I1 = {} ## Keys are k-mers, values are internal count
    for s1 in S1:
        if S1[s1]>MinStartCount1:
            I1[s1] = 0
    if RestrictToReference1 and rD1 and not(rfD1):
        I2 = {}
        for rn1 in rD1:
            rs1 = rD1[rn1]
            ars1 = arD1[rn1]
            for ri1 in range(len(rs1)-klen1+1):
                rk1 = rs1[ri1:ri1+klen1]
                ark1 = ars1[ri1:ri1+klen1]
                if rk1 in I1:
                    I2[rk1]=0
                if ark1 in I1:
                    I2[ark1]=0
        I1 = I2
    Summary1[7] = len(I1)
    if I1:
        F1.seek(0)
        if F2Open1:
            F2.seek(0)
        sD1 = {} ## keys are kmers (sense or antisense), values are a 2-element list with first being orientation ('p' or 'm') and second being a StartK object
        for i0,L0 in enumerate(chain(F1,F2)):
            if i0%LineGran1!=1: continue
            if StopAfterRead1 and (i0>StopAfterRead1*LineGran1): break
            if len(L0)>klen1:
                for i1 in range(FirstInternal1,len(L0)-klen1):
                    s1 = L0[i1:i1+klen1]
                    if s1 in I1:
                        I1[s1] += 1
            if (Threads1==1) and (i0%(LineGran1*ReportGran1)==1):
                vLog(FileMnemonic1+" Collecting Internals: "+str(i0//LineGran1))
        MyKs = sorted(list(I1.keys()), key=lambda x:(S1[x]+regularization1)/(I1[x]+regularization1), reverse=True)
        n1 = 0
        for s1 in MyKs:
            if "N" in s1: continue
            Tetritis1 = False
            for i2 in range(len(s1)-tetritisK1+1):
                tk1 = s1[i2:i2+tetritisK1]
                if tk1 in T1:
                    Tetritis1 = True
                    break
            if Tetritis1: continue
            n1+=1
            if (S1[s1] < I1[s1]*MinStartToInternalRatio1) and n1>MinReported1: continue
            if n1>=MaxReported1: break
            if not(s1 in sD1):
                sD1[s1] = ['p',StartK(s1)]
                sD1[sD1[s1][1].a] = ['m',sD1[s1][1]]
        Summary1[8] = len(sD1)
        if sD1:
            F1.seek(0)        
            if F2Open1:
                F2.seek(0)
            for i0,(L0,M0) in enumerate(zip_longest(F1,F2,fillvalue='')):
                L0 = L0.strip()
                M0 = M0.strip()
                if i0%LineGran1!=1: continue
                if StopAfterRead1 and (i0>StopAfterRead1*LineGran1): break
                if len(L0)>=klen1:
                    for i1 in range(len(L0)-klen1+1):
                        s1 = L0[i1:i1+klen1]
                        if s1 in sD1:
                            o1,K1 = sD1[s1]
                            PairedIDTag1 = M0[PairedIDStart1:PairedIDEnd1]
                            if not(PairedIDTag1):
                                PairedIDTag1 = L0[-PairedIDEnd1:-PairedIDStart1]
                            if o1=='p':
                                if not (i1+1 in K1.readpos):
                                    K1.readpos[i1+1] = Counter()
                                K1.readpos[i1+1][(PairedIDTag1,1)] += 1
                                for j1 in range(-1,-i1-1,-1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    
                                    K1.readcontext[j1][L0[i1+j1]] += 1
                                for j1 in range(len(L0)-i1-klen1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    K1.readcontext[j1][L0[j1+i1+klen1]] += 1
                                if RecordAllInstances1:
                                    K1.instances[L0[:i1].lower()+s1+L0[i1+klen1:].lower()] += 1
                            elif o1 == 'm':
                                i11 = -(len(L0)-klen1-i1)-1
                                if not (i11 in K1.readpos):
                                    K1.readpos[i11] = Counter()
                                K1.readpos[i11][(PairedIDTag1,1)] += 1
                                aL0 = antisense(L0)
                                h1 = len(aL0)-klen1-i1
                                for j1 in range(-1,-h1-1,-1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()                            
                                    K1.readcontext[j1][aL0[h1+j1]] += 1
                                for j1 in range(len(aL0)-h1-klen1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    K1.readcontext[j1][aL0[j1+h1+klen1]] += 1
                                if RecordAllInstances1:
                                    K1.instances[aL0[:i1]+K1.a.lower()+aL0[i1+klen1:]] += 1
                if len(M0)>=klen1:
                    for i1 in range(len(M0)-klen1+1):
                        s1 = M0[i1:i1+klen1]
                        if s1 in sD1:
                            o1,K1 = sD1[s1]
                            PairedIDTag1 = L0[PairedIDStart1:PairedIDEnd1]
                            if o1=='p':
                                if not (i1+1 in K1.readpos):
                                    K1.readpos[i1+1] = Counter()
                                K1.readpos[i1+1][(PairedIDTag1,2)] += 1
                                for j1 in range(-1,-i1-1,-1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    
                                    K1.readcontext[j1][M0[i1+j1]] += 1
                                for j1 in range(len(M0)-i1-klen1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    K1.readcontext[j1][M0[j1+i1+klen1]] += 1
                                if RecordAllInstances1:
                                    K1.instances[M0[:i1].lower()+s1+M0[i1+klen1:].lower()] += 1
                            elif o1 == 'm':
                                i11 = -(len(M0)-klen1-i1)-1
                                if not (i11 in K1.readpos):
                                    K1.readpos[i11] = Counter()
                                K1.readpos[i11][(PairedIDTag1,2)] += 1
                                aM0 = antisense(M0)
                                h1 = len(aM0)-klen1-i1
                                for j1 in range(-1,-h1-1,-1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()                            
                                    K1.readcontext[j1][aM0[h1+j1]] += 1
                                for j1 in range(len(aM0)-h1-klen1):
                                    if not(j1 in K1.readcontext):
                                        K1.readcontext[j1] = Counter()
                                    K1.readcontext[j1][aM0[j1+h1+klen1]] += 1
                                if RecordAllInstances1:
                                    K1.instances[aM0[:i1]+K1.a.lower()+aM0[i1+klen1:]] += 1
                if (Threads1==1) and (i0%(LineGran1*ReportGran1)==1):
                    vLog(FileMnemonic1+" Collecting Instances: "+str(i0//LineGran1))
            if rD1:
                for rn1 in rD1:
                    rs1 = rD1[rn1]
                    ars1 = arD1[rn1]
                    for ri1 in range(len(rs1)-klen1+1):
                        rk1 = rs1[ri1:ri1+klen1]
                        if rk1 in sD1:
                            o1,K1 = sD1[rk1]
                            if o1=='p':
                                K1.refpos.append((rn1,ri1+1))
                                j1 = -1
                                p1 = ri1-1
                                while j1>=-2*klen1:
                                    if p1<0:
                                        if ReferenceCircular1:
                                            p1 += rL1[rn1]
                                        else:
                                            break
                                    if not j1 in K1.refcontext:
                                        K1.refcontext[j1] = Counter()
                                    K1.refcontext[j1][rs1[p1]] += 1
                                    j1 -= 1
                                    p1 -= 1
                                j1 = 0
                                p1 = ri1+klen1
                                while j1<2*klen1:
                                    if p1>=rL1[rn1]:
                                        if ReferenceCircular1:
                                            p1 -= rL1[rn1]
                                        else:
                                            break
                                    if not j1 in K1.refcontext:
                                        K1.refcontext[j1] = Counter()
                                    K1.refcontext[j1][rs1[p1]] += 1
                                    j1 += 1
                                    p1 += 1
                            elif o1=='m':
                                K1.refpos.append((rn1,-(ri1+1)))
                                j1 = -1
                                p1 = ri1+klen1
                                while j1>=-2*klen1:
                                    if p1>=rL1[rn1]:
                                        if ReferenceCircular1:
                                            p1 -= rL1[rn1]
                                        else:
                                            break
                                    if not j1 in K1.refcontext:
                                        K1.refcontext[j1] = Counter()
                                    K1.refcontext[j1][antibase[rs1[p1]]] += 1
                                    j1 -= 1
                                    p1 += 1
                                j1 = 0
                                p1 = ri1-1
                                while j1<2*klen1:
                                    if p1<0:
                                        if ReferenceCircular1:
                                            p1 += rL1[rn1]
                                        else:
                                            break
                                    if not j1 in K1.refcontext:
                                        K1.refcontext[j1] = Counter()
                                    K1.refcontext[j1][antibase[rs1[p1]]] += 1
                                    j1 += 1
                                    p1 -= 1
            n1=1
            for s1 in sD1:
                o1,K1 = sD1[s1]
                if o1=='m': continue
                if len(Summary1[-1])<KsInSummary1:
                    Summary1[-1].append(s1)
                reups1,rekme1,redwn1 = ContextToConsensus1(K1.readcontext,s1)
                if rD1:
                    rfups1,rfkme1,rfdwn1 = ContextToConsensus1(K1.refcontext,s1)
                DedupStats1 = DictionaryDedupStats1(K1.readpos)
                Results1 += '\t'.join(map(myStr1,(Fn11,n1,S1[s1],I1[s1],DedupStats1,reups1,rekme1,redwn1)))
                if RecordAllInstances1:
                    if K1.instances:
                        Results1 += '\t'+K1.instances
                    else:
                        Results1 += '\t '
                if rD1:
                    if K1.refpos:
                        Results1 += '\t'+myStr1(K1.refpos)
                        Results1 += '\t'+'\t'.join((rfups1,rfkme1,rfdwn1))
                    else:
                        Results1 += '\t \t \t \t '
                if RunsToSamplesD1:
                    Results1 += '\t'+Sample1
                Results1 += '\t '+delimiter1
                n1 += 1
    F1.close()
    if F2Open1:
        F2.close()
    if Threads1 == 1:
        vLog("Finished "+ FileMnemonic1)
    Summary1 = '\t'.join(map(myStr1,Summary1))+delimiter1
    if Threads1 ==1:
        return (Results1,Summary1)
    else:
        MyQueue1.put((Results1,Summary1))


def main():
    OutFile1 = open(OutFileName1,mode='a')
    OutFile2 = open(OutFileName2,mode='a')
    if Threads1 == 1:
        for xP1 in enumerate(DataFile1):
            ResSum1 = ProminentEndsProcess1(xP1,1)
            OutFile1.write(ResSum1[0])
            OutFile2.write(ResSum1[1])
    else:
        import multiprocessing as mp
        from time import sleep
        mp.set_start_method('fork')
        MyQueue1 = mp.Queue(100)
        LPL1 = len(DataFile1)
        OngoingProcesses1 = []
        CurrentTask1 = {} ## This is a dictionary associating each thread with the current process
        for NextProcess1 in range(min(Threads1,LPL1)):
            OngoingProcesses1.append(mp.Process(target=ProminentEndsProcess1, args=((NextProcess1,DataFile1[NextProcess1]),MyQueue1)))
            OngoingProcesses1[-1].start()
            CurrentTask1[NextProcess1] = os.path.basename(DataFile1[NextProcess1]).split('_1')[0]+' ('+str(NextProcess1+1)+'/'+str(LPL1)+'), Thread-'+str(NextProcess1)
            vLog('Starting '+CurrentTask1[NextProcess1])
        NextProcess1 += 1
        while OngoingProcesses1:
            for i,p1 in enumerate(OngoingProcesses1):            
                if not(p1.is_alive()):
                    OngoingProcesses1[i].join()
                    ResSum1 = MyQueue1.get()
                    OngoingProcesses1[i].close()
                    OutFile1.write(ResSum1[0])
                    OutFile2.write(ResSum1[1])
                    vLog("Finishing "+CurrentTask1[i])
                    if NextProcess1<LPL1:
                        OngoingProcesses1[i] = mp.Process(target=ProminentEndsProcess1, args=((NextProcess1,DataFile1[NextProcess1]),MyQueue1))
                        OngoingProcesses1[i].start()
                        CurrentTask1[i] = os.path.basename(DataFile1[NextProcess1]).split('_1')[0]+' ('+str(NextProcess1+1)+'/'+str(LPL1)+'), Thread-'+str(i)
                        NextProcess1 += 1
                        vLog('Starting '+CurrentTask1[i])
                    else:
                        del(OngoingProcesses1[i])
            sleep(0.01)
        OutFile1.close()
        OutFile2.close()

if __name__ == '__main__':
    main()

##
##  Copyright 2022, Andrew Fire and Stanford University
##  Revision List (partial)
##
