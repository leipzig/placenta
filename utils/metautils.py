import pandas
import os
import sys

#different readlengths for trimmomatic outputs
READLENGTHS = [36,75]

#from utils import metautils
#metautils.srpMeta('SRP141397')
#study_accession            experiment_accession            sample_accession            run_accession            experiment_title            experiment_attribute            taxon_id            library_selection            library_layout            library_strategy            library_source            library_name            bases            spots            adapter_spec            avg_read_length
#      SRP141397                      SRX3980107                  SRS3205258                SRR7049033                S78_shotgun                                                   0                       RANDOM                  PAIRED -                         WGS               METAGENOMIC             S78_shotgun        229716627          1009555                                 227.54245880610765


#read in an SPR metadata

class srpMeta():
    def __init__(self,SRP):
        self.srp=SRP
        self.st = pandas.read_csv("metadata/"+SRP+".metadata",sep="\t")
        #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz', 

        # split into 16s and shotgun sequences
        self.st["subject"] = self.st["experiment_title"].str.split("_", expand = True)[0]
        self.st["experimental_strategy"] = self.st["experiment_title"].str.split("_", expand = True)[1]
        
        self.ALL_BNAMES = [os.path.basename(f) for f in self.getFilesFromRunList(self.populationRuns('ALL'))]

        #>>> metautils.ALL_TRIMMED
        #['NA06984.1.M_111124_4_1.trimmed.fastq.gz', 'NA06984.1.M_111124_4_2.trimmed_75.fastq.gz',
        self.ALL_TRIMMED = [os.path.basename(f).replace('.fastq','.trimmed_{}.fastq'.format(tlen)) for tlen in READLENGTHS for f in self.getFilesFromRunList(self.populationRuns('ALL'))]

        #['ERR188/ERR188325/NA06984.1.M_111124_4_1.fastq.gz',
        self.ALL_PATHS = [f for f in self.getFilesFromRunList(self.populationRuns('ALL'))]
        self.ALL_LOOKUP =  dict(zip(self.ALL_BNAMES, self.ALL_PATHS))

    def getSRRs(self):
        return(self.st['run_accession'].tolist())

    def get16SFiles(self):
        files=[]
        filesofinterest=self.st[self.st['experimental_strategy']=='16S']
        for pair in ["1","2"]:
            #'study_accession','experiment_accession','run_accession'
            filesofinterest["pair"+pair]="raw/"+filesofinterest['study_accession']+"/"+filesofinterest['experiment_accession']+"/"+filesofinterest['run_accession']+"_"+pair+".fastq"
        files=list(filesofinterest['pair1'].astype(str))+list(+filesofinterest['pair2'].astype(str))
        return(files)

    def getFilesFromRunList(self,runs):
        files=[]
        
        for run in runs:
            layout=self.st.loc[self.st['run_accession']==run]['library_layout']
            if layout.to_string(index=False).lstrip().startswith('PAIRED'):
                files+=[self.srp+"/fastq/"+run+"_1.fastq.gz",self.srp+"/fastq/"+run+"_2.fastq.gz"]
            else:
                files+=[self.srp+"/fastq/"+run+"_1.fastq.gz"]
        return(files)

    
    def getAllFiles(self):
        getFilesFromRunList(self.populationRuns('ALL'))

    #>>> metautils.calculateExpected('ALL')
    #Expecting 413 fastq and 253 bams
    def calculateExpected(self,population):
        fastq = self.getFilesFromRunList(self.populationRuns(population))
        bams = self.populationRuns(population)
        print("Expecting {0} fastq and {1} bams".format(len(fastq),len(bams)))
    
    #GBR FIN ALL YRI TSI 
    #>>> metautils.populationRuns('ALL')
    #['NA12286.2.MI_120126_5', 'NA10851.4.M_120208_1',
    def populationRuns(self,population):
        return(self.getSRRs())

    #Check to see if paired-end sequencing was used in each sample so appropriate stAR parameters are called
    def isPaired(self,run):
        if self.st.loc[self.st['Assay Name']==run]['Comment[LIBRARY_LAYOUT]'].all()=='PAIRED':
            return(True)
    
    #Either run star with paired or unpaired options
    def starProgram(self,run):
        if self.isPaired(run):
            return("star_paired.sh")
        else:
            return("star_unpaired.sh")
    
    #Run trimmomatic bash script with paired or unpaired options
    def trimProgram(self,run):
        if self.isPaired(run):
            return("trimmomatic_paired.sh")
        else:
            return("trimmomatic_unpaired.sh")
            
    #pairedOrSinglefastqInput('NA06984.1.M_111124_4')
    #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz']
    #pairedOrSinglefastqInput('NA06985.1.MI_120104_3')
    #['NA06985.1.MI_120104_3_1.fastq.gz']
    
    #set file names for trimmomatic based on whether a trimmed input or output is expected
    def pairedOrSinglefastqInput(self,run,trimmed=False,readlength=None):
        if trimmed:
            extName=".trimmed_{0}.fastq.gz".format(readlength)
        else:
            extName=".fastq.gz"
        if self.isPaired(run):
            return([run+"_1"+extName,run+"_2"+extName])
        else:
            return([run+"_1"+extName])
    

    def getPathFromName(self,name):
        return(ALL_LOOKUP[name])
    
    #ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188404/ERR188404_1.fastq.gz

    
    def getType(self,runName):
        return(self.getProperty(runName, 'Platform'))
    
    
    def getProperty(self, runName, column):
        stseries = self.st.loc[self.st['Run'] == runName, column]
        if len(self.stseries) == 1:
            return(self.stseries.to_self.string(index=False))  # pandas is ridiculous
        elif len(self.stseries) > 1:
            raise ValueError("Not expecting to match multiple run names, just 1")
        else:
            raise ValueError("Can't find that run {0}".format(runName))
    
    
    def getMemory(self, runName):
        return(int(getProperty(runName, 'size_MB')))
    
    #Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed
    def cleanJobname(self, name):
        return(name.replace('.','_'))
    
    def getECS(self, runName, units, program):
        if program == 'minimap':
            mb = 32000
        elif program == 'star':
            superIntensive = ['NA12716.7.M_120219_6', 'NA12775.7.M_120219_5',
                              'NA12814.7.M_120219_3', 'NA11831.7.M_120219_2',
                              'NA11994.2.M_111216_7', 'NA11893.7.M_120219_3']
            if runName in superIntensive:
                mb = 64000
            else:
                mb = 48000
        elif program == 'IsoModule':
            mb = 192000
        elif program == 'samtoolsindex':
            mb = 16000
        elif program == 'samtoolsmerge':
            mb = 16000
        elif program == 'samtoolssubsample':
            mb = 16000
        elif program == 'fastx':
            mb = 8000
        elif program == 'trimmomatic':
            mb = 8000
        elif program == 'rmats':
            mb = 16000
        else:
            raise ValueError
        if units == 'bytes':
            return(1048576*mb)
        elif units == 'mb':
            return(mb)
    # 2  # number of sample
    # 3  # number of replicates in sample #1
    # sam1_rep1.bam
    # sam1_rep2.bam
    # sam1_rep3.bam
    # 3  # number of replicates in sample #2
    # sam2_rep1.bam
    # sam2_rep2.bam
    # sam2_rep3.bam
    
    #Used for a two way comparison between various samples
    def twoSampleComparisonManifest(self, biosamp1, biosamp2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for biosamp in [biosamp1, biosamp2]:
            run = getRunBamsFromBioSample(biosamp, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    
        #Used for a two way comparison between various samples
    def twoTreatmentComparisonManifest(self, treatment1, treatment2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for treatment in [treatment1, treatment2]:
            run = getRunBamsFromTreatment(treatment, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    # "panorama-clk-repro/SRP091981/
    
    #>>> metautils.getRunBamsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam.bai']
    #Get a List of bams originating from the same biological source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunBamsFromBioSample(self, biosamp, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromBioSample(biosamp)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    
    def getRunBamsFromTreatment(self, treatment, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromTreatment(treatment)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    #>>> metautils.getRunsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4', 'NA12778.1.M_111124_4', 'NA12778.1.MI_120104_3']
    #Get a List of runs from the same source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunsFromBioSample(self, biosample,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['Source Name'] == biosample) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))
    
    def getRunsFromTreatment(self, treatment,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['treatment'] == treatment) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))

    #metautils.populationBiosamples('ALL')
    #['NA12778', 'NA12045', 'NA12144', ...
    #Returns either all runs that pass the QC check for a specific population or all populations depending on arguments
    def populationBiosamples(self,population):
        if(population=='ALL'):
            return(list(set(self.st.loc[(self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
        else:
            return(list(set(self.st.loc[(self.st['Characteristics[population]'] == population) & (self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
    
    #This returns all sequencing runs over 70 nt (longruns) looks through the population being considered appends them to the sample List using the mergedbam file naming convention, and then appends all runs from that population with a 36nt trimmed suffix
    def getListOfLengthOutputs(self,population):
        x = 0
        biosamps = self.populationBiosamples(population)
        biosamp_lengthList = []
        longruns = list(set(self.st.loc[self.st['Comment[SEQUENCE_LENGTH]'].astype(int)<70]['Source Name'].toList()))
        for item in longruns:
            if (item in biosamps) and (item != 'NA07000'):
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        for item in biosamps:
            if item != 'NA07000':
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        return(biosamp_lengthList)

if __name__== "__main__":
    srpMeta(sys.argv[1]).get16SFiles()