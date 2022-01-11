#!/usr/bin/env python

from urllib import request, parse
from sys import argv,stderr
import subprocess
import argparse
import platform
import os
import re


# Available platforms: Windows (native CMD, WSL and Cygwin), Linux, Unix and Darwin (i.e. MacOS)

#import the config.py file containing the paths to needed programs and pHMM files
import config as cnf

version="1.1"

##### Command Line Interface (CLI) #####
parser = argparse.ArgumentParser(prog="rtk-pred.py", description="RTK-PRED: automated detection and annotation of Receptor Tyrosine Kinases (RTKs) with profile Hidden Markov Models (pHMMs)")
parser.version=version
help=parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--input", type=str, action='store', help="input file in FASTA format", required=True)
required.add_argument("-o", "--output", type=str, action='store', help="the RTK-PRED output directory prefix", required=True)
parser.add_argument("--mkdir", action='store_true', help="if set, it automatically creates the output directory specified by '-o'.")
parser.add_argument("-wp","--webphobius", action='store_true', help="use the web-server version of Phobius instead of a locally installed one (useful for Windows and MacOS machines, for which the 'decodeanhmm' binary used by Phobius is not available, as well as those users without a Phobius license), requires internet connection")
parser.add_argument("--hmmerdir", action='store', help="set the location of the HMMER binaries directory (overrides config.py)")
parser.add_argument("--phobiusdir", action='store', help="set the location of the Phobius binaries directory (overrides config.py), not applicable when -wp/--webphobius is set.  Also, not applicable in Windows.")
parser.add_argument("-v", "--version", action='version', help="display version and exit")
parser._action_groups.append(help)
args = parser.parse_args()

##### General-purpose Classes and Functions used throughout the script #####
def readFile(filehandle, list=False):
    f=open(filehandle,"r")
    if list==True:
        content=f.readlines()
    else:
        content=f.read()
    f.close()
    return content

class Fasta():
    def __init__(self):
        self.title=str()
        self.sequence=str()
    def parseFasta(self, input, fasta_headers_dict, fasta_sequences_dict):
        if input[0] == ">":
            fa = input.split(">", 1)[1]
            fa = fa.split("\n")
        else:
            fa=input.split("\n")

        # The fasta header is stored to a variable.
        # If the fasta header is not in the keys of the dictionary, a new pair is formed which has the fasta header as the key and the number 1 as its value. The fasta header remains unchanged.
        # If the fasta header is in the keys of the dictionary, the value of that key is the counter which shows how many times this fasta header has been found within the fasta sequences.
        # This counter is increased by one, is stored to that key and it is concatenated with the fasta header thus a unique fasta header is generated and stored to the title of the object.
        fasta_header = fa.pop(0)
        if fasta_header in fasta_headers_dict.keys():
            current_counter = fasta_headers_dict[fasta_header]
            new_counter = current_counter + 1
            fasta_headers_dict[fasta_header] = new_counter
            fasta_header = "Unique_" + str(new_counter) + "_" + fasta_header
        else:
            fasta_headers_dict[fasta_header] = 1
        
        if fasta_header in fasta_sequences_dict.keys():
            pass
        else:
            fasta_header_for_dict = fasta_header.split()[0]
            fasta_sequences_dict[fasta_header_for_dict]="".join(fa)
        
        self.title=fasta_header
        self.sequence="".join(fa)

        # The dictionary with the fasta headers is returned.
        fastas_properties_dict = {
          1: fasta_headers_dict,
          2: fasta_sequences_dict
        }
        return fastas_properties_dict

def readMultiFasta(input):
    # A dictionary for the headers of the fasta sequences is initialized.
    fasta_headers_dict = {}
    fasta_sequences_dict = {}

    f=readFile(input, list=False)
    mfa=f.split("\n>")
    mfasta=list()
    for i in mfa:
        if i.strip()!="":
            fasta=Fasta()

            # The dictionary with the fasta headers is stored to a variable.
            fastas_properties_dict = fasta.parseFasta(i, fasta_headers_dict, fasta_sequences_dict)
            fasta_headers_dict = fastas_properties_dict[1]
            fasta_sequences_dict = fastas_properties_dict[2]

            mfasta.append(fasta)
    fasta_dictionary = {
      1: mfasta,
      2: fasta_sequences_dict
    }
    return fasta_dictionary

def writeMFasta(mfasta, output):
    o=open(output,"w")
    for fasta in mfasta:
        o.write(">%s\n%s\n" %(fasta.title, fasta.sequence))
    o.close()





##### The RTK class #####
class RTK():
    def __init__(self):
        self.id=str()
        self.signal=list()
        self.domains=list()
        self.transmem=list()
        self.ptk=list()
        self.family=str()
        self.sequence=str()


#create dict from list of RTK objects
def rtk_dict(rtks):
    dct={}
    for rtk in rtks:
        dct[rtk.id]=rtk
    return dct


##### HMMER-related Classes and Methods #####
class HMMER():
    def __init__(self):
        self.title=str()
        self.domain=str()
        self.domain_id=str()
        self.score=float()
        self.start=int()
        self.end=int()
        self.index=int() #for indexing purposes, in proteins with multiple copies of the same domain

# create dict from list of HMMER objects
def hmm_dict(hmms):
    dct={}
    for h in hmms:
        if " " in h.title:
            dct[h.title.split()[0]]=h
        else:
            dct[h.title]=h
    return dct

def run_hmmer(exe, hmm, input, outname, gathering=True):
    cmd="%s " %exe
    if gathering==True:
        cmd+="--cut_ga "
    cmd+="%s %s >%s/%s" %(hmm, input, output, outname)
    subprocess.call(cmd, shell=True)

def parse_hmmer_result(input, runtype):
    # A dictionary for the titles of the hmmer results.
    hmmer_headers_dict = {}

    hmm_return=list()
    hmmout=readFile(input, list=False)
    if runtype=="hmmsearch":
        results=hmmout.split("\n>> ")
        header=results.pop(0)
        for r in results:
            align=r.split("\n")
            domain=align[0].split()[0]
            table=r.split("Alignments for each domain:")[0]
            aln=table.split("\n")
            for line in aln[1:]:
                stats=line.split()
                if len(stats)>0:
                    # If the current "stats" list line starts from "#" it means that a new set of results for this specific domain has been found so the boolen variable for whether its a new hmmer
                    # title for a new set of results of domains is set to True.
                    if stats[0] == "#":
                        different_hmmer_title = True

                    if stats[0] not in ["#", "---", "[No"] and stats[1]=="!":
                        # The hmmer title is stored.
                        hmmer_title=aln[0]
                        # If the hmmer title exists within the keys of the dictionary with the hmmer titles then...
                        if hmmer_title in hmmer_headers_dict.keys():
                            # If the info refers to a new set of domain results then the hmmer title is going to change if not the object must have the same hmmer title.
                            if different_hmmer_title:
                                # Same procedure as for the Fasa objects. If the hmmer title exists within the keys of the dictionary then the counter/value of the key is increased by one,
                                # it is stored again for the same key and it is concatenated with the hmmer title thus making it unique (and same with the title of one of the Fasta objects).
                                current_hmmer_counter = hmmer_headers_dict[hmmer_title]
                                new_hmmer_counter = current_hmmer_counter + 1
                                hmmer_headers_dict[hmmer_title] = new_hmmer_counter
                                hmmer_title = "Unique_" + str(new_hmmer_counter) + "_" + hmmer_title
                            # Even if the hmmer title is not for a different set of domains that hmmer title must take the last unique hmmer title as every domain corresponds to a title.
                            else:
                                current_hmmer_counter = hmmer_headers_dict[hmmer_title]
                                if current_hmmer_counter != 1:
                                    hmmer_title = "Unique_" + str(current_hmmer_counter) + "_" + hmmer_title
                        # If not the hmmer title is stored as a new key for the dictionary with the value of 1.
                        else:
                            hmmer_headers_dict[hmmer_title] = 1
                            
                        # After the first hmmer is stored to the dictionary the boolean is set to False as a new indication for a new set of domain results must be found so to change the hmmer
                        # title. If this did not happen then two domain hits for the same fast sequence (hmemr title) would be stored to hmmer objets of different titles which would make it then
                        # impossible it to be identified that they have indeed more than one domain.
                        
                        different_hmmer_title = False

                        result=HMMER()
                        result.domain=domain
                        result.title=hmmer_title
                        result.score=float(stats[2])
                        result.start=int(stats[9])
                        result.end=int(stats[10])
                        hmm_return.append(result)

    else:
        records=hmmout.split("\n//")
        for record in records:
            if record != "\n[ok]\n":
                results=record.split(">> ")
                header=results.pop(0).split("\n")
                for line in header:
                    if line!="" and line!=" ":
                        if line.split()[0]=="Query:":
                            id=line.split()[1]
                for r in results:
                    align=r.split("\n")
                    domain=align[0].split()[0]
                    table=r.split("Alignments for each domain:")[0]
                    aln=table.split("\n")
                    for line in aln[1:]:
                        stats=line.split()
                        if len(stats)>1:
                            if stats[0] not in ["#", "---", "[No"] and stats[1]=="!":
                                result=HMMER()
                                result.domain=domain
                                result.title=id
                                result.score=float(stats[2])
                                result.start=int(stats[9])
                                result.end=int(stats[10])
                                hmm_return.append(result)

    return hmm_return


##### Phobius-related Classes and Methods #####
class Segment():
    def __init__(self, type, start, end):
        self.type=type
        self.start=int(start)
        self.end=int(end)

class Phobius():
    def __init__(self, id, topology):
        self.id=id
        self.topology=topology

def run_phobius(input, output, web=False):
    if web==False:
        cmd="perl %s -long %s > %s/phobius.txt" %(phobius, input, output)
        subprocess.call(cmd, shell=True)
    else:
        #url="http://phobius.sbc.su.se/cgi-bin/predict.pl"
        url=cnf.phobius_url
        data=parse.urlencode({"protseq":input, "format":"nog"}).encode()
        req=request.Request(url, data=data)
        response=request.urlopen(req)
        result=response.read().decode()
        d=re.split("<pre>|<\/pre>",result)
        f=open("%s/phobius.txt" %output, "w")
        f.write(d[1])
        f.close()

def parse_phobius_result(input):
    f=readFile(input, list=False)
    entries=f.split("\n//")
    results=list()
    for entry in entries:
        text=entry.split("\n")
        topology=list()
        id=""
        for line in text:
            i=line.split()
            if len(i)>0: #this solves the problem of whitespace lines that remain in the list
                if i[0]=="ID":
                    id=i[1]
                if i[0]=="FT":
                    if i[1]=="SIGNAL":
                        signal=Segment("SIGNAL", i[2], i[3])
                        topology.append(signal)
                    elif i[1]=="TOPO_DOM" or i[1]=="DOMAIN":
                        if i[4]=="NON" or i[4]=="CYTOPLASMIC.":
                            if i[4]=="NON":
                                i[4]="NON CYTOPLASMIC"
                            topo=Segment(i[4].rstrip("."), i[2], i[3])
                            topology.append(topo)
                    elif i[1]=="TRANSMEM":
                        tm=Segment("TRANSMEM", i[2], i[3])
                        topology.append(tm)
                    else:
                        continue
        if id!="":
            prediction=Phobius(id, topology)
            results.append(prediction)
    return results

def get_rtks_from_phobius(input):
    rtks=[]
    for res in input:
        sp=None
        tm=[]
        extra=[]
        cyto=[]
        #print(res.id)
        for p in res.topology:
            #print(p.type)
            if p.type=="SIGNAL":
                sp=[p.start, p.end]
            if p.type=="TRANSMEM":
                tm.append([p.start, p.end])
            if p.type=="NON CYTOPLASMIC":
                extra.append([p.start, p.end])
            if p.type=="CYTOPLASMIC":
                cyto.append([p.start, p.end])
        if len(tm)==1:
            rtk=RTK()
            rtk.id=res.id
            rtk.signal=sp
            rtk.tm=tm[0]
            rtks.append(rtk)
            del rtk
    return rtks


##### Pfam-related and RTK classification stuff #####
pfam ={
	"Recep_L_domain"  : "PF01030",
	"Cadherin"        : "PF00028",
	"Furin-like"      : "PF00757",
	"ig"              : "PF00047",
	"Ig_2"            : "PF13895" ,
	"Ig_3"            : "PF13927" ,
	"fn3"             : "PF00041",
	"Sema"            : "PF01403",
	"TIG"             : "PF01833",
	"Ig_Tie2_1"       : "PF10430" ,
	"LRR_8"           : "PF13855" ,
	"F5_F8_type_C"    : "PF00754",
	"Fz"              : "PF01534",
	"Kringle"         : "PF00051",
	"Ldl_recept_a"    : "PF00057",
	"MAM"             : "PF00629",
	"Ephrin_lbd"   :"PF01404",
	"WIF": "PF02019",
	"Frizzled": "PF01534"

}
types ={
	"Furin-likeRecep_L_domain"      : "Type 1 (EGF receptor subfamily)",
	"Furin-likeRecep_L_domainfn3"   : "Type 2 (Insulin receptor subfamily)",
	"LRR_8ig"                       : "Type 7 (Trk subfamily)",
	"LRR_8"                         : "Type 7 (Trk subfamily)",
	"FzKringleig"                   : "Type 8 (ROR subfamily)",
	"FzKringle"                     : "Type 8 (ROR subfamily)",
	"Fzig"                          : "Type 9 (MuSK subfamily)",
	"SemaTIG"                       : "Type 10 (HGF receptor subfamily)",
	"Ephrin_lbdfn3"                 : "Type 13 (Eph subfamily)",
	"Cadherin"                      : "Type 14 (RET subfamily)",
	"WIF"                           : "Type 15 (RYK subfamily)",
	"F5_F8_type_C"                  : "Type 16 (DDR subfamily)",
	"fn3"                           : "Type 17 (ROS subfamily)",
	"fn3ig"                         : "Type 11 (TAM subfamily)",
	"fn3Ig_2"                       : "Type 11 (TAM subfamily)",
	"fn3Ig_3"                       : "Type 11 (TAM subfamily)",
	"MAM"                           : "Type 19 (ALK subfamily)",
	"Ldl_recept_aMAM"               : "Type 19 (ALK subfamily)",
	"Ephrin_lbd"                 : "Type 13 (Eph subfamily)"
}

types_jm = {
"JM_12" : "Type 12 (TIE receptor subfamily)",
"JM_3" : "Type 3 (PDGF receptor subfamily)",
"JM_4" : "Type 4 (VEGF receptor subfamily)",
"JM_5" : "Type 5 (FGF receptor subfamily)",
"JM_6" : "Type 6 (CCK receptor subfamily)"
}


def classify_RTKs(EC="EC.res", JM="JM.res", rtk_results=None):
    ec=parse_hmmer_result(EC, "hmmscan")
    jm=parse_hmmer_result(JM, "hmmscan")
    proteins=[]
    family=None
    for h in ec:
        if h.title not in proteins:
            proteins.append(h.title)
    for p in proteins:
        domains=[]
        for h in ec:
            if h.title==p:
                if h.domain not in domains:
                    domains.append(h.domain)
        d=''.join(sorted(domains))
        if d in types.keys():
            family = "Uncategorized"
            family=types[d]
        #else:
        for j in jm:
                if j.title==p:
                    if d.lower() =="ig":
                        family=types_jm[j.domain]
                    elif d.lower()=="fn3ig":
                        if j.domain=="JM_12":
                            family=types_jm[j.domain]
                        else:
                            family="Type 11 (TAM receptor subfamily)"
                    else:
                        continue
        for rtk in rtk_results:
            if rtk.id==p:
                rtk.family=family
                for h in ec:
                    if h.title==p and h.score>0:
                        dct={"domain":h.domain, "id":pfam[h.domain], "score":h.score, "start":h.start, "end":h.end}
                        rtk.domains.append(dct)
    return rtk_results
##### Setting up variable values according to CLI input #####

web=False


if(args.output) and args.output!=".":
    output=args.output
    #os.mkdir(output)
else:
    output="."

if args.mkdir and output!=".":
    if not os.path.isdir(output):
        os.mkdir(output)

if(args.webphobius):
    web=True
else:
    web=False

if args.hmmerdir:
    hmmerdir=args.hmmerdir
else:
    hmmerdir=cnf.hmmer_dir

if args.phobiusdir:
    phobiusdir=args.phobiusdir
else:
    phobiusdir=cnf.phobius_dir

input=args.input



##### Setting up environmental values according to platform

system_os=platform.system()

if system_os=="Windows":
    hmmerdir=hmmerdir.replace("/", "\\")
    phobiusdir=phobiusdir.replace("/","\\")
    PTKhmm = cnf.PTKhmm_dir.replace("/","\\")+"\\PTK.hmm"
    JM = cnf.JM_dir.replace("/","\\")+"\\JM.hmm"
    Pfam = cnf.Pfam_dir.replace("/","\\")+"\\EC.hmm"
    hmmsearch="%s\\hmmsearch" %hmmerdir
    hmmscan="%s\\hmmscan" %hmmerdir
    phobius="%s\\phobius.pl" %phobiusdir #this may be moot, since Phobius's decodeanhmm has no Windows version
else:
    PTKhmm = cnf.PTKhmm_dir+"/PTK.hmm"
    JM = cnf.JM_dir+"/JM.hmm"
    Pfam = cnf.Pfam_dir+"/EC.hmm"
    hmmsearch="%s/hmmsearch" %hmmerdir
    hmmscan="%s/hmmscan" %hmmerdir
    phobius="%s/phobius.pl" %phobiusdir



##### Actually running the RTK-PRED pipeline #####

### STEP 0.  Saving the Input in a MultiFasta object and a dictionary which corresponds the headers of the fasta reads to their sequences.
# This will be useful THROUGHOUT the pipeline process
fasta_dictionary=readMultiFasta(input)
mfasta = fasta_dictionary[1]
fasta_sequences_dict = fasta_dictionary[2]

### STEP 1. Identifying PTKs ###

#1.1 Running HMMSEARCH against the PTK pHMM
run_hmmer(hmmsearch, PTKhmm, input, "PTKhmm.res", gathering=True)
#1.2 Parsing results and getting the PTKs
ptk_results=parse_hmmer_result("%s/PTKhmm.res" %output, "hmmsearch")
#1.3 Isolating the PTK sequences that were identified 
PTKs=list()

for ptk in ptk_results:
    if " " in ptk.title:
        ptk_title = ptk.title.split()[0]
    else:
        ptk_title = ptk.title
    ptkfa=Fasta()
    ptkfa.title=ptk_title
    ptkfa.sequence=fasta_sequences_dict[ptk_title]
    PTKs.append(ptkfa)
writeMFasta(PTKs, "%s/PTKs.fasta" %output)

#1.4 Isolating only single-domain PTKs for the next step
domain_no=dict()
for ptk in ptk_results:
    if " " in ptk.title:
        title=ptk.title.split()[0]
    else:
        title=ptk.title
    if title not in domain_no.keys():
        domain_no[title]=1
    else:
        domain_no[title]+=1
singlePTKs=list()
for fasta in mfasta:
    if " " in fasta.title:
        title=fasta.title.split()[0]
    else:
        title=fasta.title
    if title in domain_no.keys():
        if domain_no[title]==1:
            singlePTKs.append(fasta)
writeMFasta(singlePTKs, "%s/singlePTKs.fasta" %output)

### STEP 2.  Identifying RTKs ###

# This is for not getting an error by Phobius in case of sending empty files for analysis.
single_ptks_file_path = output + "/singlePTKs.fasta"
if os.stat(single_ptks_file_path).st_size == 0:
    print("No PTKs with only one tyrosine kinase domain found. No RTKs found.\n")
    rtk_results_classified = []
else:
    #2.1 Running Phobius
    #2.1.1 checking whether local or web-based Phobius is used and formatting the input
    if web==True:
        #phobius_input=readFile("%s/PTKs.fasta" %output, list=False)
        phobius_input=readFile("%s/singlePTKs.fasta" %output, list=False)
    else:
        #phobius_input="%s/PTKs.fasta" %output
        phobius_input="%s/singlePTKs.fasta" %output
    #2.2.2 Actually running Phobius
    run_phobius(phobius_input, output, web=web)
    #2.2 Parsing results and isolating RTKs
    phobius_results=parse_phobius_result("%s/phobius.txt" %output)
    rtk_results=get_rtks_from_phobius(phobius_results)
    #2.3 creating a subset of the MultiFasta input with only RTKs for the next step
    RTKs=list()

    for rtk in rtk_results:
        if " " in rtk.id:
            rtk_id = rtk.id.split()[0]
        else:
            rtk_id = rtk.id
        rtkfa=Fasta()
        rtkfa.title=rtk_id
        rtkfa.sequence=fasta_sequences_dict[rtk_id]
        RTKs.append(rtkfa)
    writeMFasta(RTKs, "%s/RTKs.fasta" %output)

    ### STEP 3. Functional annotation / classification of RTKs ###

    #3.1 Run HMMSCAN on RTKs against pHMM libraries
    #3.1.1 Running against the Extracellular domains library, unless the file with the RTKs is empty.
    rtks_file_path = output + "/RTKs.fasta"
    if os.stat(rtks_file_path).st_size == 0:
        print("No RTKs found.\n")
        rtk_results_classified = []
    else:
        run_hmmer(hmmscan, Pfam, "%s/RTKs.fasta" %output, "EC.res", gathering=False)
        #3.1.2 Running against the Juxtamembrane domains libtrary, for additional classification
        # of more 'difficult' sub-families
        run_hmmer(hmmscan, JM, "%s/RTKs.fasta" %output, "JM.res", gathering=True)
        #3.2 Parse results, add extracelular domains and annotate RTKs
        rtk_results_classified=classify_RTKs(EC="%s/EC.res" %output, JM="%s/JM.res" %output, rtk_results=rtk_results)


### STEP 4. Printing the Summary ###
ptk_dict=hmm_dict(ptk_results)
rtk_dict=rtk_dict(rtk_results_classified)
summary=open("%s/summary.txt" %output,"w")
for fasta in mfasta:
    nrtk=False
    rtk=False
    if " " in fasta.title:
        if fasta.title.split()[0] in ptk_dict.keys():
            if fasta.title.split()[0] in rtk_dict.keys():
                rtk=True
            else:
                nrtk=True
        id=fasta.title.split()[0]
    else:
        if fasta.title in ptk_dict.keys():
            if fasta.title.split()[0] in rtk_dict.keys():
                rtk=True
            else:
                nrtk=True
        id=fasta.title

    # If the phrase "Unique_" is found in the final fasta title this phrase and the following number and "_" are removed from the title.
    if id[:7] == "Unique_":
        splited = id.split("_")
        final_title = "_".join(splited[2:])
    else:
        final_title = id

    summary.write(">> Query: %s\n" %final_title)
    if rtk==True:
        summary.write("Classification: Receptor Tyrosine Kinase (RTK)\n")
        if rtk_dict[id].signal !=None:
            summary.write("SIGNALP:\t\t\t\tFrom:%d\tTo:%d\n"%(rtk_dict[id].signal[0], rtk_dict[id].signal[1]))
        domains=sorted(rtk_dict[id].domains, key = lambda i: i['start'])
        for dom in domains:
            summary.write("EC Domains (%s):\tFrom:%d\tTo:%d\t%s\t\tScore: %.1f\n" %(dom['id'], dom['start'], dom['end'], dom['domain'], dom['score']))
        summary.write("TRANSMEM:\t\t\t\tFrom:%d\tTo:%d\n"%(rtk_dict[id].tm[0], rtk_dict[id].tm[1]))
        summary.write("Kinase Domain:\t\t\t\tFrom:%d\tTo:%d\t\tScore: %.1f\n" %(ptk_dict[id].start, ptk_dict[id].end, ptk_dict[id].score))
        if rtk_dict[id].family != "":
            summary.write("RTK subfamily: %s\n" %rtk_dict[id].family)
        else:
            if rtk_dict[id].signal !=None:
                ec_reg=rtk_dict[id].tm[0] -rtk_dict[id].signal[1]+1
            else:
                ec_reg=rtk_dict[id].tm[0]
            if ec_reg <= 50:
                summary.write("RTK subfamily: Type 18 (LMR) or Type 20 (STYK1) subfamily\n")
            else:
                summary.write("RTK subfamily: Uncategorized\n")
    elif rtk==False and nrtk==True:
        summary.write("Classification: Non-Receptor Tyrosine Kinase (nRTK)\n")
        for ptk in ptk_results:
            if " " in ptk.title:
                ptitle=ptk.title.split()[0]
            else:
                ptitlre=ptk.title
            if ptitle==id:
                summary.write("Kinase Domain:\t\t\t\tFrom:%d\tTo:%d\t\tScore: %.1f\n" %(ptk.start, ptk.end, ptk.score))
    else:
        summary.write("Classification: Not Tyrosine Kinase\n")
    summary.write("//\n")
summary.close()

