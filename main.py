import Bio.Blast.Applications
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
import copy

db_name = "clankiDB"
"""
genome search
"""
def gbkGenomeSearch(fileName):
    inputSequenceObject = Bio.SeqIO.read(fileName, "genbank")
    file = open("queryTempSequence.txt", "w")
    cnt = 0
    for feature in inputSequenceObject.features:
        if feature.type != "CDS" or len(feature.extract(inputSequenceObject).seq) == len(inputSequenceObject.seq): continue
        description = "Query_" + str(cnt) + "|" + "location: " + str(feature.location)
        cnt += 1
        hasTrans = False
        for key in feature.qualifiers.keys():
            if key == "translation":
                hasTrans = True
                continue
            description += "|" + str(key) + ": " + str(feature.qualifiers[key]).replace("[", "").replace("]", "")
        description += str(feature.type)

        if hasTrans:
            seq = str(feature.qualifiers["translation"][0])
        else:
            seq = str(feature.extract(inputSequenceObject).seq.transcribe().translate(table="Bacterial")).replace("*", "")
        file.writelines(">" + description + "\n")
        file.writelines(seq)
        file.writelines("\n")
    file.close()


def genomeSearch(fileName, cutoff = 140):
    # read sequence from file
    inputSequenceObject = Bio.SeqIO.read(fileName, "fasta", IUPAC.unambiguous_dna)
    sequence = inputSequenceObject.seq
    # RNA from DNA
    sequenceRNA = sequence.transcribe()
    # protein from RNA
    sequenceProt = sequenceRNA.translate(table="Bacterial")  # stop_symbol="@"
    # write out possible proteins, blast searches with every entry
    toSearchSequence = sequenceProt.split("*")
    file = open("queryTempSequence.txt", "w")
    cnt = 0 # tracking order of writing
    for i in range(0, len(toSearchSequence)):
        item = toSearchSequence[i]
        if len(item) > cutoff:
            file.writelines(">myQuery_" + str(cnt) + "_index_" + str(i) + "\n")
            file.writelines(item)
            file.writelines("\n")
            cnt += 1
    file.close()
"""
BLAST
"""
def proteinBLAST(proteinFile, eval = 0.001, format = 5, outfile = "resultsPyP.txt"):
    blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(
        evalue = eval,
        query = proteinFile,
        db = db_name,
        outfmt = format,
        out = outfile
    )
    # print(blastp_cline)
    stdout, stderr = blastp_cline()


def nucleotideBLAST(geneFile, eval = 0.001, format = 5, outfile = "resultsPyN.txt"):
    blastn_cline = Bio.Blast.Applications.NcbiblastxCommandline(
        evalue = eval,
        query = geneFile,
        db = db_name,
        outfmt = format,
        out = outfile
    )
    # print(blastn_cline)
    stdout, stderr = blastn_cline()
"""
DATA PROCESSING
"""

def makePossiblePULs(blast_records):
    possiblePULs = set()
    for i in range(0, len(blast_records)):
        blast_records[i].originalIndex = i
        if len(blast_records[i].alignments) == 0: continue
        alignment = blast_records[i].alignments[0]  # SusCD should be best hit
        queryCover = 0
        for hsp in alignment.hsps: queryCover += hsp.query_end - hsp.query_start
        if queryCover < alignment.length / 2: continue
        if i + 1 == len(blast_records): break
        if "susc" in alignment.hit_def.split("|")[2].lower():
            for alignment2 in blast_records[i + 1].alignments:
                if "susd" in alignment2.hit_def.split("|")[2].lower(): possiblePULs.add((i, i + 1))
        elif "susd" in alignment.hit_def.split("|")[2].lower():
            for alignment2 in blast_records[i + 1].alignments:
                if "susc" in alignment2.hit_def.split("|")[2].lower(): possiblePULs.add((i, i + 1))
    return possiblePULs

def makePULBySubstrate(substrate, possiblePULs, blast_records, maxDist):
    foundPULs = []  # [(Sus, Sus), (lowI, highI)]
    for candidate in possiblePULs:
        borderLow = min(candidate)
        borderHigh = max(candidate)
        # to increase border, same substrate hit must be next to it
        while borderLow - 1 != 0:
            # if distance from SusCD is now larger than allowed
            if abs(blast_records[borderLow - 1].originalIndex - blast_records[
                min(candidate)].originalIndex) > maxDist: break
            # if we found next SusCD
            check = True
            if len(blast_records[borderLow - 1].alignments) > 0:
                for alignment in blast_records[borderLow - 1].alignments:
                    if "susc" in alignment.hit_def.split("|")[2].lower() or "susd" in alignment.hit_def.split("|")[2].lower():
                        check = False
                        break
            if check == False: break
            # if substrate is not in new hits
            substratesInAlignments = [alignment.hit_def.split("|")[-1] for alignment in blast_records[borderLow-1].alignments]
            if substrate not in substratesInAlignments and len(substratesInAlignments) > 0: break
            borderLow = borderLow - 1

        while borderHigh + 1 != len(blast_records):
            if abs(blast_records[borderHigh + 1].originalIndex - blast_records[
                max(candidate)].originalIndex) > maxDist: break
            check = True
            if len(blast_records[borderHigh + 1].alignments) > 0:
                for alignment in blast_records[borderHigh + 1].alignments:
                    if "susc" in alignment.hit_def.split("|")[2].lower() or "susd" in alignment.hit_def.split("|")[2].lower():
                        check = False
                        break
            if check == False: break
            substratesInAlignments = [alignment.hit_def.split("|")[-1] for alignment in
                                     blast_records[borderHigh + 1].alignments]
            if substrate not in substratesInAlignments and len(substratesInAlignments) > 0: break
            borderHigh = borderHigh + 1
        if borderHigh == max(candidate) and borderLow == min(candidate): continue
        if borderHigh - borderLow > 2: foundPULs.append([(candidate), (borderLow, borderHigh)])
    return foundPULs

def resultsBLASTwrite(substrate, maxDist = 10):
    # open results file
    blast_records = list(NCBIXML.parse(open("resultsPyP.txt")))
    # make PULs
    possiblePULs = makePossiblePULs(blast_records)
    # widen until SAME kind of substrate hit can be found in ANY hits
    foundPULs = []  # [(Sus, Sus), (lowI, highI)]
    for candidate in possiblePULs:
        borderLow = min(candidate)
        borderHigh = max(candidate)
        # to increase border, same substrate hit must be next to it
        while borderLow - 1 != 0:
            # if distance from SusCD is now larger than allowed
            if abs(blast_records[borderLow - 1].originalIndex - blast_records[
                min(candidate)].originalIndex) > maxDist: break
            # if we found next SusCD
            check = True
            if len(blast_records[borderLow - 1].alignments) > 0:
                for alignment in blast_records[borderLow - 1].alignments:
                    if "susc" in alignment.hit_def.split("|")[2].lower() or "susd" in alignment.hit_def.split("|")[2].lower():
                        check = False
                        break
            if check == False: break
            # if substrate is not in new hits
            substratesInAlignments = [alignment.hit_def.split("|")[-1] for alignment in blast_records[borderLow-1].alignments]
            if substrate not in substratesInAlignments and len(substratesInAlignments) > 0: break
            borderLow = borderLow - 1

        while borderHigh + 1 != len(blast_records):
            if abs(blast_records[borderHigh + 1].originalIndex - blast_records[
                max(candidate)].originalIndex) > maxDist: break
            check = True
            if len(blast_records[borderHigh + 1].alignments) > 0:
                for alignment in blast_records[borderHigh + 1].alignments:
                    if "susc" in alignment.hit_def.split("|")[2].lower() or "susd" in alignment.hit_def.split("|")[2].lower():
                        check = False
                        break
            if check == False: break
            substratesInAlignments = [alignment.hit_def.split("|")[-1] for alignment in
                                     blast_records[borderHigh + 1].alignments]
            if substrate not in substratesInAlignments and len(substratesInAlignments) > 0: break
            borderHigh = borderHigh + 1
        if borderHigh == max(candidate) and borderLow == min(candidate): continue
        foundPULs.append([(candidate), (borderLow, borderHigh)])
    print("foundPULs", foundPULs)

    PULrecords = []
    for pul in foundPULs:
        borderLow = pul[1][0]
        borderHigh = pul[1][1]
        # all hits should have the same substrate, except SusCD hits
        temp = []
        for i in range(borderLow, borderHigh + 1):
            newAls = []
            #print(i, ":", len(blast_records[i].alignments), end=" ,")
            if len(blast_records[i].alignments) == 0: continue
            if i in pul[0]: newAls.append(blast_records[i].alignments[0])
            else:
                for alignment in blast_records[i].alignments:
                    if substrate in alignment.hit_def.split("|")[-1]:
                        newAls.append(alignment)
                        break
            blast_records[i].alignments = newAls
            temp.append(blast_records[i])
        PULrecords.append(temp)
        #print("\n")
    return PULrecords

def searchPULs(maxDist = 10):
    blast_records = list(NCBIXML.parse(open("resultsPyP.txt")))
    temp = list()
    allSubstrates = set()
    cnt = 0
    # get only items with results, compute coverage
    for blast_record in blast_records:
        if len(blast_record.alignments) > 0:
            queryCover = 0
            newAlignments = []
            for alignment in blast_record.alignments:
                allSubstrates.add(alignment.hit_def.split("|")[3])
                for hsp in alignment.hsps: queryCover += hsp.query_end - hsp.query_start
                alignment.queryCover = queryCover/blast_record.query_length
                if alignment.queryCover > 0.5: newAlignments.append(alignment)
            blast_record.alignments = newAlignments
            if len(blast_record.alignments) > 0: temp.append(blast_record)
    blast_records = temp
    # get indexes of SusCD/DC pairs
    possiblePULs = makePossiblePULs(blast_records)
    # widen window around SusCD pairs until limit or next Sus
    foundPULs = {}
    for substrate in allSubstrates:
        foundPULs[substrate] = makePULBySubstrate(substrate, possiblePULs, blast_records, maxDist)

    # foundPULs = [] [(Sus, Sus), (lowI, highI)]
    # apply selection criteria to found PULs
    PULrecords = {}
    for substrate in allSubstrates:
        for pul in foundPULs[substrate]:
            borderLow = pul[1][0]
            borderHigh = pul[1][1]
            temp = [copy.deepcopy(blast_records[i]) for i in range(borderLow, borderHigh + 1)]
            # all hits should have the same substrate, except SusCD hits
            for i in range(borderLow, borderHigh + 1):
                if i in pul[0]:
                    temp[i - borderLow].alignments = [temp[i - borderLow].alignments[0]]
                else:  # substrate not from susCD hits
                    substrateAligment = None
                    for alignment in blast_records[i].alignments:
                        if substrate == alignment.hit_def.split("|")[-1]:
                            substrateAligment = alignment
                            break
                    if substrateAligment == None:
                        print("Warning! No suitable alignments found when making PUL.")
                    temp[i-borderLow].alignments = [substrateAligment]
            if substrate not in PULrecords.keys(): PULrecords[substrate] = [temp]
            else: PULrecords[substrate].append(temp)
    return PULrecords


def countSubstrateExamples():
    #file = open("izClankovZAccZimeni.fasta")
    file = open("PULDB_merged.fasta")
    cnt = 0
    examples = {}
    for line in file:
        if line[0] == ">":
            cnt += 1
            substrate = line.split("|")[3].replace("\n", "")
            if substrate in examples.keys(): examples[substrate] += 1
            else: examples[substrate] = 1
    print(examples)
    for key in examples.keys(): print(str(key) + ": " + str(examples[key]))
    print(cnt)
    print(sum([examples[key] for key in examples.keys()]))





if __name__ == '__main__':
    #resultsBLASTwrite("Starches")
    #gbkGenomeSearch("prevotele_iz_članka_tabelaPULs\\Prevotella_ruminicola_23.gbk")
    searchPULs()
    #countSubstrateExamples()
    #gbkGenomeSearch("prevotele_iz_članka_tabelaPULs\\Prevotella_ruminicola_23.gbk")
    #gbkGenomeSearch("prevotele_iz_članka_tabelaPULs\Prevotella_sp._AGR2160.gbkU")
    # genomeSearch("B_ovatusGenome.fna")
    #genomeSearch("B_ceccaeGenome.fasta")

    #proteinBLAST("queryTemp.txt")

    #proteinBLAST("queryTempSequence.txt")
    # resultsBLASTwrite("protein")



# outfmt
     # 0 = pairwise,
     # 1 = query-anchored showing identities,
     # 2 = query-anchored no identities,
     # 3 = flat query-anchored, show identities,
     # 4 = flat query-anchored, no identities,
     # 5 = XML Blast output,
     # 6 = tabular,
     # 7 = tabular with comment lines,
     # 8 = Text ASN.1,
     # 9 = Binary ASN.1,
    # 10 = Comma-separated values,
    # 11 = BLAST archive format (ASN.1) 
