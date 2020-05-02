import urllib.parse
import urllib.request
import xmlschema
import json
from Bio import Entrez, SeqIO
from Bio import ExPASy
from Bio import SwissProt


db_name = "LinCAZyp"
# AccNumTest = "AWI10243.1"
AccNumTest = "4V1S"

# Get all accession numbers from fasta database
def extactAccNumDBCAZy():
    numbers = []
    file = open("CAZyDB.07312019.fa", "r")
    for line in file:
        if line[0] == ">":
            newline = line.split("|")[0]
            numbers.append(newline[1:])
    file.close()

    file = open("AsNumCAZy.txt", "w")
    for num in numbers:
        file.writelines(num + "\n")
    file.close()

# Get all GI numbers from fasta file
def extactGINCBI():
    numbers = []
    file = open("izClankov.fasta", "r")
    cnt = 0
    for line in file:
        if line[0] == ">":
            cnt += 1
            newline = line.split("|")[1]
            numbers.append(newline)
    file.close()

    print(cnt)
    file = open("GINumNCBI.txt", "w")
    for num in numbers:
        file.writelines(num + "\n")
    file.close()
    print(numbers, len(numbers))

# Count number of entries in fasta files from NCBI
def justcnt():
    cnt = 0
    for x in range(1, 14):
        file = open("IRP BF članki Accetto/encimi raznih PULov/" + str(x) + ".fasta", "r")
        for line in file:
            if line[0] == ">":
                cnt += 1
        file.close()
    print(cnt)

# merge AccNum with NCBI fasta
def addAccNum():
    file = open("izClankovAccN.acc", "r")
    numbers = []
    for line in file:
        numbers.append(line.strip("\n"))
    file.close()

    file = open("izClankov.fasta", "r")
    file2 = open("izClankovZAcc.fasta", "w")
    cnt = 0
    for line in file:
        if line[0] == ">":
            file2.writelines(line.strip("\n") + "|" + numbers[cnt] + "\n")
            cnt += 1
        else:
            file2.writelines(line)

    file.close()
    file2.close()
    print(cnt)

# merge organism name with fasta
def addName():
    file = open("izClankovTaxonomy.txt", "r")
    taxa = json.load(file)
    names = []
    for key in taxa.keys():
        names.append(taxa[key]["name"])
    print(names)
    file.close()

    file = open("izClankovZAcc.fasta", "r")
    file2 = open("izClankovZAccZimeni.fasta", "w")
    cnt = 0
    for line in file:
        if line[0] == ">":
            nline = line.strip("\n").split("|")
            file2.writelines(">" + nline[3] + "|" + names[cnt] + "|" + nline[4] + "|" + nline[5] + "\n")
            cnt += 1
        else:
            file2.writelines(line)

    file.close()
    file2.close()
    print(cnt)

"""
EMBL
"""
Entrez.email = 'A.N.Other@example.com'

# Get complete record from EMBL
def EMBLcomplete(accNum, type = "protein"):
    handle = Entrez.efetch(db=type, id=accNum, rettype="fasta")
    record = SeqIO.read(handle, "fasta")
    print(record)
    print(record.description)
    return record

# Get only summary from EMBL
def EMBLsummary(accNum, type = "protein"):
    handle = Entrez.esummary(db=type, id=accNum)
    record = Entrez.read(handle)
    print(record)
    print(record[0]["Title"])
    return record[0]["Title"]

def EMBLGItoAccNum():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
    gis = ""
    file = open("GINumNCBI.txt", "r")
    for line in file:
        gis = gis + line.strip("\n") + ","
    file.close()
    gis = gis[:-1]
    print(gis)
    url = url + gis + "&rettype=acc"
    print(url)
"""
UniProt
"""
# Convert accesion number to UniProt ID
def getUniProtID(accNum):
    try:
        url = 'https://www.uniprot.org/uploadlists/'
        params = {
            'from': 'EMBL',
            'to': 'ID',
            'format': 'tab',
            'query': accNum
        }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        resText = response.decode('utf-8')
        GenBankID = resText.split("\n")[1].split("\t")[1]
        print(GenBankID)
        return GenBankID
    except:
        print("Error in function getUniProtID!")

# get uniProt xml and parse it
def UniProtData(GenBankID):
    url = "https://www.uniprot.org/uniprot/" + GenBankID + ".xml"
    req = urllib.request.Request(url)
    response = urllib.request.urlopen(req).read()

    schema = xmlschema.XMLSchema('https://www.uniprot.org/docs/uniprot.xsd')
    xmlText = response.decode('utf-8')
    entry_dict = schema.to_dict(xmlText)
    # print(entry_dict)
    uniProtContent = entry_dict["entry"][0]
    print(uniProtContent["accession"])
    print(uniProtContent["protein"])
    print(uniProtContent["gene"])
    print(uniProtContent["organism"])
    print(uniProtContent["comment"])
    print(uniProtContent["dbReference"])
    print(uniProtContent["proteinExistence"])
    print(uniProtContent["keyword"])

# --------------------------------------------

#upID = getUniProtID(AccNumTest)
#UniProtData(upID)

# extactGINCBI()
# EMBLGItoAccNum()

# addAccNum()
# addName()





# file = open("izClankovAccN.acc", "r")
# s = ""
# for x in file:
#     s += x.strip("\n") + ","
# print(s)