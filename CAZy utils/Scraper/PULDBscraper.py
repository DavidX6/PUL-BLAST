from requests import get
from requests.exceptions import RequestException
from contextlib import closing
from bs4 import BeautifulSoup
import re
import os

class PUL: pass
class PULmember: pass

def simple_get(url):
    """Make a HTTP GET request and return response"""
    try:
        with closing(get(url, stream=True)) as resp:
            if is_good_response(resp):
                return {"content":resp.content, "id": resp.url.split("=")[1]}
            else:
                return None

    except RequestException as e:
        print(e)
        return None

def is_good_response(resp):
    """Check if response status is 200"""
    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200
            and content_type is not None
            and content_type.find('html') > -1)

def extractPUL(url):
    """From HTML displaying PUL extract links to members, modularity and substrate"""
    response = simple_get(url)
    content = BeautifulSoup(response["content"], features="html.parser")
    table = content.select("table")[2]
    rows = table.select("tr")
    newPUL = PUL()
    newPUL.membersLinks = []
    newPUL.modularity = ""
    newPUL.id = response["id"]
    for i in range(0, len(rows)):
        if i == 0: newPUL.name = rows[i].text.split("\n")[2]
        elif i == 2: newPUL.modularity = rows[i].select("td")[1].text
        elif i == 4: newPUL.substrate = rows[i].text.split("\n")[2]
        elif i == 5:
            proteinTable = rows[i].select("td")[1].select("table")[0]
            for tr in proteinTable.select("tr"):
                links = tr.findAll("a", attrs={'href': re.compile("^index")})
                if len(links) > 0: newPUL.membersLinks.append(links[0].get("href"))
    return newPUL

def extractPULmembers(newPUL):
    """With links to PUL members make HTTP GET requests and extract
    their information"""
    url = "http://www.cazy.org/PULDB/"
    newPUL.members = []
    for link in newPUL.membersLinks:
        print("Processing " + url + link)
        response = simple_get(url + link)
        content = BeautifulSoup(response["content"], features="html.parser")
        tempDesc = ""
        tempSeq = ""
        rows = content.findAll("tr", attrs={'align': "left"})
        for row in rows:
            temp = row.text.split("\n")
            if len(temp) != 4 and temp[1] != "CAZy links" and temp[1] != "Demonstrated function":
                print(temp)
                print(row)
            elif len(temp) > 4 and temp[1] == "Demonstrated function":
                for i in range(2, len(temp)-1):
                    tempDesc += temp[i].replace("▶", "").replace("◀", "").replace("\xa0", "").replace("[EC]", "").replace("|", " ") + ", "
            elif (temp[1] == "Organism name/NCBI id" or temp[1] == "Locus tag / Protein name"):
                tempDesc += temp[2] + "|"
            elif (temp[1] == "IMG/M-HMP description" or temp[1] == "Modularity"
                  or temp[1] == "Demonstrated function"):
                tempDesc += temp[2].replace("▶", "").replace("◀", "").replace("\xa0", "").replace("[EC]", "").replace("|", " ") + ", "
            elif temp[1] == "Amino-acid sequence":
                tempSeq += temp[2]
                tempDesc = tempDesc[:-2]
                break
        if tempSeq == "":
            print("Warning!")
        newPUL.members.append({"description": tempDesc, "sequence": tempSeq})
    return newPUL


def writeToFile(myPUL, substrate):
    """Write out compiled PUL in fasta format"""
    f = open(substrate + "_" + myPUL.id + ".fasta", "w", encoding="utf-8")
    for member in myPUL.members:
        f.write(">" + member["description"] + "|" + substrate + "\n")
        f.write(member["sequence"] + "\n")
        f.write("\n")
    f.close()
    f = open("modularityDescription.txt", "a", encoding="utf-8")
    f.write(myPUL.id + " " + substrate + "\n")
    f.write(myPUL.modularity + "\n")
    f.close()

def validateFiles():
    """Check if all fasta files are in correct format"""
    for file in os.listdir():
        if file.endswith(".fasta"):
            f = open(file, "r", encoding="utf-8")
            for line in f:
                if line[0] == ">" and len(line.split("|")) != 4:
                    print(file, line)
                    return False
            f.close()
    return True

def mergeFiles():
    """Merge all fasta files in one file"""
    merged = open("PULDB_merged.fasta", "w", encoding="utf-8")
    for file in os.listdir():
        print(file)
        if file.endswith(".fasta") and file != "PULDB_merged.fasta":
            f = open(file, "r", encoding="utf-8")
            for line in f: merged.write(line)
            f.close()
    merged.close()

# validateFiles()
# mergeFiles()


# ----------------------
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=293")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Galactomannan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32782")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Arabinogalactan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=25076")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Agar")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=25078")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Agar_Porphyran")
# #
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=33613")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Xyloglucan[33613]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=364")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Xyloglucan[364]")
# #
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=33608")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Homogalacturonan_Rhamnogalacturonan-I")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32854")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Rhamnogalacturonan-I")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32849")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "High mannose mammalian N-glycan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32868")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Yeast alpha-mannan[32868]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32846")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Yeast alpha-mannan[32846]")
# #
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32869")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Rhamnogalacturonan-II[32869]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32870")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Rhamnogalacturonan-II[32870]")
# #
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=33573")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Complex xylans")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=33593")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Simple xylans")
# #
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32799")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Levan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=33551")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Beta-glucan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32843")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Starches")

# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=269")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Barley beta-glucan[269]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=324")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Wheat arabinoxylan_oat spelt xylan[324]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=351")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Homogalacturonan")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=344")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Wheat arabinoxylan_oat spelt xylan[344]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=302")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Barley beta-glucan[302]")
#
# myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32848")
# myPUL = extractPULmembers(myPUL)
# writeToFile(myPUL, "Mucin O-glycans")

