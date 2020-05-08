from requests import get
from requests.exceptions import RequestException
from contextlib import closing
from bs4 import BeautifulSoup
import re

class PUL: pass
class PULmember: pass

def simple_get(url):
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
    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200
            and content_type is not None
            and content_type.find('html') > -1)

def extractPUL(url):
    response = simple_get(url)
    content = BeautifulSoup(response["content"], features="html.parser")
    table = content.select("table")[2]
    rows = table.select("tr")
    newPUL = PUL()
    newPUL.membersLinks = []
    newPUL.id = response["id"]
    for i in range(0, len(rows)):
        if i == 0: newPUL.name = rows[i].text.split("\n")[2]
        elif i == 4: newPUL.substrate = rows[i].text.split("\n")[2]
        elif i == 5:
            proteinTable = rows[i].select("td")[1].select("table")[0]
            for tr in proteinTable.select("tr"):
                links = tr.findAll("a", attrs={'href': re.compile("^index")})
                if len(links) > 0: newPUL.membersLinks.append(links[0].get("href"))
    return newPUL

def extractPULmembers(newPUL):
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
                tempDesc = tempDesc[:-2] + "|"
                break
        if tempSeq == "":
            print("Warning!")
        newPUL.members.append({"description": tempDesc, "sequence": tempSeq})
    return newPUL


def writeToFile(myPUL, substrate):
    f = open(substrate + "_" + myPUL.id + ".fasta", "w")
    for member in myPUL.members:
        f.write(">" + member["description"] + "|" + substrate + "\n")
        f.write(member["sequence"] + "\n")
        f.write("\n")
    f.close()


myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=32782")
myPUL = extractPULmembers(myPUL)
writeToFile(myPUL, "Arabinogalactan")

myPUL = extractPUL("http://www.cazy.org/PULDB/index.php?pul=293")
myPUL = extractPULmembers(myPUL)
writeToFile(myPUL, "Galactomannan")


