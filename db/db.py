from Bio import Entrez
import xml.etree.ElementTree as ET
import time
from pprint import pprint
import math

QUERY = "((epilepsy[MeSH Major Topic]) AND (2000/1/1[ppdat] : 2014/12/31[ppdat]))"
PUBMED_DB = "pubmed"
MAX_DOWNLOADABLE_RECORDS = 10000
Entrez.email = "qtw1vvk8ut3cjrp@jetable.org"  # e-mail address only valid until 2019-01-07

# Search the DB and get the total number of records
handle = Entrez.esearch(db=PUBMED_DB, term=QUERY, retmax=0)
records = Entrez.read(handle)
RecordsCount = int(records["Count"])

# Get all identifiers
handle = Entrez.esearch(db=PUBMED_DB, term=QUERY, retmax=RecordsCount)
records = Entrez.read(handle)
identifiers = records['IdList']
print("Number of records found: {}".format(len(identifiers)))

# Workaround the max number of downloadable records
id_list = []
for i in range(math.ceil(RecordsCount / MAX_DOWNLOADABLE_RECORDS)):
    id_list.append(identifiers[i * MAX_DOWNLOADABLE_RECORDS:(
        i + 1) * MAX_DOWNLOADABLE_RECORDS])

# Get the data
for idl in id_list:
    handle = Entrez.efetch(
        db=PUBMED_DB, id=idl, retmax=RecordsCount, retmode="xml")
    xml_text = handle.read()

    # Parse the data
    root = ET.fromstring(xml_text)
    print("Number of downloaded articles data: {}".format(len(root)))
    for PubmedArticle in root.iter("PubmedArticle"):
        PMID = PubmedArticle.find("MedlineCitation").find("PMID").text
        ArticleTitle = PubmedArticle.find(".//ArticleTitle").text
        PublicationDate = PubmedArticle.find(".//PubDate").find("MedlineDate").text
        Authors = PubmedArticle.find(".//AuthorList").findall("Author")
        for auth in Authors:
            LastName = auth.find("LastName").text
            ForeName = auth.find("ForeName").text
            ForeName = auth.find("ForeName").text
            Initials = auth.find("Initials").text


handle.close()
