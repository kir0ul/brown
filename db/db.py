from Bio import Entrez
import xml.etree.ElementTree as ET
import math
from db_model import db_session, Article, Author

QUERY = "((epilepsy[MeSH Major Topic]) AND (2000/1/1[ppdat] : 2014/12/31[ppdat]))"
QUERY_DB = "pubmed"


class db():

    max_downloadable_records = 10000

    def __init__(self, query, query_db, email=None):
        self.query = query
        self.query_db = query_db
        if email:
            Entrez.email = email
        else:
            # e-mail address only valid until 2019-01-07
            Entrez.email = "qtw1vvk8ut3cjrp@jetable.org"

    def get_total_records_nb(self):
        """Search the DB and get the total number of records"""
        handle = Entrez.esearch(db=self.query_db, term=self.query, retmax=0)
        records = Entrez.read(handle)
        RecordsCount = int(records["Count"])
        handle.close()
        return RecordsCount

    def get_identifiers_list(self, RecordsCount):
        """Get all identifiers"""
        handle = Entrez.esearch(
            db=self.query_db, term=self.query, retmax=RecordsCount)
        records = Entrez.read(handle)
        identifiers = records['IdList']
        print("Number of records found: {}".format(len(identifiers)))
        handle.close()
        return identifiers

    def parse_and_commit_to_db(self, xml_text):
        """Parse the data"""
        root = ET.fromstring(xml_text)
        print("Chunk of downloaded articles data: {}".format(len(root)))

        with db_session() as session:
            for PubmedArticle in root.iter("PubmedArticle"):
                PMID = PubmedArticle.find("MedlineCitation").find("PMID").text
                ArticleTitle = PubmedArticle.find(".//ArticleTitle").text

                PubDate = PubmedArticle.find(".//PubDate")
                Year = None
                Month = None
                MonthEnd = None
                if PubDate.find("Year") is not None:
                    Year = PubDate.find("Year").text
                if PubDate.find("Month") is not None:
                    Month = PubDate.find("Month").text
                if PubDate.find("Day") is not None:
                    import ipdb; ipdb.set_trace()
                elif PubDate.find("MedlineDate") is not None:
                    MedlineDate = PubDate.find("MedlineDate").text.split(" ")
                    Year = MedlineDate[0]
                    if len(MedlineDate) == 2:
                        Range = MedlineDate[1].split("-")
                        Month = Range[0]
                        MonthEnd = Range[1]
                    if len(MedlineDate) == 3:
                        Month = MedlineDate[1]

                # Add record to DB
                article = Article(
                    PMID=PMID,
                    Title=ArticleTitle,
                    PublicationYear=Year,
                    PublicationMonth=Month,
                    PublicationMonthEnd=MonthEnd)
                session.add(article)

                if PubmedArticle.find(".//AuthorList") is not None:
                    Authors = PubmedArticle.find(".//AuthorList").findall(
                        "Author")
                    for auth in Authors:
                        LastName = None
                        ForeName = None
                        Initials = None
                        if auth.find("LastName") is not None:
                            LastName = auth.find("LastName").text
                        if auth.find("ForeName") is not None:
                            ForeName = auth.find("ForeName").text
                        if auth.find("Initials") is not None:
                            Initials = auth.find("Initials").text

                        # Add record to DB
                        author = Author(
                            LastName=LastName,
                            ForeName=ForeName,
                            Initials=Initials)
                        session.add(author)

    def populate_data(self, identifiers):

        # Workaround the max number of downloadable records
        RecordsCount = len(identifiers)
        id_list = []
        for i in range(
                math.ceil(RecordsCount / self.max_downloadable_records)):
            id_list.append(identifiers[i * self.max_downloadable_records:(
                i + 1) * self.max_downloadable_records])

        # Get the data
        for idl in id_list:
            handle = Entrez.efetch(
                db=self.query_db, id=idl, retmax=RecordsCount, retmode="xml")
            xml_text = handle.read()
            self.parse_and_commit_to_db(xml_text)
        handle.close()

    def main(self):
        RecordsCount = self.get_total_records_nb()
        identifiers = self.get_identifiers_list(RecordsCount)
        self.populate_data(identifiers)


if __name__ == '__main__':
    db = db(QUERY, QUERY_DB)
    db.main()
