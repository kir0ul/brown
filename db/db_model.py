from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager

USR = "postgres"
PASS = "postgres"
DOMAIN = "localhost"
PORT = "5432"
POSTGRES_URL = "".join(
    ["postgresql://", USR, ":", PASS, "@", DOMAIN, ":", PORT])

engine = create_engine(POSTGRES_URL, echo=True)
Base = declarative_base()


@contextmanager
def db_session():
    """Provide a transactional scope around a series of operations."""

    session = Session()
    try:
        yield session
        session.commit()
    except BaseException:
        session.rollback()
        raise
    finally:
        session.close()


class Article(Base):
    __tablename__ = 'articles'

    id = Column(Integer, primary_key=True)
    PMID = Column(String, nullable=False, unique=True)
    Title = Column(String)
    PublicationYear = Column(Integer)
    PublicationMonth = Column(Integer)
    MedlineDate = Column(String)

    # parent = relationship("Article", back_populates="children")

    def __init__(self, PMID, Title, PublicationYear, PublicationMonth,
                 MedlineDate):
        self.PMID = PMID,
        self.Title = Title,
        self.PublicationYear = PublicationYear,
        self.PublicationMonth = MapStrMonthToInt(PublicationMonth)
        self.MedlineDate = MedlineDate

    def __repr__(self):
        return "<Article(PMID={}, Title={}, PublicationDate={})>".format(
            self.PMID, self.Title,
            "/".join(self.PublicationYear, self.PublicationMonth)
            or self.MedlineDate)


class Author(Base):
    __tablename__ = 'authors'

    id = Column(Integer, primary_key=True)
    LastName = Column(String)
    ForeName = Column(String)
    Initials = Column(String)
    ArticlePMID = Column(
        String)  #, ForeignKey('articles.PMID'), nullable=False)

    # children = relationship("Article", back_populates="parent")

    def __repr__(self):
        return "<Author(LastName={}, ForeName={}, PMID={})>".format(
            self.LastName, self.ForeName, self.ArticlePMID)


def MapStrMonthToInt(str):
    Mapping = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
        "Nov", "Dec"
    ]
    res = None
    for idx, val in enumerate(Mapping):
        if str == val:
            res = idx + 1
            break

    return res


Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)

# # check register table exist to Modal
# for _t in Modal.metadata.tables:
#     print(_t)

# # check all table in database
# meta = MetaData(engine, reflect=True)
# for _t in meta.tables:
#     print(_t)

# # check table names exists via inspect
# ins = inspect(engine)
# for _t in ins.get_table_names():
#     print(_t)

for _t in Base.metadata.tables:
    print("Table: ", _t)
