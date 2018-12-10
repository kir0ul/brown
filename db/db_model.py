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
Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)


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
    parent = relationship("Article", back_populates="children")

    def __init__(self, PublicationMonth):
        self.PublicationMonth = MapStrMonthToInt(PublicationMonth)

    def __repr__(self):
        return "<Article(PMID={}, Title={}, PublicationDate={})>".format(
            self.PMID, self.Title, "/".join(self.PublicationYear,
                                            self.PublicationMonth))


class Author(Base):
    __tablename__ = 'authors'

    id = Column(Integer, primary_key=True)
    LastName = Column(String)
    ForeName = Column(String)
    Initials = Column(String)
    ArticlePMID = Column(Integer, ForeignKey('articles.PMID'), nullable=False)
    children = relationship("Article", back_populates="parent")

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
