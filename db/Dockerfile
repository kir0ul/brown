FROM python:3.5-slim-stretch

RUN mkdir -p /data/
WORKDIR /data/

RUN apt update && apt install -y wait-for-it

COPY . .
RUN pip install -U pip pipenv
RUN pipenv install --system --deploy --ignore-pipfile

CMD ["python", "db.py"]
