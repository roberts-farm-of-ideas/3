#
FROM python:3.7.8-slim-buster AS base

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt
#RUN pip install \
#    numpy       \
#    sympy       \
#    pulp        

RUN pip freeze > pipfreeze 

COPY . .

CMD [ "python", "./test-cycle-finding-algorithm.py" ]