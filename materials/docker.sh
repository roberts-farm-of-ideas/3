#!/bin/bash
#docker build -t asdf . && docker run --rm -it asdf
docker build -t asdf . && docker run --rm -it asdf python ./experimentation-wedding-problem.py
#docker build -t asdf . && docker run --rm -it asdf python ./wedding_initial.py