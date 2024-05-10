# pull official python base image
FROM python:3.7
# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED 1
# Copy code
RUN mkdir /code
WORKDIR /code
COPY requirements.txt /code/
RUN pip install -r requirements.txt
# COPY . /code/