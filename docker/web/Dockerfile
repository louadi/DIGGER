# pull official python base image
FROM python:3.7
# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get -y install cron
# Copy code
RUN mkdir /code
WORKDIR /code
COPY requirements.txt /code/
COPY auto_clean.sh /code/
RUN chmod +x auto_clean.sh
# cron job to run auto_clean.sh
RUN (crontab -l ; echo "0 0 * * * /code/auto_clean.sh") | crontab -

RUN pip install -r requirements.txt
# COPY . /code/