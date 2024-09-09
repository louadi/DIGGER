# pull official python base image
FROM python:3.7

# Unbuffered stdin / out for faster duping of logs
ENV PYTHONUNBUFFERED=1

# Install cron
RUN apt-get update && apt-get -y install cron

# Set up work directory
RUN mkdir /code
WORKDIR /code

# Copy necessary files
COPY requirements.txt /code/
COPY auto_clean.sh /usr/local/bin/

# Make auto_clean.sh executable
RUN chmod +x /usr/local/bin/auto_clean.sh

# Add cron job
RUN (crontab -l ; echo "0 0 * * * /usr/local/bin/auto_clean.sh >> /var/log/auto_clean.log 2>&1") | crontab -

# Install Python dependencies
RUN pip install -r requirements.txt
# COPY . /code/

# Copy the startup script
COPY start.sh /usr/local/bin/start.sh

# Make the startup script executable
RUN chmod +x /usr/local/bin/start.sh

# Use the custom startup script as the entry point
ENTRYPOINT ["/usr/local/bin/start.sh"]
