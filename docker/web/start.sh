#!/bin/bash

# Start the cron service
service cron start

# Start Gunicorn with your application settings
# --preload imports the app (incl. the in-memory pickle graphs in load_data.py) once in
# the master process; workers are forked after, so the graphs are copy-on-write shared
# instead of loaded 5x. Cuts per-worker RAM and cold-start time.
exec gunicorn DomainExplorer.wsgi:application --bind 0.0.0.0:8000 --preload --timeout 360 --workers=5 --threads=2 --max-requests=1000 --max-requests-jitter=50
