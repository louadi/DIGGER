#!/bin/bash

# Start the cron service
service cron start

# Start Gunicorn with your application settings
exec gunicorn DomainExplorer.wsgi:application --bind 0.0.0.0:8000 --timeout 360 --workers=5 --threads=2 --max-requests=1000 --max-requests-jitter=50
