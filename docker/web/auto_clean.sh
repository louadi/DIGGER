#!/bin/bash
# This script is used to clean the files in the /code/nease_events, /code/run/media directories + any subdirectory
# after 7 days since the last access time.
# The script is scheduled to run every day at 00:00 using a cron job

# check if theres an environment variable for the time to delete files
DELETE_AFTER=${DELETE_AFTER:-7}


# Find and delete files in /code/nease_events directory
echo "Found $(find /code/nease_events -type f -atime +$DELETE_AFTER | wc -l) files to delete"
find /code/nease_events -type f -atime +$DELETE_AFTER -delete

# Find and delete files in /code/run/media directory
echo "Found $(find /code/run/media -type f -atime +$DELETE_AFTER | wc -l) files to delete"
find /code/run/media -type f -atime +$DELETE_AFTER -delete
