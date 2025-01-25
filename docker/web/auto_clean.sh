#!/bin/bash
# This script is used to clean the files in the subdirectories of /code/nease_events, /code/run/media directories + any subdirectory
# The script is scheduled to run every day at 00:00 using a cron job

# define the different folders that should be cleaned
FOLDER_ONE="zero_days"
FOLDER_TWO="seven_days"
FOLDER_THREE="thirtyone_days"
FOLDER_FOUR="onehundredeightysix_days"
FOLDERS=($FOLDER_ONE $FOLDER_TWO $FOLDER_THREE $FOLDER_FOUR)

# define the number of days after which the files should be deleted
ONE_DELETE_AFTER=0
TWO_DELETE_AFTER=7
THREE_DELETE_AFTER=31
FOUR_DELETE_AFTER=186
DELETION_TIMES=($ONE_DELETE_AFTER $TWO_DELETE_AFTER $THREE_DELETE_AFTER $FOUR_DELETE_AFTER)

# Loop over the folders and deletion times
for i in "${!FOLDERS[@]}"; do
    # Find and delete files in the corresponding /code/nease_events directory
    echo "Found $(find /code/nease_events/${FOLDERS[$i]} -type f -mtime +${DELETION_TIMES[$i]} | wc -l) files to delete, as they are older than ${DELETION_TIMES[$i]} days."
    find /code/nease_events/${FOLDERS[$i]} -type f -mtime +${DELETION_TIMES[$i]} -delete
done

# Find and delete files in /code/run/media directory
echo "Found $(find /code/run/media -type f -mtime +FOUR_DELETE_AFTER | wc -l) files to delete"
find /code/run/media -type f -mtime +FOUR_DELETE_AFTER -delete
