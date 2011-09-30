#!/bin/bash
# NAME: biology_tests.sh

# make sure job runs in current directory
#$ -cwd
# make stdout and stderr go to same file
#$ -j y
# Make use computers that are in the blue room.
#$ -l hostname=blue*
# Run multiple instances of the script (maximum please)
#$ -t 1-50

# load system profile (because job scripts do not do this by default)
. /etc/profile

# Feed script into matlab
echo "proccessnetworks" | matlab -nodisplay
