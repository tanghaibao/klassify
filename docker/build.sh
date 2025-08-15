#!/bin/bash

# Build the image (replace the tag if you pinned a different version)
docker build -t pbccs .

# Show the ccs help text to confirm it works
docker run --rm pbccs --help

# Typical usage example (mount local data):
docker run --rm -v $(pwd):/data pbccs \
        /data/subreads.bam /data/output.ccs.bam \
        --minPasses 3

