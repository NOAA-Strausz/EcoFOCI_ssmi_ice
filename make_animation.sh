#!/bin/bash

#making this to remember the ffmpeg settings for making the ice animations
#11/16/2026
#First check that 2 arguments were provided


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 \"<input_glob>\" <output_file>"
    echo "Example: $0 \"*.png\" output.mp4"
    exit 1
fi

INPUT_PATTERN="$1"
OUTPUT_FILE="$2"

#run ffmpeg including some compatibility settings
ffmpeg -framerate 8 -f image2 -pattern_type glob -i "$INPUT_PATTERN" \
    -c:v libx264 -r 30 -pix_fmt yuv420p -movflags +faststart \
    "$OUTPUT_FILE"
