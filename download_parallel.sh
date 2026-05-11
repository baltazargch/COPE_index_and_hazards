#!/bin/bash

# Path to your list of download links
LINKS_FILE="outputs/monthly_climates/links_download.txt"

# Output folder
OUTDIR="outputs/monthly_climates/NAclim"
mkdir -p "$OUTDIR"

# Function to download a single file into the NAclim folder
download_file() {
  url="$1"
  filename=$(basename "$url")
  echo "Downloading $filename..."
  wget -q --show-progress -O "$OUTDIR/$filename" "$url"
}

# Export the function so GNU parallel can use it
export -f download_file

# Run parallel downloads (adjust -j 4 to your CPU/network capacity)
cat "$LINKS_FILE" | parallel -j 4 download_file {}

echo "✅ All files downloaded into: $OUTDIR/"

