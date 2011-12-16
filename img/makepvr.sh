#!/bin/bash

if [ -z "$1" ]; then
    echo "usage: $0 [2|4]"
    echo "Please specify whether you'd like to create a 2bpp or 4bpp compressed texture"
    exit
fi

/Developer/Platforms/iPhoneOS.platform/Developer/usr/bin/texturetool -m -e PVRTC --alpha-is-opacity --bits-per-pixel-$1 -p firefox-preview-$1bpp.png -o firefox-$1bpp.pvr -f Raw firefox-original.png