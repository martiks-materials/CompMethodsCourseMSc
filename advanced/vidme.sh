#!/bin/bash

ffmpeg -r 60 -f image2 -s 1920x1080 -i images/im_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4
