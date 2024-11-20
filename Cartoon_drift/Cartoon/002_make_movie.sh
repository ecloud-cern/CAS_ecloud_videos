#!/bin/bash

ffmpeg -i frames/frame_%d.png -c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=1.*PTS" -profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22  -codec:a aac pinch.mp4
