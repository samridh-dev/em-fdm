#!/bin/sh
FRAMERATE=60
ffmpeg -y \
  -framerate ${FRAMERATE} -i out/t%d.png \
  -threads $(nproc)\
  -c:v libx264 -pix_fmt yuv420p vid.mp4
