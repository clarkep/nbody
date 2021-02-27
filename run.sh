# Constants
N=1000
INIT_FILE=products/init.csv
OUTPUT_FILE=products/output.csv
OUTPUT_MOVIE=products/nbody.mp4
TIME_STEP=0.1
DURATION=10000.0
OUTPUT_EVERY=10     # time steps between writing to the output file
G=1.0               # gravitational constant
EPS=2.0             # Smoothing constant to avoid singularities in the forces
THETA=0.0           # Barnes-Hut constant to decide when to combine far-away bodies.
                    #   Lower is more accurate and slower.


python initial.py $N

gcc -Ofast -Wall -std=c99\
             -D N=$N \
             -D TIME_STEP=$TIME_STEP \
             -D DURATION=$DURATION \
             -D OUTPUT_EVERY=$OUTPUT_EVERY \
             -D G=$G \
             -D EPS=$EPS \
             -D THETA=$THETA \
             nbody.c -o nbody.o -lm

./nbody.o $INIT_FILE > $OUTPUT_FILE

rm products/frames/*
python animator_matplotlib.py $OUTPUT_FILE

ffmpeg -r 30 -f image2 -s 1920x1080 -i products/frames/%d.png \
       -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 25  \
       -pix_fmt yuv420p -y $OUTPUT_MOVIE

open $OUTPUT_MOVIE
