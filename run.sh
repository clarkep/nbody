# Constants
N=10000
INIT_FILE=products/init.csv
OUTPUT_FILE=products/output10k.csv
OUTPUT_MOVIE=products/nbody.mp4
TIME_STEP=0.1     # seconds
DURATION=10000.0  # seconds
OUTPUT_EVERY=10     # time steps between writing to the output file
G=0.2               # gravitational constant
EPS=2.0             # Smoothing constant to avoid singularities in the forces
THETA=0.5           # Barnes-Hut constant to decide when to combine far-away bodies.
                    #   Lower is more accurate and slower.
BRUTE_FORCE=0     # Whether to use the O(n^2) brute force algorithm


python initial.py $N

gcc-7 -Ofast -Wall -std=c99 -fopenmp \
             -D N=$N \
             -D TIME_STEP=$TIME_STEP \
             -D DURATION=$DURATION \
             -D OUTPUT_EVERY=$OUTPUT_EVERY \
             -D G=$G \
             -D EPS=$EPS \
             -D THETA=$THETA \
             -D BRUTE_FORCE=$BRUTE_FORCE \
             nbody.c -o nbody.o -lm

nice ./nbody.o $INIT_FILE > $OUTPUT_FILE

# rm products/frames/*

python animator_vpython.py $OUTPUT_FILE

# ffmpeg -r 30 -f image2 -s 1920x1080 -i products/frames/%d.png \
#        -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 25  \
#        -pix_fmt yuv420p -y $OUTPUT_MOVIE


# open $OUTPUT_MOVIE
