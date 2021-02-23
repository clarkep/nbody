python initial.py
gcc -O3 -Wall nbody.c -o nbody.o
./nbody.o products/init.csv > products/output.csv
python animator.py
