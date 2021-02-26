python initial.py
gcc -g -Wall nbody.c -o nbody.o
./nbody.o products/init.csv > products/output.csv
python animator.py
