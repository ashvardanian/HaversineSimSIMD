echo "100 thousand iterations with SkLearn"
time python 0.py

echo "10 million iterations in Python"
time python 1.py
time python 2.py
time python 3.py

echo "1 billion iterations in C"
gcc -O3 3.c -o 3 -lm && time ./3
gcc -O3 4.c -o 4 -lm && time ./4

echo "Google Benchmark"
cmake -B ./build_release
cmake --build ./build_release --config Release
./build_release/main
