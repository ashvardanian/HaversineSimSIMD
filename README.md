# Haversine Distance implementations for SimSIMD

- 0.py: SkLearn version. 2.9 secs for 1e5 iterations. __34.5 KOps/sec__.
- 1.py: Raw Python version. 8.5 secs for 1e7 iterations. __1.18 MOps/sec__.
- 2.py: NumBa version. 3.2 secs for 1e7 iterations. __3.12 MOps/sec__.
- 3.c: C version. 1.1 secs for 1e8 iterations. __90.9 MOps/sec__.
- 4.c: SIMD version. 2.5 secs for 1e9 iterations. __400 MOps/sec__.
