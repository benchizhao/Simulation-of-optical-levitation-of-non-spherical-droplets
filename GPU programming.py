# -*- coding: utf-8 -*-
"""
Time: 2020/Sep/4

@author: Benchi Zhao
Python version: 3.7
Functionality: This module is a code to test how fast the GPU programming could fast the code.
"""

from numba import jit, cuda
import numpy as np
# to measure exec time
from timeit import default_timer as timer

# normal function to run on cpu
def func(a):
    for i in range(10000000):
        a[i] += 1

    # function optimized to run on gpu

@jit
def func2(a):
    for i in range(10000000):
        a[i] += 1

if __name__ == "__main__":
    n = 10000000
    a = np.ones(n, dtype=np.float64)
    b = np.ones(n, dtype=np.float32)

    start = timer()
    func(a)
    print("without GPU:", timer() - start)

    start = timer()
    func2(a)
    print("with GPU:", timer() - start)