# -*- coding: utf-8 -*-
"""
Time: 2020/Sep/4

@author: Benchi Zhao
Python version: 3.7
Functionality: This module is a code to test how fast the Mulitiprocess could fast the code.
"""

from multiprocessing import Process,Pool
import os, time, random

def fun1(name):
    print('Run task %s (%s)...' % (name, os.getpid()))
    start = time.time()
    time.sleep(2)
    end = time.time()
    print('Task %s runs %0.2f seconds.' % (name, (end - start)))


if __name__=='__main__':
    start = time.time()

    pool = Pool(5) # create a pool with 5 processes

    for i in range(10):
        pool.apply_async(func=fun1, args=(i,))

    pool.close()
    pool.join()
    print("with multiprocesses:", time.time() - start)

    start = time.time()
    for i in range(10):
        time.sleep(2)
    print("without multiprocesses:", time.time() - start)
    print('end test')