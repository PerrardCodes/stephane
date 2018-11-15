#!/usr/bin/env python
#Sort of like a worker in multiprocessing, but allows you not to have to preload all the data when doing something like a map.

from multiprocessing import Pipe, Process, cpu_count
import time
from math import sqrt
import sys

def worker_func(conn, func):
    while True:
        input = conn.recv()
        if input is None:
            break
        
        else:
            try:
                result = func(input)
            except Exception, e:
                result = e
            
            conn.send(result)
            
def fib(n):
    a, b = 0, 1
    for i in range(n):
        a, b = b, a+b
    return (n, b)
    
#Here's a function that should take some time...
def find_factors(n):
    nr = n
    factors = []
    while True:
        for i in range(2, int(sqrt(nr))+1):
            if not nr % i:
                nr //= i
                factors.append(i)
                break
        else:
            if nr != 1:
                factors.append(nr)
            break
    return n, factors

class WorkerPool(object):
    def __init__(self, func, num_proc=None, cycle_time=1E-2, print_result_number=True, done_mesg='-> Done in %.1f s.'):
        self.num_proc = num_proc if num_proc else cpu_count()
        self.cycle_time = cycle_time
        self.worker_conn = []
        self.worker = []
        self.ready = []
        self.print_result_number = print_result_number
        self.done_mesg = done_mesg
        self.start = None
        
        for i in range(self.num_proc):
            conn1, conn2 = Pipe(duplex=True)
            self.worker.append(Process(target=worker_func, args=(conn2, func)))
            self.worker_conn.append(conn1)
            self.worker[-1].start()
            self.ready.append(i)
            
        self.results = []
            
    def get_results(self):
        if not self.worker:
            raise ValueError('No workers, is this pool closed?')
        
        for i, conn in enumerate(self.worker_conn):
            if conn.poll():
                self.results.append(conn.recv())
                if isinstance(self.results[-1], Exception):
                    self.force_close()
                    print "Error in worker %d" % i
                    raise self.results[-1]
                
                self.ready.append(i)
                
                if self.print_result_number:
                    print '%3d' % len(self.results),
                    sys.stdout.flush()
    
    def send_job(self, args):
        if self.start is None: self.start = time.time()
        
        #Blocks until someone is ready
        while not len(self.ready):
            try:
                time.sleep(self.cycle_time)
                self.get_results()
            except KeyboardInterrupt:
                print "\n!!! KEYBOARD INTERRUPT !!!\nForcing workers closed!"
                self.force_close()
                sys.exit()
            
        i = self.ready.pop()
        #print "This is the args file", args
        self.worker_conn[i].send(args)
        
    def return_results(self):
        while len(self.ready) != self.num_proc:
            time.sleep(self.cycle_time)    
            results = self.get_results()
            
        if self.print_result_number:
            if '%' in self.done_mesg:
                print self.done_mesg % (time.time() - self.start)
            else:
                print self.done_mesg
            sys.stdout.flush()
               
        results = self.results     
        self.results = []
        self.start = None
        return results
            
    def close(self):
        for conn in self.worker_conn: conn.send(None)
        self.ready = []
        self.worker = []
        self.worker_conn = []
        
    def force_close(self):
        for worker in self.worker:
            worker.terminate()
        self.ready = []
        self.worker = []
        self.worker_conn = []
        
    def __del__(self):
        self.close
        
        
if __name__ == '__main__':
    pool = WorkerPool(find_factors, print_result_number=True)
    
    i0 = 10**14
    for i in range(30): pool.send_job(i + i0)
    
    
    results = pool.return_results()
    results.sort()
    pool.close()
    
    for i, results in results: print i, results
    