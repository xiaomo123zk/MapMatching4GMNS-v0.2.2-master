""" MapMatching4GMNS

Based on input network and given GPS trajectory data, the map-matching program
of Matching2Route aims to find most likely route in terms of node sequence in
the underlying network, with the following data flow chart.

The code is adopted and modified from
https://github.com/asu-trans-ai-lab/MapMatching4GMNS
"""

import time
import ctypes
import collections
import heapq
import os.path
from sys import platform
import shutil

print("To avoid complex data folder settings, please always first put the input data on the current directory.")
print('call MapMatching4GMNS  dynamic library')
if platform.startswith('win32'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.dll')
elif platform.startswith('linux'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.so')
elif platform.startswith('darwin'):
    _dll_file = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS.dylib')
else:
    raise Exception('Please build the shared library compatible to your OS\
                    using source files in engine_cpp!')

_cdll = ctypes.cdll.LoadLibrary(_dll_file)
_cdll.MapMatching4GMNS.argtypes = [ctypes.c_int]


def _optimal_MapMatching4GMNS_CAPI(mode):
    '''Two modes:
        modes 0: reading trace.csv, then automaticically generate input_agent.csv and output agent.csv
        model 1: reading input_agent.csv, then output agent.csv and link_performance.csv
    '''
    print('\MapMatching4GMNS run starts')
    _cdll.MapMatching4GMNS(mode)
    print('\MapMatching4GMNS run completes')


if __name__ == "__main__":
    # put the input data to current path

    src = './dataset/node.csv'
    dst = './'
    shutil.copy(src, dst)
    src = './dataset/link.csv'
    dst = './'
    shutil.copy(src, dst)
    src = './dataset/trace.csv'
    dst = './'
    shutil.copy(src, dst)
    src = './dataset/input_agent.csv'
    dst = './'
    shutil.copy(src, dst)

    # choose mode
    mode = 0

    start = time.time()
    _optimal_MapMatching4GMNS_CAPI(mode)
    end = time.time()

    print('MapMatching4GMNS time cost: %.6f seconds' % (end - start))
    print("The output data is generated!")
