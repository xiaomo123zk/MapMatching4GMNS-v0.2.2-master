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

if platform.startswith('win32'):
    _dll_file_noOMP = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS_noOMP.dll')
elif platform.startswith('linux'):
    _dll_file_noOMP = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS_noOMP.so')
elif platform.startswith('darwin'):
    _dll_file_noOMP = os.path.join(os.path.dirname(
        __file__), 'bin/MapMatching4GMNS_noOMP.dylib')
else:
    raise Exception('Please build the shared library compatible to your OS\
                    using source files in engine_cpp!')

_cdll = ctypes.cdll.LoadLibrary(_dll_file)
_cdll.MapMatching4GMNS.argtypes = [ctypes.c_int]

_cdll_noOMP = ctypes.cdll.LoadLibrary(_dll_file_noOMP)
_cdll_noOMP.MapMatching4GMNS.argtypes = [ctypes.c_int]


def mmg_CAPI(parallel_mode, input_data_mode):
    '''Two choice:
        parallel_mode=0: Using non-parallel version, which only uses the single thread to compute.
        parallel_mode=0: Using parallel version, which applies openMP for multithreading acceleration

        input_data_mode=0: reading trace.csv, then automaticically generate input_agent.csv and output agent.csv
        input_data_mode=1: reading input_agent.csv, then output agent.csv and link_performance.csv
    '''
    if(parallel_mode == 0):
        print('\MapMatching4GMNS run starts')
        _cdll_noOMP.MapMatching4GMNS(input_data_mode)
        print('\MapMatching4GMNS run completes')
    elif (parallel_mode == 1):
        print('\MapMatching4GMNS run starts')
        _cdll.MapMatching4GMNS(input_data_mode)
        print('\MapMatching4GMNS run completes')


def Data_cleaning():
    '''
    input_agent.csv
    agent.csv
    link_performance.csv
    '''
    # input_agent.csv
    '''
    There will be an extra comma(',') behind the geometry value in the output input_agent.csv 
    file due to different data specifications  by different programing languages. Therefore, we have to clean the data.
    '''
    if(os.path.exists("input_agent.csv")):
        # utf-8 按照 电脑环境；2）按照二进制方式读取，提前读取一遍 读取灵活
        f = open('input_agent.csv', encoding="utf-8")
        content = f.read()
        f.close()
        t = content.replace(",)", ")")  # delete the comma(',')
        with open("input_agent.CSV", "w", encoding='utf-8') as f1:
            f1.write(t)

    '''
    In the colab environment, if the output file already exists in the current directory, 
    an error will be reported when the newly generated output file is needed.
    '''
    # agent.csv
    if(os.path.exists("agent.csv")):
        os.remove("agent.csv")
    # link_performance.csv
    if(os.path.exists("link_performance.csv")):
        os.remove("link_performance.csv")


if __name__ == "__main__":
    # first, please to put the input data to current path

    # second, data cleaning
    Data_cleaning()

    # third choose mode
    parallel_mode = 1
    input_data_mode = 0

    # fourth,Call the mapmatching4gmns library to calculate and output the result in the current directory.
    start = time.time()
    mmg_CAPI(parallel_mode, input_data_mode)
    end = time.time()
    print('MapMatching4GMNS time cost: %.6f seconds' % (end - start))
    print("The output data is generated!")

    # fifth, Visualization which is donging by Zanyang Cui
