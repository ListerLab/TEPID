#! /usr/local/python2.7
"""
mean_size.py
Created by Tim Stuart
"""

import numpy as np


def get_data(inp):
    lengths = []
    for line in inp:
        if line.startswith('@'):
            pass
        else:
            line = line.rsplit()
            length = int(line[8])
            if length > 0:
                lengths.append(length)
            else:
                pass
    return lengths


def reject_outliers(data, m=2.):
    """
    rejects outliers more than 2
    standard deviations from the median
    """
    median = np.median(data)
    std = np.std(data)
    for item in data:
        if abs(item - median) > m * std:
            data.remove(item)
        else:
            pass


def calc_size(data):
    mn = int(np.mean(data))
    std = int(np.std(data))
    return mn, std


if __name__ == "__main__":
    import sys
    lengths = get_data(sys.stdin)
    reject_outliers(lengths)
    mn, std = calc_size(lengths)
    print mn, std
