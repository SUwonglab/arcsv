import os
import numpy as np
from math import ceil, floor

# uses dictionary instead of array
class SparseSignalTrack(object):
    valid_types = ['int', 'array']

    def __init__(self, chrom_name, signal_type='int'):
        self.chrom_name = chrom_name
        if signal_type not in self.valid_types:
            raise ValueError('signal_type is invalid, choose from ' + str(self.valid_types))
        self.signal_type = signal_type
        self.signal = {}

    def __len__(self):
        return len(list(self.signal.keys()))

    def __add__(self, other):
        new_track = SparseSignalTrack(self.chrom_name, self.signal_type)
        all_keys = set(self.signal.keys()).union(set(other.signal.keys()))
        new_signal = {}
        d = 0 if new_track.signal_type == 'int' else []
        new_signal = {k: self.signal.get(k, d) + other.signal.get(k, d) for k in all_keys}
        new_track.signal = new_signal
        return new_track

    def __radd__(self, other):  # need this to do use sum()
        if other == 0:
            return self
        else:
            raise TypeError('addition of SparseSignalTrack and {0} not supported'.format(type(other)))

    def add(self, location, value = 1):
        if self.signal_type == 'int':
            self.signal[location] = self.signal.get(location, 0) + value
        elif self.signal_type == 'array':
            self.signal[location] = self.signal.get(location, []) + [value]

    def add_all(self, locations):
        for loc in locations:
            self.add(location = loc, value = 1)

    def windowed_signal(self, location, window):
        # loc_high - loc_low + 1 = window
        loc_low = location - int(ceil((window - 1) / 2.0))
        loc_high = location + int(floor((window - 1) / 2.0))
        if self.signal_type == 'int':
            sig = 0
        elif self.signal_type == 'array':
            sig = []
        for i in range(loc_low, loc_high + 1):
            if self.signal_type == 'int':
                sig += self.signal.get(i, 0)
            elif self.signal_type == 'array':
                sig.extend(self.signal.get(i, []))
        return sig

    def write_bed(self, fileprefix, type = 'count', every = 1, window = 1,
                  mu = None, sigma = None):
        if every > window:
            print('Warning: sliding window size less than sliding amount')
        file = open(fileprefix + '.bed', 'w')
        if len(self) == 0:
            file.write('\n')
            file.close()
            return
        minloc = min(list(self.signal.keys()))
        maxloc = max(list(self.signal.keys()))
        for loc in range(minloc, maxloc + 1, every):
            windowed = self.windowed_signal(loc, window)
            if windowed == [] or windowed == 0:
                continue
            if every > 1:
                loc_low = loc - int(floor(every/2.0))
                loc_high = loc + int(ceil(every/2.0))
            else:
                loc_low = loc
                loc_high = loc + 1
            if self.signal_type == 'int' and type != 'count':
                raise ValueError('invalid type argument.')
            elif self.signal_type == 'int':
                value = windowed/float(window)
            elif type == 'count': # signal_type == 'array'
                value = len(windowed)/float(window)
            elif type == 'mean':
                value = np.mean(windowed)
            elif type == 'zscore':
                value = zscore(windowed, mu, sigma)
            file.write('{chr}\t{start}\t{end}\t{val}\n'.format(chr = self.chrom_name,
                                                              start = loc_low,
                                                              end = loc_high,
                                                              val = value))
        file.close()

    def write_bigwig(self, fileprefix, type = 'count', every = 1, window = 1,
                     mu = None, sigma = None):
        self.write_bed(fileprefix, type, every, window, mu, sigma)
        os.system('bedGraphToBigWig {file}.bed /scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {file}.bigwig'.format(file=fileprefix))

# not used
class SignalTrack(object):
    def __init__(self, chrom_name, start, end):
        self.chrom_name = chrom_name
        self.start = start
        self.end = end
        self.len = end - start
        self.signal = [None] * self.len

    def __len__(self):
        return self.len

    # def add(self, location, value):
    #     idx = self.get_index(location)
    #     if self.signal[idx] is None:
    #         self.signal[idx] = value
    #     else:
    #         self.signal[idx] += value

    def append(self, location, value = 1):
        idx = self.get_index(location)
        if self.signal[idx] is None:
            self.signal[idx] = [value]
        else:
            self.signal[idx].append(value)

    # convert from real location to array index
    def get_index(self, location):
        idx = location - self.start
        if idx < 0 or idx >= self.len:
            raise IndexError
        return idx

    # convert from array index to real location
    def get_location(self, idx):
        if idx < 0 or idx >= self.len:
            raise IndexError
        return self.start + idx

    def windowed_signal(self, idx, window):
        if idx < 0 or idx >= self.len:
            raise IndexError
        # idx_high - idx_low + 1 = window
        idx_low = idx - int(ceil((window - 1) / 2.0))
        idx_low = max(0, idx_low)
        idx_high = idx + int(floor((window - 1) / 2.0))
        idx_high = min(self.len - 1, idx_high)
        sig = []
        for i in range(idx_low, idx_high + 1):
            if self.signal[i] is not None:
                sig.extend(self.signal[i])
        return sig

    # type 0:raw values
    # type 1:mean
    # type 2:z scores
    # currently every = 1 is required
    def write_file(self, filename, type = 'count',
                   every = 1, window = 1, mu = None, sigma = None):
        file = open(filename, 'w')
        for idx in range(0, self.len, every):
            windowed = self.windowed_signal(idx, window)
            if windowed == []:
                continue
            if type == 'count':
                value = len(windowed)/float(window)
            if type == 'mean':
                value = np.mean(windowed)
            if type == 'zscore':
                value = zscore(windowed, mu, sigma)
            location = self.get_location(idx)
            file.write('{chr}\t{loc}\t{locp}\t{val}\n'.format(chr = self.chrom_name,
                                                              loc = location,
                                                              locp = location + 1,
                                                              val = value))
        file.close()

def zscore(L, mu, sigma):
    return (np.mean(L) - float(mu))/(float(sigma)/sqrt(len(L)))

def write_trackdb(file, libname, trackname, extension, tracktype,
                  itemRgb = False, heightPixels = None, color = None,
                  visibility = None, viewMin = None, viewMax = None):
    out = "track {0}-{1}\n" + "bigDataUrl {0}-{1}.{2}\n" + "shortLabel {0}-{1}\n" + \
        "longLabel {0}-{1}\n" + "type {3}\n"
    out = out.format(libname, trackname, extension, tracktype)
    if visibility:
        out += "visibility {0}\n".format(visibility)
    if itemRgb:
        out += "itemRgb on\n"
    if heightPixels:
        out += "maxHeightPixels 100:{0}:8\n".format(heightPixels)
    elif extension == "bigwig":
        out += "maxHeightPixels 100:32:8\n"
    if viewMin is not None and viewMax is not None:
        out += "viewLimits {0}:{1}\n".format(viewMin, viewMax)
        out += "viewLimitsMax {0}:{1}\n".format(min(viewMin - 20, 0),
                                                10 * viewMax)
    if color == 'orange':
        out += "color 240,162,29\n"
    elif color == 'magenta':
        out += "color 240,29,222\n"
    out += "\n"
    file.write(out)

def write_array_bed(arr, chrom_name, fileprefix, start = None, end = None):
    file = open(fileprefix + '.bed', 'w')
    ucsc_chrom = get_ucsc_name(chrom_name)
    if start is None:
        start = 0
    if end is None:
        # TODO need to specify reference
        end = get_chrom_size(chrom_name) - 1
    for i in range(start, end):
        file.write('{chr}\t{start}\t{end}\t{val}\n'.format(chr = ucsc_chrom, start = i, end = i+1,
                                                           val = arr[i]))
    file.close()

def write_array_bigwig(arr, chrom_name, fileprefix, start = None, end = None, delete_bed = False):
    write_array_bed(arr, chrom_name, fileprefix, start, end)
    os.system('bedGraphToBigWig {file}.bed /scratch/PI/whwong/svproject/reference/hg19.chrom.sizes {file}.bigwig'.format(file=fileprefix))
    if delete_bed:
        os.system('rm %s' % (fileprefix + '.bed'))

def testsignal():
    s = SignalTrack('chr1', 0, 100)
    s.append(0,10)
    s.append(0,15)
    s.append(0,14)

    s.write_file('zscore.txt', type='zscore', mu=11, sigma = 4)
    s.write_file('mean.txt', type='mean')
    s.write_file('val.txt')

def testsparsesignal():
    s = SparseSignalTrack(chrom_name='chr1', signal_type='int')
    s.add(0)
    s.add(1)
    s.add(100)
    s.write_file('sparse11.txt', every = 1, window = 1)
    s.write_file('sparse13.txt', every = 1, window = 3)
    s.write_file('sparse1-10.txt', every=1, window=10)

