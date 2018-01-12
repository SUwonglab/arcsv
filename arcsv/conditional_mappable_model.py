import os
import pysam
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')           # required if X11 display is not present
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import Counter, defaultdict
from sklearn import linear_model
from sklearn.metrics import log_loss
from patsy import dmatrix, dmatrices, build_design_matrices

from arcsv.constants import *

# 1 = mapped
# 0 = unmapped
#-1 = distant
PAIR_CLASSES = [(1,1), (1,0), (0,1),
                (1,-1), (-1,1)]
PAIR_CLASS_DICT = {PAIR_CLASSES[i]: i for i in range(len(PAIR_CLASSES))}

# INSERTIONS max distance undefined further down

# aln2 may be None
def process_aggregate_mapstats(pair, mapstats, min_mapq, max_distance):
    aln1, aln2 = pair
    aln1_pass = aln1.mapq >= min_mapq
    aln1_un = aln1.is_unmapped
    label = None
    add_mirror = False          # need to double count intrachromosomal 'distant' reads
                                # since they'll contribute to 2 regions
    if aln2 is not None:
        aln2_pass = aln2.mapq >= min_mapq
        aln2_un = aln2.is_unmapped
        is_distant = (aln1.rname != aln2.rname) or abs(aln1.pos - aln2.pos) > max_distance
        is_intra = (aln1.rname == aln2.rname)
        if aln1_pass and aln2_pass and not (aln1_un or aln2_un or is_distant):
            label = PAIR_CLASS_DICT[(1,1)]
        elif is_distant and not (aln1_un or aln2_un):
            if aln1_pass:
                label = PAIR_CLASS_DICT[(1,-1)]
            elif aln2_pass:
                label = PAIR_CLASS_DICT[(-1,1)]
            if is_intra:
                add_mirror = True
        elif aln1_pass and not aln1_un and aln2_un:
            label = PAIR_CLASS_DICT[(1,0)]
        elif aln2_pass and not aln2_un and aln1_un:
            label = PAIR_CLASS_DICT[(0,1)]
    elif aln1_pass and not aln1_un: # aln1 needs to be passing and mapped
        if aln1.mate_is_unmapped:
            label = PAIR_CLASS_DICT[(1,0)]
        elif (aln1.rname != aln1.mrnm) or \
             abs(aln1.pos - aln1.mpos) > max_distance: # 2nd read distant
            label = PAIR_CLASS_DICT[(1,-1)]
    if label is None:
        return

    mapstats[label] += 1
    if add_mirror:
        mirrored = PAIR_CLASS_DICT[(tuple(reversed(PAIR_CLASSES[label])))]
        mapstats[mirrored] += 1

# add mirrored observations to enforce symmetry
def add_dummy_obs_mapstats(mapstats):
    # print('[add dummy obs] before: {0}'.format(mapstats))
    num_label = len(mapstats)
    to_add = {}
    for cl in PAIR_CLASSES:       # (1,1), (1,0), etc.
        cl_mirrored = tuple(reversed(cl))
        label = PAIR_CLASS_DICT[cl]
        label_mirrored = PAIR_CLASS_DICT[cl_mirrored]
        to_add[label_mirrored] = mapstats[label]
    for label, amt in to_add.items():
        mapstats[label] += amt
    # print('[add dummy obs] after: {0}'.format(mapstats))

# input: list of mappability stats objects [defaultdict(int)]
# output: mappable_model, class_probs, rlen_stats
def model_from_mapstats(mapstats):
    # add dummy obs
    add_dummy_obs_mapstats(mapstats)

    n = sum(mapstats.values())
    pairs = list(mapstats.items())
    pairs.sort()
    class_prob = [p[1]/n for p in pairs]
    predicted_prob = lambda qmean1, rlen1, qmean2, rlen2, cp = tuple(class_prob): cp

    return predicted_prob, class_prob

def process_mappable_stats(aln1, aln2, mappable_stats, min_mapq, use_rlen = True):
    qmean1, qmean2 = np.mean(aln1.query_qualities), np.mean(aln2.query_qualities)
    rlen1, rlen2 = aln1.query_length, aln2.query_length
    aln1_pass = (aln1.mapq >= min_mapq)
    aln2_pass = (aln2.mapq >= min_mapq)
    aln1_un = aln1.is_unmapped
    aln2_un = aln2.is_unmapped

    # read pairs mapping farther than max_dist are considered the same as interchromosomal
    is_distant = (aln1.rname != aln2.rname) or abs(aln1.pos - aln2.pos) > max_distance

    label = None
    if aln1_pass and aln2_pass and not aln1_un and not aln2_un and not is_distant:
        label = PAIR_CLASS_DICT[(1,1)]
    elif is_distant and not aln1_un and not aln2_un:
        if aln1_pass:
            label = PAIR_CLASS_DICT[(1,-1)]
        elif aln2_pass:
            label = PAIR_CLASS_DICT[(-1,1)]
    elif aln1_pass and not aln1_un and aln2_un:
        label = PAIR_CLASS_DICT[(1,0)]
    elif aln2_pass and not aln2_un and aln1_un:
        label = PAIR_CLASS_DICT[(0,1)]
    if label == None:
        return

    mappable_stats['label'].append(label)
    mappable_stats['qmean1'].append(qmean1)
    mappable_stats['qmean2'].append(qmean2)
    # interchromosomal reads may need to be added twice
    if (label == PAIR_CLASS_DICT[(1,-1)] and aln2_pass) or (label == PAIR_CLASS_DICT[(-1,1)] and aln1_pass):
        mappable_stats['label'].append(PAIR_CLASS_DICT[tuple(reversed(PAIR_CLASSES[label]))])
        mappable_stats['qmean1'].append(qmean2)
        mappable_stats['qmean2'].append(qmean1)
    if use_rlen:
        mappable_stats['rlen1'].append(rlen1)
        mappable_stats['rlen2'].append(rlen2)
        if (label == PAIR_CLASS_DICT[(1,-1)] and aln2_pass) or (label == PAIR_CLASS_DICT[(-1,1)] and aln1_pass):
            mappable_stats['rlen1'].append(rlen2)
            mappable_stats['rlen2'].append(rlen1)

# create a model using aggregate stats only
# DEPRECATED?
def create_aggregate_mappable_model(bam_name, meta, outdir, min_mapq,
                                    sampling_rate = 1):
    bam = pysam.AlignmentFile(bam_name, 'rb')
    lib_patterns, lib_stats = parse_library_stats(meta)
    lib_dict = {}
    nlib = len(lib_stats)

    seen = {}
    rejected = set()
    mappable_stats = [{'label': [],
                       'rlen1': [],
                       'rlen2': [],
                       'qmean1': [],
                       'qmean2': []} for i in range(nlib)]
    nreads = 0
    reads = bam.fetch()
    for aln in reads:
        nreads += 1
        if nreads % 1000000 == 0:
            if aln.rname >= 0:
                print('chromosome {0} fetched {1} reads; pos {2}; len(seen) = {3}; len(rejected) = {4}'.format(bam.getrname(aln.rname), nreads, aln.pos, len(seen), len(rejected)))
        if aln.is_secondary or aln.is_supplementary:
            continue
        if aln.qname in rejected:
            rejected.remove(aln.qname)
            continue
        if aln.is_duplicate or aln.is_unmapped:
            continue
        if aln.qname in seen:
            mate = seen[aln.qname]
            del seen[aln.qname]
            lib_idx = get_lib_idx(aln.get_tag('RG'), lib_dict, lib_patterns)
            process_mappable_stats(aln, mate, mappable_stats[lib_idx], min_mapq, use_rlen = True)
        elif aln.mate_is_unmapped and \
                (sampling_rate == 1 or np.random.rand() <= sampling_rate):
            mate = pysam.AlignedSegment()
            mate.is_unmapped = True
            mate.pos = aln.pos
            mate.rname = aln.rname
            mate.seq = aln.seq
            mate.query_qualities = aln.query_qualities
            lib_idx = get_lib_idx(aln.get_tag('RG'), lib_dict, lib_patterns)
            process_mappable_stats(aln, mate, mappable_stats[lib_idx], min_mapq, use_rlen = True)
        else:
            if sampling_rate == 1 or np.random.rand() <= sampling_rate:
                seen[aln.qname] = aln
            else:
                rejected.add(aln.qname)

    for l in range(nlib):
        use_rlen = lib_stats[l]['readlen'] > 0
        add_dummy_obs(mappable_stats[l], use_rlen)

    # class_prob
    for l in range(nlib):
        n = len(mappable_stats[l]['label'])
        cts = Counter(mappable_stats[l]['label'])
        print('lib {0}: {1} pairs'.format(l, n))
        print('\t' + str(cts))
        pairs = list(cts.items())
        pairs.sort()
        freq = [p[1]/n for p in pairs]
        print('\tfrequencies {0}'.format(freq))
        outname = os.path.join(outdir, 'mapstats_{0}_{1}.pkl'.format(l, os.path.basename(bam_name)))
        outfile = open(outname, 'wb')
        pickle.dump(freq, outfile)
        outfile.close()

def create_mappable_model(bam_name, meta, outdir, min_mapq,
                          sampling_rate = 1):
    bam = pysam.AlignmentFile(bam_name, 'rb')
    lib_patterns, lib_stats = parse_library_stats(meta)
    lib_dict = {}
    nlib = len(lib_stats)

    seen = {}
    rejected = set()
    mappable_stats = [{'label': [],
                       'rlen1': [],
                       'rlen2': [],
                       'qmean1': [],
                       'qmean2': []} for i in range(nlib)]
    nreads = 0
    reads = bam.fetch(until_eof = True) # fetch unmapped pairs as well
    for aln in reads:
        nreads += 1
        if nreads % 1000000 == 0:
            if aln.rname >= 0:
                print('chromosome {0} fetched {1} reads; pos {2}; len(seen) = {3}; len(rejected) = {4}'.format(bam.getrname(aln.rname), nreads, aln.pos, len(seen), len(rejected)))
            # print('\nreads missing pairs are on these chromosomes:')
            # print(Counter([bam.getrname(a.rname) for a in seen.values()]))
            # print('\nreads missing pairs have mates on these chromosomes:')
            # print(Counter([bam.getrname(a.mrnm) for a in seen.values()]))
            # print('')

        # MEMORY MINOR if dups removed, unmapped reads with dup mates stick around forever...
        if aln.is_secondary or aln.is_supplementary: # split -> don't remove from rejected yet
            continue
        if aln.qname in rejected: # ignore b/c not sampled or mate is duplicate
            rejected.remove(aln.qname)
            continue
        if aln.is_duplicate:
            if aln.qname in seen: # happens if mate is unmapped --> not marked as dup
                del seen[aln.qname]
            else:
                rejected.add(aln.qname)
            continue
        if aln.qname in seen:
            mate = seen[aln.qname]
            del seen[aln.qname]
            lib_idx = get_lib_idx(aln.get_tag('RG'), lib_dict, lib_patterns)
            process_mappable_stats(aln, mate, mappable_stats[lib_idx], min_mapq, use_rlen = True)
        else:
            if sampling_rate == 1 or np.random.rand() <= sampling_rate:
                seen[aln.qname] = aln
            else:
                rejected.add(aln.qname)

    for l in range(nlib):
        use_rlen = lib_stats[l]['readlen'] > 0
        add_dummy_obs(mappable_stats[l], use_rlen)

    # save aggregate statistics
    # class_prob
    for l in range(nlib):
        n = len(mappable_stats[l]['label'])
        cts = Counter(mappable_stats[l]['label'])
        print('lib {0}: {1} pairs'.format(l, n))
        print('\t' + str(cts))
        pairs = list(cts.items())
        pairs.sort()
        freq = [p[1]/n for p in pairs]
        print('\tfrequencies {0}'.format(freq))
        outname = os.path.join(outdir, 'mapstats_{0}_{1}.pkl'.format(l, os.path.basename(bam_name)))
        outfile = open(outname, 'wb')
        pickle.dump(freq, outfile)
        outfile.close()
    # and rlen
    for l in range(nlib):
        rlen_short = [min(r1, r2) for (r1,r2) in zip(mappable_stats[l]['rlen1'],
                                                     mappable_stats[l]['rlen2'])]
        rlen_long = [max(r1, r2) for (r1,r2) in zip(mappable_stats[l]['rlen1'],
                                                    mappable_stats[l]['rlen2'])]
        rlen_stats = (int(np.mean(rlen_long)), int(np.mean(rlen_short)))
        outname = os.path.join(outdir, 'rlen_{0}_{1}.pkl'.format(l, os.path.basename(bam_name)))
        outfile = open(outname, 'wb')
        pickle.dump(rlen_stats, outfile)
        outfile.close()

    # fit model to each library and save results
    for l in range(nlib):
        use_rlen = lib_stats[l]['readlen'] > 0
        model, builder, qr = fit_mappable_model(mappable_stats[l], use_rlen)
        prefix = os.path.join(outdir, 'fit_{0}_{1}'.format(l, lib_stats[l]['name']))
        do_plot_fits(prefix, model, qr, builder, use_rlen, mappable_stats[l])
        if use_rlen:
            qgrid, rgrid, fitted_prob = predict_grid(model, builder,
                                                     mappable_stats[l],
                                                     use_rlen)
            outname = os.path.join(outdir, 'pmappable_{0}_{1}.pkl'.format(l, os.path.basename(bam_name)))
            save_model(outname, fitted_prob, qgrid, rgrid)
        else:
            qgrid, fitted_prob = predict_grid(model, builder,
                                              mappable_stats[l],
                                              use_rlen)
            outname = os.path.join(outdir, 'pmappable_{0}_{1}.pkl'.format(l, os.path.basename(bam_name)))
            save_model(outname, fitted_prob, qgrid)

def add_dummy_obs(mappable_stats, use_rlen):
    n = len(mappable_stats['label'])
    print(n)
    quantiles = [np.percentile(m, (10,50, 90)) for m in mappable_stats.values() if len(m) > 0]
    print('\n'.join([str(q) for q in quantiles]))
    ncl = []
    nq1 = []
    nq2 = []
    if use_rlen:
        nr1 = []
        nr2 = []
    for i in range(n):
        cl = mappable_stats['label'][i]
        q1 = mappable_stats['qmean2'][i]
        q2 = mappable_stats['qmean1'][i]
        if use_rlen:
            r1 = mappable_stats['rlen2'][i]
            r2 = mappable_stats['rlen1'][i]
        cl = PAIR_CLASS_DICT[tuple(reversed(PAIR_CLASSES[cl]))]
        ncl.append(cl), nq1.append(q1), nq2.append(q2)
        if use_rlen:
            nr1.append(r1), nr2.append(r2)
    mappable_stats['label'].extend(ncl)
    mappable_stats['qmean1'].extend(nq1)
    mappable_stats['qmean2'].extend(nq2)
    if use_rlen:
        mappable_stats['rlen1'].extend(nr1)
        mappable_stats['rlen2'].extend(nr2)
    print(len(mappable_stats['label']))
    quantiles = [np.percentile(m, (10,50, 90)) for m in mappable_stats.values() if len(m) > 0]
    print('\n'.join([str(q) for q in quantiles]))

def fit_mappable_model(mappable_stats, use_rlen, max_obs = 5000000):
    # fit model
    if use_rlen:
        design = 'label ~ (qmean1 + qmean2 + rlen1 + rlen2)*(qmean1 + qmean2 + rlen1 + rlen2)'
    else:
        design = 'label ~ qmean1 * qmean2'
    y, X = dmatrices(design, mappable_stats)
    builder = X.design_info.builder
    X = np.asarray(X)
    y = np.ravel(y)
    if max_obs < len(y):
        factor = max_obs/len(y)
        subsample = [i for i in range(len(y)) if np.random.rand() <= factor]
        X = X[subsample, :]
        y = y[subsample]
    print(X.shape)
    print(y.shape)
    model = linear_model.LogisticRegression(multi_class = 'multinomial', solver = 'newton-cg', fit_intercept = False, C = 1e6, max_iter = 50000)
    model.fit(X, y)
    print('loss on training data: {0}'.format(log_loss(y, predict_prob(model, X))))

    qmin, qmax = min(mappable_stats['qmean1']), max(mappable_stats['qmean2'])

    return model, builder, (qmin, qmax)

def predict_prob(model, X):
    tmp = X.dot(model.coef_.T)
    M = np.max(tmp, 1, keepdims = True)
    L = np.exp(tmp - M)
    return L / np.sum(L, 1, keepdims=True)

def predict_grid(model, builder, mappable_stats, use_rlen, qmean_res = 1, rlen_res = 10):
    qmin = min(mappable_stats['qmean1'] + mappable_stats['qmean2'])
    qmax = max(mappable_stats['qmean1'] + mappable_stats['qmean2'])
    qmin_rounded = int( np.floor(qmin / qmean_res) * qmean_res )
    qmax_rounded = int( np.ceil(qmax / qmean_res) * qmean_res )
    qrange = np.arange(qmin_rounded, qmax_rounded + qmean_res, qmean_res)
    if use_rlen:
        rmin = min(mappable_stats['rlen1'] + mappable_stats['rlen2'])
        rmax = max(mappable_stats['rlen1'] + mappable_stats['rlen2'])
        rmin_rounded = int( np.floor(rmin / rlen_res) * rlen_res )
        rmax_rounded = int( np.ceil(rmax / rlen_res) * rlen_res )
        rrange = np.arange(rmin_rounded, rmax_rounded + rlen_res, rlen_res)
    if use_rlen:
        q1g, q2g, r1g, r2g = np.meshgrid(qrange, qrange, rrange, rrange)
        q1ravel, q2ravel = np.ravel(q1g), np.ravel(q2g)
        r1ravel, r2ravel = np.ravel(r1g), np.ravel(r2g)
        m = len(q1ravel)
        x_new = build_design_matrices([builder],
                                      {'qmean1': np.ravel(q1g),
                                       'qmean2': np.ravel(q2g),
                                       'rlen1': np.ravel(r1g),
                                       'rlen2': np.ravel(r2g)})[0]
        x_new = np.asarray(x_new)
        pr = predict_prob(model, x_new)
        grid_dict = {(q1ravel[i], r1ravel[i], q2ravel[i], r2ravel[i]): tuple(pr[i]) for i in range(m)}
        return (qmin_rounded, qmax_rounded, qmean_res), (rmin_rounded, rmax_rounded, rlen_res), grid_dict
    else:
        q1g, q2g = np.meshgrid(qrange, qrange)
        q1ravel, q2ravel = np.ravel(q1g), np.ravel(q2g)
        m = len(q1ravel)
        x_new = build_design_matrices([builder],
                                      {'qmean1': np.ravel(q1g),
                                       'qmean2': np.ravel(q2g)})[0]
        x_new = np.asarray(x_new)
        pr = predict_prob(model, x_new)
        grid_dict = {(q1ravel[i], q2ravel[i]): tuple(pr[i]) for i in range(m)}
        return (qmin_rounded, qmax_rounded, qmean_res), grid_dict

def save_model(filename, fitted_prob, qgrid, rgrid = None):
    outfile = open(filename, 'wb')
    pickle.dump(fitted_prob, outfile)
    pickle.dump(qgrid, outfile)
    pickle.dump(rgrid, outfile)
    outfile.close()

def load_aggregate_model(model_dir, bam_name, lib_stats):
    nlib = len(lib_stats)
    predicted_prob = [None] * nlib
    class_prob = [None] * nlib
    rlen_stats = [None] * nlib
    for l in range(nlib):
        stats_name = '{0}mapstats_{1}_{2}.pkl'.format(model_dir, l, os.path.basename(bam_name))
        with open(stats_name, 'rb') as stats_file:
            class_prob[l] = pickle.load(stats_file)
        rlen_stats[l] = (0,0)

        # create constant functions for mappability model
        predicted_prob[l] = lambda qmean1, rlen1, qmean2, rlen2, cp = tuple(class_prob[l]): cp
    return predicted_prob, class_prob, rlen_stats

def load_model(model_dir, bam_name, lib_stats):
    nlib = len(lib_stats)
    predicted_prob = [None] * nlib
    class_prob = [None] * nlib
    rlen_stats = [None] * nlib
    for l in range(nlib):
        stats_name = '{0}mapstats_{1}_{2}.pkl'.format(model_dir, l, os.path.basename(bam_name))
        with open(stats_name, 'rb') as stats_file:
            class_prob[l] = pickle.load(stats_file)
        rlen_name = '{0}rlen_{1}_{2}.pkl'.format(model_dir, l, os.path.basename(bam_name))
        with open(rlen_name, 'rb') as rlen_file:
            rlen_stats[l] = pickle.load(rlen_file)
        model_name = '{0}pmappable_{1}_{2}.pkl'.format(model_dir, l, os.path.basename(bam_name))
        with open(model_name, 'rb') as model_file:
            use_rlen = lib_stats[l]['readlen'] > 0
            pred_dict = pickle.load(model_file)
            (qmin, qmax, q_res) = pickle.load(model_file)
            if use_rlen:
                (rmin, rmax, r_res) = pickle.load(model_file)
                predicted_prob[l] = lambda qmean1, rlen1, qmean2, rlen2, \
                    qm=qmin, qM=qmax, qr=q_res, rm=rmin, rM=rmax, rr=r_res, pd=pred_dict: \
                    pd[round_to_grid(qmean1, qm, qM, qr), \
                       round_to_grid(rlen1, rm, rM, rr), \
                       round_to_grid(qmean2, qm, qM, qr), \
                       round_to_grid(rlen2, rm, rM, rr)]
            else:
                predicted_prob[l] = lambda qmean1, rlen1, qmean2, rlen2, \
                    qm=qmin, qM=qmax, qr=q_res, pd=pred_dict: \
                    pd[round_to_grid(qmean1, qm, qM, qr), \
                       round_to_grid(qmean2, qm, qM, qr)]
    return predicted_prob, class_prob, rlen_stats

def round_to_grid(x, xmin, xmax, x_res):
    return np.round(min(max(x, xmin), xmax) / x_res) * x_res

def plot_fit(model, filename, qrange, builder, rlen1 = None, rlen2 = None, mappable_stats = None, max_scatter_points = 5000, rlen_tol = 10):
    qmin, qmax = qrange
    q1g, q2g = np.meshgrid(np.arange(qmin, qmax + .1, .1),
                           np.arange(qmin, qmax + .1, .1))
    use_rlen = (rlen1 is not None and rlen2 is not None)
    if use_rlen:
        n = np.prod(q1g.shape)
        x_new = build_design_matrices([builder],
                                      {'qmean1': np.ravel(q1g),
                                       'qmean2': np.ravel(q2g),
                                       'rlen1': np.repeat(rlen1, n),
                                       'rlen2': np.repeat(rlen2, n)})[0]
    else:
        x_new = build_design_matrices([builder],
                                      {'qmean1': np.ravel(q1g),
                                       'qmean2': np.ravel(q2g)})[0]
    proba = predict_prob(model, np.asarray(x_new))
    pr = [np.reshape(col, q1g.shape) for col in proba.T]

    if mappable_stats is not None:
        if use_rlen:
            rn = range(len(mappable_stats['label']))
            which_scatter = [i for i in rn if \
                                 abs(mappable_stats['rlen1'][i] - rlen1) < rlen_tol and \
                                 abs(mappable_stats['rlen2'][i] - rlen2) < rlen_tol]
        else:
            which_scatter = list(range(len(mappable_stats['label'])))
        if len(which_scatter) > max_scatter_points:
            ratio = max_scatter_points / len(which_scatter)
            which_scatter = [i for i in which_scatter if np.random.rand() <= ratio]

    plt.close('all')
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex = 'col', sharey = 'row')
    f.set_size_inches(10, 12)
    col = 'k'
    if use_rlen:
        titles = ('both mapped (r1 = {0}, r2 = {1})'.format(rlen1, rlen2), '1 mapped, 2 unmapped', '1 unmapped, 2 mapped', '1 mapped, 2 distant', '1 distant, 2 mapped')
    else:
        titles = ('both mapped', '1 mapped, 2 unmapped', '1 unmapped, 2 mapped', '1 mapped, 2 distant', '1 distant, 2 mapped')
    i = 0
    for ax in (ax1, ax2, ax3, ax4):
        if i >= len(pr):
            i += 1
            continue
        # ax.colorbar(im, orientation='horizontal', shrink = 0.8)
        con = ax.contour(q1g, q2g, pr[i], colors = col, linewidths=1.5)
        # ex = (qrange[0], qrange[1], qrange[0], qrange[1])
        ax.pcolormesh(q1g, q2g, pr[i], cmap = cm.coolwarm)
        # im = ax.imshow(pr[i], interpolation='bilinear', origin='lower', cmap=cm.gray, extent = ex)
        ax.clabel(con, inline=1, fontsize=10)

        if mappable_stats is not None:
            which_class = [j for j in which_scatter if mappable_stats['label'][j] == i]
            for j in which_class:
                ax.scatter(mappable_stats['qmean1'][j], mappable_stats['qmean2'][j], color = 'k', alpha = .3)

        ax.set_title(titles[i])
        ax.set_xlabel('read 1 quality')
        ax.set_ylabel('read 2 quality')
        i += 1
    plt.savefig(filename)

def do_plot_fits(prefix, model, qr, builder, use_rlen, mappable_stats):
    if use_rlen:
        r1q = np.percentile(mappable_stats['rlen1'], (90, 50, 10))
        r1q = [int(j) for j in r1q]
        r2q = np.percentile(mappable_stats['rlen2'], (90, 50, 10))
        r2q = [int(j) for j in r2q]
        rlen1 = (r1q[0], r1q[0], r1q[1], r1q[2])
        rlen2 = (r2q[2], r2q[1], r2q[0], r2q[0])
        for i in range(4):
            plot_fit(model, '{0}-{1}-{2}.png'.format(prefix, rlen1[i], rlen2[i]), qr, builder, rlen1[i], rlen2[i], mappable_stats = mappable_stats)
    else:
        plot_fit(model, '{0}.png'.format(prefix), qr, builder, mappable_stats = mappable_stats)
