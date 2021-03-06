import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')           # required if X11 display is not present
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn import linear_model
from sklearn.metrics import log_loss


# 1 = mapped
# 0 = unmapped
# -1 = distant
PAIR_CLASSES = [(1, 1), (1, 0), (0, 1),
                (1, -1), (-1, 1)]
PAIR_CLASS_DICT = {PAIR_CLASSES[i]: i for i in range(len(PAIR_CLASSES))}
# INSERTIONS max distance undefined further down


# aln2 may be None
def process_aggregate_mapstats(pair, mapstats, min_mapq, max_distance):
    aln1, aln2 = pair
    aln1_pass = aln1.mapq >= min_mapq
    aln1_un = aln1.is_unmapped
    label = None
    # need to double count intrachromosomal 'distant' reads since
    # they'll contribute to 2 regions
    add_mirror = False

    if aln2 is not None:
        aln2_pass = aln2.mapq >= min_mapq
        aln2_un = aln2.is_unmapped
        is_distant = (aln1.rname != aln2.rname) or abs(aln1.pos - aln2.pos) > max_distance
        is_intra = (aln1.rname == aln2.rname)
        if aln1_pass and aln2_pass and not (aln1_un or aln2_un or is_distant):
            label = PAIR_CLASS_DICT[(1, 1)]
        elif is_distant and not (aln1_un or aln2_un):
            if aln1_pass:
                label = PAIR_CLASS_DICT[(1, -1)]
            elif aln2_pass:
                label = PAIR_CLASS_DICT[(-1, 1)]
            if is_intra:
                add_mirror = True
        elif aln1_pass and not aln1_un and aln2_un:
            label = PAIR_CLASS_DICT[(1, 0)]
        elif aln2_pass and not aln2_un and aln1_un:
            label = PAIR_CLASS_DICT[(0, 1)]
    elif aln1_pass and not aln1_un:  # aln1 needs to be passing and mapped
        if aln1.mate_is_unmapped:
            label = PAIR_CLASS_DICT[(1, 0)]
        # 2nd read distant
        elif (aln1.rname != aln1.mrnm) or abs(aln1.pos - aln1.mpos) > max_distance:
            label = PAIR_CLASS_DICT[(1, -1)]
    if label is None:
        return

    mapstats[label] += 1
    if add_mirror:
        mirrored = PAIR_CLASS_DICT[(tuple(reversed(PAIR_CLASSES[label])))]
        mapstats[mirrored] += 1


# add mirrored observations to enforce symmetry
def add_dummy_obs_mapstats(mapstats):
    # print('[add dummy obs] before: {0}'.format(mapstats))
    to_add = {}
    for cl in PAIR_CLASSES:       # (1, 1), (1, 0), etc.
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


def add_dummy_obs(mappable_stats, use_rlen):
    n = len(mappable_stats['label'])
    print(n)
    quantiles = [np.percentile(m, (10, 50, 90)) for m in mappable_stats.values() if len(m) > 0]
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
    quantiles = [np.percentile(m, (10, 50, 90)) for m in mappable_stats.values() if len(m) > 0]
    print('\n'.join([str(q) for q in quantiles]))


def load_aggregate_model(model_dir, bam_name, lib_stats):
    nlib = len(lib_stats)
    predicted_prob = [None] * nlib
    class_prob = [None] * nlib
    rlen_stats = [None] * nlib
    for l in range(nlib):
        stats_name = '{0}mapstats_{1}_{2}.pkl'.format(model_dir, l, os.path.basename(bam_name))
        with open(stats_name, 'rb') as stats_file:
            class_prob[l] = pickle.load(stats_file)
        rlen_stats[l] = (0, 0)

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
                    pd[round_to_grid(qmean1, qm, qM, qr),
                       round_to_grid(rlen1, rm, rM, rr),
                       round_to_grid(qmean2, qm, qM, qr),
                       round_to_grid(rlen2, rm, rM, rr)]
            else:
                predicted_prob[l] = lambda qmean1, rlen1, qmean2, rlen2, \
                    qm=qmin, qM=qmax, qr=q_res, pd=pred_dict: \
                    pd[round_to_grid(qmean1, qm, qM, qr), round_to_grid(qmean2, qm, qM, qr)]
    return predicted_prob, class_prob, rlen_stats


def round_to_grid(x, xmin, xmax, x_res):
    return np.round(min(max(x, xmin), xmax) / x_res) * x_res


def plot_fit(model, filename, qrange, builder, rlen1=None, rlen2=None, mappable_stats=None, max_scatter_points=5000, rlen_tol=10):
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
            which_scatter = [i for i in rn if
                             abs(mappable_stats['rlen1'][i] - rlen1) < rlen_tol and
                             abs(mappable_stats['rlen2'][i] - rlen2) < rlen_tol]
        else:
            which_scatter = list(range(len(mappable_stats['label'])))
        if len(which_scatter) > max_scatter_points:
            ratio = max_scatter_points / len(which_scatter)
            which_scatter = [i for i in which_scatter if np.random.rand() <= ratio]

    plt.close('all')
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
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
        con = ax.contour(q1g, q2g, pr[i], colors=col, linewidths=1.5)
        # ex = (qrange[0], qrange[1], qrange[0], qrange[1])
        ax.pcolormesh(q1g, q2g, pr[i], cmap=cm.coolwarm)
        # im = ax.imshow(pr[i], interpolation='bilinear', origin='lower', cmap=cm.gray, extent = ex)
        ax.clabel(con, inline=1, fontsize=10)

        if mappable_stats is not None:
            which_class = [j for j in which_scatter if mappable_stats['label'][j] == i]
            for j in which_class:
                ax.scatter(mappable_stats['qmean1'][j], mappable_stats['qmean2'][j], color='k', alpha=.3)

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
            plot_fit(model, '{0}-{1}-{2}.png'.format(prefix, rlen1[i], rlen2[i]), qr, builder, rlen1[i], rlen2[i], mappable_stats=mappable_stats)
    else:
        plot_fit(model, '{0}.png'.format(prefix), qr, builder, mappable_stats=mappable_stats)
