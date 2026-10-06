import os
import pickle
import numpy as np


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
        stats_name = f'{model_dir}mapstats_{l}_{os.path.basename(bam_name)}.pkl'
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
        stats_name = f'{model_dir}mapstats_{l}_{os.path.basename(bam_name)}.pkl'
        with open(stats_name, 'rb') as stats_file:
            class_prob[l] = pickle.load(stats_file)
        rlen_name = f'{model_dir}rlen_{l}_{os.path.basename(bam_name)}.pkl'
        with open(rlen_name, 'rb') as rlen_file:
            rlen_stats[l] = pickle.load(rlen_file)
        model_name = f'{model_dir}pmappable_{l}_{os.path.basename(bam_name)}.pkl'
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
