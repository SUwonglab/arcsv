def apply_filters(sv_list, rmsk_track=None, segdup_track=None):
    for sv in sv_list:
        if sv.type == 'INS' or (sv.type == 'BND' and sv.bnd_ins > 0):
            sv.filters.add('INSERTION')


def get_filter_string(sv, filter_criteria):
    intersection = set(sv.filters).intersection(set(filter_criteria))
    if len(intersection) > 0:
        return ','.join(sorted(intersection))
    else:
        return 'PASS'
