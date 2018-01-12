import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import os
from math import floor
import matplotlib
matplotlib.use('Agg')           # required if X11 display is not present
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import jet, rainbow, gist_rainbow
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


from arcsv.helper import GenomeInterval

PHI = 1.618

def arrow_polygon(is_inverted, width, left_pos, bottom_pos):
    if width < 1:
        print('[arrow_polygon]: width < 1')
        width = 1
    template = [[0,0], [width - 1/2,0],
                [width, PHI/2], [width-1/2, PHI],
                [0,PHI]]
    for pt in template:
        if is_inverted:
            pt[0] = width - pt[0]
        pt[0] += left_pos
        pt[1] += bottom_pos

    if is_inverted:
        text_pos = [left_pos + width/2, bottom_pos + PHI/2 - 1/8]
    else:
        text_pos = [left_pos + width/2 - 1/4, bottom_pos + PHI/2 - 1/8]

    return template, text_pos

def get_arrow_patches(blocks, path, bottom_pos,
                      start_block, end_block, left_pos = 0):
    block_seq = [int(floor(path[i]/2)) for i in range(0, len(path), 2)]
    inverted_seq = [path[i] % 2 == 1 for i in range(0, len(path), 2)]

    patches = []
    patches_ins = []
    text_coords = []
    for i in range(len(block_seq)):
        block_idx = block_seq[i]
        width = np.log10(9 + len(blocks[block_idx]))
        arrow, text_pos = arrow_polygon(inverted_seq[i], width, left_pos, bottom_pos)
        poly = Polygon(arrow, linewidth=4)#, fc = colors[block_idx], ec = 'black', alpha = .5) #
        if blocks[block_idx].is_insertion():
            patches_ins.append(poly)
        else:
            patches.append(poly)
        text_coords.append(text_pos)
        left_pos += width

    # regular blocks
    nblocks = len(patches)
    scaled_block_seq = [round((b - start_block) / (end_block - start_block) * 255) for b in block_seq if not blocks[b].is_insertion()]
    fc = [jet(x) for x in scaled_block_seq]
    ec = ['black'] * nblocks
    p = PatchCollection(patches, facecolor=fc, edgecolor=ec)
    # insertions
    nins = len(patches_ins)
    fc_ins = ['gray'] * nins
    ec_ins = ['black'] * nins
    p_ins = PatchCollection(patches_ins, facecolor = fc_ins, edgecolor = ec_ins,
                            hatch = '/')

    right_pos = left_pos

    return p, p_ins, text_coords, right_pos

def write_block_labels(text_coords, path, blocks, start):
    block_seq = [int(floor(path[i]/2)) for i in range(0, len(path), 2)]
    inverted_seq = [path[i] % 2 == 1 for i in range(0, len(path), 2)]
    for i in range(len(text_coords)):
        (x,y) = text_coords[i]
        block = block_seq[i]
        if blocks[block].is_insertion():
            letter = '='
        else:
            letter = chr(65 + block - start)
        txt = plt.text(x, y, letter,
                       color = 'white', size = 25, ha='center', va='center')
        txt.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'),
                              path_effects.Normal()])

def write_path_label(bottom_pos, label):
    plt.text(-1.5, bottom_pos + PHI/2 - 1/8 + .065, label, color = 'black', size = 30, ha = 'right', va = 'center')

def plot_rearrangement(filename, blocks, start_block, end_block,
                       path1, path2 = None, show_ref = None):
    fig, ax = plt.subplots()
    plt.axis('image')
    plt.axis('off')

    right_pos = 0
    bottom_pos = 0

    # do ref
    if show_ref:
        ref_path = list(range(start_block * 2, end_block * 2 + 2))
        p_ref, _, text_coords, ref_right_pos = get_arrow_patches(blocks, ref_path, bottom_pos,
                                                                 start_block, end_block)
        ax.add_collection(p_ref)
        write_block_labels(text_coords, ref_path, blocks, start_block)
        write_path_label(bottom_pos, 'REF')
        right_pos = max(right_pos, ref_right_pos)
        bottom_pos -= 3

    paths = (path1,) if path2 is None else (path1, path2)
    pathnum = 1
    for path in paths:
        p, p_ins, text_coords, p_right_pos = get_arrow_patches(blocks, path, bottom_pos,
                                                               start_block, end_block)
        ax.add_collection(p)
        ax.add_collection(p_ins)
        write_block_labels(text_coords, path, blocks, start_block)
        if len(paths) == 1:
            write_path_label(bottom_pos, 'ALT'.format(pathnum))
        else:
            write_path_label(bottom_pos, 'ALT {0}'.format(pathnum))
        right_pos = max(right_pos, p_right_pos)
        bottom_pos -= 3
        pathnum += 1

    ax.set_ylim((bottom_pos+(3-PHI),2))
    ax.set_xlim((-6, right_pos + 1/2))
    bottom_pos = 0
    fig = plt.gcf()
    axsize = plt.axis()
    # at dpi 150, prevent plot size larger than ~32768 pixels (o/w crash)
    dpi = 150
    fig_width = min(218, .45*(axsize[1]-axsize[0]))
    fig_height = min(218, .45*(axsize[3]-axsize[2]))
    fig.set_size_inches(fig_width, fig_height)
    plt.savefig(filename, dpi = dpi)
    plt.close()

def test_plot_rearrangement():
    blocks = [GenomeInterval('1',0,1000), GenomeInterval('1',1010, 1012),#1500),
              GenomeInterval('1',1505,2000), GenomeInterval('1',2000,4000),
              GenomeInterval('1',4000,20000), GenomeInterval('1',0,10,is_de_novo = True),
              GenomeInterval('1',0,1000,is_de_novo = True)]
    outdir = '~/tmp/'

    p1 = [0,1,4,5,6,7]
    fn = 'ACD.png'
    print(fn)
    plot_rearrangement(os.path.join(outdir, fn), blocks, 0, 4, p1, None, True)

    p1 = [0,1,10,11,4,5,6,7]
    fn = 'AICD.png'
    print(fn)
    plot_rearrangement(os.path.join(outdir, fn), blocks, 0, 4, p1, None, True)

    p1 = [0,1,2,3,3,2,6,7,8,9]
    p2 = [0,1,4,5,2,3,4,5,6,7,8,9]
    fn = 'ABB-DE_ACBCDE.png'
    print(fn)
    plot_rearrangement(os.path.join(outdir, fn), blocks, 0, 4, p1, p2, True)

    p1 = [0,1,2,3,3,2,2,3,3,2,2,3,2,3,6,7,7,6,8,9]
    p2 = p2
    fn = 'ABB--.png'
    print(fn)
    plot_rearrangement(os.path.join(outdir, fn), blocks, 0, 4, p1, p2, True)

    blocks = [GenomeInterval('1',i,i+100) for i in range(0,2000,100)]
    p1 = list(range(0,39))
    p2 = p2
    fn = 'ABCD---.png'
    print(fn)
    plot_rearrangement(os.path.join(outdir, fn), blocks, 0, 19, p1, p2, True)
