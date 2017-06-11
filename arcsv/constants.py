# misc. constants
ALTERED_QNAME_MAX_LEN = 256
MAX_MAPQ = 60

SPLIT_FIRST_PLUS = 0x1
SPLIT_SECOND_PLUS = 0x2
SPLIT_OVERLAP = 0x4
SPLIT_LEFT_FIRST = 0x8
SPLIT_TYPES = {(SPLIT_FIRST_PLUS | SPLIT_SECOND_PLUS): 'Del+',
               0x0: 'Del-',
               (SPLIT_FIRST_PLUS | SPLIT_SECOND_PLUS | SPLIT_OVERLAP): 'Dup+',
               SPLIT_OVERLAP: 'Dup-',
               (SPLIT_FIRST_PLUS | SPLIT_LEFT_FIRST): 'InvL+',
               SPLIT_FIRST_PLUS: 'InvL-',
               (SPLIT_SECOND_PLUS | SPLIT_LEFT_FIRST): 'InvR+',
               SPLIT_SECOND_PLUS: 'InvR-'}

# Junction structure constants
LEFT = 0
RIGHT = 1
SEQ = 0
QUAL = 1
REFSEQ = 2
ORIENT = 3
BPLOC = 4
NCLIP = 5
NUNIQ = 6
NDUP = 7
MAPQ = 8
NSUPP = 9
LIBS = 10
