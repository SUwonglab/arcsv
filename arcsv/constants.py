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

# softclip orientation constants
LEFT = 0
RIGHT = 1

# Low base quality symbols
LOWQUAL_CHARS = ('#', '!', '"')
