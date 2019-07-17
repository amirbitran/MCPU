"""
Counts contacts in PDB file specified by first argument
"""

import numpy as np
import Analyze_structures
import sys

native_file = sys.argv[0]

coords, resis=Analyze_structures.read_PDB(native_file, 'CA')
native_contacts=Analyze_structures.compute_contacts_matrix(coords, thresh=d_cutoff, min_seq_separation=min_seq_separation)
print(np.sum(native_contacts))

