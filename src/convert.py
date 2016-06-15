from __future__ import print_function
import re
import glob
import fileinput

pats = {
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.vars\[([\w\-\+\s*]+)\]': '\\1.var(\\2, \\3, \\4)',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.fluxes\[([\w\-\+\s*]+)\]': '\\1.flux(\\2, \\3, \\4)',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.homo_fluxes\[([\w\-\+\s*]+)\]': '\\1.homo_flux(\\2, \\3, \\4)',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.hetero_fluxes\[([\w\-\+\s*]+)\]': '\\1.hetero_flux(\\2, \\3, \\4)',
	r'SMC_cell\s*\*\*': 'SMC_type& ',
	r'EC_cell\s*\*\*': 'EC_type& ',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.JPLC': '\\1.JPLC(\\2, \\3)',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.NO': '\\1.NO(\\2, \\3)',
	r'(\w+)\[([\w\-\+\s*]+)\]\[([\w\-\+\s*]+)\]\.NE': '\\1.NE(\\2, \\3)',
}


for line in fileinput.input():
    newline = line[:]
    for old, new in pats.items():
        newline = re.sub(old, new, newline)
    print(newline, end='')
