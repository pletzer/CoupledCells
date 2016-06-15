import re
import glob
import fileinput

pats = {
	r'(\w+)\[(\w+)\]\[(\w+)\]\.vars\[([\w\-\s*]+)\]': '\\1.var(\\2, \\3, \\4)',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.fluxes\[([\w\-\s*]+)\]': '\\1.flux(\\2, \\3, \\4)',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.homo_fluxes\[([\w\-\s*]+)\]': '\\1.homo_flux(\\2, \\3, \\4)',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.hetero_fluxes\[([\w\-\s*]+)\]': '\\1.hetero_flux(\\2, \\3, \\4)',
	r'SMC_cell\s*\*\*': 'SMC_type& ',
	r'EC_cell\s*\*\*': 'EC_type& ',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.JPLC': '\\1.JPLC(\\2, \\3)',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.NO': '\\1.NO(\\2, \\3)',
	r'(\w+)\[(\w+)\]\[(\w+)\]\.NE': '\\1.NE(\\2, \\3)',
}


for line in fileinput.input():
    newline = line[:]
    for old, new in pats.items():
        newline = re.sub(old, new, newline)
    print(newline, end='')
