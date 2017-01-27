import sys
import glob

letters = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
           'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
           'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
           'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
           'TYR': 'Y', 'VAL': 'V'}

pdb_files = sorted(glob.glob('../data/*.pdb'))
with open('active_sites_combined.fasta', 'w') as out_file:
    for pdb in pdb_files:
        with open(pdb, 'r') as f:
            name = pdb[8:-4]
            # print(name)
            out_file.write('>' + name + '\n')
            prev = '-1'
            prev_chain = '-1'
            loop_count = 0
            for line in f:
                toks = line.split()
                curr_chain = line[21:22]
                if loop_count == 0:
                    prev_chain = curr_chain
                if len(toks) < 1:
                    continue
                if line[:4] != 'ATOM':
                    continue
                if line[22:26] != prev and (prev_chain == curr_chain):
                    print(line[21:22], prev_chain, curr_chain, loop_count)
                    out_file.write('%c' % letters[line[17:20]])
                    prev_chain = line[21:22]
                    print('aas written')
                prev = line[22:26]
                loop_count += 1
                print(loop_count)
            out_file.write('\n')
