ApaI = 'GGGCCC'
BamHI = 'GGATCC'
KpnI = 'GGTACC'
NotI = 'GCGGCCGC'
SacII = 'CCGCGG'
SmaI = 'CCCGGG'
SpeI = 'ACTAGT'
XmaI = 'CCCGGG'



# ---------- CUSTOMIZE HERE ----------

# Number of introns to insert. Should be a positive integer.
intronsToInsert = 3

# Whether to insert the enhancer. Should be true or false.
insertEnhancer = False

# List of rescrition sites to prevent, as defined above.
restrictionSitesToPrevent = [ApaI, BamHI, KpnI, NotI, SpeI, XmaI]  

# Sequence to optimize. Can be DNA, RNA, or protein. Non-letters are ignored.
sequence = '''
        METDTLLLWVLLLWVPGSTG
DRSANDTVVVGSIIFTEGII
VANMVAEMIEAHTDLKVVRK
LNLGGVNVNFEAIKRGGANN
GIDIYVEYTGHGLVDILGFP
EPNVYITADKQKNGIKANFK
IRHNVEDGSVQLADHYQQNT
PIGDGPVLLPDNHYLSTQSV
LSKDPNEKRDHMVLLEFVTA
AGITLGMDELYKGGTGGSMS
KGEELFTGVVPILVELDGDV
NGHKFSVRGEGEGDATNGKL
TLKFICTTGKLPVPWPTLVT
TLTYGVQCFSRYPDHMKQHD
FFKSAMPEGYVQERTISFKD
DGTYKTRAEVKFEGDTLVNR
IELKGIDFKEDGNILGHKLE
YNFPPPATTDPEGAYETVKK
EYKRKWNIVWLKPLGFNNTY
TLTVKDELAKQYNLKTFSDL
AKISDKLILGATMFFLEGPD
GYPGLQKPYNFKFKHTKSMD
MGIRYTAIDNNEVQVIDAWA
TDGLLVSHKLKILEDDKAFF
PPYYAAPIIRQDVLDKHPEL
KDVLNKLANQISLEEMQKLN
YKVDGEGQDPAKVAKEFLKE
KGLILQVDEQKLISEEDLNA
VGQDTQEVIVVPHSLPFKVV
VISAILALVVLTIISLIILI
MLWQKKPR*
'''


# ------------------------------------




import re

# Define codons.
codons = {'A': ['GCT', 'GCC', 'GCA', 'GCG'],
          'C': ['TGC', 'TGT'],
          'D': ['GAC', 'GAT'],
          'E': ['GAG', 'GAA'],
          'F': ['TTC', 'TTT'],
          'G': ['GGA', 'GGT', 'GGC', 'GGG'],
          'H': ['CAC', 'CAT'],
          'I': ['ATC', 'ATT', 'ATA'],
          'K': ['AAG', 'AAA'],
          'L': ['CTC', 'CTT', 'TTG', 'CTG', 'CTA', 'TTA'],
          'M': ['ATG'],
          'N': ['AAC', 'AAT'],
          'P': ['CCA', 'CCG', 'CCT', 'CCC'],
          'Q': ['CAG', 'CAA'],
          'R': ['CGT', 'CGC', 'AGA', 'CGA', 'AGG', 'CGG'],
          'S': ['TCC', 'TCT', 'TCA', 'TCG', 'AGC', 'AGT'],
          'T': ['ACC', 'ACT', 'ACA', 'ACG'],
          'V': ['GTC', 'GTT', 'GTG', 'GTA'],
          'W': ['TGG'],
          'Y': ['TAC', 'TAT'],
          '*': ['TAA', 'TAG', 'TGA']}

intron = 'gtaagtttaaacatgattttactaactaactaatctgatttaaattttcag'
enhancer = 'CCCGGGATTGGCCAAAGGACCCAAAGgtatgtttcgaatgatactaacataacatagaacattttcagGAGGACCCTTGGCTAGCGTCGACGgatccaaaaa'

clean = lambda x: ''.join(s for s in x.upper() if s.isalpha())
restrictionSitesToPrevent = map(clean, restrictionSitesToPrevent)
sequence = clean(sequence)


# Translate sequence first if DNA or RNA.
def translate(seq):
    reverseCodons = {v: k for k, vs in codons.items() for v in vs}
    return ''.join([reverseCodons[seq[i*3:i*3+3]] for i in range(int(len(seq)/3))])
if not set(sequence) - set('ACGTU'):
    print('---\nWARNING: Sequence is ' + ('RNA' if 'U' in sequence else 'DNA') + '. Translating..')
    sequence = translate(sequence.replace('U', 'T'))


# Codon optimization.
print('---\n 1) Codon optimization:')
import random
sequenceWorm = ''
for i, aa in enumerate(sequence):
    k = random.randint(0, 1)
    if len(codons[aa]) < 2:
        k = 0
    sequenceWorm += codons[aa][k]
    # Check for undesirable restriction sites.
    for site in restrictionSitesToPrevent:
        if site in sequenceWorm:
            for backIdx in range(4):
                if i-backIdx < 0:
                    continue
                for aa2 in codons[sequence[i-backIdx]][1:]:
                    tempSeq = sequenceWorm[:-3*(backIdx+1)] + aa2 + (sequenceWorm[-3*(backIdx+1)+3:] if backIdx else '')
                    if not site in tempSeq:
                        sequenceWorm = tempSeq
                        printArgs = [i-backIdx+1, sequence[i-backIdx], aa2, codons[sequence[i-backIdx]][0], site]
                        print('At residue {0} ({1}): {2} used instead {3} to prevent {4}.'.format(*printArgs))
                        break
                else:
                    continue
                break
            else:
                print(':::\n::: WARNING: Could not prevent formation of the site:', site, '\n:::')
print('Length of transcript:', len(sequenceWorm))

# Insert introns.
if intronsToInsert:
    print('---\n 2) Inserting introns:')
    if any([s in intron.upper() for s in restrictionSitesToPrevent]):
        print(':::\n::: WARNING: Could not add introns, as the restriction site', site, 'is the intron.\n:::')
    else:
        ggPositions = [m.start() for m in re.finditer('(?=GG)', sequenceWorm) if not (m.start()+1) % 3]
        ggs = len(ggPositions)
        exonLength = len(sequenceWorm)/(intronsToInsert+1)
        tempSeq, last, insertedAt = '', 0, []
        for i in range(intronsToInsert):
            insertAt = (i+1)*exonLength
            while len(ggPositions):
                closest = min(ggPositions, key=lambda x:abs(x-insertAt))
                ggPositions.remove(closest)
                testSeq = tempSeq + sequenceWorm[last:closest+1] + intron + sequenceWorm[closest+1:]
                # Check for undesirable restriction sites.
                for site in restrictionSitesToPrevent:
                    if site in testSeq.upper() + sequenceWorm[last:]:
                        break
                else:
                    tempSeq += sequenceWorm[last:closest+1] + intron
                    break
            else:
                print(':::\n::: WARNING: Did not find enough available GGs to insert introns.\n:::')
                break
            
            last = closest+1
            insertedAt.append(last)
        tempSeq += sequenceWorm[last:]
        sequenceWorm = tempSeq
    
        print('Found', ggs, 'GGs in frame.')
        print('Inserted intron{0} at: {1} (of {2})'.format('s' if intronsToInsert > 1 else '',
                                                           ', '.join([str(a) for a in insertedAt]),
                                                           len(sequenceWorm)))

# Insert enhancer.
if insertEnhancer:
    print('---\n 3) Inserting enhancer:')
    if any([s in enhancer.upper() for s in restrictionSitesToPrevent]):
        print(':::\n::: WARNING: Could not add enhancer, as the restriction site', site, 'is the enhancer.\n:::')
    else:
        insertAt = sequenceWorm.find('ATG')
        tempSeq = sequenceWorm[:insertAt] + enhancer + sequenceWorm[insertAt:]
        # Check if insertion caused undesirable restriction sites.
        if any([s in tempSeq.upper() for s in restrictionSitesToPrevent]) and not any([s in sequenceWorm.upper() for s in restrictionSitesToPrevent]):
            print(':::\n::: WARNING: Adding the enhancer causes the restriction site', site, 'to occur.\n:::')
        else:
            sequenceWorm = tempSeq
            print('Inserted at:', insertAt, '\nComplete.')

# Check that protein sequence did not change.
if sequence != translate(sequenceWorm.replace(intron, '').replace(enhancer, '')):
    print('ERROR: protein sequence changed! Something is wrong with the script!')
    print(re.sub(r'(.{10})', r'\1 ', sequenceWorm))
    print('ERROR: protein sequence changed! Something is wrong with the script!')
else:
    print('---\nFinal sequence:')
    print(re.sub(r'(.{10})', r'\1 ', sequenceWorm))




