#!/usr/bin/env python
import os
import sys
import timeit

usage = '\t --------\n' \
        '\t| usage  : python sort_plant_hit_vs_plant_refseq_prot.py f1 f2\n' \
        '\t| input  : f1 = blast.tsv (vs plant_ref_seq)\n' \
        '\t| input  : f2 = seqs.fasta (blastx queries)\n' \
        '\t| indir  : prot.accession2taxid.gz (to map saccver to staxids)\n' \
        '\t| output : "f2"_1.fasta (plant_hit)\n' \
        '\t| output : "f2"_2.fasta (non_plant_hit = a control here, should be empty)\n' \
        '\t| output : "f2"_3.fasta (no_hit)\n' \
        '\t --------'

if len(sys.argv) != 3:
    print(usage)
    sys.exit()

##############
### Step 1 ###
##############
print('\n\tStep 1) Map saccver to staxids')
t0 = timeit.default_timer()

cpt = 0
qseqid_set, saccver_set = set(), set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        cpt += 1
        columns = row.split('\t')
        qseqid, saccver = columns[0].rstrip(), columns[14].rstrip()
        # Get rid of accession number version
        saccver = saccver.split('.')[0]
        if not qseqid in qseqid_set:
            saccver_set.add(saccver)
        qseqid_set.add(qseqid)

# Use saccver_set as local query against prot.accession2taxid
# create saccver.temp as input for grep
with open('saccver.temp', 'w') as temp:
    for item in saccver_set:
        temp.write(item.rstrip() + "\n")

# map saccver to staxids
cmd = 'zfgrep -f saccver.temp ../Data/all_prot.accession2taxid.gz'
cmd_result = os.popen(cmd).read()
# First attempt
# cmd = ('zgrep -f saccver.temp prot.accession2taxid.gz | awk \'{print $2\"\t\"$3}\' $_')
# grep on unzipped file is faster (3.42 vs 5.28), but file is too big to please me
# cmd = 'fgrep -f saccver.temp all_prot.accession2taxid'

# store grep result in a temp file (easier to parse)
with open('grep.temp', 'w') as temp:
    temp.write(cmd_result)

# Finally, map saccver to staxids in a dic
saccver_to_staxids_dic = {}
with open('grep.temp', 'r') as grep:
    for line in grep:
        field = line.split('\t')
        saccver_to_staxids_dic.setdefault(field[0].rstrip(), field[2].rstrip())

print('\t\t- TSV have ' + str(cpt) + ' lines, containing ' +
      str(len(saccver_set)) + ' unique saccver')
print('\t\t- grep result contains ' + str(len(saccver_set)) + ' lines')
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 2 ###
##############
print('\tStep 2) Retrieve taxonomic infos for each staxids with efetch')
t0 = timeit.default_timer()

# Use staxids_set as query with efetch & store result
# Don't give more than let's say 500 entries at a time to avoid timeout
staxids_set = set()
for k, v in saccver_to_staxids_dic.items():
    staxids_set.add(v)

staxids_li = list(staxids_set)
staxids_sub_li = [staxids_li[x:x + 500]
                  for x in range(0, len(staxids_li), 500)]
efetch_li = []
for item in staxids_sub_li:
    staxids_input = ','.join(str(z) for z in item)
    # Details about "cmd" : https://www.biostars.org/p/163595/#271497
    cmd = ('efetch -db taxonomy -id ' + staxids_input + ' -format xml | xtract '
           '-pattern Taxon -sep \'@\' -element TaxId,ScientificName -division '
           'LineageEx -group Taxon -if Rank -equals superkingdom -or Rank '
           '-equals kingdom -or Rank -equals phylum -or Rank -equals class'
           ' -or Rank -equals order -or Rank -equals family -or Rank -equals'
           ' genus -sep \'@\' -element Rank,ScientificName')
    cmd_result = os.popen(cmd).read()
    cmd_result_split = cmd_result.split('\n')
    for i in cmd_result_split:
        efetch_li.append(i)
# 257314@Lactobacillus johnsonii NCC
# 533\tsuperkingdom@Bacteria\tphylum@Firmicutes\tclass@Bacilli\torder@Lactobacillales\tfamily@Lactobacillaceae\tgenus@Lactobacillus'

# Create a dict associating key=staxid with value=list=tax_infos
taxonomy_dic = {}
for line in efetch_li:
    field = line.split('\t')
    tax_ids = field[0].split('@')
    # Sometimes more than one staxid is associated to an entry
    # e.g. "170850@3666@Cucurbita hybrid cultivar"
    for i in tax_ids[:-1]:
        taxonomy_dic.setdefault(i, [None, None, None, None, None, None, None])
        for item in field:
            if 'superkingdom@' in item:
                taxonomy_dic[i][0] = item.split('@')[-1]
            elif 'kingdom@' in item:
                taxonomy_dic[i][1] = item.split('@')[-1]
            elif 'phylum@' in item:
                taxonomy_dic[i][2] = item.split('@')[-1]
            elif 'class@' in item:
                taxonomy_dic[i][3] = item.split('@')[-1]
            elif 'order@' in item:
                taxonomy_dic[i][4] = item.split('@')[-1]
            elif 'family@' in item:
                taxonomy_dic[i][5] = item.split('@')[-1]
            elif 'genus@' in item:
                taxonomy_dic[i][6] = item.split('@')[-1]
# '41840': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Sphagnopsida', 'Sphagnales', 'Sphagnaceae', 'Sphagnum']

print('\t\t- taxonomy_dic len is ' + str(len(taxonomy_dic)))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 3 ###
##############
print('\tStep 3) Assign contigs best hit to plant or non-plant')
t0 = timeit.default_timer()

# Assign contigs best hits to plant or non-plant based on taxonomy_dic infos
qseqid_set, viridi_hit_set, non_viridi_hit_set = set(), set(), set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        columns = row.split('\t')
        qseqid, saccver = columns[0].rstrip(), columns[14].rstrip()
        # Get rid of accession number version
        saccver = saccver.split('.')[0]
        # Check if we encounter qseqid for the first time <=> best hit
        if not qseqid in qseqid_set:
            if taxonomy_dic[saccver_to_staxids_dic[saccver]][1] == 'Viridiplantae':
                viridi_hit_set.add(qseqid)
            else:
                non_viridi_hit_set.add(qseqid)
        qseqid_set.add(qseqid)

print('\t\t- qseqid_set len is ' + str(len(qseqid_set)))
print('\t\t- viridi_hit_set len is ' + str(len(viridi_hit_set)))
print('\t\t- non_viridi_hit_set len is ' + str(len(non_viridi_hit_set)))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 4 ###
##############
print('\tStep 4) Find contigs with no hits')
t0 = timeit.default_timer()

# Read initial FASTA (f2) & check intersection with viridi_hit_set & non_viridi_hit_set
# We can deduce contigs with no hit from this intersection
cpt = 0
no_hit_set = set()
with open(sys.argv[2], 'r') as fa:
    for line in fa:
        if line.startswith('>'):
            cpt += 1
            line = line.lstrip('>')
            # if line in no_hit_set:
            #    print(line)
            # no_hit_set.add(line)
            fields = line.split()
            no_hit_set.add(fields[0])

before_union = len(no_hit_set)
no_hit_set = no_hit_set - viridi_hit_set
no_hit_set = no_hit_set - non_viridi_hit_set

print('\t\t- number of seqs (>) in FASTA : ' + str(cpt))
print('\t\t- number of headers added to initial set is ' + str(before_union))
print('\t\t- number of res in no_hit_set : ' + str(len(no_hit_set)))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')


##############
### Step 5 ###
##############
print('\tStep 5) Create output files')
t0 = timeit.default_timer()

# Create input files (sequence IDs list) for seqtk
file_2 = sys.argv[2].split('/')
sample = file_2[-1].split('.')[0]
with open(sample + '_plant_hit.temp', 'w') as out:
    for item in viridi_hit_set:
        out.write(item + "\n")
with open(sample + '_non_plant_hit.temp', 'w') as out:
    for item in non_viridi_hit_set:
        out.write(item + "\n")
with open(sample + '_no_hit.temp', 'w') as out:
    for item in no_hit_set:
        field = item.split()
        # out.write(item)
        out.write(field[0] + '\n')

# Create output files (plant_hit = 1, non-plant_hit = 2, no-hit = 3) with seqtk
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_plant_hit.temp > ' + sample + '_1.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_non_plant_hit.temp > ' + sample + '_2.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_no_hit.temp > ' + sample + '_3.fasta')

sum_seqs = int(len(viridi_hit_set)) + \
    int(len(non_viridi_hit_set)) + int(len(no_hit_set))
print('\t\t- number of seqs with plant_hit : ' + str(len(viridi_hit_set)))
print('\t\t- number of seqs in non_plant_hit : ' + str(len(non_viridi_hit_set)))
print('\t\t- number of seqs with no_hit : ' + str(len(no_hit_set)))
print('\t\t- sum of the 3 above : ' + str(sum_seqs) +
      ', which should be equal to : ' + str(cpt))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds')

os.system('rm *.temp')
