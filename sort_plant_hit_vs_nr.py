#!/usr/bin/env python
import os
import sys
import timeit

usage = '\t --------\n' \
        '\t| usage  : python sort_plant_hit.py f1 f2\n' \
        '\t| input  : f1 = blast.tsv\n' \
        '\t| input  : f2 = seqs.fasta\n' \
        '\t| output : f2_plant_hit.fasta\n' \
        '\t| output : f2_non_plant_hit.fasta\n' \
        '\t| output : f2_no_hit.fasta\n' \
        '\t --------'

if len(sys.argv) != 3:
    print(usage)
    sys.exit()

##############
### Step 1 ###
##############
print('\n\tStep 1) Retrieve taxonomic infos for each staxids with efetch')
t0 = timeit.default_timer()
# For each line in TSV (f1), fill staxids_set
staxids_set = set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        columns = row.split('\t')
        # Sometimes you have more than one staxids for an entry
        staxids = columns[15].split(';')
        for i in staxids:
            staxids_set.add(i)

# Use staxids_set as query with efetch & store result
# Don't give more than let's say 500 entries at a time to avoid errors
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
# 257314@Lactobacillus johnsonii NCC 533\tsuperkingdom@Bacteria\tphylum@Firmicutes\tclass@Bacilli\torder@Lactobacillales\tfamily@Lactobacillaceae\tgenus@Lactobacillus'

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
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds')
# '41840': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Sphagnopsida', 'Sphagnales', 'Sphagnaceae', 'Sphagnum']

##############
### Step 2 ###
##############
print('\tStep 2) Assign contigs best hits to plant or non-plant')
t0 = timeit.default_timer()
# Assign contigs best hits to plant or non-plant based on taxonomy_dic infos
qseqid_set, viridi_hit_set, non_viridi_hit_set = set(), set(), set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        columns = row.split('\t')
        qseqid, staxids = columns[0], columns[15].split(';')[0]
        # Check if we encounter qseqid for the first time <=> best hit
        if not qseqid in qseqid_set:
            if taxonomy_dic[staxids][1] == 'Viridiplantae':
                viridi_hit_set.add(qseqid)
            else:
                non_viridi_hit_set.add(qseqid)
        qseqid_set.add(qseqid)
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds')

##############
### Step 3 ###
##############
print('\tStep 3) Find contigs with no hits')
t0 = timeit.default_timer()
# Read initial FASTA (f2) & check intersection with viridi_hit_set & non_viridi_hit_set
# We can deduce contig with no hit from this intersection
no_hit_set = set()
cpt = 0
with open(sys.argv[2], 'r') as fa:
    for line in fa:
        if line.startswith('>'):
            cpt += 1
            line = line.lstrip('>')
            fields = line.split()
            no_hit_set.add(fields[0])
# Some colisions to solves here (less seqid than number of lines with >, wtf)
no_hit_set = no_hit_set - viridi_hit_set
no_hit_set = no_hit_set - non_viridi_hit_set
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds')

##############
### Step 4 ###
##############
print('\tStep 4) Create output files')
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
        out.write(item + "\n")

# Create output files (plant_hit, non-plant_hit, no-hit) with seqtk
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_plant_hit.temp > ' + sample + '_1.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_non_plant_hit.temp > ' + sample + '_2.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_no_hit.temp > ' + sample + '_3.fasta')
os.system('rm *.temp')
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds')
