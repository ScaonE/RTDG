#!/usr/bin/env python
import os
import sys
import timeit

usage = '\t --------\n' \
        '\t| usage  : python sort_plant_hit_vs_nt.py f1 f2\n' \
        '\t| input  : f1 = blastn.tsv (vs nt)\n' \
        '\t| input  : f2 = seqs.fasta (blastn queries)\n' \
        '\t| output : "f2"_1.fasta (plant_hit)\n' \
        '\t| output : "f2"_2.fasta (non_plant_hit)\n' \
        '\t| output : "f2"_3.fasta (no_hit)\n' \
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
# Don't give more than let's say 500 entries at a time to avoid timeout
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

print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 2 ###
##############
print('\tStep 2) Assign contigs best hit to plant or non-plant')
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

print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 3 ###
##############
print('\tStep 3) Find contigs with no hits')
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
            fields = line.split()
            no_hit_set.add(fields[0])

before_union = len(no_hit_set)
no_hit_set = no_hit_set - viridi_hit_set
no_hit_set = no_hit_set - non_viridi_hit_set

print('\t\t- number of seqs (>) in FASTA : ' + str(cpt))
print('\t\t- number of headers added to initial set is ' + str(before_union))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

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

# Create output files (plant_hit = 1, non-plant_hit = 2, no-hit = 3) with seqtk
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_plant_hit.temp > ' + sample + '_1.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_non_plant_hit.temp > ' + sample + '_2.fasta')
os.system('seqtk subseq ' + sys.argv[2] + ' ' + sample +
          '_no_hit.temp > ' + sample + '_3.fasta')

sum_seqs = int(len(viridi_hit_set)) + int(len(non_viridi_hit_set)) + int(len(no_hit_set))
print('\t\t- number of seqs with plant_hit : ' + str(len(viridi_hit_set)))
print('\t\t- number of seqs in non_plant_hit : ' + str(len(non_viridi_hit_set)))
print('\t\t- number of seqs with no_hit : ' + str(len(no_hit_set)))
print('\t\t- sum of the 3 above : ' + str(sum_seqs))
print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

os.system('rm *.temp')
