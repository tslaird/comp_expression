#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:26:05 2019
@author: tslaird
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 15:12:11 2019
@author: tslaird
"""

import re
import textwrap
import multiprocessing as mp
import glob
import os
import sys

def gbff2transcripts(gbff_file):
    locus_re= re.compile('(?<=LOCUS\s{7})\w+')
    definition_re = re.compile('(?<=DEFINITION\s{2})[\s\S]+?(?=\nACC)')
    definition_sub = re.compile('\n|(\s{12})')
    locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?locus_tag[\S\s]+?(?="\n)')
    locus_tag_sub = re.compile('[\s\S]+"')
    old_locus_tag_re = re.compile('(?<=gene\W{12})[\S\s]+?old_locus_tag[\S\s]+?(?="\n)')
    product_re = re.compile('(?<=product=")[\S\s]+?(?=")')
    coords_re = re.compile('(?<=(gene)(\s|A){12})\S+')
    protein_id_re = re.compile('(?<=protein_id=")[\s\S]+?(?=")')
    protein_re = re.compile('(?<=translation=")[\S\s]+?(?=")')
    biosample_re = re.compile('(?<=BioSample:\s)\w+')
    assembly_re = re.compile('(?<=Assembly:\s)\S+')
    features_re = re.compile('(?<=\n)\W{5}(?=gene\W{12})')
    product_sub = re.compile('\n|\s{19}')
    protein_out_sub = re.compile('\n|(\s+)')
    separate_re = re.compile('//\n')
    pseudogene_re = re.compile('(?<=\/)pseudo(?=\n)')
    with open(gbff_file) as file:
        file_text=file.read()
    loci = separate_re.split(file_text)
    all_features = []
    for i in loci:
        if locus_re.search(i):
            locus = locus_re.search(i).group(0)
            definition = definition_re.search(i).group(0)
            definition = definition_sub.sub('',definition)
            biosample =  biosample_re.search(i)
            if biosample:
                biosample_out = biosample.group(0)
            else:
                biosample_out = 'NULL'
            if assembly_re.search(i):
                assembly = assembly_re.search(i).group(0)
            else:
                assembly = re.search('GCF_\d+?(?=\.)',gbff_file).group(0)
            features = features_re.split(i)
            dna_raw=re.search('(?<=ORIGIN)[\s\S]+',i).group(0)
            dna_clean= re.sub('\d|\W','',dna_raw).upper()
            for f in features[1:]:
                 locus_tag = locus_tag_re.search(f).group(0)
                 locus_tag = locus_tag_sub.sub('', locus_tag)
                 if old_locus_tag_re.search(f):
                     old_locus_tag = old_locus_tag_re.search(f).group(0)
                     old_locus_tag = locus_tag_sub.sub('', old_locus_tag)
                 else:
                     old_locus_tag = 'NULL'
                 if product_re.search(f):
                     product = product_re.search(f).group(0)
                     product = product_sub.sub('',product)
                 else:
                     product = 'NULL'
                 coords = coords_re.search(f).group(0)
                 start= int(re.findall('\d+',coords)[0])
                 end= int(re.findall('\d+',coords)[1])
                 if re.search('complement',coords):
                     cds_raw= dna_clean[start-1:end][::-1]
                     rev_comp=''
                     for base in cds_raw:
                         if base == 'G':
                             rev_comp += 'C'
                         if base == 'C':
                             rev_comp += 'G'
                         if base == 'A':
                             rev_comp += 'T'
                         if base == 'T':
                             rev_comp += 'A'
                     cds=rev_comp
                 else:
                     cds= dna_clean[start-1:end]
                 cds= textwrap.fill(cds, 80)
                 protein_id = protein_id_re.search(f)
                 if pseudogene_re.search(f):
                     pseudogene = "PSEUDOGENE"
                 else:
                     pseudogene = "NULL"
                 if protein_id:
                     protein_id_out = protein_id.group(0)
                 else:
                     protein_id_out = 'NULL'
                 protein = protein_re.search(f)
                 if protein:
                     protein_out = protein.group(0)
                     protein_out = protein_out_sub.sub('', protein_out)
                     protein_out = textwrap.fill(protein_out, 70)
                     almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (assembly,locus,locus_tag,old_locus_tag,coords,product, protein_id_out,pseudogene)
                     almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                     all_items = "%s\n%s\n" % (almost_all_items, cds)
                 else:
                     almost_all_items = ">%s!!%s!!%s!!%s!!%s!!%s!!%s!!%s" % (assembly,locus,locus_tag,old_locus_tag,coords,product,protein_id_out,pseudogene)
                     almost_all_items = re.sub('\n|\s{4,}','',almost_all_items)
                     all_items = "%s\n%s\n" % (almost_all_items, cds)
                 #all_proteins.add(">"+assembly+'•'+locus+'•'+locus_tag+'•'+biosample+'•'+product+'•'+coords+'•'+protein_id+'\n'+protein+'\n')
                 all_features.append(all_items)
    result= '\n'.join(all_features)+'\n'
    result = re.sub(' ','_', result)
    outname = [re.sub('.gbk|.gbff','_transcripts.fasta', gbff_file)]
    outname = ''.join(outname)
    with open(outname,'w+') as output:
        output.writelines(result)
    print("made transcriptome for: " + gbff_file)

input_file = sys.argv[1]

gbff2transcripts(input_file)
