#!/usr/bin/python

# usage:

# multiGtfHTSeq.py 
#    -o output_suffix [default=mHTSeq]
#    --gr1=<gtf_file_of_group1>,<gtf_file_of_group1> 
#    --gr2=<gtf_file_of_group2>,<gtf_file_of_group2> 
#    --gr3=<gtf_file_of_group3>,<gtf_file_of_group3> 
#    <inputSAMfile>

# example:

# multiGtfHTSeq.py --gr1=/data/BA/sdepreti/annotation/mm9/exons_mm9.gtf,/data/BA/sdepreti/annotation/mm9/tRNA_genes_mm9.gtf,/data/BA/sdepreti/annotation/mm9/miRNA_genes_mm9.gtf --gr2=/data/BA/sdepreti/annotation/mm9/introns_mm9.gtf --gr3=/data/BA/sdepreti/annotation/mm9/repmask_mm9.gtf prova.sam

# this software count reads aligned ( in sam format ) using HTSeq python package (intersection
# non-empty mode) over many different gtf files. if the gtf files belong to the same group, 
# reads that overlap  features from them are counted as ambiguos. If gtf files belong to 
# different grops, reads that overlap features from different gtf files are assigned to the 
# feature belonging to  higher ranked gtf file. This is paricularly helpful if one of the
# gtf file has features  spanning most of the genome (like genomic repeats) that would make 
# most of the reads being  ambiguous.  

import sys, getopt, os, re
# import pysam
import numpy
import HTSeq

# read arguments in
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,'ho:',["exons=","introns=","gr3="])
except getopt.GetoptError:
    print 'multiGtfHTSeq.py -o output_suffix --exons=<exons_gtf_file> --introns=<introns_gtf_file> <inputSAMfile>'
    sys.exit(2)

if len(args) != 1:
    print 'multiGtfHTSeq: Error: One input file have to be provided.'
    sys.exit(2)
else:
    inputSAMfile = args[0]
    if inputSAMfile == "-":
        read_seq = iter( HTSeq.SAM_Reader( sys.stdin ) )
        inputSAMfile = 'genericSAM'
    else:
        if not os.path.exists(inputSAMfile):
            print 'multiGtfHTSeq: Error: File ' + inputSAMfile + ' not found.'
            sys.exit(2)
        read_seq = HTSeq.SAM_Reader( inputSAMfile )

gr1_gtf_files = []
gr2_gtf_files = []
gr3_gtf_files = []

outputDIR = 'mHTSeq'

for opt, arg in opts:
    if opt == '-h':
        print 'multiGtfHTSeq.py -o output_suffix --exons=<exons_gtf_file> --introns=<introns_gtf_file> <inputSAMfile>'
        sys.exit(2)
    if opt == '-o':
        outputDIR = arg
    if opt == '--exons':
        gr1_gtf_files = arg
    if opt == '--introns':
        gr2_gtf_files = arg
        introns_mode = True
    else:
        introns_mode = False

if len(gr1_gtf_files) == 0:
    print 'multiGtfHTSeq: Error: at least one gtf file must be provided.'
    sys.exit(2)

print 'Input file: '+inputSAMfile
print 'Exons gtf file: '+gr1_gtf_files
if introns_mode:
    print 'Introns gtf file: '+gr2_gtf_files

# group 1
gr1_gtf_names = 'exons'
# for gtf_file in gr1_gtf_files: 
#     gr1_gtf_names.append( re.sub('.*/','',gtf_file) )

gr1_counts = {}
gr1_gtf_file_reader = {}

# for i in range(len(gr1_gtf_files)):
#     gr1_counts[ gr1_gtf_names[i] ] = {}
#     gr1_gtf_file_reader[ gr1_gtf_names[i] ] = HTSeq.GFF_Reader(gr1_gtf_files[i])
gr1_counts[ gr1_gtf_names ] = {}
gr1_gtf_file_reader[ gr1_gtf_names ] = HTSeq.GFF_Reader(gr1_gtf_files)

gr1_annotation = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for feature in gr1_gtf_file_reader[ gr1_gtf_names ]:
    if feature.type == "exon":
        gr1_annotation[ feature.iv ] += (feature.name, gr1_gtf_names)
        gr1_counts[ gr1_gtf_names ][ feature.name ] = 0

# group 2
gr2_gtf_names = 'introns'
# for gtf_file in gr2_gtf_files: 
#     gr2_gtf_names.append( re.sub('\.gtf$','',re.sub('.*/','',gtf_file)) )

gr2_counts = {}
gr2_gtf_file_reader = {}

gr2_counts[ gr2_gtf_names ] = {}
gr2_gtf_file_reader[ gr2_gtf_names ] = HTSeq.GFF_Reader(gr2_gtf_files)

gr2_annotation = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for feature in gr2_gtf_file_reader[ gr2_gtf_names ]:
    if feature.type == "exon":
        gr2_annotation[ feature.iv ] += (feature.name, gr2_gtf_names)
        gr2_counts[ gr2_gtf_names ][ feature.name ] = 0

# # group 3
# gr3_gtf_names = []

# gr3_counts = {}
# gr3_gtf_file_reader = {}

# for i in range(len(gr3_gtf_files)):
#     gr3_counts[ gr3_gtf_names[i] ] = {}
#     gr3_gtf_file_reader[ gr3_gtf_names[i] ] = HTSeq.GFF_Reader(gr3_gtf_files[i])

# gr3_annotation = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
# for name in gr3_gtf_names:
#     for feature in gr3_gtf_file_reader[ name ]:
#         if feature.type == "exon":
#             gr3_annotation[ feature.iv ] += (feature.name, name)
#             gr3_counts[ name ][ feature.name ] = 0

general_stats = {'feature': 0, 'no_feature': 0, 'ambiguous': 0, 'too_low_aQual': 0, 'not_aligned': 0, 'alignment_not_unique': 0}
global_counts = {}
for name in [gr1_gtf_names, gr2_gtf_names]: #+ gr3_gtf_names:
    global_counts[ name ] = 0

#import itertools
#for alnmt in itertools.islice( sam_file, 199 ):
for alnmt in read_seq:
    if alnmt.aligned:
        if alnmt.optional_field( "NH" ) > 1:
            general_stats[ 'alignment_not_unique' ] += 1
        else:
            # try group 1
            intersection_set = set([])
            for iv2, step_set in gr1_annotation[ alnmt.iv ].steps():
                if step_set: # meaning that is not considered when null (intersection non-empty)
                    if not intersection_set:
                        intersection_set = step_set.copy()
                    else:
                        intersection_set.intersection_update( step_set )
            if len( intersection_set ) == 1:
                name = list(intersection_set)[0][1]
                gr1_counts[name][ list(intersection_set)[0][0] ] += 1
            else:
                if len( intersection_set ) == 0: 
                    if len(gr2_gtf_files):
                    # try group 2
                        intersection_set = set([])
                        for iv2, step_set in gr2_annotation[ alnmt.iv ].steps():
                            if step_set: # meaning that is not considered when null (intersection non-empty)
                                if not intersection_set:
                                    intersection_set = step_set.copy()
                                else:
                                    intersection_set.intersection_update( step_set )
                        if len( intersection_set ) == 1:
                            name = list(intersection_set)[0][1]
                            gr2_counts[name][ list(intersection_set)[0][0] ] += 1
                        else:
                            if len( intersection_set ) == 0: 
                                # if len(gr3_gtf_files):
                                #     # try group 3
                                #     intersection_set = set([])
                                #     for iv2, step_set in gr3_annotation[ alnmt.iv ].steps():
                                #         if step_set: # meaning that is not considered when null (intersection non-empty)
                                #             if not intersection_set:
                                #                 intersection_set = step_set.copy()
                                #             else:
                                #                 intersection_set.intersection_update( step_set )
                                #     if len( intersection_set ) == 1:
                                #         name = list(intersection_set)[0][1]
                                #         gr3_counts[name][ list(intersection_set)[0][0] ] += 1
                                #     else:
                                #         if len( intersection_set ) == 0:
                                #             general_stats[ 'no_feature' ] += 1
                                #         else:
                                #             general_stats[ 'ambiguous' ] += 1
                                # else:
                                general_stats[ 'no_feature' ] += 1
                            else:
                                general_stats[ 'ambiguous' ] += 1
                    else:
                        general_stats[ 'no_feature' ] += 1
                else:
                    general_stats[ 'ambiguous' ] += 1
    else:
        general_stats[ 'not_aligned' ] += 1

# for name in gr1_gtf_names:
#     global_counts[ name ] = sum( gr1_counts[ name ].values() )
global_counts[ 'exons' ] = sum( gr1_counts[ 'exons' ].values() )

# for name in gr2_gtf_names:
#     global_counts[ name ] = sum( gr2_counts[ name ].values() )
global_counts[ 'introns' ] = sum( gr2_counts[ 'introns' ].values() )

# for name in gr3_gtf_names:
#     global_counts[ name ] = sum( gr3_counts[ name ].values() )

general_stats['feature'] = sum(global_counts.values())

#global_counts
#general_stats

if not os.path.exists(outputDIR): os.makedirs(outputDIR)

# write group 1 on file
f = open( outputDIR+'/'+gr1_gtf_names+'_gr1.counts', 'w' )
for key in gr1_counts[ gr1_gtf_names ].keys():
    f.write('{0}\t{1}\n'.format(key, gr1_counts[ gr1_gtf_names ][ key ]))
f.close()

# write group 2 on file
if len(gr2_gtf_files)>0:
    f = open( outputDIR+'/'+gr2_gtf_names+'_gr2.counts', 'w' )
    for key in gr2_counts[ gr2_gtf_names ].keys():
        f.write('{0}\t{1}\n'.format(key, gr2_counts[ gr2_gtf_names ][ key ]))
    f.close()

# # write group 3 on file
# if len(gr3_gtf_names)>0:
#     for name in gr3_gtf_names:
#         f = open( outputDIR+'/'+name+'_gr3.counts', 'w' )
#         for key in gr3_counts[ name ].keys():
#             f.write('{0}\t{1}\n'.format(key, gr3_counts[ name ][ key ]))
#         f.close()

# write statistics on file
f = open( outputDIR+'/general_stats.counts', 'w' )
for key in global_counts.keys():
    f.write('{0}\t{1}\n'.format(key, global_counts[ key ]))
f.write('\n')
for key in general_stats.keys():
    f.write('{0}\t{1}\n'.format(key, general_stats[ key ]))
f.close()









