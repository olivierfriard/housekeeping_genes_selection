#!/usr/bin/python3

nGenes = 2

import xlrd
import math
import statistics
import sys

workbook = xlrd.open_workbook(sys.argv[1])

worksheet = workbook.sheet_by_index(0)
num_rows = worksheet.nrows - 1
num_cells = worksheet.ncols - 1
curr_row = -1

genesList = []
expr = {}

while curr_row < num_rows:
    curr_row += 1
    row = worksheet.row(curr_row)

    curr_cell = -1
    while curr_cell < num_cells:
        curr_cell += 1

        cell_type = worksheet.cell_type(curr_row, curr_cell)
        cell_value = worksheet.cell_value(curr_row, curr_cell)
        
        if curr_row == 0 and curr_cell > 0:
            genesList.append( cell_value )
            expr[ cell_value ] = []
        
        if curr_row > 0 and curr_cell > 0:
            expr[ genesList[ curr_cell - 1 ] ].append( cell_value )  

print( 'Genes list:', genesList )
#print( expr )

while len(genesList) > nGenes:

    v = {}

    for j in genesList:
        v[j] = []

        stdev = []
        for k in genesList:
            if k == j:
                continue

            logRatio = []
            for r in range( curr_row ):
                logRatio.append( math.log2(  expr[ j ][r]  /  expr[ k ][r]) )
    
            stdev.append( statistics.stdev( logRatio ))
        v[j] = statistics.mean( stdev )
    
    for gene in sorted(v, key=v.get, reverse=False):
        print(gene, v[gene])


    
    maxGene =  [ x for x in v if v[x] == max(v.values())][0]
    print( 'maxGene', maxGene)
    genesList.remove( maxGene )
    #input('key')
    
print( genesList )



