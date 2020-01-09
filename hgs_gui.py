#!/usr/bin/python3

"""
Housekeeping gene selection

Copyright 2014 Olivier Friard

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.

"""

__version__ = 0.4



import PyQt5
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *


import xlrd
import math
import sys

def average(s):
    return sum(s) * 1.0 / len(s)

def stdev(s):
    avg = average(s)
    variance = list(map(lambda x: (x - avg)**2, s))
    return math.sqrt(average(variance))


class showResults_class(QWidget):
    '''
    class for displaying results
    '''

    def __init__(self):
        super(showResults_class, self).__init__()

        self.label = QLabel()
        self.label.setText('')
        self.te = QTextEdit()

        hbox = QVBoxLayout(self)

        hbox.addWidget(self.label)
        hbox.addWidget(self.te)

        hbox2 = QHBoxLayout(self)


        self.pbSave = QPushButton('Save results')
        hbox2.addWidget(self.pbSave)

        spacerItem = QSpacerItem(241, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        hbox2.addItem(spacerItem)


        self.pbClose = QPushButton('Close')
        hbox2.addWidget(self.pbClose)

        hbox.addLayout(hbox2)

        self.setWindowTitle('Results')


        self.pbClose.clicked.connect(self.pbClose_clicked)
        self.pbSave.clicked.connect(self.pbSave_clicked)



    def pbClose_clicked(self):
        self.close()

    def pbSave_clicked(self):
        '''
        save results in TXT format
        '''
        fd = QFileDialog(self)
        filename = fd.getSaveFileName(self,'Save results','','All files (*)')
        if filename:
            with open(filename,'w') as outfile:
                outfile.writelines(str(self.te.toPlainText()))


class main(QWidget):

    def __init__(self, parent=None):
        super(main, self).__init__(parent)

        self.setWindowTitle('Housekeeping gene selection -  v.%.1f' % __version__)

        self.genesList = []
        self.expr = {}

        self.lineEdit = QLineEdit(self)
        self.view = QTableView(self)
        self.view.setEnabled(False)
        self.label = QLabel(self)

        self.lbInfo = QLabel(self)
        self.lbInfo.setText('Olivier Friard - Life Sciences and System Biology Department - University of Torino - 2014')
        self.lbInfo.setEnabled(False)
        self.load     = QPushButton(self)
        self.load.setText('Load Ct values')

        self.submit     = QPushButton(self)
        self.submit.setText('Select genes')

        self.gridLayout = QGridLayout(self)

        self.gridLayout.addWidget(self.load, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.label, 0, 1, 1, 1)
        self.gridLayout.addWidget(self.lineEdit, 0, 2, 1, 1)
        self.gridLayout.addWidget(self.submit, 0, 3, 1, 1)
        self.gridLayout.addWidget(self.view, 1, 0, 1, 4)
        self.gridLayout.addWidget(self.lbInfo, 2, 0, 1, 4)

        self.label.setText("Number of genes to select")

        self.model = QStandardItemModel(self)

        self.proxy = QSortFilterProxyModel(self)
        self.proxy.setSourceModel(self.model)

        self.view.setModel(self.proxy)

        self.horizontalHeader = self.view.horizontalHeader()

        self.load.clicked.connect( self.load_clicked )
        self.submit.clicked.connect( self.submit_clicked )


    def load_clicked(self):

        fd = QFileDialog(self)
        fileName, _ = fd.getOpenFileName(self, 'Load a file', '', 'Excel files (*.xls);;All files (*)')

        if fileName:

            workbook = xlrd.open_workbook(fileName)

            worksheet = workbook.sheet_by_index(0)
            num_rows = worksheet.nrows - 1
            num_cells = worksheet.ncols - 1
            curr_row = -1

            self.genesList = []
            self.expr = {}

            while curr_row < num_rows:
                curr_row += 1
                row = worksheet.row(curr_row)

                curr_cell = -1
                while curr_cell < num_cells:
                    curr_cell += 1

                    cell_type = worksheet.cell_type(curr_row, curr_cell)
                    cell_value = worksheet.cell_value(curr_row, curr_cell)

                    if curr_row == 0 and curr_cell > 0:
                        self.genesList.append( cell_value )
                        self.expr[ cell_value ] = []

                    if curr_row > 0 and curr_cell > 0:
                        self.expr[ self.genesList[ curr_cell - 1 ] ].append( cell_value )

            maxCt = {}
            for gene in self.genesList:
                maxCt[gene] = max(self.expr[gene])

            print(self.expr)

            for gene in self.genesList:
                self.expr[gene] = [ 2**+( maxCt[gene] - x) for x in self.expr[gene] ]


            self.model.setHorizontalHeaderLabels(self.genesList)
            self.model.setRowCount(0)
            for row in range(len(self.expr[ list(self.expr.keys())[0]])):
                cols = []
                for gene in self.genesList:
                    cols.append( QStandardItem( str(self.expr[gene][row])) )
                self.model.appendRow( cols )





    def submit_clicked(self):
        '''
        apply algorithm
        '''

        try:
            nGenes = int(self.lineEdit.text())
        except:
            QMessageBox.critical(self, 'Housekeeping gene selection', 'Check the number of final genes!')
            return

        out = ''
        while len(self.genesList) > int(self.lineEdit.text()):

            v = {}

            for j in self.genesList:
                v[j] = []

                stdv = []
                for k in self.genesList:
                    if k == j:
                        continue

                    logRatio = []
                    for r in range(len(   self.expr[ list(self.expr.keys())[0]])):
                        logRatio.append( math.log(  self.expr[ j ][r]  /  self.expr[ k ][r], 2) )

                    stdv.append( stdev( logRatio ))
                v[j] = average( stdv )

            for gene in sorted(v, key=v.get, reverse=False):
                out += '%s:\t%f\n' % ( gene, v[gene] )
                print(gene, v[gene])



            maxGene =  [ x for x in v if v[x] == max(v.values())][0]
            print( 'Remove: %s gene' % maxGene)
            out += 'Remove: %s gene\n' % maxGene

            out += '='*20 + '\n\n'

            self.genesList.remove( maxGene )

        print( self.genesList )
        out += 'The most stable genes are:\n%s\n'  % ('\n'.join(self.genesList))

        self.showResults = showResults_class()
        self.showResults.setWindowTitle('Results')
        self.showResults.te.setText(out)
        self.showResults.show()



if __name__ == "__main__":

    app  = QApplication(sys.argv)
    main = main()
    main.show()
    main.resize(800, 600)
    sys.exit(app.exec_())
