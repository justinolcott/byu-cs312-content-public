#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	def computeEditDistance(self, seq1, seq2):
		rows = len(seq1) + 1
		cols = len(seq2) + 1
		matrix = {}
		prev = {}

		# set vals in the first row
		for i in range(rows):
			matrix[i,0] = i*INDEL
			prev[i,0] = (i-1,0)
		# set vals in the first column
		for j in range(cols):
			matrix[0,j] = j*INDEL
			prev[0,j] = (0,j-1)
		
		# fill in the rest of the matrix
		# aka for each row
		for i in range(1,rows):
			# and for each column
			for j in range(1,cols):
				# Fill in the value

				# Find the diagonal distance
				# if equal, take match
				# if not equal, take substitution
				diagonal = float('inf')
				if seq1[i-1] == seq2[j-1]:
					diagonal = matrix[i-1,j-1]+MATCH
				else:
					diagonal = matrix[i-1,j-1]+SUB

				# find the up and left distance
				up = matrix[i-1,j]+INDEL
				left = matrix[i,j-1]+INDEL

				# set it to the minimum of the three
				matrix[i,j] = min(diagonal,up,left)

				# set the previous node in the prev dictionary
				# with the left then up then diagonal if they're all the same
				if matrix[i,j] == left:
					prev[i,j] = (i,j-1)
				elif matrix[i,j] == up:
					prev[i,j] = (i-1,j)
				else:
					prev[i,j] = (i-1,j-1)		
		return matrix[rows-1,cols-1], prev
	
	def computeEditDistanceBanded(self, seq1, seq2, n=MAXINDELS):
		rows = len(seq1) + 1
		cols = len(seq2) + 1
		matrix = {}
		prev = {}

		# set vals in the first row and column
		# it won't be the full values, just the ones that are in the band
		for i in range(n+1):
			matrix[i,0] = i*INDEL
			prev[i,0] = (i-1,0)
		for j in range(n+1):
			matrix[0,j] = j*INDEL
			prev[0,j] = (0,j-1)
		
		# for all the rows
		for i in range(1,rows):
			# and for the columns within the band
			for j in range(i-n,i+n+1):
				# if it's out of bounds, skip it
				if j < 1 or j > cols-1:
					continue
				
				# diagonal will always exist
				d = matrix[i-1,j-1]

				# left and right will only exist if they're in the band, so we init to inf
				u = float('inf')
				l = float('inf')

				# if they're in the matrix, set them to the value
				if (i-1,j) in matrix:
					u = matrix[i-1,j]
				if (i,j-1) in matrix:
					l = matrix[i,j-1]
				
				# Find the diagonal distance
				diagonal = d+SUB
				# If they're equal, take the match
				if seq1[i-1] == seq2[j-1]:
					diagonal = d+MATCH
				
				# find the up and left distance
				up = u+INDEL
				left = l+INDEL

				# take the minimum of the directions
				matrix[i,j] = min(diagonal,up,left)

				# set the previous node in the prev dictionary
				# with the left then up then diagonal if they're all the same
				if matrix[i,j] == left:
					prev[i,j] = (i,j-1)
				elif matrix[i,j] == up:
					prev[i,j] = (i-1,j)
				else:
					prev[i,j] = (i-1,j-1)
		
		# if we can reach the end, return the value
		if (rows-1, cols-1) in matrix:
			return matrix[rows-1,cols-1], prev
		# else return infinity
		else:
			return float('inf'), None
		
	# Function for finding the string alignment
	def findAlign(self, prev, seq1, seq2):
		# if there's no alignment, return None in case of not making it to the end
		if prev == None:
			return "No Alignment Possible.","No Alignment Possible."
		
		i = len(seq1)
		j = len(seq2)
		align1 = ""
		align2 = ""

		# while we're not at the beginning empty cell
		while i > 0 or j > 0:
			if j == prev[i,j][1]: # if it came from the top aka same column
				align1 = seq1[i-1] + align1
				align2 = "-" + align2
				i, j = prev[i,j]
			elif i == prev[i,j][0]: # if it came from the left aka same row
				align1 = "-" + align1
				align2 = seq2[j-1] + align2
				i, j = prev[i,j]
			else: # if diagonal, not same col or row, and we can take both
				align1 = seq1[i-1] + align1
				align2 = seq2[j-1] + align2
				i, j = prev[i,j]
		return align1, align2

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		# Set the max length
		seq1 = seq1[:align_length]
		seq2 = seq2[:align_length]

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = 0
		alignment1 = ''
		alignment2 = ''
		if banded:
			score, prev = self.computeEditDistanceBanded(seq1,seq2)
			alignment1, alignment2 = self.findAlign(prev,seq1,seq2)
		else:
			score, prev = self.computeEditDistance(seq1,seq2)
			alignment1, alignment2 = self.findAlign(prev,seq1,seq2)

###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1[:100], 'seqj_first100':alignment2[:100]}