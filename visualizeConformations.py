import lzma as l
import time as t
import argparse as a
import decimal as d
import numpy as n
import math as m
import os
from tqdm import tqdm
import sys as s 
import copy as c

# from contextlib import ExitStack 
# ExitStack is a good module to load when the two files have equal number of lines. If not, the iteration stops when the shortest input file ends.

# Input dihedral angles from dihedral.dump
# Find atoms of coils, extended coils and helices
# Change the atom type for visualization, in dump file
# Create a new dump file output

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield d.Decimal (item)
		except:
			pass

def computeConformations (dihedralList):
	trans1_low = -180
	trans1_high = -117
	trans2_low = 116
	trans2_high = 180
	gauchePlus_low = 34
	gauchePlus_high = 115
	gaucheMinus_low = -116
	gaucheMinus_high = -38

	for x in range (len (dihedralList)):
		try:
			if ((dihedralList[x][1] > trans1_low and dihedralList[x][1] < trans1_high) or (dihedralList[x][1] > trans2_low and dihedralList[x][1] < trans2_high)):
				dihedralList[x].append ("T")
			elif (dihedralList[x][1] > gaucheMinus_low and dihedralList[x][1] < gaucheMinus_high):
				dihedralList[x].append ("GP")
			elif (dihedralList[x][1] > gauchePlus_low and dihedralList[x][1] < gauchePlus_high):
				dihedralList[x].append ("GM")
			else:
				dihedralList[x].append ("NA")
		except:
			print ("Error occured at {} iteration while computing conformations".format (x))

	return dihedralList

def computeChirality (dihedralList):

	for x in range (len (dihedralList) - 1):
		try:
			if ((dihedralList[x][6] == "T" and dihedralList[x + 1][6] == "T") or (dihedralList[x][6] == "T" and dihedralList[x - 1][6] == "T")):
				dihedralList[x].append ("TT")
			else:
				dihedralList[x].append ("NTT")
		except:
			print ("Exception raised in computeChirality () [TT] for {} iteration".format (x))

	for x in range (len (dihedralList) - 1):
		try:
			if ((dihedralList[x][6] == "T" and dihedralList[x + 1][6] == "GP") or (dihedralList[x][6] == "GP" and dihedralList[x + 1][6] == "T") or (dihedralList[x][6] == "T" and dihedralList[x - 1][6] == "GP") or (dihedralList[x][6] == "GP" and dihedralList[x - 1][6] == "T")):
				dihedralList[x].append ("TGP")
			else:
				dihedralList[x].append ("NTGP")
		except:
			print ("Exception raised in computeChirality () [TGP] for {} iteration".format (x))

	for x in range (len (dihedralList) - 1):
		try:
			if ((dihedralList[x][6] == "T" and dihedralList[x + 1][6] == "GM") or (dihedralList[x][6] == "GM" and dihedralList[x + 1][6] == "T") or (dihedralList[x][6] == "T" and dihedralList[x - 1][6] == "GM") or (dihedralList[x][6] == "GM" and dihedralList[x - 1][6] == "T")):
				dihedralList[x].append ("TGM")
			else:
				dihedralList[x].append ("NTGM")
		except:
			print ("Exception raised in computeChirality () [TGM] for {} iteration".format (x))

	for x in range (len (dihedralList) - 1):
		try:
			if ((("G" in dihedralList[x][6]) and ("G" in dihedralList[x + 1][6])) or (("G" in dihedralList[x][6]) and ("G" in dihedralList[x - 1][6]))):
				dihedralList[x].append ("GG")
			else:
				dihedralList[x].append ("NGG")
		except:
			print ("Exception raised in computeChirality () [GG] for {} iteration".format (x))

	return dihedralList

def removePendantDihedrals (dihedralList):
	backbone = []
	newDihedralList = []
	for line in dihedralList:
		if (line[3] not in backbone):
			backbone.append (line[3])
		if (line[4] not in backbone):
			backbone.append (line[4])

	for x in range ((len (dihedralList))):
		if ((dihedralList[x][2] not in backbone) or (dihedralList[x][5] not in backbone)):
			pass
		else:
			newDihedralList.append (dihedralList[x])

	return newDihedralList

def convertListToString (inputList):
	seperator = "   "
	outputString = ""

	for x in inputList:
		outputString = ("{}{}{}".format (outputString, seperator, x))
	
	outputString = ("{}\n".format (outputString))
	return outputString

def main (inputDump, inputDihedral, maxTimeframes, outputCoil, outputExCoil, outputRight, outputLeft, outputHelix):
	nLines = 0
	nTimeframe = 0
	dihedralList = []
	nDihedrals = 0
	dumpList = []
	nAtoms = 0
	dihLine = "notEmpty"
	dumpLine = "notEmpty"
	nBackboneAtoms = 0

	firstFrame = 1
	dihHeader = []
	dumpHeader = []

	dihedralFile = open (inputDihedral, "r")
	dumpFile = open (inputDump, "r")

	# Output section
	coilList = []
	exCoilList = []
	rightList = []
	leftList = []
	helixList = []

	coilFile = open (outputCoil, "w")
	exCoilFile = open (outputExCoil, "w")
	rightFile = open (outputRight, "w")
	leftFile = open (outputLeft, "w")
	helixFile = open (outputHelix, "w")
	statsFile = open ("stats_visual.csv", "w")

	coilFileLine = 0
	exCoilFileLine = 0
	rightFileLine = 0
	leftFileLine = 0
	helixFileLine = 0

	statsFile.write ("coil,exCoil,right,left,helix\n")

	while ((len (dihLine) + len (dumpLine)) > 0):

		nTimeframe += 1

		if (nTimeframe == maxTimeframes):
			print ("\n\nReached final timeframe (specified by the user)...")
			break

		print ("Scanning timeframe: {}...".format (nTimeframe), end = "\r")
		for x in range (4):
			dihLine = dihedralFile.readline ()
			dumpLine = dumpFile.readline ()

			if (firstFrame == 1):
				dihHeader.append (dihLine)
				dumpHeader.append (dumpLine)

			nDihedrals = list (extract_numbers (dihLine))
			nAtoms = list (extract_numbers (dumpLine))

		for x in range (5):
			dihLine = dihedralFile.readline ()
			dumpLine = dumpFile.readline ()

			if (firstFrame == 1):
				dihHeader.append (dihLine)
				dumpHeader.append (dumpLine)

		try:
			for x in range (int (nDihedrals[0])):
				dihLine = dihedralFile.readline ()
				dihedralList.append (list (extract_numbers (dihLine)))

			for x in range (int (nAtoms[0])):
				dumpLine = dumpFile.readline ()
				dumpList.append (list (extract_numbers (dumpLine)))
		except:
			s.exit ("\n\nEnd of file detected... Program will stop now...")

		dihedralList = removePendantDihedrals (dihedralList)
		dihedralList = computeConformations (dihedralList)
		dihedralList = computeChirality (dihedralList)

		for dihedral in dihedralList:
			try:
				if (dihedral[7] == "TT"):
					if (int (dihedral[3]) not in exCoilList):
						exCoilList.append (int (dihedral[3]))
					if (int (dihedral[4]) not in exCoilList):
						exCoilList.append (int (dihedral[4]))

				if (dihedral[10] == "GG"):
					if (int (dihedral[3]) not in coilList):
						coilList.append (int (dihedral[3]))
					if (int (dihedral[4]) not in coilList):
						coilList.append (int (dihedral[4]))

				if (dihedral[8] == "TGP"):
					if (int (dihedral[3]) not in rightList):
						rightList.append (int (dihedral[3]))
					if (int (dihedral[4]) not in rightList):
						rightList.append (int (dihedral[4]))

				if (dihedral[9] == "TGM"):
					if (int (dihedral[3]) not in leftList):
						leftList.append (int (dihedral[3]))
					if (int (dihedral[4]) not in leftList):
						leftList.append (int (dihedral[4]))
			except:
				pass

		helixList.append (rightList)
		helixList.append (leftList)

		if (firstFrame == 1):
			for atom in dumpList:
				if (atom[2] != 3):
					nBackboneAtoms += 1

		nBackboneAtomsString = ("{}\n".format (nBackboneAtoms))

		coilFile.writelines (dumpHeader[:3])
		exCoilFile.writelines (dumpHeader[:3])
		rightFile.writelines (dumpHeader[:3])
		leftFile.writelines (dumpHeader[:3])
		helixFile.writelines (dumpHeader[:3])

		coilFile.write (str (nBackboneAtomsString))
		exCoilFile.write (str (nBackboneAtomsString))
		rightFile.write (str (nBackboneAtomsString))
		leftFile.write (str (nBackboneAtomsString))
		helixFile.write (str (nBackboneAtomsString))

		coilFile.writelines (dumpHeader[4:])
		exCoilFile.writelines (dumpHeader[4:])
		rightFile.writelines (dumpHeader[4:])
		leftFile.writelines (dumpHeader[4:])
		helixFile.writelines (dumpHeader[4:])

		dumpListTemp = []
		dumpListTemp = c.deepcopy (dumpList)
		for atom in dumpListTemp:
			if (atom[2] != 3):
				if (int (atom[0]) in exCoilList):
					exCoilFileLine += 1
					atom[0] = exCoilFileLine
					atom[1] = 2
				else:
					exCoilFileLine += 1
					atom[0] = exCoilFileLine
					atom[1] = 1

				outputLine = convertListToString (atom)
				exCoilFile.write (outputLine)

		dumpListTemp = []
		dumpListTemp = c.deepcopy (dumpList)
		for atom in dumpListTemp:
			if (atom[2] != 3):
				if (int (atom[0]) in coilList):
					coilFileLine += 1
					atom[0] = coilFileLine
					atom[1] = 2
				else:
					coilFileLine += 1
					atom[0] = coilFileLine
					atom[1] = 1

				outputLine = convertListToString (atom)
				coilFile.write (outputLine)

		dumpListTemp = []
		dumpListTemp = c.deepcopy (dumpList)
		for atom in dumpListTemp:
			if (atom[2] != 3):
				if (int (atom[0]) in rightList):
					rightFileLine += 1
					atom[0] = rightFileLine
					atom[1] = 2
				else:
					rightFileLine += 1
					atom[0] = rightFileLine
					atom[1] = 1

				outputLine = convertListToString (atom)
				rightFile.write (outputLine)

		dumpListTemp = []
		dumpListTemp = c.deepcopy (dumpList)
		for atom in dumpListTemp:
			if (atom[2] != 3):
				if (int (atom[0]) in leftList):
					leftFileLine += 1
					atom[0] = leftFileLine
					atom[1] = 2
				else:
					leftFileLine += 1
					atom[0] = leftFileLine
					atom[1] = 1

				outputLine = convertListToString (atom)
				leftFile.write (outputLine)

		dumpListTemp = []
		dumpListTemp = c.deepcopy (dumpList)
		for atom in dumpListTemp:
			if (atom[2] != 3):
				if (int (atom[0]) in helixList):
					helixFileLine += 1
					atom[0] = helixFileLine
					atom[1] = 2
				else:
					helixFileLine += 1
					atom[0] = helixFileLine
					atom[1] = 1

				outputLine = convertListToString (atom)
				helixFile.write (outputLine)

		statsString = ("{},{},{},{},{}\n".format(len (coilList), len (exCoilList), len (rightList), len (leftList), len (helixList)))
		statsFile.write (statsString)

			# if (atom[2] != 3):
			# 	if (int (atom[0]) in exCoilList):
			# 		coilFileLine += 1
			# 		atom[0] = coilFileLine
			# 		atom[1] = 2
			# 	else:
			# 		coilFileLine += 1
			# 		atom[0] = coilFileLine

			# 	print ("{} {}".format (atom[0], atom[1]))
			# 	t.sleep(1)

		dihedralList = []
		dumpList = []

		coilFileLine = 0
		exCoilFileLine = 0
		rightFileLine = 0
		leftFileLine = 0
		helixFileLine = 0

		coilList = []
		exCoilList = []
		rightList = []
		leftList = []
		helixList = []

		firstFrame = 0


	dihedralFile.close ()
	dumpFile.close ()
	coilFile.close ()
	exCoilFile.close ()
	rightFile.close ()
	leftFile.close ()
	helixFile.close ()

			
if __name__ == '__main__':
	parser = a.ArgumentParser (description = "Compute global order parameter from input LAMMPS dump file in XZ format.")
	parser.add_argument ("--dump", "-trj", required = True, type = str, help = "Input LAMMPS dump file name in compressed XZ format [required arg]")
	parser.add_argument ("--dihedral", "-dih", required = True, type = str, help = "Input LAMMPS dihedral file name in compressed XZ format [required arg]")
	parser.add_argument ("--timeframes", "-tf", default = 0, type = int, help = "Number of timeframes to scan")
	parser.add_argument ("--outputCoil", "-oc", type = str, default = "coil.lammpstrj", help = "Output file containing marked coils")
	parser.add_argument ("--outputExCoil", "-oxc", type = str, default = "exCoil.lammpstrj", help = "Output file containing marked extended coils")
	parser.add_argument ("--outputRight", "-or", type = str, default = "right.lammpstrj", help = "Output file containing marked right handed helices")
	parser.add_argument ("--outputLeft", "-ol", type = str, default = "left.lammpstrj", help = "Output file containing marked left handed helices")
	parser.add_argument ("--outputHelix", "-oh", type = str, default = "helix.lammpstrj", help = "Output file containing marked all helices (both left and right handed helices)")
	args = parser.parse_args()

	main(args.dump, args.dihedral, args.timeframes, args.outputCoil, args.outputExCoil, args.outputRight, args.outputLeft, args.outputHelix)