#!/usr/bin/python

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2013 by Raul Mera-Adasme
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ------------------------------


##This work is dedicated to the long life of the Ven. Khempo Phuntzok Tenzin Rinpoche.##


from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
import json



#returns dictionary with information returned by the go program
def get_info(proc):
	first=False
	info=json.loads(proc.stdout.readline())
	return info


#puts coordinates returned by the go program into the model model. Coordinates of atom with names in the 
#list name are the only one replaced, or the ones not replaced depending on whether included is True.
#proc is the process object for the go program, atom is the lenght of the model.
def get_coords(proc,model,names,included,info,number):
	atomsread=0
	atom=info["AtomsPerMolecule"][number]
	rcoords=True
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline()
		if "Coords" in v and not "Molname" in v and rcoords:
#			print "yuy"
			coords=json.loads(v)
			if included:
				if model.atom[atomsread].name in names:
					model.atom[atomsread].coord=coords["Coords"]
			else:
				if not model.atom[atomsread].name in names:
					model.atom[atomsread].coord=coords["Coords"]
			atomsread=atomsread+1
			if atomsread==atoms:
				rcoords=False
				atomsread=0
				if  info["Bfactors"]:
					rbfactors=True
				if info["SS"]:
					rss=True
			continue
		#In many cases this part will not be needed
		if "Bfactors" in v and rbfactors:
#			print "YAYAYAYAYA"
			bf=json.loads(v)
			model.atom[atomsread].b=bf["Bfactors"]
			atomsread=atomsread+1
			if atomsread==atoms-1:
				atomsread=0
				rbfactors=False
				if info["SS"]:
					rss=True
			continue
		#This one should be needed only seldomly
		if "SS" in v and rss:
#			print "yey"
			SS=json.loads(v)
			model.atom[atomsread].ss=SS["SS"]
			++atomsread
			if atomsread==atoms-1:
				atomsread=0
				rss=False
			continue
		break
#		print "me fui con una deuda de 500"
#	print "ATOMS!", len(model.atom),atoms, atomsread, model.atom[-1].coord, model.atom[-2].coord, model.atom[-1] 
	return model



#similar to get_gochem but it doesnt take a reference model, hence, all the data is taken from the output of the go program.
#Proc is the gochem process object, atoms is the number of atoms that should be added.
#Notice that PyMOL 1.6 will NOT guess bonds for objects created with these functions, so lines, sticks, cartoons representations
#willl be NOT available. 
def get_model(proc, info,number):
	vmodel=Indexed()
	atoms=info["AtomsPerMolecule"][number]
	atomsread=0
	ratoms=True
	rcoords=False
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline()
		if "Molname" in v and ratom:
#			print "YAY", atoms, atomsread, v
			ad=json.loads(v)
#			print "v atom", json.loads(v)
			at=Atom()
			at.name=ad["Name"]
			at.symbol=ad["Symbol"]
			at.chain=ad["Chain"]
			at.id=ad["Id"]
			at.resi_number=ad["Molid"]
			at.resn=ad["Molname"]
			vmodel.atom.append(at)
			atomsread=atomsread+1
			if atomsread==atoms-1:
				ratoms=False
				rcoords=True
				atomsread=0
			continue
		if "Coords" in v and not "Molname" in v and rcoords:
#			print "yuy"
			coords=json.loads(v)
#			print "coords!!", coords["Coords"], atomsread
			vmodel.atom[atomsread].coord=coords["Coords"]
#			print "coords in mol!", vmodel.atom[atomsread]
			atomsread=atomsread+1
#			print "reeeaddd", atomsread
			if atomsread==atoms:
				rcoords=False
				atomsread=0
				if  info["Bfactors"]:
					rbfactors=True
				if info["SS"]:
					rss=True
			continue
		#In many cases this part will not be needed
		if "Bfactors" in v:
#			print "YAYAYAYAYA"
			bf=json.loads(v)
			vmodel.atom[atomsread].b=bf["Bfactors"]
			atomsread=atomsread+1
			if atomsread==atoms-1:
				atomsread=0
				rbfactors=False
				if info["SS"]:
					rss=True
			continue
		#This one should be needed only seldomly
		if "SS" in v:
#			print "yey"
			SS=json.loads(v)
			vmodel.atom[atomsread].ss=SS["SS"]
			++atomsread
			if atomsread==atoms-1:
				atomsread=0
				rss=False
			continue
#		print "me fui con una deuda de 500"
		break
#	print "ATOMS!", len(vmodel.atom),atoms, atomsread, vmodel.atom[-1].coord, vmodel.atom[-2].coord, vmodel.atom[-1] 
	return vmodel





def Atom2gcRef(i):
	Ref=json.dumps({"Name":i.name,"Id":i.id,"Molname":i.resn,"Symbol":i.symbol,"Molid":int(i.resi_number),"Chain":i.chain})
	coords=json.dumps({"Coords":i.coord})
	return Ref,coords
	

