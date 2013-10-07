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

#Replaces the coordinates, secondary structure and b-factor of some selected atoms in model for the ones obtained from a Gochem program whose process object is proc. 
#Proc is the Gochem process object,
#model is the original model on which the gochem program operated. Names are the names of atoms whose coordinates will or will not be changed for the
#Gochem program's output. include determines whether the names in names are or not replaced in the original object. In addition to the modified object, it returns a
#dictionary with information returned by the go program
def get_gochem(proc,model,names,included):
	atoms=len(model.atom)
	atomsread=0
	ratoms=False
	rcoords=False
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline()
#		print "VVV", v
#		print "LULA LULA LULA", first, ratoms, rcoords
		if first:
			print "first!"
			first=False
			info=json.loads(v)
			atoms=info["AtomsPerMolecule"][0]
			ratoms=True
#			print "LALA"
			continue
		if "Coords" in v and not "Molname" in v:
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
		if "Bfactors" in v:
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
		if "SS" in v:
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
	return model, info



#similar to get_gochem but it doesnt take a reference model, hence, all the data is taken from the output of the go program.
#Proc is the gochem process object.
#Notice that at this point, PyMOL will NOT guess bonds for objects created with this functions, so lines, sticks, cartoons representations
#willl be NOT available. 
def get_gochem_newmodel(proc):
	vmodel=Indexed()
	atoms=0
	atomsread=0
	ratoms=False
	rcoords=False
	rbfactors=False
	rss=False
	first=True
	while(True):
		v=proc.stdout.readline()
		if first:
			first=False
			info=json.loads(v)
			atoms=info["AtomsPerMolecule"][0]
			ratoms=True
			continue
		if "Molname" in v:
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
		if "Coords" in v and not "Molname" in v:
			print "veeeeee", v
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
	return vmodel, info





def Atom2gcRef(i):
	Ref=json.dumps({"Name":i.name,"Id":i.id,"Molname":i.resn,"Symbol":i.symbol,"Molid":int(i.resi_number),"Chain":i.chain})
	coords=json.dumps({"Coords":i.coord})
	return Ref,coords
	

