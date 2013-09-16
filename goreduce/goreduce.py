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
import tkSimpleDialog
import tkMessageBox
import json
from subprocess import Popen, PIPE

		

	
#This reads 2 structures, calculates RMSD for each atom in backbone and assigns that value to b-factor of the atom. It sets b-factors for all other atoms to 0
def jsoner(sel):
	lens=[]
	states=[]
	q1=cmd.get_model(sel)
	lens.append(len(q1.atom))
	states.append(1)
	proc = Popen("goreduce", shell=True, stdin=PIPE) #, stdout=PIPE)
	options=json.dumps({"SelNames":[sel],"AtomsPerSel":lens,"StatesPerSel":states})  #, "IntOptions":[[5, 11]] })
	proc.stdin.write(options+"\n")
	for i in q1.atom:
		atom,coords=Atom2gcRef(i)
		proc.stdin.write(atom+"\n")
		proc.stdin.write(coords+"\n")
	proc.stdin.close()
#	if  proc.wait() != 1:
#		print "There were some errors"
#	for j in proc.stderr:
#		print(json.loads(j))
#	proc.stderr.close()	
#	model,info=get_go_output(proc.stdout)
#	cmd.load_model(model,sel+"_H")





#parses the json output from go. Just copypaste it in your plugin.
def get_go_output(stdout):
	model=Indexed()
	atoms=0
	atomsread=0
	readingatoms=False
	readingcoords=False
	readingbfactors=False
	readingss=False
	first=True
	for v in stdout:
		if first:
			first=False
			print "v", v
			info=json.loads(v)
			atoms=info["AtomsPerMolecule"]
			if atoms!=len(q1):
				print "wrong number of atoms!", atoms, len(q1)
				break
			readingatoms=True
			continue
		if readingatoms:
			ad=json.loads(v)
			at=Atom()
			at.name=ad["Name"]
			at.symbol=ad["Symbol"]
			at.chain=ad["Chain"]
			at.id=ad["Id"]
			at.resi_number=ad["Molid"]
			at.resn=ad["Molname"]
			model.atom.append(at)
			++atomsread
			if atomsread==atoms:
				readingatoms=False
				readingcoords=True
				atomsread=0
				continue
			if readingcoords:
				coords=json.loads(v)
				model.atom[atomsread].coord=coords
				++atomsread
				if atomsread==atoms:
					readingcoords=False
					atomsread=0
					if  info["Bfactors"]:
						readingbfactors=True
						continue
					if info["SS"]:
						readingss=True
						continue
			#In many cases this part will not be needed
			if readingbfactors:
				bf=json.loads(v)
				model.atom[atomsread].b=bf
				++atomsread
				if atomsread==atoms:
					readingbfactors=False
					if info["SS"]:
						readingss=True
						continue
			#This one should be needed only seldomly
			if readingss:
				SS=json.loads(v)
				model.atom[atomsread].ss=SS
				++atomsread
				if atomsread==atoms:
					readingss=False #not really needed
			return model, info










			
	#some errors maybe just warnings
#	for i in proc.stderr:
#		try:
#			print json.loads(i)
#		except:
#			print i
#
		

def Atom2gcRef(i):
	Ref=json.dumps({"Name":i.name,"Id":i.id,"Molname":i.resn,"Symbol":i.symbol,"Molid":int(i.resi_number),"Chain":i.chain})
	coords=json.dumps({"Coords":i.coord})
	return Ref,coords
	


#The Tk interface (I know that  the interface is awful, but I have no time to study Tk at present)
#Please feel free to improve it if you can/want, and please e-mail me your version.
#maybe it would be better to make this app a pymol script, so we don't need this interface        
def goReduceDialog(app): 

	sel = tkSimpleDialog.askstring("goReduce",
                                       'Enter the selection to protonate',
                                       parent=app.root)
	jsoner(sel)
 
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
                             'goReduce',
                             label = 'goReduce',
	command = lambda s=self : goReduceDialog(s))


