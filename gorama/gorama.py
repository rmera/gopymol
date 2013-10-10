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

		

	
#this reads 2 structures, calculates RMSD for each atom in backbone and assigns that value to b-factor of the atom. It sets b-factors for all other atoms to 0
def jsoner(sel1):
	q1=[]
	lens=[]
	states=[]
	for j in sel1:
		q1.append(cmd.get_model(j))
		print "current sel: ", j, len(q1[-1].atom)
		lens.append(len(q1[-1].atom))
		states.append(1)
	proc = Popen("gorama", shell=True, stdin=PIPE)  #, stderr=PIPE)
	options=json.dumps({"SelNames":sel1,"AtomsPerSel":lens,"StatesPerSel":states,"StringOptions":[["GLY"]]})  #, "IntOptions":[[5, 11]] })
	proc.stdin.write(options+"\n")
	for k in q1:
		for i in k.atom:
			atom,coords=Atom2gcRef(i)
			proc.stdin.write(atom+"\n")
			proc.stdin.write(coords+"\n")
	proc.stdin.close()
	if  proc.wait() != 0:
		print "There were some errors"
	#some errors maybe just warnings
#	for i in proc.stderr:
#		try:
#			print json.loads(i)
#		except:
#			print i
#
		

def Atom2gcRef(i):
	Ref=json.dumps({"Name":i.name,"Id":i.id,"Molname":i.resn,"Symbol":i.symbol,"Molid":int(i.resi_number),"SS":i.ss,"Chain":i.chain})
	coords=json.dumps({"Coords":i.coord})
	return Ref,coords
	


#The Tk interface (I know that  the interface is awful, but I have no time to study Tk at present)
#Please feel free to improve it if you can/want, and please e-mail me your version.
#maybe it would be better to make this app a pymol script, so we don't need this interface        
def goRamaDialog(app): 

	sel = tkSimpleDialog.askstring("goRama",
                                       'Enter selection names separated by commas',
                                       parent=app.root)
	selsp=sel.split(",")
	jsoner(selsp)
 
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
                             'goRama',
                             label = 'goRama',
	command = lambda s=self : goRamaDialog(s))


