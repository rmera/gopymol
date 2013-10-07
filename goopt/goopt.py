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
from subprocess import Popen, PIPE
import array
import gochem
import json
		


#Reads a selection, cuts the Ca--CO and N--Ca bonds, replaces CO and N with nitrogens, 
def jsoner(sel):
	lens=[]
	states=[]
	q1=cmd.get_model(sel)
	lens.append(len(q1.atom))
	states.append(1)
	proc = Popen("goopt", shell=True, stdin=PIPE,stdout=PIPE)
	charge=-1       #in the future these 2 can be read from the plugin meny
	multiplicity=1
	options=json.dumps({"SelNames":[sel],"AtomsPerSel":lens,"StatesPerSel":states,"IntOptions":[[charge,multiplicity]]}) 
	proc.stdin.write(options+"\n")
	for i in q1.atom:
		atom,coords=gochem.Atom2gcRef(i)
		proc.stdin.write(atom+"\n")
		proc.stdin.write(coords+"\n")
	proc.stdin.close()
	if  proc.wait() != 1:
		print "There were some errors"
	mod, info=gochem.get_gochem(proc,q1,["CA","HA","HA2","HA3", "O","N","H","C"],False)
	cmd.load_model(mod,sel+"_H",discrete=1,zoom=1)





#The Tk interface (I know that  the interface is awful, but I have no time to study Tk at present)
#Please feel free to improve it if you can/want, and please e-mail me your version.
#maybe it would be better to make this app a pymol script, so we don't need this interface        
def goOptDialog(app): 

	sel = tkSimpleDialog.askstring("goOpt",
                                       'Enter the selection to optimize',
                                       parent=app.root)
	jsoner(sel)
 
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
                             'goOpt',
                             label = 'goOpt',
	command = lambda s=self : goOptDialog(s))


