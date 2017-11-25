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

from __future__ import print_function
try:
    # for Python2
    import Tkinter.simpledialog   ## notice capitalized T in Tkinter 
    import Tkinter.messagebox

except ImportError:
    # for Python3
    import tkinter.simpledialog
    import tkinter.messagebox



from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
import json
from subprocess import Popen, PIPE
import array
import gochem

		

	
def jsoner(sel):
	lens=[]
	states=[]
	q1=cmd.get_model(sel)
	lens.append(len(q1.atom))
	states.append(1)
	proc = Popen("goreduce", shell=True, stdin=PIPE, stdout=PIPE)
	options=json.dumps({"SelNames":[sel],"AtomsPerSel":lens,"StatesPerSel":states})  #, "IntOptions":[[5, 11]] })
	proc.stdin.write((options+"\n").encode(encoding='UTF-8'))
	for i in q1.atom:
		atom,coords=gochem.Atom2gcRef(i)
		proc.stdin.write((atom+"\n").encode(encoding='UTF-8'))
		proc.stdin.write((coords+"\n").encode(encoding='UTF-8'))
	proc.stdin.close()
#	if  proc.wait() != 1:
#		print "There were some errors"
#	for j in proc.stderr:
#		print(json.loads(j))
#	proc.stderr.close()	
#	for i in proc.stdout:
#		print json.loads(i)
	info = gochem.get_info(proc)
	mod=gochem.get_model(proc,info,0)

	cmd.load_model(mod,sel+"_Htmp",discrete=1,zoom=1)
	modR=cmd.get_pdbstr(sel+"_Htmp")
	cmd.read_pdbstr(modR,sel+"_H")
	cmd.delete(sel+"_Htmp")
	print("Jumalauta y wea")

	


#The Tk interface (I know that  the interface is awful, but I have no time to study Tk at present)
#Please feel free to improve it if you can/want, and please e-mail me your version.
#maybe it would be better to make this app a pymol script, so we don't need this interface        
def goReduceDialog(app): 

	sel = tkinter.simpledialog.askstring("goReduce",
                                       'Enter the selection to protonate',
                                       parent=app.root)
	jsoner(sel)
 
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
                             'goReduce',
                             label = 'goReduce',
	command = lambda s=self : goReduceDialog(s))


