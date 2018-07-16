from __future__ import print_function

'''
The GUI for this program is based on that of Optimize, by Osvaldo Martin.  (http://www.pymolwiki.org/index.php/optimize)
License: GNU General Public License
'''

try:
    # for Python2
    from Tkinter import *   ## notice capitalized T in Tkinter 
except ImportError:
    # for Python3
    from tkinter import * 


import pickle
from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
from subprocess import Popen, PIPE
import array
import gochem
import json
import Pmw


#Reads a selection, cuts the Ca--CO and N--Ca bonds, replaces CO and N with nitrogens, 
def goQM(selside="",selbb="",selwats="",dielectric="-1", charge=0,multiplicity=1,dryrun=False):
	lens=[]
	states=[]
	q1=[]
	side=False
        wats=False
	if selwats and not " " in selwats:
		wats=True
		m=cmd.get_model(selwats)
		q1.append(m)
		lens.append(len(q1[0].atom))
		states.append(1)

	if selside and not " " in selside:
		side=True
		m=cmd.get_model(selside)
		q1.append(m)
		lens.append(len(q1[-1].atom))
		states.append(1)
	bb=selbb.split(",")
	for i in bb:
		if i=="":
			continue
		q1.append(cmd.get_model(i))
		lens.append(len(q1[-1].atom))
		states.append(1)
	if side:
		bb.insert(0,selside)
        if wats:
                bb.insert(0,selwats)
        print("HERE WE GO!!!!!!!!!!!!")
	proc = Popen("goqmmm", shell=True, stdin=PIPE) #,stdout=PIPE)
	options=json.dumps({"SelNames":bb,"AtomsPerSel":lens,"StatesPerSel":states,"IntOptions":[[int(charge),int(multiplicity)]],"FloatOptions":[[]],"StringOptions":[[]],"BoolOptions":[[side,wats,dryrun]]})
	print("side", side)
	proc.stdin.write(options.encode(encoding='UTF-8'))
	proc.stdin.write("\n".encode(encoding='UTF-8'))
	for j in q1:
		for i in j.atom:
			atom,coords=gochem.Atom2gcRef(i)
			proc.stdin.write((atom+"\n").encode(encoding='UTF-8'))
			proc.stdin.write((coords+"\n").encode(encoding='UTF-8'))
	proc.stdin.close()
	if  proc.wait() != 0:
		print("There were some errors")
	if dryrun:
		return  #in a dry run there is nothing to receive. The input files should be in the current directory.
        #Here should be the parts that do the rest of the work.
        return 



#Reads a selection, cuts the Ca--CO and N--Ca bonds, replaces CO and N with nitrogens, 
def goQMpdyn(selfixed="", selside="sele",selbb="",selwats="",dielectric="-1", charge=0,multiplicity=1,dryrun=False):
	lens=[]
	states=[]
	q1=[]
        q1ats=[]
        fixed=[]
        fixedats=[]
	side=False
        wats=False
        if selfixed and not " " in selfixed:
		m=cmd.get_model(selfixed)
		fixed.append(m)
                fixedats=m.atom
                print("%d Fixed atoms!"%len(fixedats))
	if selwats and not " " in selwats:
		wats=True
		m=cmd.get_model(selwats)
		q1.append(m)
		lens.append(len(q1[0].atom))
                q1ats=q1ats+m.atom
		states.append(1)
	if selside and not " " in selside:
		side=True
		m=cmd.get_model(selside)
		q1.append(m)
                q1ats=q1ats+m.atom
		lens.append(len(q1[-1].atom))
		states.append(1)
	bb=selbb.split(",")
	for i in bb:
		if i=="":
			continue
		q1.append(cmd.get_model(i))
                q1ats=q1ats+q1[-1].atom
		lens.append(len(q1[-1].atom))
		states.append(1)
	if side:
		bb.insert(0,selside)
        if wats:
                bb.insert(0,selwats)
        print("Snap!")
        print("QM atoms")
        QM2pickle=[]
        tot=len(q1ats)
        print("[ ",end="")
        for i,v in enumerate(q1ats):
            print("%d"%(v.id-1),end="")
            if i<(tot-1):
                print(", ",end="")
            QM2pickle.append(v.id-1) #quick&dirty AND inefficient. Sorry. I'm counting on this to be pretty quick in all cases anyway
            
        print("]")
        print("\n")
        print("Fixed atoms")
        Fixed2pickle=[]
        tot=len(fixedats)
        print("[ ",end="")
        for i,v in enumerate(fixedats):
            if tot<=5000: #we don't print if there are too many, only pickle
                print("%d"%(v.id-1),end="")
                if i<(tot-1):
                    print(", ",end="")
            Fixed2pickle.append(v.id-1) #quick&dirty AND inefficient. Sorry. I'm counting on this to be pretty quick in all cases anyway

        print("]")
        print("\n")
        pickle.dump(QM2pickle, open('QMAtoms.pkl', 'wb'))
        pickle.dump(Fixed2pickle, open('FixedAtoms.pkl', 'wb'))
	if dryrun:
		return  #in a dry run there is nothing to receive. The input files should be in the current directory.
        #Here should be the parts that do the rest of the work.
        return 


def mainDialog():
    """ Creates the GUI """
    def set_goQM():
        sidesel = selside.get()
        bbsel = selbb.get()
        watsel = selwats.get()
        fixsel = selfix.get()
        cut_alpha=False
        dodryrun=False
        if dryrun.get()==1:
           dodryrun=True
        if gromacs.get()==0:
            print("FIXSEL", fixsel) ##########
            goQMpdyn(fixsel,sidesel, bbsel, watsel, charge.get(), multiplicity.get(),dodryrun)
        else:
            goQM(sidesel, bbsel, watsel, charge.get(), multiplicity.get(),dodryrun)
    master = Tk()
    master.title(' goQMMM ')
    w = Label(master, text="goQMMM: Help preparing Gromacs/pDynamo QMMM calculations\n")#,
                               # background = 'black',
                               # foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=320)
    p1 = nb.add(' QM Selection ')
    p2 = nb.add(' Help ')
    p3 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ Calculation tab #################################
    group = Pmw.Group(p1,tag_text='QM Selection')
    group.pack(fill='both', expand=1, padx=5, pady=5)
# Side Chain Cutting scheme
    Label(group.interior(), text='Charge').grid(row=5, column=0)
    charge = StringVar(master=group.interior())
    charge.set('0')
    entry_charge = Entry(group.interior(),textvariable=charge, width=15)
    entry_charge.grid(row=5, column=1)
    entry_charge.configure(state='normal')
    entry_charge.update()

    Label(group.interior(), text='Multiplicity').grid(row=7, column=0)
    multiplicity = StringVar(master=group.interior())
    multiplicity.set('1')
    entry_multiplicity = Entry(group.interior(),textvariable=multiplicity, width=15)
    entry_multiplicity.grid(row=7, column=1)
    entry_multiplicity.configure(state='normal')
    entry_multiplicity.update()

    Label(group.interior(), text='Side-chain selection').grid(row=9, column=0)
    selside = StringVar(master=group.interior())
    names = cmd.get_names('all')
    if len(names) > 0:
        selside.set('')
    else:
        selside.set('')
    entry_selside = Entry(group.interior(),textvariable=selside, width=15)
    entry_selside.grid(row=9, column=1)
    entry_selside.configure(state='normal')
    entry_selside.update()

    Label(group.interior(), text='Backbone selections').grid(row=11, column=0)
    selbb = StringVar(master=group.interior())
    selbb.set('')
    entry_selbb = Entry(group.interior(),textvariable=selbb, width=15)
    entry_selbb.grid(row=11, column=1)
    entry_selbb.configure(state='normal')
    entry_selbb.update()

    Label(group.interior(), text='Waters selections').grid(row=13, column=0)
    selwats = StringVar(master=group.interior())
    selwats.set('')
    entry_selwats = Entry(group.interior(),textvariable=selwats, width=15)
    entry_selwats.grid(row=13, column=1)
    entry_selwats.configure(state='normal')
    entry_selwats.update()

    Label(group.interior(), text='Fixed selections').grid(row=15, column=0)
    selfix = StringVar(master=group.interior())
    selfix.set('')
    entry_selfix = Entry(group.interior(),textvariable=selfix, width=15)
    entry_selfix.grid(row=15, column=1)
    entry_selfix.configure(state='normal')
    entry_selfix.update()



    dryrun = IntVar(master=group.interior())
    dryrun.set(0)
    C2 = Checkbutton(group.interior(), text = "Dry run", variable = dryrun, \
                 onvalue = 1, offvalue = 0, height=2, \
                 ).grid(row=17,column=1)
    gromacs = IntVar(master=group.interior())
    gromacs.set(0)
    C3 = Checkbutton(group.interior(), text = "Use Gromacs", variable = gromacs, \
                 onvalue = 1, offvalue = 0, height=2, \
                 ).grid(row=18,column=1)

  #  C1.pack()

# Run
    Button(p1, text="Run QM calculation!", command=set_goQM).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text ="""
SideChain Selections: Give only one selection
that can contain one or more, non-contiguous residues.

Backbone Selections: Give the selections separated
by commas (,). Each Backbone selection must contain
only contiguous residues.


In case of error, check the following:

-Is the total charge correct?

-Do the selections have hydrogens properly placed?

""",justify=LEFT).pack()
#    Button(p2, text="Reset View", command=colorize).pack(side=BOTTOM)
############################ About TAB ########################################
    Label(p3, text = """
goQMMM uses the goChem libraries to automatize (part of)
the preparation of a Gromacs or pDynamo QMMM calculations.

There is no publication on Gochem yet. If you use one 
this plugin, or any program based on Gochem before such
a publication is available, please support Gochem by 
citing it as:

Mera-Adasme, R. Savacsi, G., Pesonen, J. 
"goChem: a library for computational chemistry". 
https://www.gochem.org

Bug reports, feature requests, etc to:
devel@gochem.org


""",justify=CENTER).pack()
    master.mainloop()



def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'goQMMM',
                            label = 'goQMMM',
                            command = lambda : mainDialog())
