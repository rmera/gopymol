'''
The GUI for this program is based on that of Optimize, by Osvaldo Martin.  (http://www.pymolwiki.org/index.php/optimize)
License: GNU General Public License
'''

from chempy.models import Indexed
from chempy import Bond, Atom
from pymol import cmd
from subprocess import Popen, PIPE
import array
import gochem
import json
from Tkinter import *
import Pmw


#Reads a selection, cuts the Ca--CO and N--Ca bonds, replaces CO and N with nitrogens, 
def goQM(selside="sele",selbb="",qmprogram="MOPAC2012",method="Cheap", calctype="Optimization",dielectric="-1", charge=0,multiplicity=1):
	lens=[]
	states=[]
	q1=[]
	print selside ##################3
	m=cmd.get_model(selside)
	q1.append(m)
	lens.append(len(q1[0].atom))
	states.append(1)
	bb=selbb.split(",")
	for i in bb:
		if i=="":
			continue
		q1.append(cmd.get_model(i))
		lens.append(len(q1[-1].atoms))
		states.append(1)
	bb.insert(0,selside)
	proc = Popen("goopt", shell=True, stdin=PIPE,stdout=PIPE)
	options=json.dumps({"SelNames":bb,"AtomsPerSel":lens,"StatesPerSel":states,"IntOptions":[[int(charge),int(multiplicity)]],"FloatOptions":[[float(dielectric)]],"StringOptions":[[qmprogram,method, calctype]]}) 
	proc.stdin.write(options+"\n")
	for j in q1:
		for i in j.atom:
			atom,coords=gochem.Atom2gcRef(i)
			proc.stdin.write(atom+"\n")
			proc.stdin.write(coords+"\n")
	proc.stdin.close()
	if  proc.wait() != 1:
		print "There were some errors"
	print "WTF"
	info=gochem.get_info(proc)
	print "WTF2"
	mod=gochem.get_coords(proc,q1[0],["CA","HA","HA2","HA3", "O","N","H","C"],False,info,0)
	cmd.load_model(mod,selside+"_H",discrete=1,zoom=1)
	energy=info["Energies"][0]
	print "Final energy: ", energy, " kcal/mol"


def mainDialog():
    """ Creates the GUI """
    def set_goQM():
        qmprogram = qmprog_value.get()
        method = method_value.get()
        calctype = calc_value.get()
        sidesel = selside.get()
        bbsel = selbb.get()
        goQM(sidesel, bbsel, qmprogram, method, calctype, dielectric.get(), charge.get(), multiplicity.get())
    master = Tk()
    master.title(' goOptGUI ')
    w = Label(master, text="goQM:  Quick QM calculations\n",
                                background = 'black',
                                foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=320)
    p1 = nb.add(' Calculation ')
    p2 = nb.add(' Help ')
    p3 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ Minimization TAB #################################
    group = Pmw.Group(p1,tag_text='QM calculation options')
    group.pack(fill='both', expand=1, padx=5, pady=5)
# Force Field options
    qmprog_value = StringVar(master=group.interior())
    qmprog_value.set('MOPAC2012')
    qmprog_menu = Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = 'QM Program',
                menubutton_textvariable = qmprog_value,
                items = ['MOPAC2012', 'ORCA', 'TURBOMOLE'],
                menubutton_width = 15,
        ).grid(row=0, columnspan=2)
# Method
    method_value = StringVar(master=group.interior())
    method_value.set('Cheap')
    Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = '   Method    ',
                menubutton_textvariable = method_value,
                items = ['Cheap', 'Expensive'],
                menubutton_width = 15,
        ).grid(row=1, columnspan=2)
# Type of Calculations
    calc_value = StringVar(master=group.interior())
    calc_value.set('Optimization')
    Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = '   Calc. type  ',
                menubutton_textvariable = calc_value,
                items = ['Optimization', 'SinglePoint'],
                menubutton_width = 15,
        ).grid(row=2, columnspan=2)
#
    Label(group.interior(), text='Dielectric').grid(row=3, column=0)
    dielectric = StringVar(master=group.interior())
    dielectric.set("-1")
    entry_dielectric = Entry(group.interior(),textvariable=dielectric, width=15)
    entry_dielectric.grid(row=3, column=1)
    entry_dielectric.update()


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
        selside.set('sele')
    else:
        selside.set('all')
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

# Run
    Button(p1, text="Run QM calculation!", command=set_goQM).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Dielectric: Anything that is not a 
real number equal or greater than zero counts as no 
dielectric.

SideChain Selections: Give only one selection
that can contain several, non-contiguous residues.

Backbone Selections: Give the selections separated
by commas (,). Each Backbone selection must contain
only contiguous residues.

Note: 

-What "Cheap" or "Expensive" means
depends on the program used, and also on whether
optimized geometries or only energies are obtained.

-The selections supplied must be properly protonated.


""",justify=LEFT).pack()
#    Button(p2, text="Reset View", command=colorize).pack(side=BOTTOM)
############################ About TAB ########################################
    Label(p3, text = """
goQM uses the Gochem libraries to run a QM calculation
for one or more PyMOL selections.

There is no publication on Gochem yet. If you use one 
this plugin, or any program based on Gochem before such
a publication is available, please support Gochem by 
citing it as:

Mera-Adasme, R., Ochsenfeld, C., Pesonen, J. "Gochem: 
a library for computational chemistry". 
https://www.github.com/rmera/gochem

""",justify=CENTER).pack()
    master.mainloop()



def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'goOpt',
                            label = 'goOpt',
                            command = lambda : mainDialog())
