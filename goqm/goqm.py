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
def goQM(selside="sele",selbb="",qmprogram="MOPAC2012",method="Cheap", calctype="Optimization",dielectric="-1", charge=0,multiplicity=1, alphacut=True):
	lens=[]
	states=[]
	q1=[]
	side=False
	print selside ##################3
	if selside and not " " in selside:
		side=True
		m=cmd.get_model(selside)
		q1.append(m)
		lens.append(len(q1[0].atom))
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
	proc = Popen("goqm", shell=True, stdin=PIPE,stdout=PIPE)
	options=json.dumps({"SelNames":bb,"AtomsPerSel":lens,"StatesPerSel":states,"IntOptions":[[int(charge),int(multiplicity)]],"FloatOptions":[[float(dielectric)]],"StringOptions":[[qmprogram,method, calctype]],"BoolOptions":[[side,alphacut]]})
	print "side", side
	proc.stdin.write(options+"\n")
	for j in q1:
		for i in j.atom:
			atom,coords=gochem.Atom2gcRef(i)
			proc.stdin.write(atom+"\n")
			proc.stdin.write(coords+"\n")
	proc.stdin.close()
	if  proc.wait() != 0:
		print "There were some errors"
	info=gochem.get_info(proc)
	energy=info["Energies"][0]
	print "Final energy: ", energy, " kcal/mol"
	if calctype=="Optimization":
		for k,v in enumerate(q1):
			if k==0 and side:
				exclude=["CA","HA","HA2","HA3", "O","N","H","C"]
				idexclude=[]
			else:
				idexclude=[v.atom[0].resi,v.atom[-1].resi]
				exclude=["CTZ","NTZ","HCZ","HNZ"]
			mod=gochem.get_coords(proc,v,exclude,idexclude,False,info,k)
			cmd.load_model(mod,bb[k]+"_H",discrete=1,zoom=1)


def mainDialog():
    """ Creates the GUI """
    def set_goQM():
        qmprogram = qmprog_value.get()
        method = method_value.get()
        calctype = calc_value.get()
        sidesel = selside.get()
        bbsel = selbb.get()
        cut_alpha=True
        if alpha.get()==0:
           cut_alpha=False
        goQM(sidesel, bbsel, qmprogram, method, calctype, dielectric.get(), charge.get(), multiplicity.get(),cut_alpha)
    master = Tk()
    master.title(' goQM ')
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
############################ Calculation tab #################################
    group = Pmw.Group(p1,tag_text='QM calculation options')
    group.pack(fill='both', expand=1, padx=5, pady=5)
# QM program options
    qmprog_value = StringVar(master=group.interior())
    qmprog_value.set('MOPAC2012')
    qmprog_menu = Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = 'QM Program',
                menubutton_textvariable = qmprog_value,
                items = ['MOPAC2012', 'ORCA', 'TURBOMOLE','NWCHEM'],
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
# Side Chain Cutting scheme
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
    alpha = IntVar(master=group.interior())
    alpha.set(1)
    C1= Checkbutton(group.interior(), text = "alpha cut", variable = alpha , \
                 onvalue = 1, offvalue = 0, height=2, \
                 ).grid(row=13,column=0)
#	C2 = Checkbutton(top, text = "Actually run calculation", variable = run, \
#                 onvalue = True, offvalue = False, height=5, \
#                 width = 20)
  #  C1.pack()

# Run
    Button(p1, text="Run QM calculation!", command=set_goQM).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Dielectric: A number. If less than zero, no dielectirc
is used.

SideChain Selections: Give only one selection
that can contain one or more, non-contiguous residues.

Backbone Selections: Give the selections separated
by commas (,). Each Backbone selection must contain
only contiguous residues.

What "Cheap" or "Expensive" means
depends on the program used, and the kind of calculation.

In case of error, check the following:

-Is the total charge correct?

-Do the selections have hydrogens properly placed?

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

Mera-Adasme, R. Savacsi, G., Ochsenfeld, C., Pesonen, J. 
"goChem: a library for computational chemistry". 
https://www.gochem.org

Bug reports, feature requests, etc to:
devel@gochem.org


""",justify=CENTER).pack()
    master.mainloop()



def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'goQM',
                            label = 'goQM',
                            command = lambda : mainDialog())
