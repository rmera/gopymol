'''
The GUI for this program is based on that of Optimize, by Osvaldo Martin.  (http://www.pymolwiki.org/index.php/optimize)
License: GNU General Public License
'''
import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd

def mainDialog():
    """ Creates the GUI """
    global entry_vdw, entry_elec
    def set_minimize():
        qmprogram = qmprog_value.get()
        method = method_value.get()
        selection = sel_value.get()
        minimize(selection, forcefield, method, nsteps, conv, cutoff, cut_vdw, cut_elec)

    master = Tk()
    master.title(' goOptGUI ')
    w = Tkinter.Label(master, text="goQM:  Quick QM calculations\n",
                                background = 'black',
                                foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=320)
    p1 = nb.add(' Local optimization ')
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
                items = ['MOPAC2012', 'ORCA', 'Turbomole'],
                menubutton_width = 15,
        ).grid(row=0, columnspan=2)
# Method
    method_value = StringVar(master=group.interior())
    method_value.set('Cheap')
    Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = '   Method   ',
                menubutton_textvariable = method_value,
                items = ['Cheap', 'Expensive'],
                menubutton_width = 15,
        ).grid(row=1, columnspan=2)
# Type of Calculations
    method_value = StringVar(master=group.interior())
    method_value.set('Optimization')
    Pmw.OptionMenu(group.interior(),
                labelpos = 'w',
                label_text = '   Calc. type  ',
                menubutton_textvariable = method_value,
                items = ['Optimization', 'Single Point'],
                menubutton_width = 15,
        ).grid(row=2, columnspan=2)
#
    Label(group.interior(), text='Dielectric').grid(row=3, column=0)
    dielectric = StringVar(master=group.interior())
    dielectric.set("No dielectric")
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
    selside_value = StringVar(master=group.interior())
    names = cmd.get_names('all')
    if len(names) > 0:
        selside_value.set(names[0])
    else:
        selside_value.set('all')
    entry_selside_value = Entry(group.interior(),textvariable=selside_value, width=15)
    entry_selside_value.grid(row=9, column=1)
    entry_selside_value.configure(state='normal')
    entry_selside_value.update()
    Label(group.interior(), text='Backbone selections').grid(row=11, column=0)
    selbb_value = StringVar(master=group.interior())
    selbb_value.set('')
    entry_selbb_value = Entry(group.interior(),textvariable=selbb_value, width=15)
    entry_selbb_value.grid(row=11, column=1)
    entry_selbb_value.configure(state='normal')
    entry_selbb_value.update()

# Run
    Button(p1, text="Run QM calculation!", command=goQM).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Dielectric: Anything that is not a 
positive integer counts as no dielectric.

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
                            'goOptGUI',
                            label = 'goOptGUI',
                            command = lambda : mainDialog())
