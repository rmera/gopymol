package main

// Copyright Notice
// ================
//
// The PyMOL Plugin source code in this file is copyrighted, but you can
// freely use and copy it as long as you don't change or remove any of
// the copyright notices.
//
// ----------------------------------------------------------------------
// This PyMOL Plugin is Copyright (C) 2013 by Raul Mera-Adasme
//
//                        All Rights Reserved
//
// Permission to use, copy, modify, distribute, and distribute modified
// versions of this software and its documentation for any purpose and
// without fee is hereby granted, provided that the above copyright
// notice appear in all copies and that both the copyright notice and
// this permission notice appear in supporting documentation, and that
// the name(s) of the author(s) not be used in advertising or publicity
// pertaining to distribution of the software without specific, written
// prior permission.
//
// THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
// INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
// NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
// CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
// USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
// OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
// PERFORMANCE OF THIS SOFTWARE.
// ------------------------------

import (
	"bufio"
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/scu"
	"log"
	"os"
	"runtime"
)

//For most plugins you dont really need to report errors like this, just panicking or using log.fatal should be enough.
//Here I use JSON errors as a test. They would be useful if you want to implement some recovery behavior.

func main() {
	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	stdin := bufio.NewReader(os.Stdin)
	options, err := chem.DecodeJSONOptions(stdin)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	dielectric := options.FloatOptions[0][0]
	charge := options.IntOptions[0][0]
	multi := options.IntOptions[0][1]
	qmprogram:=options.StringOptions[0][0]
	method:=options.StringOptions[0][1]
	calctype:=options.StringOptions[0][2]
	mol, coordarray, err := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords := coordarray[0]
	resid, chains := GetResidueIds(mol)
	//	fmt.Println("resid, chains", resid, chains)
	list := chem.CutAlphaRef(mol, chains, resid)
	fmt.Fprintln(os.Stderr,"list", list)
	optcoords := chem.ZeroVecs(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ := chem.NewTopology(nil, charge, multi) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	chem.ScaleBonds(optcoords,optatoms,"CA","HA2",chem.CHDist)
	chem.ScaleBonds(optcoords,optatoms,"CA","HA3",chem.CHDist)
	//	fmt.Println("lens!", optatoms.Len(), optcoords.NumVec())
	frozen := make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		if curr.Name == "HA" || curr.Name == "CA" || curr.Name == "CB" {
			frozen = append(frozen, i)
		}
	}
	chem.PDBWrite("OPT.pdb", optatoms,optcoords,nil) /////////////////////////////////////
	//Now we set the calculation. Lets just start with mopac.
	calc := new(chem.QMCalc)
	if calctype=="Optimization"{
		calc.Optimize=true
	}
	calc.RI = true //some options, including this one, are meaningless for MOPAC
	calc.CConstraints = frozen
	calc.Dielectric = dielectric
	calc.SCFTightness = 1
	calc.Disperssion = "D3"
	calc.Method="TPSS"
	if method == "Cheap" {
		calc.BSSE="gcp"
		if qmprogram=="ORCA"{
			calc.Method="HF-3c"
			calc.RI = false
		}else{
			calc.Basis = "def2-SVP" 
		}
	} else {
		calc.Basis = "def2-TZVP"
	}
	//We will use the default methods and basis sets of each program. In the case of MOPAC, that is currently PM6-D3H4.
	var QM chem.QMRunner
	switch qmprogram {
	case "ORCA":
		QM= chem.QMRunner(chem.MakeOrcaRunner())
		QM.SetnCPU(runtime.NumCPU())
	case "TURBOMOLE":
		QM=chem.QMRunner(chem.MakeTMRunner())
	default:
		QM = chem.QMRunner(chem.MakeMopacRunner())
	}

	QM.SetName(options.SelNames[0])
	QM.BuildInput(optatoms, optcoords, calc)
	if err2 := QM.Run(true); err != nil {
		log.Fatal(err2.Error())
	}
	var newcoords *chem.VecMatrix
	var err2 error
	if calc.Optimize{
		newcoords, err2 = QM.GetGeometry(optatoms)
		if err2 != nil {
			log.Fatal(err2.Error())
		}
	}else{
		newcoords=nil
	}

	fmt.Fprint(os.Stderr, "noooooaaaaaaaaaaaaaaaaaaaaaaaaaooooooo")
	energy,err2:=QM.GetEnergy()
	if err2!=nil{
		log.Fatal(err2.Error())
	}
	fmt.Fprint(os.Stderr, "noooooooooooooooooooooooooooooo")
	chem.XYZWrite("opti.xyz", optatoms, newcoords)
	coords.SetVecs(newcoords, list)
	//Start transfering data back
	info := new(chem.JSONInfo)
	info.Molecules = 1
	info.FramesPerMolecule = []int{1}
	info.AtomsPerMolecule = []int{mol.Len()}
	info.Energies=[]float64{energy}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")
	if err := chem.TransmitMoleculeJSON(nil, []*chem.VecMatrix{coords}, nil, nil, os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err)
		log.Fatal(err)
	}

	fmt.Fprint(os.Stderr, "todo se derrumbo dentro de mi dentro de mi")

}

//This is a very inefficient way to get the residue IDs. I am kind of hoping that
//the QM optimization will allways take longer than this function.
func GetResidueIds(mol chem.Atomer) ([]int, []string) {
	residues := make([]int, 0, int(mol.Len()/10))
	chains := make([]string, 0, int(mol.Len()/10))
	for i := 0; i < mol.Len(); i++ {
		curr := mol.Atom(i)
		if !scu.IsInInt(curr.Molid, residues) {
			residues = append(residues, curr.Molid)
			chains = append(chains, curr.Chain)
		}

	}
	return residues, chains
}
