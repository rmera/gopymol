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
	mol, coordarray, err := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords := coordarray[0]
	resid, chains := GetResidueIds(mol)
	fmt.Println("resid, chains", resid, chains)
	list := chem.CutAlphaRef(mol, chains, resid)
	fmt.Println("list")
	optcoords := chem.ZeroVecs(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ := chem.MakeTopology(nil, options.IntOptions[0][0], options.IntOptions[0][1]-1) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	fmt.Println("lens!", optatoms.Len(), optcoords.NumVec())
	frozen := make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		if curr.Name == "HA" || curr.Name == "CA" {
			frozen = append(frozen, i)
		}
	}
	//Now we set the calculation. Lets just start with mopac.
	qmprogram := "mopac" //this is so we can easily add more program options later
	calc := new(chem.QMCalc)
	calc.RI = true //some options, including this one, are meaningless for MOPAC
	calc.CConstraints = frozen
	calc.Optimize = true
	calc.Dielectric = 4
	calc.SCFTightness = 2
	calc.Disperssion = "D3"
	calc.Method = "XXXX"
	//We will use the default methods and basis sets of each program. In the case of MOPAC, that is currently PM6-D3H4.
	var QM chem.QMRunner
	switch qmprogram {
	default:
		QM = chem.QMRunner(chem.MakeMopacRunner())
	}
	QM.SetName(options.SelNames[0])
	QM.BuildInput(optatoms, optcoords, calc)
	QM.Run(true) //wait for the result
	newcoords, err2 := QM.GetGeometry(optatoms)
	if err2 != nil {
		log.Fatal(err2.Error())
	}
	chem.XYZWrite("opti.xyz", optatoms, newcoords)
	coords.SetVecs(optcoords, list)
	//Start transfering data back
	info := new(chem.JSONInfo)
	info.Molecules = 1
	info.FramesPerMolecule = []int{1}
	info.AtomsPerMolecule = []int{mol.Len()}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")
	if err := chem.TransmitMoleculeJSON(nil, []*chem.CoordMatrix{coords}, nil, nil, os.Stdout); err2 != nil {
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
