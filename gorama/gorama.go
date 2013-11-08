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
	"github.com/rmera/gochem/chemplot"
	"github.com/rmera/scu"
	"log"
	"os"
	"strings"
)

//For most plugins you dont really need to report errors like this, just panicking or using log.fatal should be enough.
//Here I use JSON errors as a test. They would be useful if you want to implement some recovery behavior.

func main() {
	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	stdin := bufio.NewReader(os.Stdin)
	options, errj := chem.DecodeJSONOptions(stdin)
	if errj != nil {
		fmt.Fprint(os.Stderr, errj.Marshal())
		log.Fatal(errj)
	}
	mols := make([]*chem.Topology, 0, len(options.SelNames))
	coordset := make([]*chem.VecMatrix, 0, len(options.SelNames))
	for k, _ := range options.SelNames {
		mol, coords, errj := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[k],1)
		if errj != nil {
			fmt.Fprint(os.Stderr, errj.Marshal())
			log.Fatal(errj)
		}
		mols = append(mols, mol)
		coordset = append(coordset, coords[0])
	}
	//The program itself
	ramadata := make([][][]float64, 0, 0)
	var HLS [][]int
	var HL []int
	for k, mol := range mols {

		fmt.Println("len in go", mol.Len(), coordset[k].NVecs()) //////
		HL = []int{}
		oldres1 := mol.Atom(0).Molid + 1 //the residues should be contiguous!!!
		chem.FixNumbering(mol)
		ramalist, errj := chemplot.RamaList(mol, "ABC DEFGHI", []int{0, -1}) ////
		if errj != nil {
			log.Fatal(errj)
		}
		ramalist2, index := chemplot.RamaResidueFilter(ramalist, options.StringOptions[0], false)
		rama, errj := chemplot.RamaCalc(coordset[k], ramalist2)
		if errj != nil {
			log.Fatal(errj)
		}
		ramadata = append(ramadata, rama)
		var i int
		if options.IntOptions != nil && options.IntOptions[0] != nil {
			for i = 0; i < len(ramalist); i++ {
				fmt.Println(i, i+oldres1, index[i], options.IntOptions[0])
				if index[i] != -1 && scu.IsInInt(i+oldres1, options.IntOptions[0]) {
					HL = append(HL, index[i])
					fmt.Println("NAME:", ramalist[index[i]].Molname)
				}
			}
			HLS = append(HLS, HL)
		}
	}
	name := append(options.SelNames, "Rama")
	var err error
	if len(ramadata) == 1 {
		err = chemplot.RamaPlot(ramadata[0], HL, "Ramachandran plot", strings.Join(name, "_"))
	} else {
		err = chemplot.RamaPlotParts(ramadata, HLS, "Ramachandran plot", strings.Join(name, "_"))
	}
	if err != nil {
		fmt.Fprint(os.Stderr, chem.MakeJSONError("process", "main", err).Marshal())
	}
}
