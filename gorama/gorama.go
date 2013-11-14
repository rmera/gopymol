package main

/*
 * goRama.
 *
 * Copyright 2013 Raul Mera <devel@gochem.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

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
