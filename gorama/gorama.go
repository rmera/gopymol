// +build plot
package main

import (
	"bufio"
	"github.com/rmera/gochem"
	"os"
)

func main() {
	mol, coords, err := chem.DecodePyMolStream(bufio.NewReader(os.Stdin))
	if err != nil {
		panic(err)
	}
	chem.FixNumbering(mol)
	ramalist, err := chem.RamaList(mol, "ABC D", []int{0, -1}) ////
	if err != nil {
		panic(err)
	}
	//	ramalist2, index := chem.RamaResidueFilter(ramalist, []string{"GLY"}, true)
	rama, err := chem.RamaCalc(coords, ramalist)
	if err != nil {
		panic(err)
	}
	//	var i int
	//	for i = 0; i < len(ramalist); i++ {
	//		if index[i] != -1 {
	//			break
	//		}
	//	}
	err = chem.RamaPlot(rama, []int{}, "Ramachandran", "Ramachandran plot")
	if err != nil {
		panic(err)
	}

}
