// +build plot
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
