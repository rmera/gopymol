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
//	"github.com/rmera/scu"
	"log"
	"os"
	"strings"
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
	mol, coordarray, err := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[0],1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords := coordarray[0]
	var rep *os.File
	//The program itself
	for {
		reportname:=strings.Join([]string{"reduce_report",options.SelNames[0],"log"},".")
		rep,err2:=os.Create(reportname)
			if err2!=nil{
				break
			}
		defer rep.Close()
		break
	}
	var build int
	if len(options.IntOptions)==0 || len(options.IntOptions[0])==0{
		build=2
	}else{
		build=options.IntOptions[0][0]
	}
	newmol,err2:=chem.Reduce(mol,coords,build,rep)
	if err2!=nil{
		fmt.Fprint(os.Stderr,chem.MakeJSONError("process","chem.Reduce",err2)) //Not always fatal.
		if  newmol==nil{ // !strings.Contains(err2.Error(),"invalid argument"){
			log.Fatal(err2)
			}
	}
	//Start transfering data back
	info:=new(chem.JSONInfo)
	info.Molecules=1
	info.FramesPerMolecule=[]int{1}
	info.AtomsPerMolecule=[]int{newmol.Len()}
	if err2:=info.Send(os.Stdout); err2!=nil{
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")
	if err2:=chem.TransmitMoleculeJSON(newmol,newmol.Coords,nil,nil,os.Stdout); err2!=nil{
		fmt.Fprint(os.Stderr,err2)
		log.Fatal(err2)
	}

	fmt.Fprint(os.Stderr,"todo se derrumbo dentro de mi dentro de mi")
}


