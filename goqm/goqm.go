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
	var osidemol  *chem.Topology
	var osidecoords, sidecoords *chem.VecMatrix
	var sidelist, sidefrozen []int
	selindex:=0
	selections:=len(options.AtomsPerSel)
	total:=0
	if options.BoolOptions[0][0]{  //sidechain selections exist
		sidecoords,osidecoords,osidemol,sidelist, sidefrozen=SideChains(stdin,options)
		selections--
		total+=osidemol.Len()
		selindex++
	}

	obbmol:=make([]*chem.Topology,selections,selections)
	obbcoords:=make([]*chem.VecMatrix,selections,selections)
	bbcoords:=make([]*chem.VecMatrix,selections,selections)
	bblist:=make([][]int,selections,selections)
	bbfrozen:=make([][]int,selections,selections)
	for i:=0;i<selections;i++{
		bbcoords[i],obbcoords[i],obbmol[i],bblist[i], bbfrozen[i]=BackBone(stdin,options,selindex)
		total+=obbmol[i].Len()
		selindex++
	}
//Now we put the juit together
	bigC:=chem.ZeroVecs(total)
	var bigA *chem.Topology
	bigFroz:=make([]int,0,total)
	setoffset:=0
	if options.BoolOptions[0][0]{
		bigC.SetMatrix(0,0,osidecoords)
		setoffset+=osidecoords.NVecs()
		bigA=osidemol
		bigFroz=append(bigFroz,sidefrozen...)
	}
	for k,v:=range(obbcoords){
		bigC.SetMatrix(setoffset,0,v)
		bigA,_=chem.MergeAtomers(bigA,obbmol[k])
		tmpfroz:=SliceOffset(bbfrozen[k],setoffset)
		bigFroz=append(bigFroz,tmpfroz...)
		setoffset+=v.NVecs()

	}
	bigA.SetCharge(charge)
	bigA.SetMulti(multi)
	chem.PDBWrite("OPT.pdb", bigA,bigC,nil) /////////////////////////////////////
//Ok, we have now one big matrix and one big atom set, now the optimization

	calc := new(chem.QMCalc)
	if calctype=="Optimization"{
		calc.Optimize=true
	}
	calc.RI = true //some options, including this one, are meaningless for MOPAC
	calc.CConstraints = bigFroz
	calc.Dielectric = dielectric
	calc.SCFTightness = 1
	calc.Disperssion = "D3"
	calc.Method="TPSS"
	if method == "Cheap" {
		calc.BSSE="gcp"
		if qmprogram=="ORCA"{
			calc.Method="HF-3c"
			calc.RI = false
		}else if qmprogram=="MOPAC2012"{
			calc.Method="PM6-D3H4 MOZYME"
		}else {
			calc.Basis = "def2-SVP"
		}
	} else {
		calc.Basis = "def2-TZVP"
	}
	//We will use the default methods and basis sets of each program. In the case of MOPAC, that is currently PM6-D3H4.
	var QM chem.QMRunner
	switch qmprogram {
	case "ORCA":
		QM= chem.QMRunner(chem.NewOrcaRunner())
		QM.SetnCPU(runtime.NumCPU())
	case "TURBOMOLE":
		QM=chem.QMRunner(chem.NewTMRunner())
	default:
		QM = chem.QMRunner(chem.NewMopacRunner())
	}

	QM.SetName(options.SelNames[0])
	QM.BuildInput(bigA, bigC, calc)
	if err2 := QM.Run(true); err != nil {
		log.Fatal(err2.Error())
	}
//Now we ran the calculation, we must retrive the geometry and divide the coordinates among the original selections.

	var newBigC *chem.VecMatrix
	info := new(chem.JSONInfo) //Contains the return info
	var err2 error
	if calc.Optimize{
		tmp:=chem.EmptyVecs()
		newBigC, err2 = QM.GetGeometry(bigA)
		if err2 != nil {
			log.Fatal(err2.Error())
		}
		info.Molecules = len(options.AtomsPerSel)
		geooffset:=0
		if options.BoolOptions[0][0]{
			tmp.View2(newBigC,geooffset,0,len(sidelist),3)  //This is likely to change when we agree on a change for the gonum API!!!!
			sidecoords.SetVecs(tmp,sidelist)
			info.FramesPerMolecule = []int{1}
			info.AtomsPerMolecule = []int{sidecoords.NVecs()}
			geooffset+=len(sidelist)
		}
		for k,v:=range(bbcoords){
			tmp.View2(newBigC,geooffset,0,len(bblist[k]),3)  //This is likely to change when we agree on a change for the gonum API!!!!
			v.SetVecs(tmp,bblist[k])
			info.FramesPerMolecule=append(info.FramesPerMolecule,1)
			info.AtomsPerMolecule=append(info.AtomsPerMolecule,v.NVecs())

		}
//	for k,v:=range(bbcoords){
//		chem.XYZWrite(fmt.Sprintf("opti%d.xyz",k), , newcoords)
//	}
	}else{
		//nothing here, the else part will get deleted after tests
	}
	energy,err2:=QM.GetEnergy()
	if err2!=nil{
		log.Fatal(err2.Error())
	}
	//Start transfering data back

	info.Energies=[]float64{energy}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")

	// A loop again to transmit the info.
	for _,v:=range(bbcoords){
		if err := chem.TransmitMoleculeJSON(nil, []*chem.VecMatrix{v}, nil, nil, os.Stdout); err2 != nil {
			fmt.Fprint(os.Stderr, err)
			log.Fatal(err)
		}
	}

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


func SideChains(stdin *bufio.Reader, options *chem.JSONOptions)(coords, optcoords *chem.VecMatrix, optatoms *chem.Topology, list, frozen []int) {
	mol, coordarray, err := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords = coordarray[0]
	resid, chains := GetResidueIds(mol)
//	fmt.Fprintln(os.Stderr,"SIDE! resid, chains", resid, chains)
	list = chem.CutAlphaRef(mol, chains, resid)
	optcoords = chem.ZeroVecs(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ = chem.NewTopology(nil, 0, 0) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	chem.ScaleBonds(optcoords,optatoms,"CA","HA2",chem.CHDist)
	chem.ScaleBonds(optcoords,optatoms,"CA","HA3",chem.CHDist)
	frozen = make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		if curr.Name == "HA" || curr.Name == "CA" || curr.Name == "CB" {
			frozen = append(frozen, i)
		}
	}
	return coords, optcoords, optatoms, list, frozen
}



func BackBone(stdin *bufio.Reader, options *chem.JSONOptions, i int)(coords, optcoords *chem.VecMatrix, optatoms *chem.Topology, list, frozen []int) {
	mol, coordarray, err := chem.DecodeJSONMolecule(stdin, options.AtomsPerSel[i], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords = coordarray[0]

	//chem.PDBWrite("OPTpp.pdb", mol,coords,nil) /////////////////////////////////////
	resid, chain := GetResidueIds(mol)
	fmt.Fprintln(os.Stderr,"resid, chains, atomspersel, i", resid, chain, options.AtomsPerSel[i],i)
	var err2 error
	list,err2 = chem.CutBackRef(mol, []string{chain[0]}, [][]int{resid[1:len(resid)-1]}) //in each backbone selection the chain should be the same for all residues
	if err!=nil{
		panic(err2.Error())   //at least for now
	}
	optcoords = chem.ZeroVecs(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ = chem.NewTopology(nil, 0, 0) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	chem.ScaleBonds(optcoords,optatoms,"NTZ","HNZ",chem.CHDist)
	chem.ScaleBonds(optcoords,optatoms,"CTZ","HCZ",chem.CHDist)
	frozen = make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		//In the future there could be an option to see whether C and N are fixed 
		if  curr.Name == "NTZ" || curr.Name == "CTZ" || curr.Name == "C" || curr.Name == "N"  {
			frozen = append(frozen, i)
		}
	}
	return coords, optcoords, optatoms, list, frozen
}


func SliceOffset(list []int, offset int) (rlist []int) {
	rlist=make([]int,len(list))
	for k,v:=range(list){
		rlist[k]=v+offset
	}
	return rlist
}
