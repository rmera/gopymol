package main

/*
 * goQM.
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
	"github.com/rmera/gochem/qm"
	"github.com/rmera/gochem/v3"
	"github.com/rmera/gochem/chemjson"
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
	options, err := chemjson.DecodeOptions(stdin)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	mainName:=options.SelNames[0]
	if len(options.AtomsPerSel)>1{
		for _,v:=range(options.SelNames[1:]){
			mainName=mainName+"__"+v    //inefficient but there should never be THAT many selections.
		}
	}
	dielectric := options.FloatOptions[0][0]
	charge := options.IntOptions[0][0]
	multi := options.IntOptions[0][1]
	qmprogram := options.StringOptions[0][0]
	method := options.StringOptions[0][1]
	calctype := options.StringOptions[0][2]
	var osidemol *chem.Topology
	var osidecoords, sidecoords *v3.Matrix
	var sidelist, sidefrozen []int
	selindex := 0
	total := 0
	selections := len(options.AtomsPerSel)
	if options.BoolOptions[0][0] { //sidechain selections exist
		sidecoords, osidecoords, osidemol, sidelist, sidefrozen = SideChains(stdin, options)
		selections--
		total += osidemol.Len()
		selindex++
	}
	fmt.Fprint(os.Stderr, selections)
	obbmol := make([]*chem.Topology, selections, selections)
	obbcoords := make([]*v3.Matrix, selections, selections)
	bbcoords := make([]*v3.Matrix, selections, selections)
	bblist := make([][]int, selections, selections)
	bbfrozen := make([][]int, selections, selections)
	for i := 0; i < selections; i++ {
		bbcoords[i], obbcoords[i], obbmol[i], bblist[i], bbfrozen[i] = BackBone(stdin, options, selindex)
		total += obbmol[i].Len()
		selindex++
		fmt.Fprint(os.Stderr, "chetumanga")
	}
	//Now we put the juit together
	bigC := v3.Zeros(total)
	bigA,_:=chem.NewTopology([]*chem.Atom{},0,0)
	bigFroz := make([]int, 0, total)
	setoffset := 0
	if options.BoolOptions[0][0] {
		bigC.SetMatrix(0, 0, osidecoords)
		setoffset += osidecoords.NVecs()
		bigA,_ = chem.MergeAtomers(bigA,osidemol)
	//	bigA = osidemol
		bigFroz = append(bigFroz, sidefrozen...)
	}
	for k, v := range obbcoords {
		bigC.SetMatrix(setoffset, 0, v)
		bigA, _ = chem.MergeAtomers(bigA, obbmol[k])
		tmpfroz := SliceOffset(bbfrozen[k], setoffset)
		bigFroz = append(bigFroz, tmpfroz...)
		setoffset += v.NVecs()

	}
	bigA.SetCharge(charge)
	bigA.SetMulti(multi)
	chem.PDBWrite(mainName+"toOPT.pdb",bigC, bigA, nil) /////////////////////////////////////
	chem.XYZWrite(mainName+"toOPT.xyz",bigC, bigA) /////////////////////////////////////
	//Ok, we have now one big matrix and one big atom set, now the optimization

	calc := new(qm.Calc)
	if calctype == "Optimization" {
		calc.Optimize = true
	}
	calc.RI = true //some options, including this one, are meaningless for MOPAC
	calc.CConstraints = bigFroz
	calc.Dielectric = dielectric
	calc.SCFTightness = 1
	calc.Dispersion = "D3"
	calc.Method = "TPSS"
	if method == "Cheap" {
		calc.BSSE = "gcp"
		if qmprogram == "ORCA" {
			calc.Method = "HF-3c"
			calc.RI = false
		} else if qmprogram == "MOPAC2012" {
			calc.Method = "PM6-D3H4 NOMM MOZYME"
		} else {
			calc.Basis = "def2-SVP"
		}
	} else {
		calc.Basis = "def2-TZVP"
	}
	//We will use the default methods and basis sets of each program. In the case of MOPAC, that is currently PM6-D3H4.
	var QM qm.Handle
	switch qmprogram {
	case "ORCA":
		orca := qm.NewOrcaHandle()
		orca.SetnCPU(runtime.NumCPU())
		QM=qm.Handle(orca)
	case "TURBOMOLE":
		QM = qm.Handle(qm.NewTMHandle())
	case "NWCHEM":
		QM = qm.Handle(qm.NewNWChemHandle())
		calc.SCFConvHelp=1
	default:
		QM = qm.Handle(qm.NewMopacHandle())
	}

	QM.SetName(mainName)
	QM.BuildInput(bigC,bigA, calc)
	fmt.Fprint(os.Stderr,options.BoolOptions)
	if options.BoolOptions[0][2]{
		return //Dry run
	}
	if err2 := QM.Run(true); err != nil {
		log.Fatal(err2.Error())
	}
	//Now we ran the calculation, we must retrive the geometry and divide the coordinates among the original selections.
	var newBigC *v3.Matrix
	info := new(chemjson.Info) //Contains the return info
	var err2 error
	if calc.Optimize {
		newBigC, err2 = QM.OptimizedGeometry(bigA)
		if err2 != nil {
			log.Fatal(err2.Error())
		}
		if qmprogram=="NWCHEM"{ //NWchem translates/rotates the system before optimizing so we need to superimpose with the original geometry in order for them to match.
			newBigC,err2=chem.Super(newBigC,bigC,bigFroz,bigFroz)
			if err2!=nil{
				log.Fatal(err2.Error())
			}
		}
		info.Molecules = len(options.AtomsPerSel)
		geooffset := 0
		if options.BoolOptions[0][0] {
			tmp:=newBigC.View(geooffset, 0, len(sidelist), 3) //This is likely to change when we agree on a change for the gonum API!!!!
			sidecoords.SetVecs(tmp, sidelist)
			info.FramesPerMolecule = []int{1}
			info.AtomsPerMolecule = []int{sidecoords.NVecs()}
			//I DO NOT understand why the next line is += len(sidelist)-1 instead of len(sidelist), but it works. If a bug appears
			//take a look at this line, and the equivalent in the for loop that follows.
			geooffset += (len(sidelist)-1)
		}
		for k, v := range bbcoords {
			//Take a look here in case of bugs.
			tmp:=newBigC.View(geooffset, 0, len(bblist[k]), 3) //This is likely to change when we agree on a change for the gonum API!!!!
			v.SetVecs(tmp, bblist[k])
			info.FramesPerMolecule = append(info.FramesPerMolecule, 1)
			info.AtomsPerMolecule = append(info.AtomsPerMolecule, v.NVecs())
			geooffset+=(len(bblist[k])-1)

		}
		//	for k,v:=range(bbcoords){
		//		chem.XYZWrite(fmt.Sprintf("opti%d.xyz",k), , newcoords)
		//	}
	} else {
		//nothing here, the else part will get deleted after tests
	}
	energy, err2 := QM.Energy()
	if err2 != nil {
		log.Fatal(err2.Error())
	}
	//Start transfering data back

	info.Energies = []float64{energy}
	if err2 := info.Send(os.Stdout); err2 != nil {
		fmt.Fprint(os.Stderr, err2)
		log.Fatal(err2)
	}
	//	fmt.Fprint(os.Stdout,mar)
	//	fmt.Fprint(os.Stdout,"\n")

	// A loop again to transmit the info.

	if options.BoolOptions[0][0] {
		if err := chemjson.SendMolecule(nil, []*v3.Matrix{sidecoords}, nil, nil, os.Stdout); err2 != nil {
			fmt.Fprint(os.Stderr, err)
			log.Fatal(err)
		}
	}
	for _, v := range bbcoords {
		fmt.Fprintln(os.Stderr,"BB transmit!", v.NVecs())
		if err := chemjson.SendMolecule(nil, []*v3.Matrix{v}, nil, nil, os.Stdout); err2 != nil {
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
		if !scu.IsInInt(curr.MolID, residues) {
			residues = append(residues, curr.MolID)
			chains = append(chains, curr.Chain)
		}

	}
	return residues, chains
}

func SideChains(stdin *bufio.Reader, options *chemjson.Options) (coords, optcoords *v3.Matrix, optatoms *chem.Topology, list, frozen []int) {
	mol, coordarray, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords = coordarray[0]
	resid, chains := GetResidueIds(mol)
	//	fmt.Fprintln(os.Stderr,"SIDE! resid, chains", resid, chains)
	toscale:=[]string{"CA","HA2","HA3"}
	if options.BoolOptions[0][1]{
		list = chem.CutAlphaRef(mol, chains, resid)
	}else{	
		list=chem.CutBetaRef(mol,chains,resid)
		toscale=[]string{"CB","HB4","HB4"} //Yes, I am doing this twice for no reason other to have 3 elements in this slice.
	}
	optcoords = v3.Zeros(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ = chem.NewTopology(nil, 0, 0) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	chem.ScaleBonds(optcoords, optatoms, toscale[0], toscale[1], chem.CHDist)
	chem.ScaleBonds(optcoords, optatoms, toscale[0], toscale[2], chem.CHDist)
	frozen = make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		if curr.Name == "HA" || curr.Name == "CA" || curr.Name == "CB" {
			frozen = append(frozen, i)
		}
	}
	return coords, optcoords, optatoms, list, frozen
}

func BackBone(stdin *bufio.Reader, options *chemjson.Options, i int) (coords, optcoords *v3.Matrix, optatoms *chem.Topology, list, frozen []int) {
	mol, coordarray, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[i], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	coords = coordarray[0]

	//chem.PDBWrite("OPTpp.pdb", mol,coords,nil) /////////////////////////////////////
	resid, chain := GetResidueIds(mol)
	fmt.Fprintln(os.Stderr, "resid, chains, atomspersel, i", resid, chain, options.AtomsPerSel[i], i)
	var err2 error
	list, err2 = chem.CutBackRef(mol, []string{chain[0]}, [][]int{resid[1 : len(resid)-1]}) //in each backbone selection the chain should be the same for all residues
	if err != nil {
		panic(err2.Error()) //at least for now
	}
	optcoords = v3.Zeros(len(list))
	optcoords.SomeVecs(coords, list)
	optatoms, _ = chem.NewTopology(nil, 0, 0) //the last 2 options are charge and multiplicity
	optatoms.SomeAtoms(mol, list)
	chem.ScaleBonds(optcoords, optatoms, "NTZ", "HNZ", chem.CHDist)
	chem.ScaleBonds(optcoords, optatoms, "CTZ", "HCZ", chem.CHDist)
	frozen = make([]int, 0, 2*len(list))
	for i := 0; i < optatoms.Len(); i++ {
		curr := optatoms.Atom(i)
		//In the future there could be an option to see whether C and N are fixed
		if curr.Name == "NTZ" || curr.Name == "CTZ" || curr.Name == "C" || curr.Name == "N" {
			frozen = append(frozen, i)
		}
	}
	return coords, optcoords, optatoms, list, frozen
}

func SliceOffset(list []int, offset int) (rlist []int) {
	rlist = make([]int, len(list))
	for k, v := range list {
		rlist[k] = v + offset
	}
	return rlist
}
