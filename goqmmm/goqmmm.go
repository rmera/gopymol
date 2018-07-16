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
	"github.com/rmera/gochem/chemjson"
	"github.com/rmera/scu"
	"log"
	"os"
)

const (
	CN=0.133 //aprox C-N distances in peptidic bond
	CcarbH=0.110 //C(carbonyl)-H distance in acetaldehyde.
	CC=0.153
	CH=0.1
)

//For most plugins you dont really need to report errors like this, just panicking or using log.fatal should be enough.
//Here I use JSON errors as a test. They would be useful if you want to implement some recovery behavior.

func main() {
	fmt.Println("NO WEI!") ///////
	//This is the part that collects all the data from PyMOL, with all  the proper error checking.
	//There are lots of things here from goQM that we don't need (for now!) in QMMM, like everything related to
	//QM programs and calculation details. For now I'll keep the data collection things, as I don't want to risk breaking the whole function.
	//I'll just not use the data later. Sorry about that :-).
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

	var sidelist []int
	var waterlist []int
	selindex := 0
//	total := 0
	selections := len(options.AtomsPerSel)
//	fmt.Println("And the selections are...", selections, options.AtomsPerSel) //////////////////////
	var sideFirstMM []*QMMMPairs
	if options.BoolOptions[0][1] { //waters selections exist
		waterlist  = Waters(stdin, options) //no link atoms for waters
		selections--
	//	total += osidemol.Len()
		selindex++
	}

	if options.BoolOptions[0][0] { //sidechain selections exist
		sidelist, sideFirstMM = SideChains(stdin, options,selindex)
		selections--
	//	total += osidemol.Len()
		selindex++
	}
	fmt.Fprint(os.Stderr, selections)
	bblistsN:=0
	bblist := make([][]int, selections, selections)
	bbFirstMM := make([][]*QMMMPairs, selections, selections)
	for i := 0; i < selections; i++ {
		bblist[i], bbFirstMM[i] = BackBone(stdin, options, selindex)
		bblistsN=+len(bblist)
		fmt.Fprint(os.Stderr, "chetumanga")
	}
	//This is way easier than for goQM, as we dont' need to worry about capping and so on.

	//Here I just ask for more memory. If just the list with the QM atoms is too much memory for you, you are probably in trouble. This list is very unlikely to have more than O(1000) elements.

	bigFirstMM:=make([]*QMMMPairs,0,len(sideFirstMM)+len(bbFirstMM)*2) //The 2 is because each sub-sequence has 2 link-atoms, one for N-ter and another for C-ter. This may not be true for terminal-subsequences. Of course the memory loss is not great.
	bigFirstMM=append(bigFirstMM,sideFirstMM...)
	biglist:=make([]int,0,len(sidelist)+len(waterlist)+bblistsN)
	biglist=append(biglist,sidelist...)
	for i,v:=range(bblist){
		biglist=append(biglist,v...)
		bigFirstMM=append(bigFirstMM,bbFirstMM[i]...) // len(bbFirstMM) really should equal len(bblist)
	}
	biglist=append(biglist,waterlist...)
	//we go for the topology file first.
	fmt.Println("Topology file!\n")
	fmt.Print("[ virtual_sites2 ]\n")
	for _,v:=range(bigFirstMM){
			var ratio=CH/CC
			if v.t=="CN"{
				ratio=CcarbH/CN
			}
		    fmt.Printf("LA %d %d 1 %5.3f\n",v.QM,v.MM,ratio)
	}
	fmt.Print("[ constraints ]\n")
	for _,v:=range(bigFirstMM){
			var l=CC
			if v.t=="CN"{
				l=CN
			}
		    fmt.Printf("%d %d 1 %5.3f\n",v.QM,v.MM,l)
	}
	//At least for now I am _not_ adding the QM-QM bonds to the topology file because that would be a nightmare. I suspect you only need those bonds in ONIOM calculations
	//which will not be supported here. Also, I suspect that I'll need to remove the QM-QM and QM-MM bonds from the bond, angle and dihedral parts of the topology. I suppose I'll write a python script to do that based on the info generated here. QUITE annoying.


	//Now for the index file. We'll simply add a "QMatoms" section with the indices of everything QM. Unfortunately, Gromacs also wants you to include the link atoms in the index, which is redundant. Here, and, at least for now, I won't do anything about it. The user will have to add the link atoms to the PDB/gro file (probably they all can have 0.0 0.0 0.0 coordinates), then to the itp file (in the end of both, I'd suggest, and then add the indices to to the ndx file. It really shouldn't be that much work, as the link atoms should never be _that_ many. Again, I can always have a Python script do it if it's too annoying. 

	fmt.Print("[ QMAtoms ]\n")
	for _,v:=range(biglist){
		fmt.Printf("%d ",v)
	}
	fmt.Print("\n")






}

/*****Modified goChem versions for QMMM. I probably should modify the ones in the library to make them more general, perhaps introducing a variadic "QMMM" flag?. For now I'm simply implementing custom QMMM versions here ******/

//CutLateralRef will return a list with the atom indexes of the lateral chains of the residues in list
//for each of these residues it will change the alpha carbon to oxygen and change the residue number of the rest
//of the backbone to -1.
func CutBetaRef(r chem.Atomer, chain []string, list []int) ([]int, []*QMMMPairs) {   //The return values are: a list with the QM atoms and a list with the binding MM atoms.
	//	pairs := make([][]int,1,10)
	//	pairs[0]=make([]int,0,2)
	exclude:=make([]int,0,10)
	tlist:=make([]int,0,len(list)*5)
	FirstMM:=make([]*QMMMPairs,0,1)
	for i := 0; i < r.Len(); i++ {
		curr := r.Atom(i)
		if scu.IsInInt(curr.MolID,list) && scu.IsInString(curr.Chain,chain) {
			tlist=append(tlist,curr.ID)
			if curr.Name == "CA" {
				FirstMM=append(FirstMM,&QMMMPairs{MM:curr.ID,t:"CC"})
			}else if curr.Name == "CB"{
				FirstMM[len(FirstMM)-1].QM=curr.ID
			}
			if scu.IsInString(curr.Name, []string{"CA","C", "H", "HA", "O", "N"}) {
				exclude=append(exclude,curr.ID)
			}
		}
	}
	newlist:=make([]int,0,len(tlist)-len(exclude))
	for _,v:=range(tlist){
		if !scu.IsInInt(v,exclude){
			newlist=append(newlist,v)
		}
	}
	return newlist,FirstMM
}


type QMMMPairs struct{
	MM int
	QM int
	t string //the bond type, CC or CN 
}


//CutBackRef takes a list of lists of residues and deletes from r
//all atoms not in the list or not belonging to the chain chain.
//It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!. It returns the list of kept atoms
//and a list with the "first" MM atoms and their corresponding QM atoms, useful to build the link atoms later.
//I suspect this will not work if you want to include the N- or C-terminal residue in the QM part.It is not our case so I'll leave it as it is, but probably
//will have to fix it at some point.
func CutBackRef(r chem.Atomer, chains []string, list [][]int) ([]int, []*QMMMPairs,error) {
	//i:=r.Len()
	if len(chains) != len(list) {
		return nil, nil,fmt.Errorf("Mismatched chains (%d) and list (%d) slices", len(chains), len(list))
	}
	ret:=make([]int,0,80) //This will be filled with the atoms that are kept, and will be returned.
	FirstMM:=make([]*QMMMPairs,0,2)
	for k, v := range list {
		nter := v[0]
		cter := v[len(v)-1]
		nresname := "" //This may not be needed!
		for j := 0; j < r.Len(); j++ {
			if r.Atom(j).MolID == nter && r.Atom(j).Chain == chains[k] {
				nresname = r.Atom(j).Molname
				break
			}
		}
		if nresname == "" {
			//we will protest if the Nter is not found. If Cter is not found we will just
			//cut at the real Cter
			return nil,nil, fmt.Errorf("list %d contains residue numbers out of boundaries", k)

		}
		fmt.Println("list",list) ////////////////////////////
		//This works under the (reasonable, I think) assumption that you will find first the N-terminal-1 residue, then the N-terminal, then the C-terminal and then the c-terminal+1
		for j := 0; j < r.Len(); j++ {
			curr := r.Atom(j)
			if curr.Chain != chains[k] {
				continue
			}
			if scu.IsInInt(curr.MolID,v){
				ret=append(ret,curr.ID)
			}

			if (curr.MolID == nter-1 && curr.Name == "C") {
				FirstMM=append(FirstMM,&QMMMPairs{MM:curr.ID,t:"CN"})
			}
			if (curr.MolID==nter && curr.Name == "N"){
				FirstMM[len(FirstMM)-1].QM=curr.ID
			}
			if (curr.MolID == cter && curr.Name=="C") {
				FirstMM=append(FirstMM,&QMMMPairs{QM:curr.ID,t:"CN"})
			}

			if curr.MolID == cter+1 && curr.Name=="N"{
				FirstMM[len(FirstMM)-1].MM=curr.ID
			}
		}
	}
//	for _, i := range list {
//		t := chem.Molecules2Atoms(r, i, chains)
//		ret = append(ret, t...)
//	}
	fmt.Println("ret",ret) ///////////////////////////////////
	return ret, FirstMM,nil
}




/************************END Cut functions ***************************/






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

func Waters(stdin *bufio.Reader, options *chemjson.Options) (list []int) {
	mol,_, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[0], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	//a, _:= GetResidueIds(mol)
	//	fmt.Fprintln(os.Stderr,"SIDE! resid, chains", resid, chains)

	list=make([]int,len(list)*5)
	for i := 0; i < mol.Len(); i++ {
			list=append(list,mol.Atom(i).ID)
	}

	return list
}



func SideChains(stdin *bufio.Reader, options *chemjson.Options,i int) (list []int, FirstMM []*QMMMPairs) {
	mol,_, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[i], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	resid, chains := GetResidueIds(mol)
	//	fmt.Fprintln(os.Stderr,"SIDE! resid, chains", resid, chains)
	list,FirstMM=CutBetaRef(mol,chains,resid)
	return list,FirstMM
}

func BackBone(stdin *bufio.Reader, options *chemjson.Options, i int) (list []int, FirstMM []*QMMMPairs) {
	mol, _, err := chemjson.DecodeMolecule(stdin, options.AtomsPerSel[i], 1)
	if err != nil {
		fmt.Fprint(os.Stderr, err.Marshal())
		log.Fatal(err)
	}
	//chem.PDBWrite("OPTpp.pdb", mol,coords,nil) /////////////////////////////////////
	resid, chain := GetResidueIds(mol)
	fmt.Fprintln(os.Stderr, "resid, chains, atomspersel, i", resid, chain, options.AtomsPerSel[i], i)
	var err2 error
	list, FirstMM, err2 = CutBackRef(mol, []string{chain[0]}, [][]int{resid[1 : len(resid)-1]}) //in each backbone selection the chain should be the same for all residues
	if err != nil {
		panic(err2.Error()) //at least for now
	}
	return list, FirstMM
}

func SliceOffset(list []int, offset int) (rlist []int) {
	rlist = make([]int, len(list))
	for k, v := range list {
		rlist[k] = v + offset
	}
	return rlist
}
