package require bigdcd
package require pbctools

proc CalXYZ {frame} {
#	global numcp upleaflet lowleaflet
	global numcp Zcom Zlowmax Zupmin
	global upleafletindex lowleafletindex
	set upleaflet [atomselect top "index $upleafletindex"]
	set lowleaflet [atomselect top "index $lowleafletindex"]
	set Zunitcell [lindex [lindex [pbc get] 0] 2]
	set upxyz [$upleaflet get {x y z}]
	set lowxyz [$lowleaflet get {x y z}]
	set outfileup [open "DAT/upxyzf${frame}.dat" w]
	set outfilelow [open "DAT/lowxyzf${frame}.dat" w]
	set outfileM2 [open "DAT/M2xyzf${frame}.dat" w]
	puts $outfileup "X\tY\tZ" 
	puts $outfilelow "X\tY\tZ" 
	puts $outfileM2 "X\tY\tZ" 
	foreach XYZ $upxyz {
		set xi [lindex $XYZ 0] 	
		set yi [lindex $XYZ 1] 	
		set zi [lindex $XYZ 2] 	
		if {$zi < $Zupmin} {
			set zi [expr $zi + $Zunitcell]
		}
		# center the z position
		set zi [expr $zi - $Zcom]
		puts $outfileup "$xi\t$yi\t$zi" 
	}
	foreach XYZ $lowxyz {
		set xi [lindex $XYZ 0] 	
		set yi [lindex $XYZ 1] 	
		set zi [lindex $XYZ 2] 	
		if {$zi > $Zlowmax} {
			set zi [expr $zi - $Zunitcell]
		}
		# center the z position
		set zi [expr $zi - $Zcom]
		puts $outfilelow "$xi\t$yi\t$zi" 
	}
	for {set i 0} {$i < $numcp} {incr i} {
		set selM2 [atomselect top "segname M2$i and name BB"]
		set comM2 [measure center $selM2]
		set xi [lindex $comM2 0]	
		set yi [lindex $comM2 1]	
		set zi [lindex $comM2 2]	
		set zi [expr $zi - $Zcom]
		puts $outfileM2 "$xi\t$yi\t$zi" 
	}
	close $outfileup 
	close $outfilelow
	close $outfileM2
}

set dz 15
mol new minim.gro 
set allPO4 [atomselect top "name PO4"]
set COM [measure center $allPO4]
set Zcom [lindex $COM 2]
set Zlowmax [expr $Zcom + $dz]
set Zupmin  [expr $Zcom - $dz]

set upleaflet0 [atomselect top "name PO4 and z > $Zcom"]
set upleafletindex [$upleaflet0 get {index}]
set lowleaflet0 [atomselect top "name PO4 and z < $Zcom"]
set lowleafletindex [$lowleaflet0 get {index}]


atomselect macro cgprotein {resname ALA ARG ASN ASP CYS GLN GLU GLY HSD HSE HSP HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL}
set selPRO1BB [atomselect top "resname SER and resid 1 and name BB"]
set numT 396
set numE 99

set numcp [$selPRO1BB num]
for {set i 0} {$i < $numcp} {incr i} {
set A1 [expr $i*$numT + $numE*0 + 1]
set A2 [expr $i*$numT + $numE*1]
set B1 [expr $i*$numT + $numE*1 + 1]
set B2 [expr $i*$numT + $numE*2]
set C1 [expr $i*$numT + $numE*2 +1]
set C2 [expr $i*$numT + $numE*3]
set D1 [expr $i*$numT + $numE*3 +1]
set D2 [expr $i*$numT + $numE*4]
[atomselect top "serial $A1 to $A2"] set chain A$i 
[atomselect top "serial $B1 to $B2"] set chain B$i 
[atomselect top "serial $C1 to $C2"] set chain C$i 
[atomselect top "serial $D1 to $D2"] set chain D$i
[atomselect top "serial $A1 to $D2"] set segname M2$i
} 

#set upleafletResID [$upleaflet get {resid}]
#set lowleafletResID [$lowleaflet get {resid}]

#set selupleaflet [atomselect top "resid $upleafletResID and name PO4"]
#set sellowleaflet [atomselect top "resid $lowleafletResID and name PO4"]


bigdcd CalXYZ eq4_no_jump_center.xtc 
bigdcd_wait
exit

