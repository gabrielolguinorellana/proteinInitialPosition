########################################################################
# initialPosition_v1.0.tcl                                             #
# ---------------------------------------------------------------------#
# Este  script en TCL fija la  posición inicial de una  proteína  para # 
# que interactúe con una superficie.                                   #
#                                                                      #
# Opera  pera  realizando una  exploración exhausiva  de las  posibles #
# orientaciones  iniciales de  la  proteína (con  una  profundidad  de #
# búsqueda  determinada)  y en cada una  de ellas cuenta y nombra  los # 
# residuos   que  favorecen  y  desfavorecen  la  interacción  con  la #
# superficie  desde  el  aminoácido  más cercano  a  ella,  hasta  una #
# distancia dada por el usuario. Luego, fija como posición inicial más #
# optima       aquella       en        la       que          operación # 
# "residuosQueFavorecen-residuosQueDesfavorecen"  obtiene el valor más # 
# alto.                                                                #
#                                                                      # 
# Finalmente, genera los archivos de salida.                           #
#                                                                      # 
#                                                                      #    
# Entradas:                                                            #
# 		- nombre del archivo .pdb de la proteína                       #
#		- nombre del archivo .pdb de la superficie                     #
#		- distancia a la cual estará la proteína de la superficie      #
#		- ángulo de variación de la rotación                           # 
#		- distancia  hasta   la   cual se contarán residuos (desde  la #  
#		  superficie)                                                  #
#		- aminoácidos que desfavorecen la interacción                  #                
#         entre la proteína y la superficie                            #
#		- aminoácidos que favorecen la interacción                     #                
#         entre la proteína y la superficie                            #
#                                                                      #
# Salidas:                                                             #
#		- archivo .pdb de la proteína en su nueva orientación          #
#                                                                      #
########################################################################

# Entradas: 

set proteinName "XXX"     ;# Nombre del .pdb de la proteina 
set surfaceName "YYY";     # el nombre del .pdb de la superficie
set angle 15     ;# Ángulo de variación de la rotación de la proteína
set cutoff 12     ;# Distancia hasta la cual se contarán residuos
set surfProtDist 2     ;# Distancia a la cual se quiere la proteína de la superficie
set disfavorResidues "resname AAA BBB CCC"     ;# Aminoácidos que mantienen la estructura de la proteína
set favorResidues "resname DDD EEE FFF"     ;# Aminoácidos que favorecen la interacción proteina/superficie 
set disfavorAtoms "name CA"	; # don't touch it
set favorAtoms "name ND1 ND2 NE NE1 NE2 NH1 NH2 NZ OD1 OD2 OE1 OE2 OG OG1 OH"	; # don't touch it
#----------

mol delete all
mol load pdb $proteinName.pdb

set fo [open "IP_$proteinName\_angle$angle\_3.log" w]

# Carga la proteína sobre una variable y calcula su centro
set protein [atomselect top all]
set centerProt [measure center $protein]

# Determina cual es la afinidad que tiene la proteína con la superficie en su conformación original...
set mMProt [measure minmax $protein]
set interactDist [expr [lindex $mMProt 0 2]+$cutoff-$surfProtDist]

set disfavorSelection [atomselect top "protein and z<$interactDist and ($disfavorResidues) and $disfavorAtoms"] 
# (Calculando las propiedades del segmento. En orden: 1. cantidad de átomos; 2. nombre de los átomos) 
set listDisfavor [$disfavorSelection get index]

set favorSelection [atomselect top "protein and z<$interactDist and ($favorResidues) and ($favorAtoms)"]
set listFavor [$favorSelection get index]

#set maxAfinity [expr $numFavor-$numDisfavor]
set maxAfinity [expr [llength $listFavor]-[llength $listDisfavor]]
set minZValue [expr [lindex $mMProt 1 2]-[lindex $mMProt 0 2]]


# ...y la fija como la máxima.
puts $fo "Los residuos que favorecen la interacción proteína-superficie son [llength $listFavor]:" 
for {set index 0} {$index<[llength $listFavor]} {incr index} {
	  set actualAtom [atomselect top "index [lindex $listFavor $index]"]
	  set atomProp [$actualAtom get {z resname resid name}]
	  set zDist [expr [lindex $atomProp 0 0]-[lindex $mMProt 0 2]]
	  puts $fo "\t\t[lindex $atomProp 0 1][lindex $atomProp 0 2]:[lindex $atomProp 0 3] que estará a $zDist de la superficie"
} 
puts $fo "\nLos residuos que desfavorecen la interacción proteína-superficie son [llength $listDisfavor]:" 
for {set index 0} {$index<[llength $listDisfavor]} {incr index} {
	  set actualAtom [atomselect top "index [lindex $listDisfavor $index]"]
	  set atomProp [$actualAtom get {z resname resid name}]
	  set zDist [expr [lindex $atomProp 0 0]-[lindex $mMProt 0 2]]
	  puts $fo "\t\t[lindex $atomProp 0 1][lindex $atomProp 0 2]:[lindex $atomProp 0 3] que estará a $zDist de la superficie"
}
puts $fo "\nAfinidad máxima actual: $maxAfinity en la rotación x:0, y:0, z:0 y dimensión en z=$minZValue\n\n"
set maxLoops [expr 360/$angle]

# Luego, comienza el ciclo de exploración. Para cada $angle grados que se gira la proteína en x...
for {set i 1} {$i<=$maxLoops} {incr i} {
	$protein move [trans center $centerProt axis x $angle]
	# ...y en y... 
	for {set j 1} {$j<=$maxLoops} {incr j} {
		$protein move [trans center $centerProt axis y $angle]
		# ...se gira también en z, de $angle en $angle grados hasta alcanzar los 360°.
		for {set k 1} {$k<=$maxLoops} {incr k} {
			$protein move [trans center $centerProt axis z $angle]
			set mMProt [measure minmax $protein]
			set zValue [expr [lindex $mMProt 1 2]-[lindex $mMProt 0 2]]
			# Se determina nuevamente la afinidad en la nueva rotación...
			set interactDist [expr [lindex $mMProt 0 2]+$cutoff-$surfProtDist]
			
			set disfavorSelection [atomselect top "protein and z<$interactDist and ($disfavorResidues) and $disfavorAtoms"]
			set numDisfavorNext [$disfavorSelection num]
			lreplace $listDisfavor 0 [llength $listDisfavor]
			set listDisfavor [$disfavorSelection list]
			
			set favorSelection [atomselect top "protein and z<$interactDist and ($favorResidues) and ($favorAtoms)"] 
			set numFavorNext [$favorSelection num]
			lreplace $listDisfavor 0 [llength $listFavor]
			set listFavor [$favorSelection list]	
			set afinityNext [expr $numFavorNext-$numDisfavorNext]
			# ...y se evalúa si es mayor que la afinidad anterior. Sí es más afín, entonces se guardan los ángulos de esta rotación. El procedimiento se repite hasta que se alcanza la rotación completa en todos los ejes.
			if {$afinityNext>$maxAfinity || ($afinityNext==$maxAfinity && $zValue<$minZValue)} {
				set maxAfinity $afinityNext
				set minZValue $zValue
				puts $fo "Los residuos que favorecen la interacción proteína-superficie son $numFavorNext:" 
				for {set index 0} {$index<[llength $listFavor]} {incr index} {
					set actualAtom [atomselect top "index [lindex $listFavor $index]"]
					set atomProp [$actualAtom get {z resname resid name}]
					set zDist [expr [lindex $atomProp 0 0]-[lindex $mMProt 0 2]]
					puts $fo "\t\t[lindex $atomProp 0 1][lindex $atomProp 0 2]:[lindex $atomProp 0 3] que estará a $zDist de la superficie"
				}
				puts $fo "\nLos residuos que desfavorecen la interacción proteína-superficie son $numDisfavorNext:" 
				for {set index 0} {$index<[llength $listDisfavor]} {incr index} {
					set actualAtom [atomselect top "index [lindex $listDisfavor $index]"]
					set atomProp [$actualAtom get {z resname resid name}]
					set zDist [expr [lindex $atomProp 0 0]-[lindex $mMProt 0 2]]
					puts $fo "\t\t[lindex $atomProp 0 1][lindex $atomProp 0 2]:[lindex $atomProp 0 3] que estará a $zDist de la superficie"
				}
				puts $fo "\nAfinidad máxima actual: $maxAfinity en rotación x: [expr $angle*$i], rotación y: [expr $angle*$j], rotación z: [expr $angle*$k] y dimensión en z=$minZValue\n\n" 
				set xTurn [expr $angle*$i] 
				set yTurn [expr $angle*$j]
				set zTurn [expr $angle*$k]
			} 
		}
	}
}


# A continuación se aplica la rotación,
mol delete all
mol load pdb $proteinName.pdb
set protein [atomselect top all]
$protein move [trans center $centerProt axis x $xTurn]
$protein move [trans center $centerProt axis y $yTurn]
$protein move [trans center $centerProt axis z $zTurn]


set mMProt [measure minmax $protein]
set centerProt [measure center $protein]

mol load pdb $surfaceName.pdb
set surface [atomselect top all]
set mMSurf [measure minmax $surface]
set centerSurf [measure center $surface]
mol delete top

# Se calcula el vector que dará a la proteína su nueva posición
set xDifference [expr [lindex $centerProt 0]-[lindex $centerSurf 0]]
set yDifference [expr [lindex $centerProt 1]-[lindex $centerSurf 1]]
set zDifference [expr [lindex $mMProt 0 2]-[lindex $mMSurf 1 2]] 
set xDifference [expr $xDifference*-1]
set yDifference [expr $yDifference*-1]
set zDifferenceFinal [expr $zDifference*-1+$surfProtDist]
# y se aplica la translación
$protein moveby "$xDifference $yDifference $zDifferenceFinal"

set mMProt [measure minmax $protein]
puts $fo "\n########################################################"
puts $fo "Se aplicó la rotación y translación en la proteína. Los ángulos de rotación fueron x=$xTurn, y=$yTurn y z=$zTurn, el vector de traslación usado fue {$xDifference $yDifference $zDifferenceFinal} y la dimensión final en z=$minZValue"
puts $fo "########################################################"


$protein writepdb $proteinName\_IP.pdb
mol delete all
close $fo
exit
