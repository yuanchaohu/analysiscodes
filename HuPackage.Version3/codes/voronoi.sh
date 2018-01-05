#!/bin/bash
#voronoi.sh
#Carry out Voronoi analysis using voro++
#http://math.lbl.gov/voro++/about.html

for j in CuZr
do
tmp=Voronoi      #temporary directory for files
filename=$j.x     #common string in file names
dumpfile=$filename.lammpstrj    #dump file
dataatom=voro.$filename.global.dat       #single particle's value
dataneighbor=voro.$filename.neighbor.dat     #neighbor list
dataarea=voro.$filename.facearea.dat            #area of each face in voronoi cell
dataindex=voro.$filename.index.dat 
 
numatom=`sed -n '4p' ../$dumpfile`    #atom number
numconfig=`grep -c 'ITEM: TIMESTEP' ../$dumpfile` #configuration number in a dump file
((numrow=$numatom+9))
echo "Atom Number: $numatom"
echo "Configurations Nubmer: $numconfig"

mkdir -p ./$tmp
for i in `seq 1 $numconfig`  #1,2,...,numconfig
do 
  echo "id     cn     volume      area" >>./$dataatom
  echo "id     cn     neighborlist" >>./$dataneighbor
  echo "id     cn     facearea-list" >>./$dataarea
  echo "id     cnindex3-6" >>./$dataindex

  ((sstart=($i-1)*$numrow+1))  #start line of each configuration (i-1)*n+1
  ((end=$i*$numrow))          #end line of each configuration i*n
  ((next=$end+1))             #start line of the next configuration

  sed -n "$sstart,$end p; $next q" ../$dumpfile >./$tmp/dump  #extract a confuguration from dump
  cd ./$tmp

  #timeconfig=`sed -n '2p' ./dump`
  #((timeused=$timeconfig*$timestep))   #timescale of the current configuration, unit: ps
  x1=`sed -n '6p' ./dump|awk '{print $1}'`  #get the box size
  x2=`sed -n '6p' ./dump|awk '{print $2}'`
  y1=`sed -n '7p' ./dump|awk '{print $1}'`
  y2=`sed -n '7p' ./dump|awk '{print $2}'`
  z1=`sed -n '8p' ./dump|awk '{print $1}'`
  z2=`sed -n '8p' ./dump|awk '{print $2}'`

  sed -i '1,9d' ./dump #delete the head
  sort -g ./dump|awk '{print $1,$3,$4,$5}' >./dumpused  #file format for voro++

#-----------------------------------------------------------------------------------------------------------
         voro++ -p -c "%i %s %v %F @%i %A @%i %s %n @%i %s %f" $x1 $x2 $y1 $y2 $z1 $z2 ./dumpused
#*****************-p: periodic in three directions; -px: along x; see the voro++ manual*********************
#***************************http://math.lbl.gov/voro++/doc/cmd.html*****************************************
#-----------------------------------------------------------------------------------------------------------

  awk -F@ '{print $1}' ./dumpused.vol|sort -g >> ../$dataatom
  awk -F@ '{print $3}' ./dumpused.vol|sort -g >> ../$dataneighbor
  awk -F@ '{print $4}' ./dumpused.vol|sort -g >> ../$dataarea
  awk -F@ '{print $2}' ./dumpused.vol|sort -g > ./datatmp
  awk '{print $1,$5,$6,$7,$8}' ./datatmp >> ../$dataindex

  echo $i
  cd ..
done
rm -rf ./$tmp
done

echo '*****Job Done******'
