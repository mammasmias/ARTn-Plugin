#!/bin/sh

#
#  Convert xyz position to lmp input format
#   We add the box 
#

EXE="$CODE/Code_En_Vrac/xyz2lmp.x"

echo "EXECUTABLE :" $EXE

if [ ${#} -eq 0 ]
then
  echo "usage: ./xyz2lmp.sh /path_to_xyz_file/"
  exit 1
fi
FOLDER=$1
FILES=$(ls $1)

echo "LIST OF INPUT : " $FILES

for file in $FILES
do
  #flmp=$(echo $file | sed "s/.xyz/.lmp/p")
  #echo "$EXE $FOLDER/$file > $(echo "$file" | sed -n "s/.xyz/.lmp/p")"
  $EXE < $FOLDER/$file > $(echo "$file" | sed -n "s/.xyz/.lmp/p")
done
