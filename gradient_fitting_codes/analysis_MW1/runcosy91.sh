#!/bin/bash
# This script allows the .fox file in the current directory to run with cosy 
#   executable files that reside in a different directory.
#  The generated output will be put in the current directory.
#
start_time=$(date) ;
#-----define the directory where the cosy executable resides-------
fox_file=$1 ; # set variable $fox_file to fist input after command
if [ $HOSTTYPE = 'x86_64' ]; 
then
# cosypath="/cygdrive/c/Users/Mauricio/Desktop/_myfiles/cosy2014/srcfordev/" ;
 cosypath="/projects/a2400/cosy91/" ;
else
 cosypath="" ;
fi
echo $cosypath

#Use the following section only if .bin file will be used from source directory
#-------------prepare soft link to .bin file
if [ ! -f 'lnCOSY91.bin' ]; then
  echo "Creating symbolic link lnCOSY91.bin to COSY91.bin file" ;
  ln -s /projects/a2400/cosy91/COSY91.bin lnCOSY91.bin
fi

#-------------prepare foxyinp.dat file for cosy to read-------------
#### Call cosy ######################################################
#append: ##append input line into   ../foxyinp.dat ;
#append: cp $cosypath'/foxyinp.dat' . ;
#append: foxyinp='foxyinp.dat' ;
#append: echo 'update... $foxyinp' ;
#append: sed  '1i '"$no_fox"'' $foxyinp > temp ;
#append: mv temp $foxyinp ;
#create new: #Create a new foxyinp.dat file with call to .fox file
no_fox=$( echo $fox_file | sed 's/.fox//' ) ;
foxyinp='foxyinp.dat' ; rm -f $foxyinp ; echo 'new' >> $foxyinp ;
echo 'update... $foxyinp' ;
sed  '1i '"$no_fox"'' $foxyinp > temp ;
mv temp $foxyinp ;

#Run command to execute cosy with .fox file as input
echo $cosypath/cosy ;
#nohup nice $cosypath/cosy &
#nice $cosypath/cosy
$cosypath/cosy

#----------------------post cosy processing-------------------------------
make_splots='splots' ;    # =1 process SLOG output file and plot using gnuplot
#if [ $2 = 'splots' ] ; then
if [ $make_splots = 'splots' ] ; then
   infile='SLOG.TXT' ; outfile='SLOG.DAT'
   if [ -f $infile ]; then
     echo 'Process '$infile' and generate '$outfile
     sed -e 's/endline/new_line#/' $infile > tmp1
       #get rid of all new line characters
     cat tmp1 | tr -d '\n' > tmp2
     sed 's/new_line#/\n/g' tmp2 > $outfile
     rm -rf tmp1 tmp2 
       # Generate plots
     #gnuplot pic_slog.plt
     #Use PORTRAIT mode below instead for ROOT generated plots
       gnuplot scripts/pic_slog_PORTRAIT.plt
   fi
fi

stop_time=$(date)
echo 'start: '$start_time
echo 'stop:  '$stop_time
