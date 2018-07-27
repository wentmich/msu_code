#chmod u+x clean.sh; clean.sh

# delete files that are generated automatically
rm -f TMP_*.txt RKLOG.DAT *.lis foxyinp.dat TMP_Enge_params_for_COSY.txt *~ 

# Move some of the relavent generated files to some directory
#rmdir output; mkdir output
#mv COSY*.txt BA*.txt *.png output/

# Optional: Delete the abive files, instead
#rm -f COSY*.txt BA*.txt *.png