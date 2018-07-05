{
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*  Examples of read and write to file.
//*-*  Note that this file has sub sections of codes and does not
//*-*  run by itself.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

 //-------- last used to resolve a recurring problem with fgets -----------
 //  http://www.cplusplus.com/reference/cstdio/fgets/
 // Method to use to avoid "double counting" the last line when 
 //  file has an extra \n at the end
  FILE *fp; char buff[200]; int i, cols;
  fp = fopen("TMP_COSY_out.txt", "r");
  if( fp==NULL ) { printf("No file.\n"); exit(0); }
  i = 0;
  while( !feof(fp) )
  {
  if( fgets( buff, 200, fp)!=NULL )
  {//--------------
    i++;
    printf("%d: %s", i, buff);
  }//--------------
  }
  printf(" Number of lines read: %d ", i);
  fclose(fp);


  //check if file exists without fopen and fclose
  char fstr[100] = "output.txt";
  ifstream thefile(fstr);
  if( thefile ) printf(" File exists: %s \n", fstr);

  //ask if file name is ok
  char newstr[100], fstr[100] = "output.txt";
  printf(" File to read: %s \n",fstr);
  printf(" ENTER y if ok, else enter new name:\n");
  scanf("%s", &newstr);
  if( strncmp(newstr,"y",1)!=0 ) sprintf(fstr,"%s",newstr);
  ifstream thefile(fstr);
  if( !thefile ){ printf(" File does not exist. Aborting.\n"); exit(0); }

  
//read lines from file until string 'str' is found in 
//  first characters of line of length strlen(str)
void stop_file_read_at_string(char *str, FILE *fp)
{
  char buffer[200];
  bool TheEnd=0;
  int N=0;
  while( !feof(fp) && !TheEnd )
  {
    fgets(buffer, 200, fp); cout << buffer;
    N++;
    if( strncmp(buffer,str,strlen(str))==0 )
    {
      TheEnd=1;
      printf("\n Encounterd string=%s, Lines read, N= %d \n", str, N);
      return;
      continue;
    }
  }
  printf("\n End of file reached. Last line=\n %s\n",buffer);
  return;
}
void file_read_TESTS()
{
  FILE *fp;
  fp = fopen("test_data.csv", "r");
  if( fp==NULL ) { printf("File does not exist. Aborting.\n"); exit(-1); }

  //Simple subroutine that reads each line
  if( stop_file_read_at_string("END", fp) >0 )
  
  fclose(fp);
  return;
}


  //Read format types
  //http://www.cplusplus.com/reference/cstdio/scanf/

//Most recently used and simplified form.
// Skip known number of header lines and then read until out of data or 
//  string "END" is encountered 
  FILE *fp; int NL=0, TheEnd;
  char buffer1[200]="M5_r10cm_z_Br_mapper72.table";
  fp = fopen(buffer1, "r");
  printf("read file: %s\n",buffer1);
  if (fp==NULL) { printf("no file. Aborting.\n"); return; }
  fgets(buffer1,200,fp); printf("%s",buffer1); //title line
  while(!TheEnd)
  {
    fgets(buffer1,200,fp);
    if( feof(fp) || strncmp(buffer1,"END",3)==0 )
    {
      TheEnd=1; printf(" _______Stop read. Encounted END as expected.\n");
      continue ;
    } else {
      cout << buffer1;
      NL++;
    }
  }
  fclose(fp);
  printf(" valid lines NL= %i\n",NL);


//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
//Read Multiple files from a list of names from a file
  //Read central trajectory orbit data. Ignore all but the last point for 
  //  a set of files listed in flist.txt
  system("ls output/dim_n*.orb > flist_orb.txt");
  char flist[50]="flist_orb.txt";
  char targdir[200]="";
  char fname[250];
  char buffer[200];
  FILE *list;
  FILE *fp; int TheEnd, Nf ;
  Char_t namestr[20] ;
  //name of file with list of file names
  sprintf(flist,"%s%s",targdir,flist);
  cout << flist << endl;
  list = fopen(flist,"r");
  //-----------------------
  if (list==NULL){
    printf("File error: %s: Aborting\n",flist); exit(-1);
  }
  while(!TheEnd)
  {  //LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    fgets(buffer,200,list); 
    sscanf(buffer,"%s/n", &namestr); //name of file

    if( feof(list) ) { TheEnd=1; continue; }
    sprintf(fname,"%s%s",targdir,buffer);
    fname[strlen(fname)-1]='\0';
    //cout << fname << endl;
    Nf++;
    //RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR read from current file
    fp = fopen(fname,"r");
    if (fp==NULL){
     printf("File error: %s. Skipping.\n",fname);
    } else {
     //read all lines until out of data
     while( !feof(fp) ) {
       fgets(buffer,200,fp);
     }
     printf("%i, %s: %s", Nf, fname, buffer);
     fclose(fp);
    } //RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR end current file
  }  //LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL end list
  fclose(list);
  printf("  Nf= %d, files read.\n", Nf);
//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM end


  //Example of writing to file using 
  FILE *file0 ;
  char buffer[400] ;
  Int_t cnt;
  char fname[200] = "test.txt" ;
  // "w" --- Write output to file
  file0 = fopen(fname,"w");
   if (file0!=NULL) {
     fputs ("#File written \n",file0);
     fputs ("100   \n",file0);
     fputs ("  -0.1234  \n",file0);
     fputs ("\"write sring with quotes\"\n",file0);
     fprintf(file0,"Using formating, %d \n",11.5);
     //Write column formated output
     fprintf(file0,"%2s %8s \n", "i", "xi");
     fprintf(file0,"%2d %8.2f \n", 1, 3.5);
   }
   fclose(file0) ;

  // "a" --- How to append to existing file
  // See also: http://www.cplusplus.com/reference/cstdio/fopen/
  //  r=read w=write a=append   r+ w+ a+ ...
  file0 = fopen(fname,"a");
   if (file0!=NULL) {
     fputs ("Appended line1\n",file0);
     fputs ("Appended line2\n",file0);
   }
   fclose(file0) ;


  // "r" --- Read from file
  Int_t N ;
  Float_t x ;
  TString strng ;
  file0 = fopen(fname,"r");
   if (file0!=NULL) {
     fgets(buffer,sizeof(buffer),file0); cout << endl << buffer<<endl;
     //fgets(buffer,400,file0); cout << endl << buffer<<endl;
     fgets(buffer,sizeof(buffer),file0); cout << buffer<<endl;
         sscanf(buffer,"%d ",&N);
         cout << N << endl ;
     fgets(buffer,sizeof(buffer),file0); cout << buffer<< endl;
         sscanf(buffer," %f ",&x);
         cout << x << endl ;
  /* not functional yet
     fgets(buffer,sizeof(buffer),file0); cout << buffer<< endl;
         sscanf(buffer,"%s",&strng);
         cout << strng << endl ; */
   }
   fclose(file0) ;
  

  char Gbuffer[200];
  sprintf(fname,"./");
  strcat(fname,"_in.txt"); cout << fname << endl ;
  //Read a string that is followed by commas (e.g. data from .cvs)
  //without having to use strtok
  char buffer1[200]="string1,,,,";
  sscanf(Gbuffer,"%[^','],",fname);
  //more than one string
  char buffer1[200]="string1,str2,str3,str4";
  sscanf(Gbuffer,"%[^','],%[^','],%[^','],%s",fname);
  //Note: need to be read as strings and then converted to values 
  // after (if needed) using atof() atoi() etc...

  // "r" --- Read method that can skip avoid unwanted lines and stop 
  // when no more lines exist. Used with TOSCA output files.
  Float_t x,z,y,By ;
  sprintf(fname,"TOSCA.txt");
  file0 = fopen(fname,"r");
  if (file0==NULL) {cout << "ERROR, file not found: " << fname << endl;
     cout << "aborting." << endl; return ;}
  cnt = 0;
  while(!feof(file0)){
    fgets(buffer,400,file0);
    if(feof(file0)) continue ; //needed to avoid reading last line again
    if(sscanf(buffer,"%f %f %f %f ",&x,&z,&y,&By) == 4) {
      cnt++;
      cout << cnt << " " << buffer;
    }
  }
  fclose(file0) ;
  cout << "Valid lines read= " << cnt << endl;
  //----similar method but taken from 
  //http://www.phys.ufl.edu/~acosta/root_guide.pdf
  {
    FILE *fp = fopen("/user/brun/root/basic.dat","r");
    Float_t x,y,z;
    Int_t ncols;
    Int_t nlines = 0;
    while (1) {
      ncols = fscanf(fp,"%f %f %f",&x, &y, &z);
      if (ncols < 0) break;
      nlines++;
    }
    fclose(fp);
  }
return;


  //Special: read and rewrite grouped data into separate files
  //Plot at each angle and save in separate files for inspection of
  // grouped (g) data
  const int gNN=10000;
  Float_t gzC[gNN], gxC[gNN], gxB[gNN], gPhi[gNN], gBy[gNN];
  int gN, gi, gscan, gSets; //gSets of data separated by sscanf counter events
  FILE *gfp; char prenm[10]="By_set";
  char buff1[600], buff2[600];
  //
  fp2 = fopen("By_vs_xB_all.txt","r");
  fgets( buff2, 600, fp2);
  printf(buff2); //first line to be written as first line of all groups
  gSets=0;
  while( !feof(fp2) ){
    if( feof(fp2) ) continue; //ends the entire read process
    fgets( buff1, 600, fp2);
    gN=0;
    gscan =sscanf(buff1," %f %f %f %f %f \n", &gzC[gN], &gxC[gN], &gxB[gN], &gPhi[gN], &gBy[gN]);
    while( gscan==5 ){
      fgets( buff1, 600, fp2);
      gN++;
      gscan =sscanf(buff1," %f %f %f %f %f \n", &gzC[gN], &gxC[gN], &gxB[gN], &gPhi[gN], &gBy[gN]);
      if( feof(fp2) ) break;
    }
    gSets++;
    printf(buff1); //inspect the break point(s) or last line.
    //--- done reading group of data. Now document -->
    sprintf( buff1, "%s_%i.txt",prenm, gSets);
    gfp = fopen(buff1,"w");
    fprintf(gfp, buff2);
    for( gi=0; gi<gN; gi++){
      sprintf(buff1," %f %f %f %f %f \n", gzC[gi], gxC[gi], gxB[gi], gPhi[gi], gBy[gi]);
      fprintf(gfp, buff1);
    }
    fclose(gfp);
    //-------------------------------------------- END
  }//while( !feof(fp2) )  
  fclose(fp2);
return;



//sscanf: method for reading after some column number
  file=fopen("myfile.txt","r");
  fgets(line,200,file);
  sscanf(line+14,"%f", &(vec[0]) ); /* Read starting at col 14 float into 
   matrix element vec[0]  */
  fclose(file);



  //Read until some endstring is encountered
  FILE* fp; char buffer[200];
  char fname[100]="SLOG.TXT";
  char endstr[20]="endline";
  int TheEnd=0 ;
  fp = fopen(fname, "r");
  if (fp==NULL) { printf("File not found: %s\n",fname); exit(0); }
  while (!feof(fp) && !TheEnd){
  	fgets(buffer,sizeof(buffer),fp);
  	cout << buffer ;
  	if( strncmp(buffer,endstr,7)==0 ) TheEnd=1 ;
  }
  fclose(fp);

  //Same as above but extended to more complex file
  //  with many columns and lines
  Int_t N; //number of lines read
  Int_t NC; //number of columns read
  //Read until some endstring is encountered
  FILE* fp; char buffer[200];
  char fname[100]="SLOG.TXT";
  char endstr[20]="endline";
  int TheEnd=0 ;
  fp = fopen(fname, "r");
  if (fp==NULL) { printf("File not found: %s\n",fname); exit(0); }
  //Read column names
  while (!feof(fp) && !TheEnd){
  	NC++;
  	fgets(buffer,sizeof(buffer),fp);
  	cout << buffer ;
  	if( strncmp(buffer,endstr,7)==0 ) TheEnd=1 ;
  }
  //Read data lines
  while (!feof(fp)){
   TheEnd=0;
   N++;
   while (!TheEnd){
    fgets(buffer,sizeof(buffer),fp);
    cout << buffer ;
    if( strncmp(buffer,endstr,7)==0 ) TheEnd=1 ;
   }
  }
  fclose(fp);
  printf("  Lines read N= %i \n", N);
  printf("  Columns per line NC= %i \n", NC);

//-----------------------------end-----------------------
}

