int
fVol2MGH(fvol *vol, char *fn, DATATYPE dtype) {
  
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  
  ext = strrchr(fn, '.') ;
  if(ext) {
    char command[512];
    ++ext;
    
    // mgz ==> compressed. Route stdout to file.
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) {
      gzipped = 1;
      strcpy(command, "gzip -f -c > ");
      strcat(command, fn);
      errno = 0;
      fid = popen(command, "w"); // popen(): read and write to a unix pipe. E.g., fid = popen("ls -l", "r");
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
      if(errno) {
	pclose(fid);
	errno = 0;
	printf("\n %s: gzip error writing file %s\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
    else if( !strcmp(ext, "mgh") ) {
      fid = fopen(fn, "wb") ;
      if(!fid) {
	errno = 0;
	printf("\n %s: cannot open file %s for writing\n", __FILE__, fn);
	fflush(stdout);
	return EOF;
      }
    }
  }
  
  if (fid) {
    VolInfo2MGH(&(vol->info), fid, dtype);      
    int nb = GetNumberOfBytes(dtype);
    int tsize = vol->info.vxlcount;
    char *data_b, *data_l;
    data_b = (char *) malloc(tsize*nb);
    data_l = (char *) malloc(tsize*nb);
    CastfData(tsize, vol->data, data_l, dtype, 0);
    ByteSwap(data_l, data_b, tsize, dtype);
    free(data_b);
    free(data_l);
  }
  else {
    printf("\n cannot open file: %s\n", fn);
    fflush(stdout);
  }
  
  if(gzipped)
    return pclose(fid);
  else
    return fclose(fid);
}



fvol* fVolFromMGH(char *fn, int bLoadData) {
  
  FILE *fid = 0;
  char *ext;
  int  gzipped = 0;
  char command[512];
  
  ext = strrchr(fn, '.');
  if (ext) {
    ++ext;
    
    if( !strcmp(ext, "mgz") || strstr(fn, "mgh.gz") ) {
      gzipped = 1;
      strcpy(command, "zcat ");  // OR strcpy(command, "gunzip -c ");
      strcat(command, fn);
      
      errno = 0;
      fid = popen(command, "r"); // popen(): read and write to a unix pipe. E.g., fid = popen("ls -l", "r");
      if(!fid) {
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
      if(errno) {
	pclose(fid);
	errno = 0;
        printf("\n%s command: %s\n", __FILE__, command);
      }
    }
    else if(!strcmp(ext, "mgh")) {
      fid = fopen(fn, "rb");
      if(!fid) {
	errno = 0;
      }
    }
  }
  
  fvol *vol = NULL;
  if (fid) {
    volinfo info;
    VolInfoFromMGH(&info, fid);
    vol = fVolNew(info.dim, &(info.M_vxl2lph), &(info.M_pat2grad), NULL);
    if (vol && bLoadData) {
      int nb = GetNumberOfBytes(info.dtype);
      int tsize = vol->info.vxlcount;
      char *data_b, *data_l;
      data_b = (char *) malloc(tsize*nb);
      data_l = (char *) malloc(tsize*nb);
      fread((char *)data_b, nb, tsize, fid);
      ByteSwap(data_b, data_l, tsize, info.dtype);
      CastfData(tsize, vol->data, data_l, info.dtype, 1);
      free(data_b);
      free(data_l);
    }
    if(gzipped)
      pclose(fid);
    else
      fclose(fid);
  }
  else
    printf("\n cannot open file: %s", fn);
  return vol;
}
