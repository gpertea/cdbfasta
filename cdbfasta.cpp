#include <ctype.h>
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include "GBase.h"
#include "GArgs.h"
#include "GHash.hh"
#include "gcdb.h"
#ifdef ENABLE_COMPRESSION
#include "gcdbz.h"
#endif
#ifdef __WIN32__
#define VERSION "cdbfasta version 1.00 win32"
#else
#define VERSION "cdbfasta version 1.00"
#endif

#define USAGE "Usage:\n\
  cdbfasta <fastafile> [-o <index_file>] [-r <record_delimiter>]\n\
   [-z <compressed_db>] [-i] [-m|-n <numkeys>|-f<LIST>]|-c|-C]\n\
    [-w <stopwords_list>] [-s <stripendchars>] [{-Q|-G}]  [-v]\n\
   \n\
   Creates an index file for records from a multi-fasta file.\n\
   By default (without -m/-n/-c/-C option), only the first \n\
   space-delimited token from the defline is used as a key.\n\
  \n\
   <fastafile> is the multi-fasta file to index; \n\
   -o the index file will be named <index_file>; if not given,\n\
      the index filename is database name plus the suffix '.cidx'\n\
   -r <record_delimiter> a string of characters at the beginning of line\n\
      marking the start of a record (default: '>')\n\
   -Q treat input as fastq format, i.e. with '@' as record delimiter\n\
      and with records expected to have at least 4 lines\n\
   -z database is compressed into the file <compressed_db>\n\
      before indexing (<fastafile> can be \"-\" or \"stdin\" \n\
      in order to get the input records from stdin)\n\
   -s strip extraneous characters from *around* the space delimited\n\
      tokens, for the multikey options below (-m,-n,-f);\n\
      Default <stripendchars> set is: '\",`.(){}/[]!:;~|><+-\n\
   -m (\"multi-key\" option) create hash entries pointing to \n\
      the same record for all tokens found in\n\
      the defline\n\
   -n <numkeys> same as -m, but only takes the first <numkeys>\n\
      tokens from the defline; when used with -a option (see below),\n\
      only collects the first <numkeys> accessions from each defline\n\
   -f indexes *space* delimited tokens (fields) in the defline as given\n\
      by LIST of fields or fields ranges (the same syntax as UNIX 'cut')\n\
   -w <stopwordslist> exclude from indexing all the words found\n\
      in the file <stopwordslist> (for options -m, -n and -k)\n\
   -i do case insensitive indexing (i.e. create additional keys for \n\
      all-lowercase tokens used for indexing from the defline \n\
   -c for deflines in the format: db1|accession1|db2|accession2|...,\n\
      only the first db-accession pair ('db1|accession1') is taken as key\n\
   -C like -c, but also subsequent db|accession constructs are indexed,\n\
      along with the full (default) token; additionally,\n\
      all nrdb concatenated accessions found in the defline \n\
      are parsed and stored (assuming 0x01 or '^|^' as separators)\n\
   -a accession mode: like -C but indexes only the 'accession' part for all\n\
      'db|accession' constructs found, plus the default first tokens\n\
   -A like -a and -C together (both accessions and 'db|accession'\n\
      constructs are used as keys\n\
   -D index each pipe ('|') delimited token found in the record identifier\n\
      (e.g. >key1|key2|key3|.. )\n\
   -d same as -D but using a custom key delimiter <kdelim> instead of the pipe\n\
      character '|'\n\
   -G FASTA records are treated as large genomic sequences (e.g. full \n\
      chromosomes/contigs) and their formatting is checked for suitability\n\
      for fast range queries (i.e. uniform line length within each record)\n\
   -v show program version and exit\n"

/*
*/

#define ERR_TOOMANYFIELDS "Error: too many fields for -f option\n"
//16K initial defline buffer
#define KBUFSIZE 0x4000
#ifndef O_BINARY
 #define O_BINARY 0x0000
#endif
#define MAX_KEYLEN 1024
//64K input buffer
#define GREADBUF_SIZE 0x010000


typedef void (*addFuncType)(char*, off_t, uint32);

char ftmp[365];
char fztmp[365];
char record_marker[127]; //record delimiter
int  record_marker_len=1;
int num_recs;
int num_keys;

off_t last_cdbfpos=0;

int compact_plus; //shortcut key and
bool acc_mode=false;
bool acc_only=false;
bool do_compress=false; // compression used
bool fastq=false;
bool gFastaSeq=false;
char keyDelim=0;

FILE* zf=NULL; //compressed file handle

//store just offset and record length
const char* defWordJunk="'\",`.(){}/[]!:;~|><+-";
char* wordJunk=NULL;
bool caseInsensitive=false; //case insensitive storage
bool useStopWords=false;
unsigned int numFields=0;
// we have fields[numFields-1]=MAX_UINT as defined in gcdb.h -- 
// as an indicator of taking every single token in the defline, 
// or for open ended ranges (e.g. -f5- )
unsigned int fields[255]; //array of numFields field indices (1-based)
GHash<int> stopList;
//static int datalen=sizeof(uint32)+sizeof(off_t);

char lastKey[MAX_KEYLEN]; //keep a copy of the last valid written key

GCdbWrite* cdbidx;
addFuncType addKeyFunc;

#define ERR_W_DBSTAT "Error writing the database statististics!\n"

void die_write(const char* fname) {
  GError("Error: cdbhash was unable to write into file %s\n",
      fname);

 }

void die_read(const char* infile) {
  GError("Error: cdbhash was unable to read the input file %s\n",infile);
}


void die_readformat(const char* infile) {
  GError("Error: bad format for file %s; is it a fastA file?\n",
       infile);
}

void die_gseqformat(const char* seqname) {
  GError("Error: invalid FASTA sequence format for %s (not uniform line length)\n",
       seqname);
}

void die_fastqformat(const char* seqid, int seqlen, int qvlen) {
  GError("Error: invalid FASTQ sequence format for %s (seqlen=%d, qvlen=%d)\n",
       seqid);
}


bool add_cdbkey(char* key, off_t fpos, uint32 reclen, int16_t linelen=0, byte elen=0) {
 
 unsigned int klen=strlen(key);
 if (fpos==last_cdbfpos && strcmp(key, lastKey)==0) return true;
 if (klen<1) {
    /* GMessage("Warning: zero length key found following key '%s'\n",
              lastKey);
    */
    return false;
    }
  //------------ adding record -----------------
 num_keys++;
 strncpy(lastKey, key, MAX_KEYLEN-1);
 lastKey[MAX_KEYLEN-1]='\0';
 if ((uint64)fpos>(uint64)MAX_UINT) { //64 bit file offset
  uint64 v= (uint64) fpos; //needed for Solaris' off_t issues with gcc/32
  /*
  if (linelen>0) {
     //record with line len
    CIdxSeqData recdata;
    recdata.fpos=gcvt_offt(&v);
    recdata.reclen=gcvt_uint(&reclen);  
    recdata.linelen=gcvt_int16(&linelen);
    recdata.elen=elen;
    if (cdbidx->addrec(key,klen,(char*)&recdata,IdxSeqDataSIZE)==-1)
      GError("Error adding cdb record with key '%s'\n",key);
     }
   
   else { //plain record */
    CIdxData recdata;
    recdata.fpos=gcvt_offt(&v);
    recdata.reclen=gcvt_uint(&reclen);  
    if (cdbidx->addrec(key,klen,(char*)&recdata,IdxDataSIZE)==-1)
      GError("Error adding cdb record with key '%s'\n",key);
    //}
  }
 else {//32 bit file offset is enough
   /*
  if (linelen>0) {
    CIdxSeqData32 recdata;
    uint32 v=(uint32) fpos;
    recdata.fpos=gcvt_uint(&v);
    recdata.reclen=gcvt_uint(&reclen);
    recdata.linelen=gcvt_int16(&linelen);
    recdata.elen=elen;
    if (cdbidx->addrec(key,klen,(char*)&recdata,IdxDataSIZE32)==-1)
      GError("Error adding cdb record with key '%s'\n",key);
    }
   else { */
    CIdxData32 recdata;
    uint32 v=(uint32) fpos;
    recdata.fpos=gcvt_uint(&v);
    recdata.reclen=gcvt_uint(&reclen);
    if (cdbidx->addrec(key,klen,(char*)&recdata, IdxDataSIZE32)==-1)
      GError("Error adding cdb record with key '%s'\n",key);
    //}
  }
 last_cdbfpos=fpos;
 return true;
}

//default indexing: key directly passed -- 
// as the first space delimited token
void addKey(char* key, off_t fpos,
             uint32 reclen) {
 num_recs++;
 add_cdbkey(key, fpos, reclen);
 if (caseInsensitive) {
   char* lckey=loCase(key);
   if (strcmp(lckey, key)!=0) 
      add_cdbkey(lckey, fpos, reclen);
   GFREE(lckey);
   }
}

//the whole defline is passed 
void addKeyMulti(char* defline,
                 off_t fpos, uint32 reclen) {
 char* p=defline;
 unsigned int fieldno=0;
 char* pn;
 num_recs++;
 bool stillParsing=true;
 unsigned int fidx=0; //index in fields[] array
 while (stillParsing) {
       while ((*p)==' ' || (*p)=='\t') p++;
       if (*p == '\0') break;
       //skip any extraneous characters at the beginning of the token
       while (chrInStr(*p, wordJunk)) p++;
       //skip any padding spaces or other extraneous characters 
       //at the beginning of the word
       if (*p == '\0') break;
       pn=p;
       while (*pn!='\0' && !isspace(*pn)) pn++; 
       //found next space or end of string
       fieldno++;
       while (fields[fidx]<fieldno && fidx<numFields-1) fidx++;
       //GMessage("> p=>'%s' [%d, %d, %d] (next='%s')\n",p, numFields, 
        //      fieldno, fields[numFields-1], pn+1);
       stillParsing = (((*pn)!='\0') && (fieldno+1<=fields[numFields-1]));
       char* pend = pn-1; //pend is on the last non-space in the current token
       //--strip the ending junk, if any
       while (chrInStr(*pend, wordJunk) && pend>p) pend--;
       if (pend<pn-1) *(pend+1)='\0';
                 else *pn='\0';
       
       if (strlen(p)>0) {
         if (fields[fidx]==MAX_UINT || fields[fidx]==fieldno) {
           if (useStopWords && stopList.hasKey(p)) {
             p=pn+1;
             continue;
             }
           //--- store this key with the same current record data:
           add_cdbkey(p, fpos, reclen);
           //---storage code ends here
           if (caseInsensitive) {
              char* lcp=loCase(p);
              if (strcmp(lcp,p)!=0)
                  add_cdbkey(lcp, fpos, reclen);
              GFREE(lcp);
              }
           }
         //if (isEnd) break; //if all the token were stored
         }
       p=pn+1;//p is now on the next token's start
       } //token parsing loop
}


int qcmpInt(const void *p1, const void *p2) {
 //int n1=*((int*)p1);
 //int n2=*((int*)p2);
 const unsigned int *a = (unsigned int *) p1;
 const unsigned int *b = (unsigned int *) p2;
 if (*a < *b) return -1;
      else if (*a > *b) return 1;
                    else return 0;
 }


char* parse_dbacc(char* pstart, char*& end_acc, char*& accst) {
 if (pstart==NULL || *pstart==0) return NULL;
 bool hasDigits=false;
 char* pend=pstart;
 end_acc=NULL; //end of first accession
 accst=NULL;
 while (*pstart=='|') pstart++;
 for(char* p=pstart;;p++) {
   if (hasDigits==false && *p>='0' && *p<='9') hasDigits=true;
   /* if (*p==0) { //end of seq_name
      pend=p; //doesn't matter if it's accession or "db"-like
      break;
      }*/
   if (*p=='|' || *p==0) {
      int curlen=p-pstart;
      if (*p==0 || (hasDigits && curlen>3) || curlen>7 || accst!=NULL) {//previous token was an accession
        pend=p; //advance pend
        if (end_acc==NULL) end_acc=p;
        if (accst==NULL) accst=pstart;
        break;
        }
      else { //first pipe char or no digits
        accst=p+1;
        }
      hasDigits=false;//reset flag
      } // | separator
   }
 if (pend!=pstart) return pend;
              else return NULL;
}



char* nextKeyDelim(char* pstart, char* pmax) {
 if (pstart==NULL || pstart>=pmax) return NULL;
 char* p=pstart;
 while (*p!=keyDelim && p!=pmax && *p!=0) p++;
 return p;
}

char* endSpToken(char* str) {
 if (str==NULL) return NULL;
 char* p=str;
 //while (*p!=' ' && *p!='\t' && *p!='\v' && *p!=0) p++;
 while (!isspace(*p) && *p!=0) p++;
 *p=0;
 return p;
}

#define NRDB_CHARSEP "\1\2\3\4"
#define NRDB_STRSEP "^|^"

//nrdbp is positioned at the end of current nrdb concatenated
//defline, or NULL if there is no more
inline void NRDB_Rec(char* &nrdbp, char* defline) {
 nrdbp=strpbrk(defline, NRDB_CHARSEP);
 if (nrdbp==NULL) {
   nrdbp=strstr(defline, NRDB_STRSEP);
   if (nrdbp!=NULL) {
     *nrdbp='\0';
     nrdbp+=2;
     }
   }
  else *nrdbp='\0';
}

//-c/-C/-a/-A indexing: key up to the first space or second '|'
//receives the full defline
void addKeyCompact(char* defline,
              off_t fpos, uint32 reclen) {
 //we got the first token found on the defline
 num_recs++;
 char* nrdb_end;
 //breaks defline at the next nrdb concatenation point
 NRDB_Rec(nrdb_end, defline);
 //isolate the first token in this nrdb concatenation
 char* token_end=endSpToken(defline); //this puts 0 at the first space encountered
 if (!compact_plus) { //shortcut key
   //only the first db|accession construct will be indexed, if found
   char* end_acc1=NULL; //end of first accession
   char* acc1st=NULL;
   char* dbacc_end=parse_dbacc(defline, end_acc1, acc1st);
   if (end_acc1!=NULL) { //has acceptable shortcut
     *end_acc1=0;
     add_cdbkey(defline, fpos, reclen);
     return;
     }
   if (dbacc_end!=NULL) {
      *dbacc_end=0;
      add_cdbkey(defline, fpos, reclen);
      return;
      }
   //store this whole non-space token as key:
   add_cdbkey(defline, fpos, reclen);
   return;
   }
 //from now on only -C/-a/-A treatment:
 int max_accs=255;
 if (numFields>0) max_accs=numFields;
 for(;;) {
    //defline is on the first token    
    if (strlen(defline)>0) //add whole non-space token as the "full key"
       add_cdbkey(defline, fpos, reclen);
    //add the db|accession constructs as keys
    char* dbacc_start=defline;
    char* firstacc_end=NULL;
    char* accst=NULL;
    char* dbacc_end=parse_dbacc(dbacc_start, firstacc_end, accst);
    int acc_keyed=0;
    while (dbacc_end!=NULL) {
      if (firstacc_end!=NULL && firstacc_end<dbacc_end) {
        char c=*firstacc_end;
        *firstacc_end=0;
        if (!acc_only)
          add_cdbkey(dbacc_start, fpos, reclen);
        if (acc_mode && accst && acc_keyed<max_accs) {
             add_cdbkey(accst, fpos, reclen);
             ++acc_keyed;
             }
        *firstacc_end=c;
        }
      if (dbacc_start==defline && dbacc_end==token_end) {
           if (acc_mode && accst!=NULL && accst!=dbacc_start)
              add_cdbkey(accst, fpos, reclen);
           break; //the whole seq_name was only one db entry
           }
      *dbacc_end=0; //end key here
      if (!acc_only)
        add_cdbkey(dbacc_start, fpos, reclen);
      if (acc_mode && accst && acc_keyed<max_accs) {
        add_cdbkey(accst, fpos, reclen);
        ++acc_keyed;
        }
      if (dbacc_end==token_end)
           break; //reached the end of this whole seq_name (non-space token)
      dbacc_start=dbacc_end+1;
      firstacc_end=NULL;
      dbacc_end=parse_dbacc(dbacc_start, firstacc_end, accst);
      } //parse inside a defline
    // -- get to next concatenated defline, if any:
    if (compact_plus && nrdb_end!=NULL) { 
      defline=nrdb_end+1; //look for the next nrdb concatenation
      NRDB_Rec(nrdb_end, defline);
      //isolate the first token in this nrdb record
      token_end=endSpToken(defline);
      }
     else break;
   } //for
 }

void addKeyDelim(char* defline,
              off_t fpos, uint32 reclen) {
 //we got the first token found on the defline
 num_recs++;
 char* nrdb_end;
 //breaks defline at the next nrdb concatenation point
 NRDB_Rec(nrdb_end, defline);
 //isolate the first space-delimited token in this nrdb concatenation
 char* token_end=endSpToken(defline);
 for(;;) {
    //defline is on the first token
    if (strlen(defline)==0) break; 
    //add the db|accession constructs as keys
    char* k_start=defline;
    char* k_end=NULL;
    while ((k_end=nextKeyDelim(k_start, token_end))!=NULL) {
        *k_end=0;
        add_cdbkey(k_start, fpos, reclen);
        k_start=k_end+1;
        }
    // -- get to next concatenated defline, if any:
    if (nrdb_end!=NULL) {
      defline=nrdb_end+1; //look for the next nrdb concatenation
      NRDB_Rec(nrdb_end, defline);
      //isolate the first space-delimited token in this nrdb record
      token_end=endSpToken(defline);
      }
     else break;
   } //for
 }

int readWords(FILE* f, GHash<int>& xhash) {
  int c;
  int count=0;
  char name[256];
  int len=0;
  while ((c=getc(f))!=EOF) {
    if (isspace(c) || c==',' || c==';') {
      if (len>0) {
        name[len]='\0';
        xhash.Add(name, new int(1));
        count++;
        len=0;
        }
      continue;
      }
    //a non-space
    name[len]=(char) c;
    len++;
    if (len>255) {
      name[len-1]='\0';
      GError("Error reading words file: token too long ('%s') !\n",name);
      }
    }
 if (len>0) {
   name[len]='\0';
   xhash.Add(name, new int(1));
   count++;
   }
 return count;
}



//========================== MAIN ===============================
int main(int argc, char **argv) {
  FILE* f_read=NULL;
  off_t fdbsize;
  int ch;
  char* zfilename;
  char* fname;
  char* marker; //record marker
  int maxkeys=0;
  int multikey=0;
  record_marker[0]='>';
  record_marker[1]=0;
  GArgs args(argc, argv, "icvDQCaAmn:o:r:z:w:f:s:d:");
  int e=args.isError();
  if  (e>0)
     GError("%s Invalid argument: %s\n", USAGE, argv[e] );
  if (args.getOpt('v')!=NULL) {
    printf("%s\n",VERSION);
    return 0;
    }
  fastq = (args.getOpt('Q')!=NULL);
  gFastaSeq=(args.getOpt('G')!=NULL);
  if (fastq && gFastaSeq)
    GError("Error: options -Q and -G are mutually exclusive.\n");
  if (args.getOpt('D')) keyDelim='|';
     else if (args.getOpt('d')!=NULL) keyDelim=args.getOpt('d')[0];
  multikey = (args.getOpt('m')!=NULL);
  if (multikey) {
    fields[numFields]=1;
    numFields++;
    fields[numFields]=MAX_UINT;
    numFields++;  
    }
  caseInsensitive = (args.getOpt('i')!=NULL);
  acc_only=(args.getOpt('a')!=NULL);
  acc_mode=(acc_only || args.getOpt('A')!=NULL);
  compact_plus=(args.getOpt('C')!=NULL || acc_mode);
  wordJunk = (char *)args.getOpt('s');
  if (wordJunk==NULL) wordJunk=(char*)defWordJunk;
  int compact=(args.getOpt('c')!=NULL || compact_plus);
  if ((compact && multikey) || (multikey && keyDelim)) {
    GError("%s Error: invalid flag combination.\n", USAGE);
    }
  char* argfields = (char*)args.getOpt('f');
  char* s = (char*)args.getOpt('n');
  if (s!=NULL) {
    maxkeys = atoi(s);
    if (maxkeys>1 && (multikey || keyDelim || argfields!=NULL))
        GError("%s Error: invalid options (-m, -c/C, -n, -D/-d and -f are exclusive)\n", USAGE);
    numFields=maxkeys;
    if (numFields<1) GError("Error: invalid -n option (must be a 1..255 value)\n");
    if (numFields>254) GError(ERR_TOOMANYFIELDS);
    if (!acc_only) {
         multikey=1;
         for (unsigned int i=1;i<=numFields;i++) 
              fields[i-1]=i;
         }
    }
  
  if (argfields!=NULL) { //parse all the field #s
    if (multikey || keyDelim || compact) 
        GError("%s Error: invalid options (-m, -c/C, -n, -D/-d and -f are exclusive)\n", USAGE);
    char* pbrk;
    int prevnum=0;
    char prevsep='\0';
    numFields=0;
    char sep;
    char *p=argfields;
    do {
     pbrk=strpbrk(p,",-");
     if (pbrk==NULL) {
         sep='\0';
         pbrk=p+strlen(p);
         if (prevsep == '-' && *p=='\0' && prevnum>0) {
           //open ended range -- ending with '-'
           fields[numFields]=prevnum;
           numFields++;
           if (numFields>253) GError(ERR_TOOMANYFIELDS);
           fields[numFields]=MAX_UINT;
           numFields++;
           //GMessage("--- stored %d, %d\n",prevnum, MAX_UINT);
           break;
           }// ending with '-'
         } // '\0'
      else { sep=*pbrk; *pbrk = '\0'; }
     int num = atoi(p);
     if (num<=0 || num>254 )
              GError("%s Error: invalid syntax for -f option.\n", USAGE);
     if (prevsep == '-') { //store a range
       for (int i=prevnum;i<=num;i++) {
          fields[numFields]=i;
          numFields++;
          if (numFields>254) GError(ERR_TOOMANYFIELDS);
          }
       }
      else if (sep!='-') {
         fields[numFields]=num;
         numFields++;
         if (numFields>254) GError(ERR_TOOMANYFIELDS);
         }
      
     prevsep=sep;
     prevnum=num;
     p=pbrk+1;
     } while (sep != '\0'); //range parsing loop
    if (numFields<=0 || numFields>254 )
              GError("%s Error at parsing -f option.\n", USAGE);
    //GMessage("[%d] Fields parsed (%d values):\n", sizeof(fields[0]), numFields);
    qsort(fields, numFields, sizeof(fields[0]), &qcmpInt);    
    multikey=1;
    /*-- --------debug:
    for (unsigned int i=0;i<numFields-1;i++) {
      GMessage("%d,", fields[i]);
      }
    GMessage("%d\n",fields[numFields-1]);
    exit(0); */ 
    } //fields
  if (fastq) {
    record_marker[0]='@';
    record_marker_len=1;
    } 
  if (args.getOpt('r')!=NULL) {//non-FASTA record delimiter?
   if (fastq) {
     GMessage("Option -r ignored because -Q was given (-Q sets delimiter to '@')\n");
     }
   else {
   marker=(char*)args.getOpt('r'); //
   int v=0;
   if (strlen(marker)>126) 
      GError("Error: the specified record delimiter is too long. "
        "Maximum accepted is 126\n");
   //special case: hex (0xXX) and octal codes (\XXX) are accepted, only if by themselves
   if (strlen(marker)==4 && (marker[0]=='\\' || (marker[0]=='0' && toupper(marker[1])=='X') )) {
       if (marker[0]=='\\') {
           marker++;
           v=strtol(marker, NULL, 8);
           }
          else v=strtol(marker, NULL, 16);
       if (v==0 || v>255)
         GError("Invalid record delimiter: should be only one character,\n"
                "'\\NNN' (octal value), '0xNN' (hexadecimal value)");
       record_marker[0]=v;
       record_marker_len=1;
       }
     else {
      strcpy(record_marker, marker);
      record_marker_len=strlen(record_marker);
      }
    }
   }
 char* stopwords=(char*)args.getOpt('w'); //stop words filename given?
 if (stopwords!=NULL) {
  FILE* fstopwords=NULL;
  if ((fstopwords=fopen(stopwords, "r"))==NULL)
       GError("Cannot open stop words file '%s'!\n", stopwords);
  int c=readWords(fstopwords, stopList);
  GMessage("Loaded %d stop words.\n", c);
  fclose(fstopwords);
  useStopWords=(c>0);
  }
  if ((zfilename=(char*)args.getOpt('z')) !=NULL) {
    do_compress=true;
    #ifndef ENABLE_COMPRESSION
      GError("Error: compression requested but not enabled when cdbfasta was compiled\n")
    #endif    
    strcpy(fztmp,zfilename);
    strcat(fztmp,"_ztmp");
    zf=fopen(fztmp,"wb");
    if (zf==NULL)
      GError("Error creating file '%s'\n'", fztmp);
    }
  char* outfile=(char*) args.getOpt('o');
  int numfiles = args.startNonOpt();
  if (numfiles==0)
    GError("%sError: no fasta file given.\n", USAGE);
  fname=(char*) args.nextNonOpt(); //first fasta file given
  if (do_compress)  { //-------- compression case -------------------
     if (strcmp(fname, "-")==0 || strcmp(fname, "stdin")==0)
           f_read=stdin;
       else f_read= fopen(fname, "rb");
     if (f_read == NULL) die_read(fname);
     fname=zfilename; //forget the input file name, keep the output
     }
  else {//
    int fdread= open(fname, O_RDONLY|O_BINARY);
    if (fdread == -1) die_read(fname);
    struct stat dbstat;
    fstat(fdread, &dbstat);
    fdbsize=dbstat.st_size;
    close(fdread);
    f_read= fopen(fname, "rb");
    if (f_read == NULL) die_read(fname);
    }

  char idxfile[365];
  if (outfile==NULL) {
    if (do_compress) {
      strcpy(ftmp, zfilename);
      strcat(ftmp, ".cidx");
      strcpy(idxfile, ftmp);
      strcat(ftmp, "_tmp");
      }
    else {
      strcpy(ftmp, fname);
      strcat(ftmp, ".cidx");
      strcpy(idxfile, ftmp);
      strcat(ftmp, "_tmp");
      }
    //should add the process ID, time and user to make this unique?
    }
  else {
    strcpy(ftmp, outfile);
    strcpy(idxfile, outfile);
    strcat(ftmp, "_tmp");
    }

  cdbidx=new GCdbWrite(ftmp); //test if this was successful?
  if (keyDelim>0) {
     addKeyFunc=&addKeyDelim;
     }
   else {  
    if (compact)
         addKeyFunc=&addKeyCompact;
      else if (multikey)
               addKeyFunc = &addKeyMulti;
          else addKeyFunc = &addKey;
    }
  off_t recpos=0;
  off_t r=0;
  unsigned int recsize=0;
  char* key=NULL;
  bool fullDefline=(multikey || compact_plus);
  GReadBuf *readbuf = new GReadBuf(f_read, GREADBUF_SIZE);
  if (do_compress) { //---------------- compression case -------------
     if (fastq) GError("Error: sorry, compression is not supported with fastq format\n");
     //TODO: use a clever block compression/indexing scheme with virtual offsets, like bgzf
     //      -- this should take care of the fastq compression
     fdbsize=0;
     GCdbz cdbz(zf); // zlib interface
     recpos=cdbz.getZRecPos();
     while ((key=cdbz.compress(readbuf,record_marker))!=NULL) {
       recsize=cdbz.getZRecSize();
       if (!fullDefline) {
         //find first space after the record_marker and place a '\0' there
         for (int i=record_marker_len; key[i]!='\0';i++) {
           if (isspace(key[i])) { key[i]='\0';break; }
           }
         }
       addKeyFunc(key, recpos, recsize);
       recpos = cdbz.getZRecPos();
       }
     remove(zfilename);
     cdbz.compress_end();
     fclose(zf);
     //determine the size of this file:
     int ftmp= open(fztmp, O_RDONLY|O_BINARY);
     if (ftmp == -1) die_read(fztmp);
     struct stat dbstat;
     fstat(ftmp, &dbstat);
     fdbsize=dbstat.st_size;
     //rename it to the intended file name
     if (rename(fztmp,zfilename) != 0) {
       GMessage("Error: unable to rename '%s' to '%s'\n",fztmp,zfilename);
       perror("rename");
       }
    }  //compression requested  
  else { // not compressed -- plain (buffered) file access
     bool defline=false;
     bool seen_defline=false;
     int kbufsize=KBUFSIZE;
     if (fullDefline) { GMALLOC(key, KBUFSIZE); }//large defline storage buffer, just in case
                 else { GMALLOC(key, 1024); }
     key[0]=0; //no keys have been parsed yet
     int kidx=-1;
     num_recs=0;
     num_keys=0;
     char prevch=0;
     //first iteration -- for the beginning of file
     if (readbuf->peekCmp(record_marker, record_marker_len)==0) {
           //new record start found (defline)
           recpos=readbuf->getPos(); //new record pos
           defline=true; //we're in defline
           seen_defline=true;
           readbuf->skip(record_marker_len);
           kidx=0;
           }//new record start
     int linecounter=0; //number of non-header lines for current record
     int fq_recloc=(fastq) ? 0 : 4; //0=seq header line, 1=sequence string, 2=q-header line, 3=qvstring, 4=not fastq
     int fq_lendata[4]={0,0,0,0}; //keep track of fastq seq len, qv len
     int last_linelen=0; //length of last non-header line in a FASTA record
     int first_linelen=0; //length of the first non-header line in a FASTA record
     int last_eol_len=0;
     bool mustbeLastLine=false;
     bool wasEOL=false;
     while ((ch=readbuf->getch())>0) {
        bool isEOL=false;
        if (ch=='\n' || ch=='\r') {
           isEOL=true;
           last_eol_len=1;
           }
        if (wasEOL && isEOL) { 
           //skipped double-char EoL (DOS or empty lines)
           if (prevch=='\n' && ch=='\r') last_eol_len=2;
           prevch=ch;
           wasEOL=isEOL;
           continue;
           }
        if (defline) {
           // ===> in the header line (including its EOL)
           if (isEOL) {
              // --> end of header line
              if (kidx>=0) { key[kidx]=0;kidx=-1; }
              defline=false;
              last_linelen=0;
              first_linelen=0;
              linecounter=0;
              fq_recloc= (fastq) ? 0 : 4;
              mustbeLastLine=false;
              memset((void*)fq_lendata, 0, 4*sizeof(int));
              prevch=ch;
              wasEOL=isEOL;
              continue;
              }
             else {
              // --> within header, before EOL
              if (kidx>=0) { //still parsing keys
                if ((isspace(ch) || ch<31) 
                           && fullDefline==false) {
                       key[kidx]=0;
                       kidx=-1;
                       wasEOL=false;
                       prevch=ch;
                       continue;
                       }
                key[kidx]=(char)ch;
                kidx++;
                if (kidx>=kbufsize) {
                  kbufsize+=KBUFSIZE;
                  GREALLOC(key, kbufsize);
                  }
                } //still parsing keys
               } // <-- before EoL
           } // <=== in the header line
          else {
           // ===> not in header line (could be "sequence block", or before first defline
           if (isEOL) {
              // <-- at EoL in "sequence" block (record body)
              if (seen_defline) {
                 linecounter++; //counting sequence lines
                 if (gFastaSeq && linecounter>1) {
                        if (last_linelen>first_linelen) 
                                 die_gseqformat(key);
                        if (last_linelen<first_linelen) 
                                 mustbeLastLine=true;
                        }
                 }
              } // --> at EoL
             else {
              // <-- before EoL
              // could be a new header, should we check?
              if (wasEOL) {
                 // - start of a new line
                 bool newRecStart=false;
                 bool checkNewRec = (!fastq || fq_lendata[1]<=fq_lendata[3]);
                 if (checkNewRec && ch==record_marker[0]) {
                     newRecStart=(record_marker_len>1) ?
                           (readbuf->peekCmp(&record_marker[1], record_marker_len-1)==0) :
                           true;
                     if (newRecStart) {
                       // new record start (new header line coming up)
                       seen_defline=true;
                       recsize = readbuf->getPos()-recpos-1-last_eol_len; //previous recsize
                       if (recsize>(off_t)(record_marker_len+1) && key[0]!='\0') {
                           //add previous record
                           //TODO: validate gFastaSeq or fastq here?
                           if (fastq && fq_lendata[1]!=fq_lendata[3])
                                    die_fastqformat(key, fq_lendata[1], fq_lendata[3]);
                           addKeyFunc(key, recpos, recsize);
                           }
                       recpos=readbuf->getPos()-1; //new record pos (after reading this EOL)
                       if (record_marker_len>1)
                          readbuf->skip(record_marker_len-1); //skip past the record marker
                       defline=true; //we'll be in the header line in the next iteration
                       linecounter=0;
                       mustbeLastLine=false;
                       kidx=0; //reset keys
                       prevch=record_marker[record_marker_len-1];
                       wasEOL=(prevch=='\n' || prevch=='\r');
                       continue;
                       } // new record start
                   } //checking for new record start (record_marker match)
                 if (!seen_defline) {
                   prevch=ch;
                   wasEOL=isEOL;
                   continue;
                   }
                 if (mustbeLastLine && !newRecStart)
                     die_gseqformat(key);
                 // start of a new line in the "sequence" block
                 last_linelen=0;
                 mustbeLastLine=false;
                 if (fq_recloc<4) {
                   //take care of multi-line fastq parsing
                   if (fq_recloc==1 && ch=='+') {
                      fq_recloc=2; //q-header line
                      }
                     else if (fq_recloc==0 || fq_recloc==2) fq_recloc++;
                   } //fastq parsing
                 } // - start of a new non-header line (wasEoL)
              if (!seen_defline) {
                   prevch=ch;
                   wasEOL=isEOL;
                   continue;
                   }
              if (fq_recloc<4) {
                 fq_lendata[fq_recloc]++;
                 }
              last_linelen++;
              if (linecounter==0) first_linelen++;
              } // --> record body, before EoL
           } // <=== not in header line
           
        prevch=ch;
        wasEOL=isEOL;
        }//while getch
     recsize=readbuf->getPos()-recpos;
     if (recsize>0) {//add last record, if there
         if (prevch=='\n' || prevch=='\r') recsize-=last_eol_len;
         if (fullDefline && kidx>0) {//close the defline string
               if (prevch!='\n' && prevch!='\r') kidx++;
               key[kidx-1]='\0';
               }
         addKeyFunc(key, recpos, recsize);
         linecounter=0;
         //GMessage("adding key=%s\n",key);
         }
   delete readbuf;
   }
  if (f_read!=stdin) fclose(f_read);
  if (cdbidx->finish() == -1) die_write("");

  // === add some statistics at the end of the cdb index file!
  r=lseek(cdbidx->getfd(), 0, SEEK_END);
  cdbInfo info;
  memcpy((void*)info.tag, (void*)"CDBX", 4);
  info.idxflags=0;
  if (multikey) info.idxflags |= CDBMSK_OPT_MULTI;
  if (do_compress) {
      info.idxflags |= CDBMSK_OPT_COMPRESS;
      GMessage("Input data were compressed into file '%s'\n",fname);
      }
  if (compact) {
      if (compact_plus)
            info.idxflags |= CDBMSK_OPT_CADD;
         else
            info.idxflags |= CDBMSK_OPT_C;
     }
  if (gFastaSeq) {
     info.idxflags |= CDBMSK_OPT_GSEQ;
     }
  info.num_records=gcvt_uint(&num_recs);
  info.num_keys=gcvt_uint(&num_keys);
  info.dbsize=gcvt_offt(&fdbsize);
  info.idxflags=gcvt_uint(&info.idxflags);
  int nlen=strlen(fname);
  info.dbnamelen=gcvt_uint(&nlen);
  r=write(cdbidx->getfd(), fname, nlen);
  if (r!=nlen)
        GError(ERR_W_DBSTAT);
  r=write(cdbidx->getfd(), &info, cdbInfoSIZE);
  if (r!=cdbInfoSIZE)
        GError(ERR_W_DBSTAT);
  delete cdbidx;
  GFREE(key);
  remove(idxfile);
  if (rename(ftmp,idxfile) == -1)
    GError("Error: unable to rename %s to %s",ftmp,idxfile);
  GMessage("%d entries from file %s were indexed in file %s\n",
      num_recs, fname, idxfile);
  return 0;
}
