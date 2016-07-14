/* (c) 2015 Giorgio Gonnella & Stefan Dang, ZBH, Uni-Hamburg */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <libgen.h>

#define VERSION "0.1"
#define MAX_READLEN 1000UL
#define LOG10_MAX_READLEN 4UL
#define LINEBUFSIZE (MAX_READLEN+3UL)

// %-PURPOSE-%
#define PURPOSE "Divide FastQ file into length-based buckets.\n"\
  "Each bucket contains all reads with length in the range [a,b]\n"\
  "defined by the value e such that floor(e * a) = floor(e * b)\n"\
  "and the interval is maximal."

FILE* xfopen(const char *fname, const char *mode)
{
  FILE* file;
  assert(fname != NULL && mode != NULL);
  file = fopen(fname, mode);
  if (file == NULL)
  {
    fprintf(stderr, "Error by opening '%s', %s\n",
            fname, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

void xfclose(FILE** file)
{
  assert(file != NULL);
  assert(*file != NULL);
  if (fclose(*file)) {
    perror("Error by closing file");
    exit(EXIT_FAILURE);
  }
  *file = NULL;
}

#define XMALLOC(PTR,BUFSIZE) \
  (PTR) = malloc(BUFSIZE); \
  if ((PTR) == NULL) { \
    fprintf(stderr, "Error allocating %lu bytes", (BUFSIZE)); \
    exit(EXIT_FAILURE); \
  }

typedef struct {
  char *filename;
  FILE *file;
  unsigned long nlines, from, to;
} Bucket;

void bucket_init(Bucket *bucket, const char *dirname, const char *basename,
    unsigned long from, unsigned long to)
{
  bucket->from = from;
  bucket->to = to;
  bucket->nlines = 0;
  XMALLOC(bucket->filename, strlen(dirname) + strlen(basename)
      + (2 * LOG10_MAX_READLEN) + 4UL);
  sprintf(bucket->filename, "%s/%lu.%lu.%s", dirname, from, to, basename);
}

void bucket_puts(char *buffer, Bucket *bucket)
{
  assert(bucket->filename != NULL);
  if (bucket->file == NULL)
    bucket->file = xfopen(bucket->filename, "w");
  fputs(buffer, bucket->file);
  bucket->nlines++;
}

void bucket_delete(Bucket *bucket)
{
  if (bucket->file != NULL)
  {
    xfclose(&(bucket->file));
    fprintf(stderr, "%lu\t%lu\t%lu\n", bucket->from, bucket->to,
        bucket->nlines / 4);
  }
  free(bucket->filename);
}

typedef enum {
  sdesc,
  seq,
  qdesc,
  qual
} fqState;

void parse_fastq(FILE *fq, Bucket *buckets, unsigned long *bnums)
{
  char *buf, *current_sdesc;
  fqState nextstate = sdesc;
  unsigned long readlen, lineno = 0;
  XMALLOC(buf,LINEBUFSIZE+(size_t)1);
  XMALLOC(current_sdesc, LINEBUFSIZE+(size_t)1);
  while ((fgets(buf, LINEBUFSIZE, fq)) != NULL)
  {
    lineno++;
    if (strlen(buf) >= LINEBUFSIZE)
    {
      fprintf(stderr, "line %lu: %lu long, more than %lu,"
          " increase MAX_READLEN\n",
          lineno, LINEBUFSIZE - 2UL, strlen(buf));
      exit(1);
    }
    switch (nextstate) {
      case sdesc:
        if (buf[0] == '@')
          strcpy(current_sdesc, buf);
        else
          fprintf(stderr, "line %lu: sdesc format error\n", lineno), exit(1);
        nextstate = seq;
        continue;
      case seq:
        readlen = strlen(buf) - 1UL;
        if (readlen > MAX_READLEN)
          fprintf(stderr, "line %lu: read too long, increase MAX_READLEN\n",
              lineno), exit(1);
        if (readlen == 0)
          fprintf(stderr, "line %lu: read is empty\n", lineno), exit(1);
        bucket_puts(current_sdesc, &buckets[bnums[readlen]]);
        nextstate = qdesc;
        break;
      case qdesc:
        if (buf[0] != '+')
          fprintf(stderr, "line %lu: qdesc format error\n", lineno), exit(1);
        nextstate = qual;
        break;
      case qual:
        if (strlen(buf) - 1UL != readlen)
          fprintf(stderr, "line %lu: quality too short\n", lineno), exit(1);
        nextstate = sdesc;
        break;
    }
    bucket_puts(buf, &buckets[bnums[readlen]]);
  }
  free(buf);
  free(current_sdesc);
}

void split_filename(const char *filename, char **dirname, char **basename)
{
  char *lastslash = strrchr(filename, '/');
  if (lastslash == NULL)
  {
    *basename = strdup(filename);
    *dirname = malloc(2);
    strcpy(*dirname, ".");
  }
  else
  {
    if (lastslash == (filename + strlen(filename) - 1))
    {
      fprintf(stderr, "input filename cannot be a directory (%s)\n",
          filename);
      exit(EXIT_FAILURE);
    }
    *basename = strdup(lastslash+1);
    *dirname = strdup(filename);
    (*dirname)[lastslash - filename] = '\0';
  }
}

void parse_args(int argc, char *argv[], char **if_name, double *e)
{
  size_t if_name_len;
  if ((argc == 2) && (strcmp(argv[1],"-v") == 0))
  {
    printf("%s v."VERSION"\n", argv[0]);
    exit(EXIT_SUCCESS);
  }
  if (argc != 3)
  {
    fprintf(stderr, PURPOSE"\n\n"
                    "  Usage: %s <inputfilename> <e>\n"
                    "    (show program version: %s -v)\n\n"
                    "  Output: (for each range) <a>.<b>.<inputfilename>\n",
                    argv[0], argv[0]);
    exit(EXIT_FAILURE);
  }
  *if_name = argv[1];
  *e = atof(argv[2]);
  if (*e < 0.0 || *e > 1.0)
  {
    fprintf(stderr, "e must be between 0 and 1\n");
    exit(EXIT_FAILURE);
  }
}

void prepare_buckets(double e, const char *dirname, const char *basename,
    unsigned long **bnums, unsigned long *nof_buckets, Bucket **buckets)
{
  unsigned long len, bstart, bnum;
  XMALLOC(*bnums, sizeof (**bnums) * MAX_READLEN);
  for (len = 0; len <= MAX_READLEN; len++)
  {
    (*bnums)[len] = floor(e * len);
  }
  (*nof_buckets) = floor(e * MAX_READLEN) + 1UL;
  XMALLOC(*buckets, sizeof(Bucket) * (*nof_buckets));
  bnum = 0;
  // we start at 1 as reads of length 0 are not allowed;
  // this may create a bucket 0 [1, 0] (if e=1) but it will remain empty anyway
  bstart = 1;
  for (len = 1; len <= MAX_READLEN; len++)
  {
    if ((*bnums)[len] > bnum)
    {
      bucket_init(&(*buckets)[bnum], dirname, basename, bstart, len - 1UL);
      bnum = (*bnums)[len];
      bstart = len;
    }
  }
  bnum = (*bnums)[MAX_READLEN];
  bucket_init(&(*buckets)[bnum], dirname, basename, bstart, MAX_READLEN);
}

int main(int argc, char *argv[])
{
  FILE *infile;
  Bucket *outfiles;
  double e;
  char *if_name, *if_dirname, *if_basename;
  unsigned long *bnums, nof_buckets, i;

  parse_args(argc, argv, &if_name, &e);

  split_filename(if_name, &if_dirname, &if_basename);

  prepare_buckets(e, if_dirname, if_basename, &bnums, &nof_buckets, &outfiles);

  infile = xfopen(if_name,"r");

  parse_fastq(infile, outfiles, bnums);

  xfclose(&infile);
  fprintf(stderr, "# minlen\tmaxlen\tn_reads\n");
  for (i = 0; i < nof_buckets; i++)
    bucket_delete(&outfiles[i]);
  free(outfiles);
  free(bnums);
  free(if_dirname);
  free(if_basename);

  return EXIT_SUCCESS;
}
