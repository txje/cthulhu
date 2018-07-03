#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../incl/htslib/htslib/sam.h"
#include "../incl/klib/khash.h"
#include "../incl/klib/ksort.h"
#include "../incl/klib/kseq.h"

/*
 * covstat.c
 *
 * Jeremy Wang
 * 20180628
 *
 * Compute coverage over named sequences (maybe species/chromosomes)
*/

// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(FILE*, fileread);

// creates string:[array of uint8] hash
KHASH_MAP_INIT_STR(refSeq, char*);

// creates string:uint32 hash
KHASH_MAP_INIT_STR(refLen, uint32_t);

// creates string:string hash to map accessions to names
KHASH_MAP_INIT_STR(acc2name, char*);

typedef struct refstats {
  uint32_t total_loci;
  uint32_t covered_loci;
  uint32_t total_coverage;
} refstats;

// creates string:refstats hash to names to a stat struct
KHASH_MAP_INIT_STR(name2stats, refstats);


int main(int argc, char *argv[]) {

  if(argc < 4) {
    fprintf(stderr, "Usage: covstat <BAM> <reference FASTA> <accession map>\n");
    fprintf(stderr, "Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *ref_fasta = argv[2];

  khint_t bin; // hash bin (result of kh_put)
  int l, absent;
  FILE *fp;

  fp = fopen(argv[3], "r");
  if(fp == NULL) {
    fprintf(stderr, "Error reading accession mapping file '%s'\n", argv[3]);
  }
  char line[1024]; // maximum line size is 1024 chars
  const char* delim = "\t";
  char **parts;
  char *ln;
  khash_t(acc2name) *a2n = kh_init(acc2name);
  khash_t(name2stats) *n2s = kh_init(name2stats);
  refstats rs;
  int i;
  char* token;
  while(fgets(line, sizeof line, fp) != NULL) {
    parts = malloc(sizeof(char*) * 2); // only two fields separated by a tab
    i = 0;
    token = strtok(line, delim);
    if(token != NULL)
      parts[i] = strdup(token);
    while(token != NULL) {
      i++;
      token = strtok(NULL, delim);
      if(token != NULL) {
        parts[i] = strdup(token);
        parts[i][strlen(parts[i])-1] = NULL;
      }
    }
    bin = kh_put(acc2name, a2n, parts[0], &absent);
    kh_val(a2n, bin) = parts[1];

    // initialized name2stats hash too
    bin = kh_put(name2stats, n2s, strdup(parts[1]), &absent);
    if(absent) {
      kh_val(n2s, bin).total_loci = 0;
      kh_val(n2s, bin).covered_loci = 0;
      kh_val(n2s, bin).total_coverage = 0;
    }
  }
  fclose(fp);

  // load FASTA file

  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);

  kseq_t *seq, *nextseq;

  fp = fopen(ref_fasta, "r");
  seq = kseq_init(fp);
  fprintf(stderr, "Reading fasta file: %s\n", ref_fasta);

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    // seq array
    bin = kh_put(refSeq, ref, strdup(seq->name.s), &absent);
    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    kh_val(ref, bin) = malloc(sizeof(char)*l);
    memcpy(kh_val(ref, bin), seq->seq.s, sizeof(char)*l);

    // sequence length
    bin = kh_put(refLen, rlen, strdup(seq->name.s), &absent);
    kh_val(rlen, bin) = l;
  }

  fclose(fp);
  kseq_destroy(seq);


  // load BAM file
  //
  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int ret_val;

  bam = sam_open(bam_file, "rb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", bam_file);
    return -1;
  }
  header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return -1;
  }
  // construct array from reference information so that we can look it up with read.tid
  char **ref_array = malloc(sizeof(char*) * header->n_targets);
  uint16_t **covg = malloc(sizeof(uint16_t*) * header->n_targets); // (ref_id x ref_len) array of coverages

  int *rlen_array = malloc(sizeof(uint32_t) * header->n_targets); // has to be a plain int because that's what kseq gives out
  //printf("BAM targets:\n");
  int j;
  for (i = 0; i < header->n_targets; i++) {
    covg[i] = calloc(header->target_len[i], sizeof(uint16_t));
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);
    bin = kh_get(refSeq, ref, header->target_name[i]);
    ref_array[i] = kh_value(ref, bin);
    bin = kh_get(refLen, rlen, header->target_name[i]);
    int fa_len = kh_value(rlen, bin);
    if (fa_len != header->target_len[i]) { // target_len is a uint32_t
      fprintf(stderr, "WARNING: Reference fasta length (%i) and BAM length (%u) of %s do not agree.\n", fa_len, header->target_len[i], header->target_name[i]);
    }
    rlen_array[i] = fa_len;
  }

  aln = bam_init1();

  int32_t read_len = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {

    ct++;
    if(ct % 100000 == 0) {
      //printf("%d reads processed\n", ct);
    }

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;
    int32_t* cigar = bam_get_cigar(aln); // lower 4 bits are cigar operation, upper 28 are the length
    uint8_t* qseq = bam_get_seq(aln); // 4 bits each (1: A, 2: C, 4: G, 8: T, 15: N)
    int32_t qpos = 0;
    int32_t o;
    uint8_t rev = bam_is_rev(aln);
    //printf("query of len %d aligned to target %d (%s) at pos %d (%c), cigar: ", qlen, tid, header->target_name[tid], pos, (rev ? '-' : '+'));

    for(o = 0; o < aln->core.n_cigar; o++) {
      char opc = bam_cigar_opchr(cigar[o]);
      int oplen = bam_cigar_oplen(cigar[o]);
      //printf("%c%d, ", opc, oplen);
      int32_t qstartpos = qpos;
      int32_t startpos = pos;
      if(opc == 'M') {
        while(qpos < qstartpos + oplen) {
          if(covg[tid][pos] < (1<<16)-1) {
            covg[tid][pos]++;
          }
          pos++;
          qpos++;
        }
      } else if(opc == 'I' || opc == 'S') { // insertion or softclip both consume the query only
        qpos = qpos + oplen;
      } else if(opc == 'D') {
        pos = pos + oplen;
      } else if(opc == 'H') { // hard clip, does nothing...
        continue;
      } else {
        fprintf(stderr, "WARNING: Unhandled CIGAR op character '%c' at target %d, position %d\n", opc, tid, pos);
        break;
      }
    }
  }

  printf("reference_name\tgenome_size\taligned_bp\tcovered_loci\tcovered_pct\tavg_coverage\tcoverage_of_covered_loci\n");
  char* name;
  for (i = 0; i < header->n_targets; i++) {
    bin = kh_get(acc2name, a2n, header->target_name[i]);
    // check if absent? could cause segfault if missing
    name = kh_val(a2n, bin);

    bin = kh_get(name2stats, n2s, name);
    kh_val(n2s, bin).total_loci += header->target_len[i];
    for (j = 0; j < header->target_len[i]; j++) {
      if(covg[i][j] > 0) {
        kh_val(n2s, bin).covered_loci++;
      }
      kh_val(n2s, bin).total_coverage += covg[i][j];
    }
  }

  for (bin = 0; bin < kh_end(n2s); ++bin) {
    if (kh_exist(n2s, bin)) {
      name = kh_key(n2s, bin);
      rs = kh_val(n2s, bin);
      if(rs.total_coverage > 0) {
        printf("%s\t%u\t%u\t%u\t%f\t%f\t%f\n", name, rs.total_loci, rs.total_coverage, rs.covered_loci, (float)rs.covered_loci/rs.total_loci, (float)rs.total_coverage/rs.total_loci, (float)rs.total_coverage/rs.covered_loci);
      }
    }
  }

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  bam_destroy1(aln);

  return 0;
}

