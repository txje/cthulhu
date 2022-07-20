#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "../incl/minimap2/minimap.h"
//#include "../incl/minimap2/mmpriv.h"
#include "../incl/klib/kseq.h"
#include "../incl/klib/kvec.h"
#include "../incl/klib/khash.h"
#include "taxonomy.h"
#include <getopt.h>
#include <string.h>
#include "paf.h"

#ifndef M_E
#define M_E 2.718281828459045
#endif

KSEQ_INIT(gzFile, gzread)

typedef struct seqcov {
  uint8_t* cov;
  uint32_t n;
} seqcov;

KHASH_MAP_INIT_STR(ref2cov, seqcov);

typedef struct taxa_hit {
  size_t taxid;
  int score;
} taxa_hit;

// creates string (read name):taxid/score struct
KHASH_MAP_INIT_STR(read2tax, taxa_hit);

typedef struct coverage {
  size_t taxid;
  uint32_t n_loci;
  uint32_t total_coverage;
  uint32_t covered_loci;
} coverage;

// creates string (assembly id):coverage hash
KHASH_MAP_INIT_STR(asm2cov, coverage*);

// create simple taxid : coverage to map taxid to the best accession coverage (shared with asm2cov)
KHASH_MAP_INIT_INT(tax2cov, coverage*)

void version() {
  printf("cthulhu version 0.1\n");
}

void usage() {
  printf("Usage: cthulhu [options]\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz] or precomputed index file\n");
  printf("  --paf: PAF alignment file, or '-' for stdin, instead of reads and reference\n");
  printf("  -d: Directory with *.dmp taxonomy files (NCBI taxonomy)\n");
  printf("  -t: Threads (default: 1)\n");
  printf("  -m: Maximum memory target (GB) (default: 32)\n");
  printf("  -s: Summary output file\n");
  printf("  -o: Read classification output\n");
  printf("  --alignment-output: Raw read alignment output\n");
  printf("  -i: Save index file\n");
  printf("  -p: Sequence type preset\n");
  printf("      map-ont: Oxford Nanopore (default)\n");
  printf("      map-pb:  Pacbio\n");
  printf("      sr:      Short reads (Illumina)\n");
  printf("  -f, --align-fraction: Portion of a read that must align properly (defaults to --align-length threshold)\n");
  printf("  -l, --align-length: Minimum aligned bp (default: 100)\n");
  printf("  -a, --align-accuracy: Minimum accuracy of aligned portion of a read (default: 0.6)\n");
  printf("  -c, --careful: Compute more exact alignments and coverage profiles, increases time and memory\n");
  printf("  -b, --best: Use only the best (or tied for the best) alignments for each read\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -h, --help: show this\n");
  printf("  --version: show version information\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "align-fraction",         required_argument, 0, 'f' },
  { "align-length",           required_argument, 0, 'l' },
  { "align-accuracy",         required_argument, 0, 'a' },
  { "alignment-output",       required_argument, 0, 0  },
  { "careful",                no_argument,       0, 'c' },
  { "verbose",                no_argument,       0, 'v' },
  { "help",                   no_argument,       0, 'h' },
  { "best",                   no_argument,       0, 'b' },
  { "paf",                    required_argument, 0, 0 },
  { "version",                no_argument,       0, 0 },
  { 0, 0, 0, 0}
};

void track_align(paf_rec_t *p, int verbose, float align_fraction, int align_length, float align_accuracy, int best, int covg_bin_size, khash_t(acc2asm) *a2a, taxonomy *tax, khash_t(ref2cov) *r2c, khash_t(read2tax) *read_taxa) {
  int i;
  size_t taxid;

  int aln_len = p->qe - p->qs;
  float aln_frac = (float)aln_len / p->ql;
  float accuracy = (float)p->ml / p->bl;
  if(verbose) {
    fprintf(stderr, "    matches ref %s (%d-%d)\n", p->tn, p->ts, p->te);
  }

  if(align_fraction != -1 ? (aln_len >= align_length && aln_frac >= align_fraction) : (aln_len >= align_length) && accuracy >= align_accuracy) {

    khint_t bin, bin2;

    bin = kh_get(acc2asm, a2a, p->tn);
    int absent = (bin == kh_end(a2a)); 
    if(absent) {
      fprintf(stderr, "Target/accession ID '%s' not found in acc2asm, will be assigned to 1 (root)\n", p->tn);
      taxid = 1;
    } else {
      taxid = kh_val(a2a, bin).taxid;
      if(verbose) {
        fprintf(stderr, "      matches taxon %d (%s)\n", taxid, tax->names[taxid]);
      }
    }

    // add ref:tax to r2c hash (if this ref has already been hit, it already exists, but we don't care)
    //fprintf(stderr, "get ref2cov for %s\n", p->tn);
    bin = kh_get(ref2cov, r2c, p->tn);
    absent = (bin == kh_end(r2c)); 
    if(absent) {
      char *key = malloc(strlen(p->tn)+1);
      strcpy(key, p->tn);
      bin = kh_put(ref2cov, r2c, key, &absent);
      kh_val(r2c, bin).n = 0;
    }

    // add to the coverage for any hit meeting our minimum alignment threshold
    if(kh_val(r2c, bin).n == 0) { // coverage array has never been initialized
      kh_val(r2c, bin).n = p->tl / covg_bin_size + 1;
      if(verbose) {
        fprintf(stderr, "Making new coverage array of length %u (%u / %u) for ref %s\n", kh_val(r2c, bin).n, p->tl, covg_bin_size, p->tn);
      }
      kh_val(r2c, bin).cov = calloc(sizeof(uint8_t), kh_val(r2c, bin).n);
    }
    //for(i = p->ts/covg_bin_size; i <= p->te/covg_bin_size; i++) {
    i = p->ts/covg_bin_size;
      if(verbose) {
        fprintf(stderr, "Incrementing coverage in bin %d (of %u)\n", i, kh_val(r2c, bin).n);
      }
      if(kh_val(r2c, bin).cov[i] < 255) {
        kh_val(r2c, bin).cov[i] += 1;
      }
    //}

    //fprintf(stderr, "adding read to tax\n");
    bin2 = kh_put(read2tax, read_taxa, p->qn, &absent);
    if(absent) {
      kh_key(read_taxa, bin2) = malloc(strlen(p->qn)+1); // string has to be duplicated because the query name will be freed or reassigned as soon as we're done processing this alignment
      strcpy(kh_key(read_taxa, bin2), p->qn);
      kh_val(read_taxa, bin2).taxid = 0;
      kh_val(read_taxa, bin2).score = 0;
    }

    if(kh_val(read_taxa, bin2).taxid == 0 || (best && p->ml > kh_val(read_taxa, bin2).score)) {
      kh_val(read_taxa, bin2).taxid = taxid;
      kh_val(read_taxa, bin2).score = p->ml; // total matched loci
    } else if(!best || p->ml == kh_val(read_taxa, bin2).score) {
      //fprintf(stderr, "lca\n");
      kh_val(read_taxa, bin2).taxid = lca(taxid, kh_val(read_taxa, bin2).taxid, tax);
    }
  }
}

int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;

  char* read_fasta = NULL;
  char* ref_fasta = NULL;
  char* paf_file = NULL;
  char* tax_dir = NULL;
  char* preset = "map-ont";
  char* out_file = NULL;
  char* summary_file = NULL;
  char* idx_file = NULL;
  char* align_file = NULL;
  float align_fraction = -1;
  float align_accuracy = 0.6;
  int align_length = 100;
  int n_threads = 1;
  int verbose = 0;
  int careful = 0;
  int best = 0;
  int target_memory_gb = 32;
  size_t covg_bin_size = 1000; // bin size in coverage arrays

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:d:t:f:a:p:s:o:m:i:l:bvh", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'q':
        read_fasta = optarg;
        break;
      case 'r':
        ref_fasta = optarg;
        break;
      case 'd':
        tax_dir = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        break;
      case 'p':
        preset = optarg;
        break;
      case 'm':
        target_memory_gb = atoi(optarg);
        break;
      case 'f':
        align_fraction = atof(optarg);
        break;
      case 'a':
        align_accuracy = atof(optarg);
        break;
      case 's':
        summary_file = optarg;
        break;
      case 'o':
        out_file = optarg;
        break;
      case 'i':
        idx_file = optarg;
        break;
      case 'l':
        align_length = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'b':
        best = 1;
        break;
      case 'c':
        careful = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'q' || optopt == 'r' || optopt == 'd' || optopt == 't' || optopt == 'f' || optopt == 'a' || optopt == 'p' || optopt == 'm' || optopt == 'o' || optopt == 's' || optopt == 'i' || optopt == 'l')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        // as long as all the long arguments have characters too, I don't think this section will be used
        if (long_idx == 0) align_fraction = atof(optarg); // --align-fraction
        else if (long_idx == 1) align_length = 100; // --align-length
        else if (long_idx == 2) align_accuracy = atof(optarg); // --align-accuracy
        else if (long_idx == 3) align_file = optarg; // --alignment-output
        else if (long_idx == 4) careful = 1; // --careful
        else if (long_idx == 5) verbose = 1; // --verbose
        else if (long_idx == 6) {usage(); return 0;} // --help
        else if (long_idx == 7) best = 1; // --best
        else if (long_idx == 8) paf_file = optarg; // --paf
        else if (long_idx == 9) {version(); return 0;} // --version
        break;
      default:
        usage();
        return 1;
    }
  }

  if(tax_dir == NULL) {
    fprintf(stderr, "-d taxonomy directory is required\n");
    return 1;
  }
  if(read_fasta == NULL && paf_file == NULL) {
    fprintf(stderr, "-q reads FASTA is required (or --paf)\n");
    return 1;
  }
  if(ref_fasta == NULL && paf_file == NULL) {
    fprintf(stderr, "-r reference FASTA is required (or --paf)\n");
    return 1;
  }
  if(out_file == NULL && summary_file == NULL) {
    fprintf(stderr, "Specify -o, -s, or both, otherwise we're doing all this for nothing\n");
    return 1;
  }
  if(paf_file != NULL && strcmp(paf_file, "-") != 0 && align_file != NULL) {
    fprintf(stderr, "--paf and --alignment-output specified, but --paf is not stdin, so we're ignoring --alignment-output\n");
    align_file = NULL;
  }

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt); // initialize with defaults
  mm_set_opt(preset, &iopt, &mopt); // then add ont presets
  // this isn't great, but it's close-ish
  /*
  if(target_memory_gb < 10) {
    iopt.batch_size = target_memory_gb * 100000000ULL;
  } else if(target_memory_gb < 50) {
    iopt.batch_size = target_memory_gb * 200000000ULL;
  } else if(target_memory_gb < 100) {
    iopt.batch_size = target_memory_gb * 250000000ULL;
  } else if(target_memory_gb < 200) {
    iopt.batch_size = target_memory_gb * 300000000ULL;
  } else {
    iopt.batch_size = target_memory_gb * 400000000ULL;
  }
  */
  fprintf(stderr, "Target memory usage: %d GB\n", target_memory_gb);
  fprintf(stderr, "Using batch size: %d Gbp [expect ~(40 / batch size) Gbp for Refseq]\n", iopt.batch_size/1000000000);
  if(careful)
    mopt.flag |= MM_F_CIGAR; // perform alignment

  mopt.flag |= MM_F_NO_PRINT_2ND; // skip all secondary alignments
  iopt.flag = 0;
  //iopt.k = 21; // this ignores the requested preset probably, since map-ont only sets k anyway, but could substantially increase speed and decrease sensitivity

  char* name_f = malloc((strlen(tax_dir)+11) * sizeof(char));
  strcpy(name_f, tax_dir);
  strcat(name_f, "/names.dmp");

  char* node_f = malloc((strlen(tax_dir)+11) * sizeof(char));
  strcpy(node_f, tax_dir);
  strcat(node_f, "/nodes.dmp");

  taxonomy* tax = read_taxonomy(name_f, node_f);

  char* acc2asm_f = malloc((strlen(tax_dir)+42) * sizeof(char));
  strcpy(acc2asm_f, tax_dir);
  //strcat(acc2asm_f, "/nucl_gb.accession2taxid.filtered");
  strcat(acc2asm_f, "/refseq/accession_map.txt");

  khash_t(acc2asm) *a2a = kh_init(acc2asm); // accession (NC_**) to assembly (GCF_**)
  int ret = parse_acc2tax(acc2asm_f, a2a);
  if(ret != 0) {
    return ret;
  }
  fprintf(stderr, "Parsed taxonomy files.\n");

  khash_t(ref2cov) *r2c = kh_init(ref2cov);

  // open query file for reading; you may use your favorite FASTA/Q parser
  gzFile f;
  kseq_t *ks;

  FILE* af;
  if(align_file != NULL) {
    af = fopen(align_file, "w");
  }

  // generic hashing variables
  khint_t bin, bin2; // hash bin (result of kh_put/get)
  int absent;

  // taxonomic accounting stuff
  khash_t(read2tax) *read_taxa = kh_init(read2tax);

  if(ref_fasta != NULL && read_fasta != NULL) {

    // open index reader
    fprintf(stderr, "Building mm2 index...\n");
    mm_idx_reader_t *r = mm_idx_reader_open(ref_fasta, &iopt, idx_file);
    if(!r) {
      fprintf(stderr, "Cannot open '%s'\n", ref_fasta);
      return 1;
    }
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
      // open (or re-open) the query file -- needs to be re-read through for each part of the index
      f = gzopen(read_fasta, "r");
      if(!f) {
        fprintf(stderr, "Cannot open '%s'\n", read_fasta);
        return 1;
      }
      ks = kseq_init(f); 

      fprintf(stderr, "Processing mm2 index (or fraction thereof)...\n");
      mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
      mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
      int n = 0;
      while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        mm_reg1_t *reg;
        int j, i, n_reg;
        if(verbose) {
          fprintf(stderr, "Processing read %d (%s, %u bp): %s\n", n, ks->name.s, ks->seq.l, ks->seq.s);
        }
        reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
        if(verbose) {
          fprintf(stderr, "  %d raw alignments\n", n_reg);
        }
        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
          mm_reg1_t *r = &reg[j];
          if(careful)
            assert(r->p); // with MM_F_CIGAR, this should not be NULL

          paf_rec_t p;
          p.qn = ks->name.s;
          p.tn = mi->seq[r->rid].name;
          p.qs = r->qs;
          p.qe = r->qe;
          p.ql = ks->seq.l;
          p.ts = r->rs;
          p.te = r->re;
          p.tl = mi->seq[r->rid].len;
          p.ml = r->mlen;
          p.bl = r->blen;
          p.rev = r->rev;
          track_align(&p, verbose, align_fraction, align_length, align_accuracy, best, covg_bin_size, a2a, tax, r2c, read_taxa);

          //printf("length %d, accuracy %f\n", aln_len, accuracy);
          if(align_file != NULL) {
            fprintf(af, "%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
            fprintf(af, "%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
            if(careful) {
              for (i = 0; i < r->p->n_cigar; ++i) { // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                fprintf(af, "%d%c", r->p->cigar[i]>>4, "MIDSHN"[r->p->cigar[i]&0xf]);
              }
            }
            fprintf(af, "\n");
          }

          //fprintf(stderr, "freeing r->p\n");
          free(r->p);
        }
        n++;
        //fprintf(stderr, "freeing reg\n");
        free(reg);
      }
      //fprintf(stderr, "%d reads processed\n", n);
      //fprintf(stderr, "mm_tbuf_destroy\n");
      mm_tbuf_destroy(tbuf);
      //fprintf(stderr, "mm_idx_destroy\n");
      mm_idx_destroy(mi);
      //fprintf(stderr, "kseq_destroy\n");
      kseq_destroy(ks); // close the query file
      //fprintf(stderr, "gzclose\n");
      gzclose(f);
    }
    //fprintf(stderr, "mm_idx_reader_close\n");
    mm_idx_reader_close(r); // close the index reader
  }
  else { // paf must have been given
    fprintf(stderr, "Reading from alignment file '%s'\n", paf_file);
    paf_file_t *p = paf_open(paf_file);
    if(!p) {
      fprintf(stderr, "Cannot open '%s'\n", paf_file);
      return 1;
    }
    paf_rec_t r;
    int ret = paf_read(p, &r);
    uint64_t n = 0;
    while(ret == 0) {
      track_align(&r, verbose, align_fraction, align_length, align_accuracy, best, covg_bin_size, a2a, tax, r2c, read_taxa);
      if(align_file != NULL) {
        fprintf(af, "%s\t%d\t%d\t%d\t%c\t", r.qn, r.ql, r.qs, r.qe, "+-"[r.rev]);
        fprintf(af, "%s\t%d\t%d\t%d\t%d\t%d\t%d", r.tn, r.tl, r.ts, r.te, r.ml, r.bl, r.mq);
        // cannot currently print CIGAR string in cg:Z tag because we can't guarantee the PAF has it and we aren't parsing it
        fprintf(af, "\n");
      }
      ret = paf_read(p, &r);
      /*
      if((++n) % 100000 == 0) {
        fprintf(stderr, "---- %u reads processed.\n", n);
      }
      */
    }
    paf_close(p);
    fprintf(stderr, "Closing file '%s'\n", paf_file);
  }

  // count reads per taxa
  taxtree *tree = new_tree();
  int no_hit = 0;
  FILE* o = out_file != NULL ? fopen(out_file, "w") : (FILE*)NULL;
  for (bin = 0; bin < kh_end(read_taxa); ++bin) {
    if (kh_exist(read_taxa, bin)) {
      if(kh_val(read_taxa, bin).taxid == 0) {
        no_hit++;
      } else {
        // output taxa result for this read
        if(out_file != NULL) {
          fprintf(o, "%s\t%d\t%s\t%d\n", kh_key(read_taxa, bin), kh_val(read_taxa, bin).taxid, tax->names[kh_val(read_taxa, bin).taxid], kh_val(read_taxa, bin).score);
        }

        if(summary_file != NULL) {
          add_to_tree(tax, tree, kh_val(read_taxa, bin).taxid);
        }
      }
    }
  }
  if(summary_file != NULL) {
    int j;
    size_t taxid;

    // collate coverage by accession
    khash_t(tax2cov) *t2c = kh_init(tax2cov); // taxon ID to *coverage (shared by a2c)
    coverage* cv;
    khash_t(asm2cov) *a2c = kh_init(asm2cov);
    for (bin = 0; bin < kh_end(r2c); ++bin) {
      if (kh_exist(r2c, bin)) {
        // get ref (accession) name
        char* ref_name = kh_key(r2c, bin);

        // get tax ID
        bin2 = kh_get(acc2asm, a2a, ref_name);
        if(bin2 == kh_end(a2a)) continue; // absent, this ref_name was not in the accession lookup table, so we won't count its coverage
        taxid = kh_val(a2a, bin2).taxid;
        char* as = kh_val(a2a, bin2).assembly; // assembly ID (GCF...)

        uint32_t n = kh_val(r2c, bin).n;
        uint8_t* c = kh_val(r2c, bin).cov;

        if(n > 0) {
          bin2 = kh_put(asm2cov, a2c, as, &absent); // reuses the value in r2c
          //fprintf(stderr, "asm2cov bin %u\n", bin2);
          if(absent) {
            //fprintf(stderr, "Making new coverage count for taxa %u\n", taxid);
            cv = calloc(1, sizeof(coverage));
            cv->taxid = taxid;
            kh_val(a2c, bin2) = cv;
          } else {
            cv = kh_val(a2c, bin2);
          }
          //fprintf(stderr, "ref %s -> tax %u has %u loci\n", ref_name, taxid, n);
          //fprintf(stderr, "before tax %u has %u loci, %u covered, %u total covg\n", taxid, cv->n_loci, cv->covered_loci, cv->total_coverage);
          // if no reads mapped *specifically* to this taxa, n is simply 0
          cv->n_loci += n;
          for(j = 0; j < n; j++) {
            if(c[j] > 0) {
              cv->total_coverage += c[j];
              cv->covered_loci += 1;
            }
          }

          // keep a pointer to the best coverage profile (ordered by covered_loci) in tax2cov (t2c)
          if(verbose) {
            fprintf(stderr, "adding accession %s to taxid %u\n", ref_name, taxid);
          }
          khint_t tbin = kh_put(tax2cov, t2c, taxid, &absent);
          if(absent || cv->covered_loci > kh_val(t2c, tbin)->covered_loci) {
            kh_val(t2c, tbin) = cv;
          }

          //fprintf(stderr, "after tax %u has %u loci, %u covered, %u total covg\n", taxid, cv->n_loci, cv->covered_loci, cv->total_coverage);
        }
      }
    }

    FILE* sf = fopen(summary_file, "w");
    fprintf(sf, "phylum_taxid\tclass_taxid\torder_taxid\tfamily_taxid\tgenus_taxid\tspecies_taxid\tsubspecies_taxid\t");
    fprintf(sf, "taxid\ttaxname\tread_count\tunique_read_count\tgenome_size\ttotal_coverage\tcovered_loci\tcovered_frac\tgenome_coverage\tcovered_coverage\texpect_covered\tcovered/expected\n");
    //fprintf(sf, "genome\ttaxid\ttaxname\tread_count\tunique_read_count\tgenome_size\ttotal_coverage\tcovered_loci\tcovered_frac\tgenome_coverage\tcovered_coverage\texpect_covered\tcovered/expected\n");
    fprintf(sf, "0\tno hit\t%d\t%d\t0\t0\t0\t0\t0\t0\t0\t0\n", no_hit, no_hit);
    //fprintf(sf, "no hit\t0\tno hit\t%d\t%d\t0\t0\t0\t0\t0\t0\t0\t0\n", no_hit, no_hit);

    // output single taxa counts and build full hierarchal tree
    for (bin = 0; bin < kh_end(tree); ++bin) {
      if (kh_exist(tree, bin)) {
        taxid = kh_key(tree, bin);
        //fprintf(sf, "%d\t%d\t%d\t%s\n", taxid, kh_val(tree, bin).count, kh_val(tree, bin).unique_count, tax->names[taxid]);

    //for (bin = 0; bin < kh_end(t2c); ++bin) {
    //  if (kh_exist(t2c, bin)) {
        //char* as = kh_key(t2c, bin);

        /*
        cv = kh_val(t2c, bin);
        taxid = cv->taxid;
        bin2 = kh_get(nodehash, tree, taxid);
        if(bin2 == kh_end(tree)) continue; // does not report coverage, although there may be some recorded, for taxa that did not have any reads ultimately assigned
        */

        size_t* tax_hierarchy = get_hierarchy(taxid, tax);
        fprintf(sf, "%u\t%u\t%u\t%u\t%u\t%u\t%u\t", tax_hierarchy[0], tax_hierarchy[1], tax_hierarchy[2], tax_hierarchy[3], tax_hierarchy[4], tax_hierarchy[5], tax_hierarchy[6]);

        bin2 = kh_get(tax2cov, t2c, taxid);
        if(bin2 == kh_end(t2c)) { // no coverage - may be a higher-level taxa
          fprintf(sf, "%u\t%s\t%u\t%u\t0\t0\t0\t0\t0\t0\t0\t0\n", taxid, tax->names[taxid], kh_val(tree, bin).count, kh_val(tree, bin).unique_count);
        } else {
          cv = kh_val(t2c, bin2);
          uint32_t expect_covered = (uint32_t)round(cv->n_loci - pow(M_E, cv->total_coverage * log(cv->n_loci - 1) - (cv->total_coverage - 1) * log(cv->n_loci)));
          fprintf(sf, "%u\t%s\t%u\t%u\t%u\t%u\t%u\t%f\t%f\t%f\t%u\t%f\n", taxid, tax->names[taxid], kh_val(tree, bin).count, kh_val(tree, bin).unique_count, cv->n_loci * covg_bin_size, cv->total_coverage * covg_bin_size, cv->covered_loci * covg_bin_size, (float)cv->covered_loci/cv->n_loci, (float)cv->total_coverage/cv->n_loci, (float)cv->total_coverage/cv->covered_loci, expect_covered * covg_bin_size, (float)cv->covered_loci / expect_covered);
        }
      }
    }
    //fprintf(sf, "\n");
    //}

    //fprintf(sf, "\nTree:\n");
    //depth_first_traverse(tax, tree, 1, 0, sf); // do a depth-first tree render starting at the root
  }

  // --- clean up memory ---

  // clean up a2a (acc2asm)
  for (bin = 0; bin < kh_end(a2a); ++bin) {
    if (kh_exist(a2a, bin)) {
      free((char*)kh_key(a2a, bin));
      free((char*)kh_val(a2a, bin).assembly);
    }
  }
  kh_destroy(acc2asm, a2a);
  kh_destroy(nodehash, tree);
  for (bin = 0; bin < kh_end(read_taxa); ++bin) {
    if (kh_exist(read_taxa, bin)) {
      free((char*)kh_key(read_taxa, bin));
    }
  }
  kh_destroy(read2tax, read_taxa);
  free(name_f);
  free(node_f);
  free(acc2asm_f);
  free_tax(tax);

  return 0;
}

/*

^(;,;)^

 */
