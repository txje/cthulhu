#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../incl/minimap2/minimap.h"
#include "../incl/klib/kseq.h"
#include "../incl/klib/kvec.h"
#include "../incl/klib/khash.h"
#include "taxonomy.h"

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt);
  mopt.flag |= MM_F_CIGAR; // perform alignment

  if (argc < 5) {
    fprintf(stderr, "Usage: cattax <query.fa> <target.fa> <directory w/*.dmp> <n_threads>\n");
    return 1;
  }

  int n_threads = atoi(argv[4]);

  char* name_f = malloc((strlen(argv[3])+10) * sizeof(char));
  strcpy(name_f, argv[3]);
  strcat(name_f, "/names.dmp");

  char* node_f = malloc((strlen(argv[3])+10) * sizeof(char));
  strcpy(node_f, argv[3]);
  strcat(node_f, "/nodes.dmp");

  taxonomy* tax = read_taxonomy(name_f, node_f);

  char* acc2tax_f = malloc((strlen(argv[3])+24) * sizeof(char));
  strcpy(acc2tax_f, argv[3]);
  strcat(acc2tax_f, "/nucl_gb.accession2taxid.filtered");

  khash_t(acc2tax) *a2tx = parse_acc2tax(acc2tax_f);
  printf("Parsed taxonomy files.\n");

  // open query file for reading; you may use your favorite FASTA/Q parser
  gzFile f;
  kseq_t *ks;

  // open index reader
  fprintf(stderr, "Building mm2 index...\n");
  mm_idx_reader_t *r = mm_idx_reader_open(argv[2], &iopt, 0);
  mm_idx_t *mi;
  khint_t bin; // hash bin (result of kh_put/get)
  int absent;
  size_t taxid;
  kvec_t(int) read_taxa;
  kv_init(read_taxa);
  while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
    // open (or re-open) the query file -- needs to be re-read through for each part of the index
    f = gzopen(argv[1], "r");
    assert(f);
    ks = kseq_init(f); 

    fprintf(stderr, "Processing mm2 index (or fraction thereof)...\n");
    mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    int n = 0;
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
      if(n == kv_size(read_taxa)) {
        kv_push(int, read_taxa, 0);
      }
      mm_reg1_t *reg;
      int j, i, n_reg;
      reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
      for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
        mm_reg1_t *r = &reg[j];
        assert(r->p); // with MM_F_CIGAR, this should not be NULL

        bin = kh_get(acc2tax, a2tx, mi->seq[r->rid].name);
        absent = (bin == kh_end(a2tx)); 
        if(absent) {
          fprintf(stderr, "Target/accession ID '%s' not found in acc2tax\n", mi->seq[r->rid].name);
        } else {
          taxid = kh_val(a2tx, bin);
          //printf("Matches taxon %d (%s)\n", taxid, tax->names[taxid]);
        }

        int aln_len = r->qe - r->qs;
        float aln_frac = (float)aln_len / ks->seq.l;
        float accuracy = (float)r->mlen / r->blen;

        // arbitrary thresholds right now - these should be parameterized
        if(aln_frac > 0.5 && accuracy > 0.8) {
          if(kv_A(read_taxa, n) == 0)
            kv_A(read_taxa, n) = taxid;
          else
            kv_A(read_taxa, n) = lca(taxid, kv_A(read_taxa, n), tax);

          //printf("length %d, accuracy %f\n", aln_len, accuracy);
          if(0 && taxid != 9606) {
            printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
            printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
            for (i = 0; i < r->p->n_cigar; ++i) { // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
              printf("%d%c", r->p->cigar[i]>>4, "MIDSHN"[r->p->cigar[i]&0xf]);
            }
            putchar('\n');
          }
        }
        free(r->p);
      }
      n++;
      free(reg);
    }
    fprintf(stderr, "%d reads processed\n", n);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
    kseq_destroy(ks); // close the query file
    gzclose(f);
  }
  mm_idx_reader_close(r); // close the index reader

  // clean up a2tx (acc2tax)
  for (bin = 0; bin < kh_end(a2tx); ++bin) {
    if (kh_exist(a2tx, bin))
      free((char*)kh_key(a2tx, bin));
  }
  kh_destroy(acc2tax, a2tx);

  // count reads per taxa
  taxtree *tree = new_tree();
  int i;
  int no_hit = 0;
  for(i = 0; i < kv_size(read_taxa); i++) {
    if(kv_A(read_taxa, i) == 0) {
      no_hit++;
    } else {
      // output taxa result for this read
      //printf("%d\t%d\t%s\n", i, kv_A(read_taxa, i), tax->names[kv_A(read_taxa, i)]);

      add_to_tree(tax, tree, kv_A(read_taxa, i));
    }
  }
  printf("0\t%d\t%d\tno hit\n", no_hit, no_hit);
  // output single taxa counts and build full hierarchal tree
  for (bin = 0; bin < kh_end(tree); ++bin) {
    if (kh_exist(tree, bin)) {
      taxid = kh_key(tree, bin);
      printf("%d\t%d\t%d\t%s\n", taxid, kh_val(tree, bin).count, kh_val(tree, bin).unique_count, tax->names[taxid]);
    }
  }

  printf("\n");
  depth_first_traverse(tax, tree, 1, 0); // do a depth-first tree render starting at the root


  // clean up memory
  kh_destroy(nodehash, tree);
  kv_destroy(read_taxa);
  free(name_f);
  free(node_f);
  free(acc2tax_f);
  free_tax(tax);

  return 0;
}

/*

 /\_/\
(=^x^=)

 */
