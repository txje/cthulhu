/*
 * taxonomy.c
 *
 * Jeremy Wang
 * 20180531
 * MIT License
 *
 * Represent NCBI taxonomy (from taxdump) as a doubly-linked tree
 * with each node directly accessible from an indexed list
 *
 * .dmp files:
 *   lines are delimited by \t|\n
 *   fields are delimited by \t|\t
 *
 * There are currently ~2.2M unique names in the taxonomy database
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "taxonomy.h"

#define MAX_NODES 3000000

taxonomy* read_taxonomy(char* name_f, char* node_f) {
  /*
   * read names.dmp
   */
  name_file_t name_dmp = name_init(name_f);
  char** names = malloc(MAX_NODES * sizeof(char*));
  name_line_t *a = name_read_line(&name_dmp);
  while(a != NULL) {
    // keep only scientific names
    if(strcmp(a->name_class, "scientific name") == 0) {
      char* n;
      if(strcmp(a->unique_name, "") != 0) {
        n = malloc(strlen(a->unique_name)+1);
        strcpy(n, a->unique_name);
      } else {
        n = malloc(strlen(a->name)+1);
        strcpy(n, a->name);
      }
      names[a->taxid] = n;
      //if(a->taxid == 9606)
      //  printf("%d '%s'\n", a->taxid, names[a->taxid]);
    }
    free(a);
    a = name_read_line(&name_dmp);
  }
  free(a);
  name_close(&name_dmp);

  /*
   * read nodes.dmp
   */
  node_file_t node_dmp = node_init(node_f);
  node* nodes = malloc(MAX_NODES * sizeof(node));
  node_line_t *b = node_read_line(&node_dmp);
  while(b != NULL) {
    //printf("%d %d %d\n", b->taxid, b->parent_taxid, b->rank);
    node n = {b->parent_taxid, b->rank};
    nodes[b->taxid] = n;
    free(b);
    b = node_read_line(&node_dmp);
  }
  free(b);
  node_close(&node_dmp);

  taxonomy *tax = malloc(sizeof(taxonomy));
  tax->names = names;
  tax->nodes = nodes;
  return tax;
}

void free_tax(taxonomy* tax) {
  int i;
  for(i = 0; i < MAX_NODES; i++) {
    free(tax->names[i]);
  }
  free(tax->names);
  free(tax->nodes);
  free(tax);
}

/*
 * names.dmp
 *
 * tax_id         -- the id of node associated with this name
 * name_txt       -- name itself
 * unique name    -- the unique variant of this name if name not unique
 * name class     -- (synonym, common name, ...)
 */

name_line_t *name_read_line(name_file_t* dmp) {
  //const char* name_format = "%d|%s|%s|%s|";

  char l[1024]; // maximum line size is 1024 chars
  const char* delim = "|";
  char **parts;
  int i;
  name_line_t *line = malloc(sizeof(name_line_t));
  if(fgets(l, sizeof l, dmp->fp) != NULL) {
    char* ln = malloc(strlen(l)+1);
    strcpy(ln, l);
    parts = malloc(sizeof(char*) * 4); // 4 fields per line
    i = 0;
    parts[i] = strtok(ln, delim);
    while(parts[i] != NULL && strcmp(parts[i], "\n") != 0) {
      i++;
      if(i >= 4) break;
      parts[i] = strtok(NULL, delim);
      if(parts[i][strlen(parts[i])-1] == '\t') {
        parts[i][strlen(parts[i])-1] = '\0';
      }
      if(parts[i][0] == '\t') {
        parts[i]++;
      }
      //printf("field %d: '%s'\n", i, parts[i]);
    }
    if(i != 4) { // line too short or small or -1 for EOF
      if(i == -1) {
        return NULL;
      } else {
        fprintf(stderr, "NAMES.DMP line %d does not appear to be in expected format (has %d fields) - STOPPED READING HERE\n", dmp->cur_row, i);
        return NULL;
      }
    }
    line->taxid = atoi(parts[0]);
    line->name = malloc(strlen(parts[1])+1);
    strcpy(line->name, parts[1]);
    line->unique_name = malloc(strlen(parts[2])+1);
    strcpy(line->unique_name, parts[2]);
    line->name_class = malloc(strlen(parts[3])+1);
    strcpy(line->name_class, parts[3]);
    //for(i = 0; i < 4; i++) free(parts[i]);
    free(parts);
  } else {
    return NULL;
  }
  
  dmp->cur_row++;

  return line;
}

name_file_t name_init(char* f) {
  // open bed file from path
  FILE *dmp_fp = fopen(f, "r");
  if( dmp_fp == NULL ) {
    fprintf(stderr, "Error reading NAMES.DMP file '%s'\n", f);
  }

  name_file_t dmp = {dmp_fp, 0};

  return dmp;
}

int name_close(name_file_t* dmp) {
  fclose(dmp->fp);
}

/*
 * nodes.dmp
 *
 * tax_id                              -- node id in GenBank taxonomy database
 * parent tax_id                       -- parent node id in GenBank taxonomy database
 * rank                                -- rank of this node (superkingdom, kingdom, ...) 
 * embl code                           -- locus-name prefix; not unique
 * division id                         -- see division.dmp file
 * inherited div flag  (1 or 0)        -- 1 if node inherits division from parent
 * genetic code id                     -- see gencode.dmp file
 * inherited GC  flag  (1 or 0)        -- 1 if node inherits genetic code from parent
 * mitochondrial genetic code id       -- see gencode.dmp file
 * inherited MGC flag  (1 or 0)        -- 1 if node inherits mitochondrial gencode from parent
 * GenBank hidden flag (1 or 0)        -- 1 if name is suppressed in GenBank entry lineage
 * hidden subtree root flag (1 or 0)   -- 1 if this subtree has no sequence data yet
 * comments                            -- free-text comments and citations
 */

node_line_t *node_read_line(node_file_t* dmp) {
  // pull out the first 3 fields we want, and ignore the rest for now
  /*
  const char* node_format = "%d\t|\t%d\t|\t%s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|\t%*s\t|";
  node_line_t *line = malloc(sizeof(node_line_t));
  char rank[100];
  int nfields = fscanf(dmp->fp, node_format, &(line->taxid), &(line->parent_taxid), rank);
  */

  char l[1024]; // maximum line size is 1024 chars
  const char* delim = "|";
  char **parts;
  int i;
  node_line_t *line = malloc(sizeof(node_line_t));
  if(fgets(l, sizeof l, dmp->fp) != NULL) {
    char* ln = malloc(strlen(l)+1);
    strcpy(ln, l);
    parts = malloc(sizeof(char*) * 3); // 3 fields per line
    i = 0;
    parts[i] = strtok(ln, delim);
    while(parts[i] != NULL && strcmp(parts[i], "\n") != 0) {
      i++;
      if(i >= 3) break;
      parts[i] = strtok(NULL, delim);
      if(parts[i][strlen(parts[i])-1] == '\t') {
        parts[i][strlen(parts[i])-1] = '\0';
      }
      if(parts[i][0] == '\t') {
        parts[i]++;
      }
    }
    if(i != 3) { // line too short or small or -1 for EOF
      if(i == -1) {
        return NULL;
      } else {
        fprintf(stderr, "NODES.DMP line %d does not appear to be in expected format (has %d fields) - STOPPED READING HERE\n", dmp->cur_row, i);
        return NULL;
      }
    }
    line->taxid = atoi(parts[0]);
    line->parent_taxid = atoi(parts[1]);
    if(strcmp(parts[2], "no rank") == 0) {
      line->rank = RANK_NORANK;
    } else if(strcmp(parts[2], "superkingdom") == 0) {
      line->rank = RANK_SUPERKINGDOM;
    } else if(strcmp(parts[2], "kingdom") == 0) {
      line->rank = RANK_KINGDOM;
    } else if(strcmp(parts[2], "phylum") == 0) {
      line->rank = RANK_PHYLUM;
    } else if(strcmp(parts[2], "class") == 0) {
      line->rank = RANK_CLASS;
    } else if(strcmp(parts[2], "order") == 0) {
      line->rank = RANK_ORDER;
    } else if(strcmp(parts[2], "family") == 0) {
      line->rank = RANK_FAMILY;
    } else if(strcmp(parts[2], "genus") == 0) {
      line->rank = RANK_GENUS;
    } else if(strcmp(parts[2], "species") == 0) {
      line->rank = RANK_SPECIES;
    } else if(strcmp(parts[2], "subspecies") == 0) {
      line->rank = RANK_SUBSPECIES;
    } else {
      line->rank = RANK_NORANK;
    }

    free(parts);

    dmp->cur_row++;

  } else {
    return NULL;
  }

  return line;
}

node_file_t node_init(char* f) {
  // open bed file from path
  FILE *dmp_fp = fopen(f, "r");
  if( dmp_fp == NULL ) {
    fprintf(stderr, "Error reading NODES.DMP file '%s'\n", f);
  }

  node_file_t dmp = {dmp_fp, 0};

  return dmp;
}

int node_close(node_file_t* dmp) {
  fclose(dmp->fp);
}


/*
 * Build partial top-down taxonomic tree with counts
 */

taxtree* new_tree() {
  taxtree *tree = kh_init(nodehash);
  return tree;
}

int add_to_tree(taxonomy *tax, taxtree *tree, size_t taxid) {
  khint_t bin; // hash bin (result of kh_put)
  int absent;
  // do leaf node
  bin = kh_put(nodehash, tree, taxid, &absent);
  if(absent) {
    kv_init(kh_val(tree, bin).children);
    kh_val(tree, bin).unique_count = 1;
    kh_val(tree, bin).count = 1;
  } else {
    kh_val(tree, bin).unique_count++;
    kh_val(tree, bin).count++;
  }
  // set the child taxid to be added to parent's children if this node is new
  size_t child_taxid = absent ? taxid :0;
  //printf("  taxonomic hierarchy: ");
  while(taxid != 1) { // until we reach the root
    taxid = tax->nodes[taxid].parent;
    if(taxid == 0) {
      taxid = 1;
    }
    bin = kh_put(nodehash, tree, taxid, &absent);
    if(absent) {
      kv_init(kh_val(tree, bin).children);
      kh_val(tree, bin).count = 1;
      kh_val(tree, bin).unique_count = 0;
    } else {
      kh_val(tree, bin).count++;
    }
    if(child_taxid != 0) {
      kv_push(size_t, kh_val(tree, bin).children, child_taxid);
    }
    // set the child taxid to be added to parent's children if this node is new
    child_taxid = absent ? taxid :0;
    //printf("%d,", taxid);
  }
  //printf("\n");
}

void depth_first_traverse(taxonomy *tax, taxtree *tree, size_t taxid, int indent, FILE* o) {
  khint_t bin = kh_get(nodehash, tree, taxid);
  int absent = (bin == kh_end(tree)); 
  int i;
  if(!absent) {
    if(tax->nodes[taxid].rank > 0 || kh_val(tree, bin).unique_count > 0) {
      fprintf(o, "%d\t%d\t%d\t%c\t", taxid, kh_val(tree, bin).count, kh_val(tree, bin).unique_count, RANK_CHARS[tax->nodes[taxid].rank]);
      //for(i = 0; i < tax->nodes[taxid].rank; i++) {
      for(i = 0; i < indent; i++) {
        fprintf(o, "  ");
      }
      fprintf(o, "%s\n", tax->names[taxid]);
      indent++;
    }
    for(i = 0; i < kv_size(kh_val(tree, bin).children); i++) {
      depth_first_traverse(tax, tree, kv_A(kh_val(tree, bin).children, i), indent, o);
    }
  } else {
    if(taxid == 1) {
      fprintf(stderr, "Taxonomic tree is empty :(\n");
    } else {
      fprintf(stderr, "TaxID %d is apparently in the hierarchy but not a node in our tree - this is a bug\n", taxid);
    }
  }
}


/*
 * Other utilities
 */


int parse_acc2tax(char* f, khash_t(acc2asm) *a2a) {
  //kh_resize(acc2asm, a2a, 150000000); // there are ~133m accessions that we'll end up indexing - I'm not sure doing this saves any time (or space)
  FILE *fp = fopen(f, "r");
  if( fp == NULL ) {
    fprintf(stderr, "Error reading acc2tax file '%s'\n", f);
  }
  char l[1024]; // maximum line size is 1024 chars
  const char* delim = " ";
  char **parts;
  int i;
  khint_t bin; // hash bin (result of kh_put)
  int absent;
  parts = malloc(sizeof(char*) * 3); // 3 fields per line (accession ID, assembly ID, taxon ID)
  while(fgets(l, sizeof l, fp) != NULL) {
    char* ln = malloc(strlen(l)+1);
    strcpy(ln, l);
    i = 0;
    parts[i] = strtok(ln, delim);
    while(parts[i] != NULL && strcmp(parts[i], "\n") != 0) {
      i++;
      if(i >= 3) break;
      parts[i] = strtok(NULL, delim);
    }
    if(i != 3) { // line too short or small or -1 for EOF
      if(i == -1) {
        return 1;
      } else {
        fprintf(stderr, "acc2tax does not appear to be in expected format (some line has %d fields) - STOPPED READING HERE\n", i);
        return 1;
      }
    }
    // accession ID (key)
    char* a = malloc(strlen(parts[0])+1);
    strcpy(a, parts[0]);
    bin = kh_put(acc2asm, a2a, a, &absent);
    if(!absent) {
      fprintf(stderr, "Accession ID '%s' occurs more than once...\n", parts[1]);
      return 1;
    }
    kh_val(a2a, bin).assembly = malloc(strlen(parts[1])+1);
    strcpy(kh_val(a2a, bin).assembly, parts[1]);
    kh_val(a2a, bin).taxid = atoi(parts[2]);
  }
  free(parts);
  fclose(fp);
  return 0;
}

int lca(int taxid0, int taxid1, taxonomy* tax) {
  while(taxid0 != 1 && taxid0 != taxid1) {
    //fprintf(stderr, "lca taxids: %d %d\n", taxid0, taxid1);
    if(tax->nodes[taxid0].rank == 0 || tax->nodes[taxid0].rank > tax->nodes[taxid1].rank) {
      taxid0 = tax->nodes[taxid0].parent;
    } else {
      taxid1 = tax->nodes[taxid1].parent;
    }
    if(taxid0 == 0 || taxid1 == 0) {
      // error
      return 1;
    }
  }
  return taxid0;
}
