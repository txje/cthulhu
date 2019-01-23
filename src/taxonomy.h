#include <stdint.h>
#include "../incl/klib/kvec.h"
#include "../incl/klib/khash.h"

#define RANK_NORANK 0 // root
#define RANK_SUPERKINGDOM 1 // (domain) bacteria, archaea, eukaryota, viruses, viroid
#define RANK_KINGDOM 2 // fungi, green plants, etc...
#define RANK_PHYLUM 3
#define RANK_CLASS 4
#define RANK_ORDER 5
#define RANK_FAMILY 6
#define RANK_GENUS 7
#define RANK_SPECIES 8
#define RANK_SUBSPECIES 9

static char* RANK_CHARS = " DKPSOFGSB";

/*
All levels in taxdump (the ones we actually use starred):

* no rank (root and others...)
superkingdom
* kingdom
subkingdom
superphylum
* phylum
subphylum
superclass
* class
subclass
infraclass
cohort
superorder
* order
suborder
infraorder
parvorder
superfamily
* family
subfamily
tribe
subtribe
* genus
subgenus
* species
* subspecies
varietas
forma
*/

typedef struct node {
  size_t parent;
  int rank;
} node;

typedef struct taxonomy {
  char** names;
  node* nodes;
} taxonomy;

taxonomy* read_taxonomy(char* name_f, char* node_f);
void free_tax(taxonomy* tax);

/*
 * names.dmp
 */
typedef struct name_file_t {
  FILE *fp;
  size_t cur_row;
} name_file_t;

typedef struct name_line_t {
  int taxid;
  char *name;
  char *unique_name;
  char *name_class;
} name_line_t;

name_line_t *name_read_line(name_file_t* dmp);
name_file_t name_init(char* f);
int name_close(name_file_t* dmp);


/*
 * nodes.dmp
 */
typedef struct node_file_t {
  FILE *fp;
  size_t cur_row;
} node_file_t;

typedef struct node_line_t {
  int taxid;
  int parent_taxid;
  int rank;
  // there exist many more fields (see header of taxonomy.c), but we're not using them right now
} node_line_t;

node_line_t *node_read_line(node_file_t* dmp);
node_file_t node_init(char* f);
int node_close(node_file_t* dmp);


/*
 * for building a partial top-down taxonomic hierarchy
 */

typedef struct tree_node {
  kvec_t(size_t) children;
  uint32_t count;
  uint32_t unique_count;
} tree_node;

// creates int:tree_node hash
// to taxonomy IDs to counts
KHASH_MAP_INIT_INT(nodehash, tree_node);

typedef khash_t(nodehash) taxtree;

taxtree* new_tree();
int add_to_tree(taxonomy *tax, taxtree *tree, size_t taxid);
void depth_first_traverse(taxonomy *tax, taxtree *tree, size_t taxid, int indent, FILE* o);


/*
 * Utils
 */

typedef struct asm_tax {
  char* assembly;
  size_t taxid;
} asm_tax;

// creates string:asm_tax hash
// to map accession (NC...) to assembly (GCF...), and taxid
KHASH_MAP_INIT_STR(acc2asm, asm_tax);

int parse_acc2tax(char* f, khash_t(acc2asm)*);

int lca(size_t taxid0, size_t taxid1, taxonomy* tax);

size_t* get_hierarchy(size_t taxid, taxonomy* tax);
