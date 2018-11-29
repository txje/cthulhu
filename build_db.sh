mkdir -p db
cd db

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
md5sum taxdump.tar.gz.md5
tar -xf taxdump.tar.gz
rm taxdump.tar.gz

mkdir -p refseq
cd refseq

if [ ! -e bacteria_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
  mv assembly_summary.txt bacteria_assembly_summary.txt
fi
if [ ! -e archaea_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
  mv assembly_summary.txt archaea_assembly_summary.txt
fi
if [ ! -e fungi_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
  mv assembly_summary.txt fungi_assembly_summary.txt
fi
if [ ! -e protozoa_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt
  mv assembly_summary.txt protozoa_assembly_summary.txt
fi
if [ ! -e viral_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
  mv assembly_summary.txt viral_assembly_summary.txt
fi
if [ ! -e human_assembly_summary.txt ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/assembly_summary.txt
  mv assembly_summary.txt human_assembly_summary.txt
fi
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' bacteria_assembly_summary.txt archaea_assembly_summary.txt fungi_assembly_summary.txt protozoa_assembly_summary.txt viral_assembly_summary.txt human_assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
mkdir -p complete_genomes
cd complete_genomes
for f in $(cat ../ftpfilepaths); do
  filename="${f##*/}"
  if [ ! -e $filename ]; then
    wget $f
  fi
done
cd ..

for f in $(ls complete_genomes/*.fna.gz); do
  filename="${f##*/}"
  gcf="${filename%%.*}"
  echo $filename $gcf
  gunzip -c $f | awk '{if(substr($1,1,1) == ">") {if(index(tolower($0), "plasmid") == 0 && index(tolower($0), "chloroplast") == 0 && index(tolower($0), "mitochondrion") == 0 && index(tolower($0), "plastid") == 0) {print substr($1,2)" '$gcf'";}}}' >> accession_map.txt
  gunzip -c $f | awk '{if(substr($1,1,1) == ">") {if(index(tolower($0), "plasmid") == 0 && index(tolower($0), "chloroplast") == 0 && index(tolower($0), "mitochondrion") == 0 && index(tolower($0), "plastid") == 0) {skip=0; print;} else {skip=1}} else {if(skip==0) {print;}}}' >> refseq_filtered.fasta
done

if [ ! -e nucl_gb.accession2taxid.gz ]; then
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
fi
gunzip -c nucl_gb.accession2taxid.gz | awk 'FNR==NR{a[$1]=$2;next}{if($2 in a){print $2, a[$2], $3}}' accession_map.txt - > merged_accession_map.txt

for f in $(ls complete_genomes); do
  printf "${f%.fna.gz}\t"
  gunzip -c complete_genomes/$f | head -n 1 | awk '{print substr($1,2)}'
done > file_accessions.txt

cd ../..
