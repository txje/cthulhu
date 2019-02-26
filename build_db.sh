mkdir -p db_20190226
cd db_20190226

echo "Getting NCBI taxonomy..."
wget -q ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xf taxdump.tar.gz
rm taxdump.tar.gz

mkdir -p refseq
cd refseq

echo "Getting refseq assembly summaries..."
if [ ! -e bacteria_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
  mv assembly_summary.txt bacteria_assembly_summary.txt
fi
if [ ! -e archaea_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
  mv assembly_summary.txt archaea_assembly_summary.txt
fi
if [ ! -e fungi_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
  mv assembly_summary.txt fungi_assembly_summary.txt
fi
if [ ! -e protozoa_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt
  mv assembly_summary.txt protozoa_assembly_summary.txt
fi
if [ ! -e viral_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
  mv assembly_summary.txt viral_assembly_summary.txt
fi
if [ ! -e human_assembly_summary.txt ]; then
  wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/assembly_summary.txt
  mv assembly_summary.txt human_assembly_summary.txt
fi

awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $6}' bacteria_assembly_summary.txt archaea_assembly_summary.txt fungi_assembly_summary.txt protozoa_assembly_summary.txt viral_assembly_summary.txt human_assembly_summary.txt > taxids
echo "$(wc -l taxids) assemblies found."

# check that all taxids exist in the taxonomy database - they should
missing=$(sort -n taxids | uniq | awk '{if(NR==FNR){a[$1]=1}else{if($1 in a){}else{print "taxonomy missing "$1}}}' ../names.dmp -)
if [ -z "${#missing}" ]; then # empty or only whitespace
  echo "---------------- WARNING: SOME TAXONOMY IDs MISSING FROM TAXONOMY DATABASE ----------------"
  echo $missing
fi

awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' bacteria_assembly_summary.txt archaea_assembly_summary.txt fungi_assembly_summary.txt protozoa_assembly_summary.txt viral_assembly_summary.txt human_assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
rm ftpdirpaths
echo "Downloading $(wc -l ftpfilepaths) files by ftp... this will take a long time."

mkdir -p complete_genomes
cd complete_genomes

st=`date +%s`

while read f <&3 && read t <&4; do
  filename="${f##*/}"
  if [ ! -e $filename ]; then
    #echo "get $f"
    wget -q $f
  fi
  gcf="${filename%%.*}"
  gunzip -c $filename | awk '{if(substr($1,1,1) == ">") {if(index(tolower($0), "plasmid") == 0 && index(tolower($0), "chloroplast") == 0 && index(tolower($0), "mitochondrion") == 0 && index(tolower($0), "plastid") == 0) {print substr($1,2)" '$gcf' '$t'";}}}' >> ../accession_map.txt
  gunzip -c $filename | awk '{if(substr($1,1,1) == ">") {if(index(tolower($0), "plasmid") == 0 && index(tolower($0), "chloroplast") == 0 && index(tolower($0), "mitochondrion") == 0 && index(tolower($0), "plastid") == 0) {skip=0; print;} else {skip=1}} else {if(skip==0) {print;}}}' >> ../refseq_filtered.fasta
done 3<../ftpfilepaths 4<../taxids
cd ..

en=`date +%s`

echo "$(wc -l accession_map.txt) accessions downloaded in $((en-st)) seconds."

awk '{print $1}' accession_map.txt > file_accessions.txt

cd ../..
