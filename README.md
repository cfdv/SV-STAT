# Installation and testing

## With docker

After you clone the repository, change to the sv-stat directory. Then build the image (~1 hour; ~3 GB):
```
docker build --rm -t cfdavis/docker-sv-stat .
```

**OR** use my pre-built image on dockerhub:
```
docker pull cfdavis/docker-sv-stat
```

Then run sv-stat on some test data: 
 docker run --rm -v /path/to/reference:/work/ref -v /path/to/svstat/repo:/work/sv-stat -it cfdavis/docker-sv-stat /work/sv-stat/src/svstat.sh -m /work/sv-stat/test/bam/hs1011.bam /work/sv-stat/test/metadata/ABori.txt

## Manually

Configure environment and install the prerequisites

SVSTAT_SRC_DIR=/path/to/svstat/src
CDBFASTA_BIN_DIR=/path/to/cdbfasta
BWA_BIN_DIR=/path/to/bwa
PICARD_JAR_DIR=/path/to/picard-tools
BEDTOOLS_BIN_DIR=/path/to/BEDTools/bin
REF_DIR=/path/to/assembly
export SVSTAT_SRC_DIR
export CDBFASTA_BIN_DIR
export BWA_BIN_DIR
export PICARD_JAR_DIR
export BEDTOOLS_BIN_DIR
export REF_DIR

other packages used, and their versions [known dependency]:

cdbfasta version 0.99
bwa 0.5.9-r16
cat (GNU coreutils) 8.4
grep GNU grep 2.6.3
GNU Awk 3.1.7
samtools Version: 0.1.18 (r982:295) [bam2fq requires 0.1.17+]
java version "1.6.0_20"
picard-tools-1.40
join (GNU coreutils) 8.4
sort (GNU coreutils) 8.4
tee (GNU coreutils) 8.4
uniq (GNU coreutils) 8.4
bedtools (v2.11.2)

bioperl packages:

Bio::DB::Sam

paths and file formats:

hypoDBgen assumes genomic sequence data is on $REF_DIR/chromosomes/chr*.ol (see above for REF_DIR environmental variable)
- .ol format => where fasta header is stripped, and entire sequence is stored on the first line of the file.

hard-coded SVSTAT parameters in hypoDBgen.pl -- algorithm parameters to limit the number of junctions considered

$significant_read_tail_length = 4;      #a junction is not considered if neither one has a read of this length or greater
$significant_stack_sum_tail_lengths = 8;        #a junction is not considered if the sum of the tail lengths in both stacks is this number or less
$significant_deletion_length = 2;       #forward and reverse stack coordinates have to be separated by a gap of this size or larger
$cj_per_page = 3000000; #3e6 candidate junctions per array, to stay within bwa's index size limit

Then analyze test data as shown in the docker example