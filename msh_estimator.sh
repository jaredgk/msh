#Script for uncompressed, small VCF files. Usage with large (>1GB) may cause issues


if [ "$#" -eq 0 ]
then
    echo "Usage: bash msh_estimator.sh (VCF name) (Tag for dataset) (optional name for genetic map file)"
    exit 0
fi

VCF=$1
SETTAG=$2
MAPARG=""

if [ "$#" -eq 3 ]
then
    MAP=$3
    MAPARG="--gen "$MAP
fi

echo $MAPARG

TAGD=`dirname $1`
TAGB=`basename $1 .vcf`
TAG=$TAGD"/"$TAGB
FTAG=$TAGD"/"$TAGB"_"$SETTAG

echo $TAG
LEFT_VCF=$TAG".vcf"
RIGHT_VCF=$TAG"_reversed.vcf"
SSL=$FTAG"_left_results.txt.gz"
SSR=$FTAG"_right_results.txt.gz"
SSE=$FTAG"_estimates.txt"

#If this is done outside of script, make sure to unzip file first
#then use 'tac (vcf)', otherwise memory issues may occur with piping
if [ ! -f $RIGHT_VCF ]
then
    echo "Making "$RIGHT_VCF
    grep -v "\#" $LEFT_VCF | tac > $RIGHT_VCF
fi

    #Options for msh_from_vcf:
    # --gen [filename]: Genetic map
    # --sub [filename]: File with indices of desired individuals
    # --squish: Used for real genetic maps, will ensure all regions
    #           have non-zero genetic rate
    # VCF is sent via standard input
cat $LEFT_VCF | python msh_from_vcf.py $MAPARG | gzip -c > $SSL
cat $RIGHT_VCF | python msh_from_vcf.py $MAPARG | tac | gzip -c > $SSR
    #Positional args: [left lengths] [right lengths] (either gzipped or not)
    #Options:
    # --start [n]: Start estimates on line n
    # --end [n]: End estimates after line n
    # --mut [f]: Set mutation rate
    # --rec [f]: Set recombination rate
    # --n [n]: Set (ancestral?) population size to n
    # --n0 [n]: Set (current?) population size to n
python2.7 aae_work.py $SSL $SSR --mut 1e-6 --rec 1e-6 > $SSE


