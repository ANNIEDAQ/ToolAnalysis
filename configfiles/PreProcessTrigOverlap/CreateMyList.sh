if [ "$#" -ne 2 ]; then
      echo "Usage: ./CreateMyList.sh DIR RUN"
      echo "Specified input variable must contain the directory where the files are stored, second variable run number"
      exit 1
fi

FILEDIR=$1
RUN=$2

NUMFILES=$(ls -1q ${FILEDIR}/RAWDataR${RUN}* | wc -l)

echo "NUMBER OF FILES IN ${FILEDIR}: ${NUMFILES}"

for p in $(seq 0 $NUMFILES)
do
	echo "${FILEDIR}/RAWDataR${RUN}S0p${p}" >> my_files.txt
done
