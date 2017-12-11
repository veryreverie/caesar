file_type=$1
directory=$2
num_cores=$3
seedname=$4

wd=$( pwd )

echo 'file type        : '$file_type
echo 'directory        : '$directory
echo 'no. cores        : '$num_cores
echo 'seedname         : '$seedname
echo 'working directory: '$wd
if [ ! -e "$directory/$seedname.cell" ]; then
  echo "Error: $directory/$seedname.cell does not exist."
fi
echo ''
