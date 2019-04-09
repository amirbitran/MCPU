

echo "Modifying FASTA file"
#sed -e :a -e '/^\n*$/{$d;N;};/\n$/ba' ./mcpu_prep/1ake.fasta

sed "$ s/\(.\{$1\}\)$//" ./mcpu_prep/1ake.fasta > ./mcpu_prep/1ake_trunc$1.fasta
#sed "$ s/.$//" ./mcpu_prep/1ake.fasta > ./mcpu_prep/1ake_trunc$1.fasta
#sed is a text editor
#$ character tells it to work with only the last line
#If we inputted s/.$//, that would remove the last character and replace it with an empty string
#but since we want to remove a variable number of characters, instead of usign ., we use \(.\{$1\}\), which effectively
#tells you to repeat that . N times


#sed '$'./mcpu_prep/1ake.fasta > ./mcpu_prep/1ake_trunc$1.fasta