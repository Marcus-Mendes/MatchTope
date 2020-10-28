for pdbfile in *.pdb ; do
sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
done
<<comment
for pdbfile in *.pdb ; do
sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
done
for pdbfile in *.pdb ; do
sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
done
for pdbfile in *.pdb ; do
sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
done
for pdbfile in *.pdb ; do
sed -ri 's/(.*)\s+[^\s]+$/\1/' $pdbfile
done


or

awk -i inplace '{$0=gensub(/\s*\S+/,"",3)}1' file (substitui o 3 pelo numero da coluna q vc quer)


