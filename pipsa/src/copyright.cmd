find . -name \*.f -exec echo "cat copyright.txt '{}' > tmp && mv tmp '{}'" \;>exec.sh 
