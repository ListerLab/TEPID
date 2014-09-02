for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *.sra);do
            rm -f $myfile
        done
        cd ..
    fi
done