for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *.sra);do
            fastq-dump --split-3 -v $myfile
            rm -f $myfile
        done
        cd ..
    fi
done
