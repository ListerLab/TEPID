from sys import argv


def checkArgs(arg1, arg2):
    """
    arg1 is short arg, eg h
    arg2 is long arg, eg host
    """
    args = argv[1:]
    if arg1 in args:
        index = args.index(arg1)+1
        variable = args[index]
        return variable
    elif arg2 in args:
        index = args.index(arg2)+1
        variable = args[index]
        return variable
    else:
        variable = raw_input("\nEnter {arg2}: ".format(arg2=arg2))
        return variable


def create_coords(bedfile, saveas):
    with open(bedfile, 'r') as infile, open(saveas, 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            chr1 = int(line[0])
            start1 = int(line[1])
            stop1 = int(line[2])
            strand1 = line[8]
            chr2 = int(line[3])
            start2 = int(line[4])
            stop2 = int(line[5])
            read = line[6]
            strand2 = line[9]
            if chr1 == chr2 and strand1 != strand2:
                if start2 >= stop1:
                    start = stop1
                    stop = start2
                else:
                    start = stop2
                    stop = start1
                gapsize = stop - start
                if gapsize < 20000:
                    outfile.write('{ch}\t{start}\t{stop}\t{read}\n'.format(ch=chr1,
                                                                           start=start,
                                                                           stop=stop,
                                                                           read=read))
                else:
                    pass
            else:
                pass

input_file = checkArgs('b', 'bed')
save_file = checkArgs('f', 'file')

create_coords(input_file, save_file)
