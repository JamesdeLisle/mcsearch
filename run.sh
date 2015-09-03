#! /bin/bash

#for ((n=0;n<20;n++))
#    do
#        cd ~/Dropbox/Work/MCSearch/numerics/1DSearch-v1.0/1DSearch-v1.0/source
#        python -B search1DImag.py
#        mv *.dat ~/Dropbox/Work/MCSearch/numerics/1DSearch-v1.0/1DSearch-v1.0/rundata
#    done

main()
{
    cd ~/Dropbox/Work/MCSearch/numerics/working/source
    python -B search1DImag.py
    mv *.dat ~/Dropbox/Work/MCSearch/numerics/rundata
    cd ..
}


main
