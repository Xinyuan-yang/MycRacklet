#!/bin/bash

executable=
name=
envi=
reference=

while :
do
    case "$1" in
	-e)
	    executable=$2
	    shift 2
	    ;;
	-E)
	    envi="$2"
	    shift 2
	    ;;
	-r)
	    reference="$2"
	    shift 2
	    ;;
	-n)
	    name="$2"
	    shift 2
	    ;;
	--) # End of all options
	    shift
	    break
	    ;;
	-*)
	    echo "Error: Unknown option: $1" >&2
	    show_help
	    exit 1
	    ;;
	*) #No more options
	    break
	    ;;
    esac
done
	
_args=$@

if [ -n "${envi}" ]; then
    source ${envi}
fi

if [ -n "${executable}" ]; then
    ./${executable}
fi

if [ -n "${reference}" ]; then
    ./${name} > ${name}.lastout 2> ${name}.errout
   diff -w ${name}.lastout ${reference}
fi
