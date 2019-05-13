#!/bin/sh
READ=$1
SIZE=$2
OUTPUT=$3

grep -v '^>' $READ | tr -d '\n' | fold -w $SIZE | nl -n ln -s '
' | sed 's/^/>/;N' > $OUTPUT
