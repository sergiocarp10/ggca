#!/bin/bash
# Arg 1: example program name
# Arg 2: number of repetitions
# Arg 3: dataset size in MB (100, 500 or 1500)
# Arg 4: number of threads

# Set env variables
PROGRAM_NAME=$1

 # Run tests
for T in $4
do
	echo "$PROGRAM_NAME: Testing T=$T"

	for ((i=1; i<=$2; i++))
	do
		cargo run --example $PROGRAM_NAME --no-default-features --release -q $T "../gem-$3mb.csv"
	done
done
