#!/bin/bash
# Arg 1: example program name
# Arg 2: number of repetitions
# Arg 3: number of threads

# Set env variables
PROGRAM_NAME=$1

 # Run tests
for size in 100 500 1500
do
	echo "$PROGRAM_NAME: Testing T=$3, size=${size}"

	for ((i=1; i<=$2; i++))
	do
		cargo run --example $PROGRAM_NAME --no-default-features --release -q $3 "../gem-${size}mb.csv"
	done
done
