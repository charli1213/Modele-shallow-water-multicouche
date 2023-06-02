#!/bin/bash

layers=(1 2 3)
for layer in ${layers[@]}; do
    # Pipe the input into our program
    python readata.py <<< $layer
done

