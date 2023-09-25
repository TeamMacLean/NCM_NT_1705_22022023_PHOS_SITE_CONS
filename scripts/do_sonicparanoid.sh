#!/usr/bin/bash

singularity exec $4/sonicparanoid.img sonicparanoid -i $1 -o $2 --project-id snakerun4 -t $3 --overwrite