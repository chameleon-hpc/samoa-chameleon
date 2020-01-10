#!/bin/bash

source "./config.sh"

MACHINE=$(basename -- "$1")
MACHINE=${MACHINE%.*}
export PROFILE="./machines/$MACHINE"

export CURRENT=${CURRENT:="./logs/$(date +%Y-%m-%d_%H-%M-%S)_$MACHINE"}

if [ ! -d "$PROFILE" ]; then
    echo "Profile $MACHINE does not exist."
    exit 1
else
    echo "Using profile $MACHINE"
fi

mkdir -p $CURRENT
