#!/bin/bash

bedClip $1 <(fetchChromSizes $2 2> /dev/null) $3
