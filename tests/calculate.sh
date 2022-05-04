#!/bin/bash
set -o nounset
set -o pipefail

cd $1

ls | sort

find -name '*.txt' | xargs md5sum | sort

