#!/bin/bash
find ../SOLERA/ -iname *.hpp -o -iname *.ccp | xargs clang-format -style=file -i
find ../examples/ -iname *.hpp -o -iname *.ccp | xargs clang-format -style=file -i
