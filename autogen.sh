#!/bin/sh
autoreconf -i || exit 1
./configure || exit 1
