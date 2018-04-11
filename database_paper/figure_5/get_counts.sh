#!/bin/bash

awk -F_ '{print $4}'  Clustering-of-Blast-All-Hit-Profiles-171214_D00390_0336_ACBG26ANXX.txt | sort  | uniq -c | wc
wc Clustering-of-Blast-All-Hit-Profiles-171214_D00390_0336_ACBG26ANXX.txt
