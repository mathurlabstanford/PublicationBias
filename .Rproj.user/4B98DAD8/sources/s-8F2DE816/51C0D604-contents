
# to be run by stitch.sbatch

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
start.num = as.numeric( args[1] )  # starting results number to stitch
stop.num = as.numeric( args[2] )  # stopping results number to stitch

path = "/home/groups/manishad/RRR"
setwd(path)
source("functions_RRR.R")

######## STITCH LONG FILES ########

library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
s = stitch_files(.results.singles.path = "/home/groups/manishad/RRR/sim_results/long",
                 .results.stitched.write.path = "/home/groups/manishad/RRR/sim_results/overall_stitched",
                 .name.prefix = "long_results",
                 .stitch.file.name="stitched.csv",
                 .start.num = start.num,
                 .stop.num = stop.num )




# stitch on Sherlock
# sbatch -p qsu,normal,owners /home/groups/manishad/RRR/stitch_sbatch_files/stitch_4.sbatch
# sacct --jobs=49474291 --format=User,JobID,account,Timelimit,elapsed,ReqMem,MaxRss,ExitCode

# # move it to Desktop
# scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/RRR/sim_results/overall_stitched/stitched.csv ~/Desktop
# Vegemite2019
