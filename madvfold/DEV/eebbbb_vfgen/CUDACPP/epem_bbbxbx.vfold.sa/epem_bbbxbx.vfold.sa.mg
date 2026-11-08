set stdout_level DEBUG
set apply_flavor_grouping False --no_save
set zerowidth_tchannel F
generate e+ e- > b b b~ b~
output standalone_cudacpp /data/avalassi/NEGWGT2026/testCG-MadGraph7/madvfold/DEV/eebbbb_vfgen/CUDACPP/epem_bbbxbx.vfold.sa --vfold
