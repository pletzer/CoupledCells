# @ shell = /bin/bash
#
# @ job_name = c256f_4x64
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 1:00:00
#
# @ account_no = nesi00269
#
# @ network.MPI = sn_all,not_shared,US
# @ task_affinity = core(1)
#
# @ output               = $(job_name).$(schedd_host).$(jobid).o
# @ error                = $(job_name).$(schedd_host).$(jobid).e
# @ notification         = never
# @ class                = General
# @ node = 4
# @ tasks_per_node = 64
#
# @ queue

exe=/home/pletzera/CoupledCells/build-sundials/coupledCellsModel

# On the compute nodes /home is /gpfs_external/filesets/nesi/home
exe=`echo $exe | perl -ne "s#home#gpfs_external/filesets/nesi/home#;print;"`

cmd="time poe $exe -args \"-f config.txt -S solution -T profiling -t 100.0 -w 1.0 -i 1.e-2\""
echo "running..."
echo "$cmd"
$cmd
