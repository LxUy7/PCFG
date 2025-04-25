#test.sh
#!/bin/bash
#PBS -N test
#PBS -l nodes=1

/usr/local/bin/pssh -h $PBS_NODEFILE mkdir -p /home/sTest 1>&2

scp master_ubss1:/home/sTest/test /home/sTest 
/usr/local/bin/pscp -h $PBS_NODEFILE /home/sTest/test /home/sTest 1>&2
/home/sTest/test 
