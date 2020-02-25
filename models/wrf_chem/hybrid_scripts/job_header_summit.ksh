#!/bin/ksh -aeux
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
rm -rf job.header
cat << EOF > job.header
#!/bin/ksh -aeux
#SBATCH --job-name $1
#SBATCH --output $1.log
#SBATCH --account $6
#SBATCH --qos $2
#SBATCH --time $3
#SBATCH --nodes $4
#SBATCH --ntasks $5
#SBATCH --partition shas
#
EOF
#
rm -rf job.ksh
cat job.header $7 > job.ksh
rm -rf job.header
#
# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/hybrid_scripts/job_script_summit.ksh $
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
# $Revision: 13133 $
# $Date: 2019-04-25 15:47:54 -0600 (Thu, 25 Apr 2019) $
