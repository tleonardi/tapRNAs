#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

mkdir -p $BASEDIR/data/fastq/single/mmu
mkdir -p $BASEDIR/data/fastq/single/hsa
mkdir -p $BASEDIR/data/fastq/paired/mmu
mkdir -p $BASEDIR/data/fastq/paired/hsa


echo TISSUE ATLAS SINGLE END, BOTH MOUSE AND HUMAN
echo 1-5
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306757.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306757.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306757/SRR306757.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306757.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306758.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306758.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306758/SRR306758.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306758.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306759.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306759.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306759/SRR306759.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306759.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306760.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306760.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306760/SRR306760.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306760.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306761.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306761.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306761/SRR306761.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306761.fastq.gz 
fi

echo 6-10
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306762.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306762.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306762/SRR306762.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306762.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306763.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306763.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306763/SRR306763.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306763.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306764.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306764.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306764/SRR306764.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306764.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306765.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306765.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306765/SRR306765.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306765.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306766.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306766.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306766/SRR306766.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306766.fastq.gz 
fi

echo 11-15
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306767.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306767.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306767/SRR306767.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306767.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306768.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306768.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306768/SRR306768.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306768.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306769.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306769.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306769/SRR306769.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306769.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306770.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306770.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306770/SRR306770.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306770.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306771.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306771.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306771/SRR306771.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306771.fastq.gz 
fi

echo 16-20
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306772.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306772.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306772/SRR306772.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306772.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306773.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306773.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306773/SRR306773.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306773.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306774.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306774.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306774/SRR306774.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306774.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306775.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306775.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306775/SRR306775.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306775.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR306776.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR306776.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306776/SRR306776.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR306776.fastq.gz 
fi

echo 16-25
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306838.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306838.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306838/SRR306838.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306838.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306839.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306839.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306839/SRR306839.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306839.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306841.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306841.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306841/SRR306841.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306841.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306843.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306843.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306843/SRR306843.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306843.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306844.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306844.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306844/SRR306844.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306844.fastq.gz 
fi

echo 25-30
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306845.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306845.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306845/SRR306845.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306845.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306846.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306846.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306846/SRR306846.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306846.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306847.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306847.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306847/SRR306847.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306847.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306848.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306848.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306848/SRR306848.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306848.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306849.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306849.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306849/SRR306849.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306849.fastq.gz 
fi

echo 31-35
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306850.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306850.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306850/SRR306850.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306850.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306851.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306851.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306851/SRR306851.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306851.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306852.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306852.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306852/SRR306852.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306852.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306853.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306853.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306853/SRR306853.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306853.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306854.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306854.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306854/SRR306854.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306854.fastq.gz 
fi

echo 36-40
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306855.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306855.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306855/SRR306855.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306855.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306856.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306856.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306856/SRR306856.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306856.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306857.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306857.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306857/SRR306857.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306857.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/single/hsa/SRR306858.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/hsa/SRR306858.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306858/SRR306858.fastq.gz > $BASEDIR/data/fastq/single/hsa/SRR306858.fastq.gz 
fi


echo Paired brain from atlas, HUMAN
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR306840_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR306840_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306840/SRR306840_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR306840_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR306840_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR306840_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306840/SRR306840_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR306840_2.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR306842_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR306842_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306842/SRR306842_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR306842_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR306842_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR306842_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306842/SRR306842_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR306842_2.fastq.gz 
fi



echo ENCODE SINGLE END, MOUSE
# Encode LICR_RnaSeq_ES-Bruce4_E0 SINGLE REP1
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR496249.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR496249.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496249/SRR496249.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR496249.fastq.gz 
fi

# Encode LICR_RnaSeq_ES-Bruce4_E0 SINGLE REP2
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR496250.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR496250.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496250/SRR496250.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR496250.fastq.gz 
fi

# Encode PSU CH12 SINGLE REP1
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549363.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549363.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549363/SRR549363.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549363.fastq.gz 
fi
# Encode PSU CH12 SINGLE REP2
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549364.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549364.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549364/SRR549364.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549364.fastq.gz 
fi

# Encode LICR MEL SINGLE REP1
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR496221.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR496221.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496221/SRR496221.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR496221.fastq.gz 
fi
# Encode LICR MEL SINGLE REP2
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR496222.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR496222.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496222/SRR496222.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR496222.fastq.gz 
fi

# Encode PSU MEL SINGLE REP1
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549339.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549339.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549339/SRR549339.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549339.fastq.gz 
fi
# Encode PSU MEL SINGLE REP2
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549340.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549340.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549340/SRR549340.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549340.fastq.gz 
fi
# Encode PSU MEL SINGLE REP1 DMSO
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549335.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549335.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549335/SRR549335.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549335.fastq.gz 
fi
# Encode PSU MEL SINGLE REP2 DMSO
if [[ ! -f $BASEDIR/data/fastq/single/mmu/SRR549336.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/single/mmu/SRR549336.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR549/SRR549336/SRR549336.fastq.gz > $BASEDIR/data/fastq/single/mmu/SRR549336.fastq.gz 
fi



echo ENCODE PAIRED END, HUMAN
echo 1-6
#SRR307919 CshlLong_RnaSeq_H1-hESC_cytosol_longPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307919_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307919_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307919/SRR307919_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307919_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307919_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307919_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307919/SRR307919_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307919_2.fastq.gz 
fi


#SRR317040 CshlLong_RnaSeq_H1-hESC_cytosol_longNonPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR317040_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR317040_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/SRR317040/SRR317040_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR317040_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR317040_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR317040_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/SRR317040/SRR317040_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR317040_2.fastq.gz 
fi

#SRR307925 CshlLong_RnaSeq_H1-hESC_nucleus_longPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307925_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307925_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307925/SRR307925_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307925_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307925_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307925_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307925/SRR307925_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307925_2.fastq.gz 
fi


echo 7-12
#SRR317039 CshlLong_RnaSeq_H1-hESC_nucleus_longNonPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR317039_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR317039_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/SRR317039/SRR317039_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR317039_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR317039_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR317039_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/SRR317039/SRR317039_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR317039_2.fastq.gz 
fi

#SRR307911 CshlLong_RnaSeq_H1-hESC_cell_longPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307911_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307911_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307911/SRR307911_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307911_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307911_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307911_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307911/SRR307911_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307911_2.fastq.gz 
fi


#SRR307912 CshlLong_RnaSeq_H1-hESC_cell_longPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307912_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307912_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307912/SRR307912_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307912_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307912_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307912_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307912/SRR307912_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307912_2.fastq.gz 
fi


echo 13-16
#SRR307923 CshlLong_RnaSeq_H1-hESC_cell_longNonPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307923_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307923_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307923/SRR307923_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307923_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307923_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307923_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307923/SRR307923_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307923_2.fastq.gz 
fi


#SRR307924 CshlLong_RnaSeq_H1-hESC_cell_longNonPolyA
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307924_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307924_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307924/SRR307924_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307924_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR307924_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR307924_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307924/SRR307924_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR307924_2.fastq.gz 
fi

# Caltech_RnaSeq_HSMM_2x75_200

# Rep1
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521516_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521516_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521516/SRR521516_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521516_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521516_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521516_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521516/SRR521516_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521516_2.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521517_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521517_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521517/SRR521517_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521517_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521517_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521517_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521517/SRR521517_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521517_2.fastq.gz 
fi

# Rep2
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521518_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521518_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521518/SRR521518_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521518_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521518_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521518_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521518/SRR521518_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521518_2.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521519_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521519_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521519/SRR521519_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521519_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521519_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521519_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521519/SRR521519_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521519_2.fastq.gz 
fi


# Caltech RNASeq GM12878

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521447_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521447_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521447/SRR521447_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521447_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521447_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521447_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521447/SRR521447_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521447_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521448_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521448_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521448/SRR521448_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521448_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521448_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521448_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521448/SRR521448_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521448_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521449_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521449_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521449/SRR521449_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521449_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521449_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521449_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521449/SRR521449_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521449_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521450_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521450_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521450/SRR521450_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521450_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521450_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521450_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521450/SRR521450_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521450_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521451_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521451_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521451/SRR521451_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521451_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521451_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521451_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521451/SRR521451_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521451_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521452_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521452_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521452/SRR521452_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521452_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521452_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521452_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521452/SRR521452_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521452_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521453_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521453_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521453/SRR521453_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521453_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521453_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521453_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521453/SRR521453_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521453_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521454_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521454_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521454/SRR521454_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521454_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521454_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521454_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521454/SRR521454_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521454_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521455_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521455_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521455/SRR521455_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521455_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521455_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521455_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521455/SRR521455_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521455_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521456_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521456_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521456/SRR521456_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521456_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521456_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521456_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521456/SRR521456_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521456_2.fastq.gz 
fi

# Caltech RNASeq K562

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521457_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521457_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521457/SRR521457_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521457_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521457_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521457_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521457/SRR521457_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521457_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521458_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521458_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521458/SRR521458_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521458_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521458_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521458_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521458/SRR521458_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521458_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521459_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521459_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521459/SRR521459_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521459_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521459_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521459_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521459/SRR521459_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521459_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521460_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521460_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521460/SRR521460_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521460_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521460_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521460_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521460/SRR521460_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521460_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521461_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521461_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521461/SRR521461_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521461_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521461_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521461_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521461/SRR521461_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521461_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521462_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521462_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521462/SRR521462_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521462_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521462_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521462_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521462/SRR521462_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521462_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521463_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521463_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521463/SRR521463_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521463_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521463_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521463_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521463/SRR521463_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521463_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521464_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521464_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521464/SRR521464_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521464_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521464_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521464_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521464/SRR521464_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521464_2.fastq.gz 
fi

if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521465_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521465_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521465/SRR521465_1.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521465_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/hsa/SRR521465_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/hsa/SRR521465_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR521/SRR521465/SRR521465_2.fastq.gz > $BASEDIR/data/fastq/paired/hsa/SRR521465_2.fastq.gz 
fi









echo ENCODE PAIRED END, MOUSE

# Stanford_RnaSeq_ES-E14 PAIRED REP1
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR530639_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR530639_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530639/SRR530639_1.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR530639_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR530639_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR530639_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530639/SRR530639_2.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR530639_2.fastq.gz 
fi

# Stanford_RnaSeq_ES-E14 PAIRED REP2
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR530640_1.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR530640_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530640/SRR530640_1.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR530640_1.fastq.gz 
fi
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR530640_2.fastq.gz  ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR530640_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530640/SRR530640_2.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR530640_2.fastq.gz 
fi

# Caltech_RnaSeq_C2C12_2x75 PAIRED
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR496442_1.fastq.gz ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR496442_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496442/SRR496442_1.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR496442_1.fastq.gz
fi
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR496442_2.fastq.gz ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR496442_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496442/SRR496442_2.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR496442_2.fastq.gz
fi

# Caltech_RnaSeq_C2C12_2x75 PAIRED EqS_2.0pct_60hr 
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR496443_1.fastq.gz ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR496443_1.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496443/SRR496443_1.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR496443_1.fastq.gz
fi
if [[ ! -f $BASEDIR/data/fastq/paired/mmu/SRR496443_2.fastq.gz ]]; then
	wget -c -o $BASEDIR/data/fastq/paired/mmu/SRR496443_2.log -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR496/SRR496443/SRR496443_2.fastq.gz > $BASEDIR/data/fastq/paired/mmu/SRR496443_2.fastq.gz
fi

# 

