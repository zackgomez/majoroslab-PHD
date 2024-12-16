# PHD
Piled Higher and Deeper


## Building
```
apt install libgsl23 libgsl-dev
```
```
git clone git@github.com:bmajoros/BOOM.git
cd BOOM
make
```

```
make
```




## Useful commands

Convert BAM to sorted SAM
```
samtools sort NA12878_chr21.bam.gz -Osam -o NA12878_chr21.sorted.sam.gz
```