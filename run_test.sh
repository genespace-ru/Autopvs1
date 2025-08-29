export PATH=$PATH:/mnt/disk1/ivan/home/tools/ensembl-vep
DEBUG="-Xdebug -Xrunjdwp:transport=dt_socket,address=0.0.0.0:8123,server=y,suspend=y"
java $DEBUG -cp target/autopvs1*.jar ru.biosoft.autopvs1.Main -g hg38 -i test-data/1/test.vcf -o  test-data/1/out/
