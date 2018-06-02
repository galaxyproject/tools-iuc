ids=c(1:50)
io1 = rep(c(0,0,0,0,1),10)
ip2 = runif(50)+0.1
ip2[which(ip2>1.0)] = 1.0
ip1 = runif(50)+0.05
ip1[which(ip1>1.0)] = 1.0
df=data.frame(id=ids,input1_observed=io1,input1_predicted=ip1,input2_predicted=ip2)
fout='test-data/nri_test1.xls'
write.table(df,file=fout, quote=FALSE, sep="\t",row.names=F)
# planemo test --job_output_files /home/rlazarus/tmp --test_output /home/rlazarus/tmp/startest/foo.html --update_test_data --galaxy_root /home/rlazarus/galaxy rg_nri.xml
