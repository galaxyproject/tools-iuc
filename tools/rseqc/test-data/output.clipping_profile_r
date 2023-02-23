pdf("output.clipping_profile.pdf")
read_pos=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)
clip_count=c(16.0,12.0,11.0,8.0,7.0,6.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0,1.0,1.0,2.0,3.0,4.0,4.0)
nonclip_count= 40 - clip_count
plot(read_pos, nonclip_count*100/(clip_count+nonclip_count),col="blue",main="clipping profile",xlab="Position of read",ylab="Non-clipped %",type="b")
dev.off()
