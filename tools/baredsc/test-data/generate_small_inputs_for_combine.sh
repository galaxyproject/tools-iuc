for nnorm in 1 2; do
    baredSC_1d --input test-data/nih3t3_generated_2d_2.txt --geneColName 0.5_0_0_0.5_x --nnorm ${nnorm} --output test-data/small_${nnorm}gauss --figure test-data/small_${nnorm}gauss.png --nx 10 --prettyBins 100 --nsampMCMC 20000 --force
done
for nnorm in 1 2; do
    baredSC_2d  --input test-data/nih3t3_generated_2d_2.txt --geneXColName '0.5_0_0_0.5_x' --geneYColName '0.5_0_0_0.5_y' --nnorm ${nnorm} --output test-data/2d_small_${nnorm}gauss --figure test-data/2d_small_${nnorm}gauss.png --nx 10 --ny 12 --prettyBinsx 50 --prettyBinsy 50 --nsampMCMC 20000 --force
done
