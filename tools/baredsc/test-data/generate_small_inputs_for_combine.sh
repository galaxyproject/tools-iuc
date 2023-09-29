for nnorm in 1 2; do
    baredSC_1d --input test-data/nih3t3_generated_2d_2.txt --geneColName 0.5_0_0_0.5_x --nnorm ${nnorm} --output test-data/small_${nnorm}gauss --figure test-data/small_${nnorm}gauss.png --nx 10 --prettyBins 100 --nsampMCMC 20000 --force
done
