for i in M N Eff_d5 do
        multiBigwigSummary bins --bwfiles Philip_ATAC/${i}_r*.bw --region chr4 -bs 10 --outRawCounts 
${i}_chr4.tab -o tmp.npz -p 5
        bedtools sort -i ${i}_chr4.tab > ${i}_chr4_sort.tab
        python averaging_bigwigs.py --infile ${i}_chr4_sort.tab --outfile ${i}_chr4 --summaryType median 
--outfileType bigWig --oriBigWig Philip_ATAC/Eff_d7_rep1.bw
        gzip ${i}*tab &
done

