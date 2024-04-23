for fname in tt_mlfm*.txt; do
    nfile=$(printf "%04d" $(echo $fname | sed 's/tt_mlfm_c//' | sed 's/\.txt//' | bc))
    out="tt_mittleff_c${nfile}.txt"

    sed 's/ \+/ /g;s/^ /+/g;s/\([0-9]\) \([0-9]\)/\1 +\2/g' $fname > $out
    rm $fname
done

