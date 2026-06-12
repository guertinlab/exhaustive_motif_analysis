# Input that needs to be changed
home=/Users/guertinlab/Desktop/iterative_meme
blacklist=$home/hs1/hs1-blacklist.v2.bed
sizes=$home/hs1/hs1.chrom.sizes.txt
genome=$home/hs1/hs1.fa
hs1_bkgrnd=$home/hs1/hs1_bkgrnd.txt
dir=$home/Motif_analysis_automation
input_summits=$home/ZNF143_ChIP_summits.bed
mock_summits=$home/ZNF143_mock_summits.bed


mkdir -p $dir/TOMTOM/

#Initialize variables
factor=ZNF143
input_peak_number=1000

tool=meme
objfun=classic
round=1
count=1
slop=50

minsites=$(echo "$input_peak_number * 0.1" | bc)

ncore=1
e_thresh=.05
tomtom_thresh=0.05
# Maximum MAST p-value threshold; use the smaller of this and the dynamic threshold
mast_max_thresh=1e-4
goal=0.05


#Generated from previous variables
all_peaks=${factor}_ChIP_summit_100window
topPeaks=${factor}_ChIP_top${input_peak_number}_summit_100window
control=${factor}_mock_summit_100window
topMockPeaks=${factor}_mock_top${input_peak_number}_summit_100window


#Where to put things
mkdir -p $dir/MEME
meme_out=$dir/MEME/${factor}_${tool}_${objfun}_${round}
hits=MAST_${factor}_${tool}_${objfun}_PSWM_in_peaks_${round}
motifs=$dir/MEME/${factor}_motifs.txt
without=${factor}_without_${tool}_${objfun}_motifs_${round}

append_first_motif_sliding_windows () {
  local input_motifs=$1
  local output_motifs=$2
  local win=${3:-6}

  # Append only sliding-window motifs derived from the first motif.
  # Do not rewrite the TOMTOM database and do not re-add the original first motif.
  awk -v win="$win" '
  BEGIN {
    in_first = 0
    collecting = 0
    done = 0
  }

  done == 0 && /^MOTIF / && in_first == 0 {
    in_first = 1
    motif_name = $2
    next
  }

  in_first == 1 && /^letter-probability matrix:/ {
    matrix_header = $0

    if (match($0, /w=[[:space:]]*[0-9]+/)) {
      width_text = substr($0, RSTART, RLENGTH)
      gsub(/[^0-9]/, "", width_text)
      width = width_text + 0
    } else {
      print "ERROR: could not find motif width in matrix header" > "/dev/stderr"
      exit 1
    }

    collecting = 1
    row_count = 0
    next
  }

  collecting == 1 && NF == 4 {
    row_count++
    matrix[row_count] = $0

    if (row_count == width) {
      for (start = 1; start <= width - win + 1; start++) {
        end = start + win - 1

        printf "MOTIF %s_win%02d_%02d\n", motif_name, start, end

        new_header = matrix_header
        sub(/w=[[:space:]]*[0-9]+/, "w= " win, new_header)
        print new_header

        for (i = start; i <= end; i++) {
          print matrix[i]
        }

        print ""
      }

      done = 1
      exit
    }

    next
  }
  ' "$input_motifs" >> "$output_motifs"
}

tomtom_motifs=$dir/MEME/${factor}_tomtom_motifs.txt


#Remove junk and get a window centered on the summit
cd $dir/Peaks

grep -v "random" ${input_summits} | \
  grep -v "chrUn" | \
  grep -v "chrEBV" | \
  grep -v "chrM" | \
  grep -v "alt" | \
  intersectBed -v -a - -b $blacklist | \
  slopBed -b 50 -i stdin -g $sizes | \
  sort -k1,1 -k2,2n > ${all_peaks}.bed

fastaFromBed -fi $genome -bed ${all_peaks}.bed -fo ${all_peaks}.fasta


#Sort out the top X peaks
sort -nrk5,5 ${all_peaks}.bed | head -n $input_peak_number > ${topPeaks}.bed
fastaFromBed -fi $genome -bed ${topPeaks}.bed -fo ${topPeaks}.fasta


#Repeat for mock peaks
grep -v "random" ${mock_summits} | \
  grep -v "chrUn" | \
  grep -v "chrEBV" | \
  grep -v "chrM" | \
  grep -v "alt" | \
  intersectBed -v -a - -b $blacklist | \
  slopBed -b 50 -i stdin -g $sizes | \
  sort -k1,1 -k2,2n > ${control}.bed

fastaFromBed -fi $genome -bed ${control}.bed -fo ${control}.fasta


#Sort out the top X mock peaks
sort -nrk5,5 ${control}.bed | head -n $input_peak_number > ${topMockPeaks}.bed
fastaFromBed -fi $genome -bed ${topMockPeaks}.bed -fo ${topMockPeaks}.fasta


#Print starting peak count for later comparison
peaks_remaining=$(wc -l $all_peaks.bed | awk '{print $1}')
peaks_starting=$peaks_remaining

echo "$peaks_remaining peaks to account for"
echo "Starting round $round, using $tool with objective function $objfun"


###############################################################################
# ROUND 1: Find initial motif
###############################################################################

meme -p $ncore \
  -oc $meme_out \
  -nmotifs 1 \
  -minsites $minsites \
  -objfun $objfun \
  -minw 8 \
  -maxw 18 \
  -searchsize 0 \
  -csites 20000 \
  -revcomp \
  -dna \
  -bfile $hs1_bkgrnd \
  $topPeaks.fasta


#Start minimal motif database file
{
  grep "MEME version" $meme_out/$tool.txt
  printf '\n'
  grep "ALPHABET" $meme_out/$tool.txt
  printf '\n'
  grep -E "strands" $meme_out/$tool.txt
  printf '\n'
  grep -A1 "Background letter" $meme_out/$tool.txt
  printf '\n'
} > $motifs

cp "$motifs" "$tomtom_motifs"

#Get motif name
motif=$(grep "position-specific probability matrix" $meme_out/$tool.txt | awk '{print $2}')


#Add motif to the growing motif database
{
  echo "MOTIF $motif"
  awk '/^letter-probability matrix/,/^-/' $meme_out/$tool.txt | sed '$d'
  printf '\n'
} | tee -a "$motifs" >> "$tomtom_motifs"

append_first_motif_sliding_windows "$motifs" "$tomtom_motifs" 6

###############################################################################
# ROUND 1: Use MAST to find/remove peaks with the initial motif
###############################################################################

meme_sites=$(
  awk '/sites =/ {
    for(i=1; i<=NF; i++)
      if ($i ~ /^[0-9]+$/)
        print int($i * 1.0)
  }' $meme_out/$tool.txt | awk 'NR==2'
)

if [[ -z "$meme_sites" || "$meme_sites" -lt 1 ]]
then
  meme_sites=1
fi

p=$(
  mast -mt 1 \
    -nostatus \
    -hit_list \
    -best \
    -bfile $hs1_bkgrnd \
    -m "$motif" \
    $motifs \
    $topPeaks.fasta | \
    awk 'NR>3 {print last} {last=$0}' | \
    sort -k 8,8g | \
    head -n $meme_sites | \
    tail -n 1 | \
    awk '{print $8}'
)

# Use the smaller, more stringent threshold: dynamic MAST threshold or mast_max_thresh
p_dynamic=$p
p=$(
  awk -v dyn="$p_dynamic" -v cap="$mast_max_thresh" 'BEGIN {
    if (dyn == "") exit
    if (dyn < cap) print dyn
    else print cap
  }'
)

echo "Round 1 dynamic MAST p-value threshold:"
echo "$p_dynamic"
echo "Round 1 capped MAST p-value threshold:"
echo "$p"

mast -mt $p \
  -nostatus \
  -hit_list \
  -best \
  -bfile $hs1_bkgrnd \
  -m "$motif" \
  $motifs \
  $all_peaks.fasta > $hits.txt

echo "Round 1 MAST hits:"
wc -l $hits.txt
head $hits.txt


#Convert MAST hits to BED
awk -v slop="$slop" 'BEGIN {OFS="\t"} NR > 2 {split($1,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $2, $7, $8}' $hits.txt > $hits.bed

echo "Round 1 BED hits:"
wc -l $hits.bed
head $hits.bed


#Remove peaks with the motif
slopBed -b -$slop -i $all_peaks.bed -g $sizes | \
  intersectBed -v -a stdin -b $hits.bed > $without.bed


peaks_remaining=$(wc -l $without.bed | awk '{print $1}')
echo "$count motif(s) found so far"
echo "$peaks_remaining peaks remaining"


#Build top remaining peaks for next round
without_top=${without}_top${input_peak_number}

sort -nrk5,5 $without.bed | \
  head -n $input_peak_number | \
  slopBed -b 50 -i stdin -g $sizes | \
  sort -k1,1 -k2,2n > ${without_top}.bed

fastaFromBed -fi $genome -bed ${without_top}.bed -fo ${without_top}.fasta


###############################################################################
# SUBSEQUENT ROUNDS
###############################################################################

tool_strikes=0

while (( tool_strikes != 3 ))
do
  for tool in meme streme glam2
  do
    if (( tool_strikes == 3 ))
    then
      break
    fi

    if [[ "$tool" == "meme" ]]
    then
      objfun_strikes=0
    elif [[ "$tool" == "streme" ]]
    then
      objfun_strikes=1
    elif [[ "$tool" == "glam2" ]]
    then
      objfun_strikes=2
    fi

    while (( objfun_strikes != 3 ))
    do
      for objfun in classic cd de NA
      do
        while (( objfun_strikes != 3 ))
        do

          #Skip combos that do not make sense
          if [[ "$objfun" == "classic" && "$tool" == "streme" || "$objfun" != "NA" && "$tool" == "glam2" || "$objfun" == "NA" && "$tool" != "glam2" ]]
          then
            break
          fi

          ((round++))

          meme_out=$dir/MEME/${factor}_${tool}_${objfun}_${round}
          tomtom_out=$dir/TOMTOM/TOMTOM_${factor}_${tool}_${objfun}_${round}
          new_motifs=$dir/MEME/${factor}_${tool}_${objfun}_motifs_${round}.txt
          hits=${factor}_${tool}_${objfun}_PSWM_in_peaks_${round}

          echo "Starting round $round, using $tool with objective function $objfun"

          if [[ "$objfun" == "cd" ]]
          then
            slop=50
            slop_for_discovery=$((slop * 3))
          else
            slop=50
            slop_for_discovery=$slop
          fi

          #Expand current remaining peaks for motif discovery
          slopBed -b $slop_for_discovery -i $without.bed -g $sizes | \
            sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed

          fastaFromBed -fi $genome -bed $without.bed -fo $without.fasta

          #Rebuild top remaining peaks using the current discovery window.
          #Important: do not slop again here because $without.bed is already expanded.
          without_top=${without}_top${input_peak_number}

          sort -nrk5,5 $without.bed | \
            head -n $input_peak_number | \
            sort -k1,1 -k2,2n > ${without_top}.bed

          fastaFromBed -fi $genome -bed ${without_top}.bed -fo ${without_top}.fasta

          echo "Current discovery slop:"
          echo "$slop_for_discovery"

          echo "Current top discovery windows:"
          awk 'BEGIN {OFS="\t"} NR <= 5 {print $1, $2, $3, $3-$2}' ${without_top}.bed


          #######################################################################
          # Discover motifs
          #######################################################################

          if [[ "$objfun" == "de" && "$tool" == "streme" ]]
          then
            streme --verbosity 1 \
              --oc $meme_out \
              --thresh $e_thresh \
              --evalue \
              --dna \
              --objfun $objfun \
              --minw 6 \
              --maxw 15 \
              --bfile $hs1_bkgrnd \
              --n $control.fasta \
              --p $without.fasta

          elif [[ "$objfun" != "de" && "$tool" == "streme" ]]
          then
            streme --verbosity 1 \
              --oc $meme_out \
              --thresh $e_thresh \
              --evalue \
              --dna \
              --objfun $objfun \
              --minw 6 \
              --maxw 15 \
              --bfile $hs1_bkgrnd \
              --p $without.fasta

          elif [[ "$objfun" == "de" && "$tool" == "meme" ]]
          then
            meme -p $ncore \
              -oc $meme_out \
              -brief 100000 \
              -evt $e_thresh \
              -minsites $minsites \
              -objfun $objfun \
              -minw 8 \
              -maxw 18 \
              -wg 0 \
              -ws 0 \
              -searchsize 0 \
	      -csites 20000 \
              -revcomp \
              -dna \
              -bfile $hs1_bkgrnd \
              -neg $control.fasta \
              $without.fasta

          elif [[ "$objfun" == "classic" && "$tool" == "meme" ]]
          then
            meme -p $ncore \
              -oc $meme_out \
              -brief 100000 \
              -evt $e_thresh \
              -minsites $minsites \
	      -searchsize 0 \
	      -csites 20000 \
              -objfun $objfun \
              -minw 8 \
              -maxw 18 \
              -wg 0 \
              -ws 0 \
              -searchsize 0 \
              -revcomp \
              -dna \
              -bfile $hs1_bkgrnd \
              ${without_top}.fasta

          elif [[ "$objfun" != "classic" && "$objfun" != "de" && "$tool" == "meme" ]]
          then
            meme -p $ncore \
              -oc $meme_out \
              -brief 100000 \
              -evt $e_thresh \
              -minsites $minsites \
              -objfun $objfun \
              -minw 8 \
              -maxw 18 \
              -wg 0 \
              -ws 0 \
              -searchsize 0 \
              -revcomp \
              -dna \
              -bfile $hs1_bkgrnd \
              ${without_top}.fasta

          elif [[ "$tool" == "glam2" ]]
          then
            glam2 -Q \
              -n 1000000 \
              -z $minsites \
              -a 8 \
              -b 30 \
              -w 16 \
              -O $meme_out \
              -2 n \
              $without.fasta
          fi


          #######################################################################
          # Break if no motifs meet discovery threshold
          #######################################################################

          if [[ "$tool" == "streme" ]]
          then
            if ! grep -q "STREME-4" $meme_out/$tool.txt
            then
              echo "No more motifs with e-value < $e_thresh"

              ((objfun_strikes++))
              if (( objfun_strikes == 3 ))
              then
                ((tool_strikes++))
              fi

              slopBed -b -$slop_for_discovery -i $without.bed -g $sizes | \
                sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed

              break
            fi

          elif [[ "$tool" == "meme" ]]
          then
            if ! grep -q "MEME-1" $meme_out/$tool.txt
            then
              echo "No more motifs with e-value < $e_thresh"

              ((objfun_strikes++))
              if (( objfun_strikes == 3 ))
              then
                ((tool_strikes++))
              fi

              slopBed -b -$slop_for_discovery -i $without.bed -g $sizes | \
                sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed

              break
            fi
          fi


          #######################################################################
          # For STREME, remove last 3 motifs because they are above threshold
          #######################################################################

          if [[ "$tool" == "streme" ]]
          then
            total=$(grep 'STREME-' $meme_out/$tool.txt | tail -n 1 | awk -F'-' '{ print $3 }')
            total=$(( total - 2 ))

            awk -v pattern="STREME-$total" '
            {
              if (state == 0) {
                if ($0 ~ pattern) {
                  state = 1
                  blank_count = 0
                  next
                }
                print
              } else if (state == 1) {
                if ($0 == "") {
                  blank_count++
                  if (blank_count >= 3) {
                    state = 0
                  }
                }
              }
            }' $meme_out/$tool.txt > tmp.txt && mv tmp.txt $meme_out/$tool.txt
          fi


          #######################################################################
          # TOMTOM: match discovered motifs to current motif database
          #######################################################################

          if [[ "$tool" == "glam2" ]]
          then
            tomtom -verbosity 1 \
              -m 1 \
              -eps \
              -thresh $tomtom_thresh \
              -oc $tomtom_out \
              $meme_out/glam2.meme \
              $tomtom_motifs
          else
            tomtom -verbosity 1 \
              -eps \
              -thresh $tomtom_thresh \
              -oc $tomtom_out \
              $meme_out/$tool.txt \
              $tomtom_motifs
          fi


          #######################################################################
          # Select only one motif: lowest TOMTOM E-value among passing q-values
          #######################################################################

          best_tomtom_query=$(
            awk -F'\t' -v qthresh="$tomtom_thresh" '
              NR == 1 {
                for (i = 1; i <= NF; i++) {
                  if ($i == "Query_ID") qid = i
                  if ($i == "E-value") eval = i
                  if ($i == "q-value") qval = i
                }
                next
              }
              /^#/ { next }
              qid == "" || eval == "" || qval == "" { next }
              $qid == "" || $eval == "" || $qval == "" { next }
              $qval <= qthresh {
                if (best == "" || $eval < best_eval) {
                  best = $qid
                  best_eval = $eval
                }
              }
              END {
                if (best != "") print best
              }
            ' "$tomtom_out/tomtom.tsv"
          )

          if [[ -z "$best_tomtom_query" ]]
          then
            echo "No more matches to motif database with q-value < $tomtom_thresh"

            ((objfun_strikes++))
            if (( objfun_strikes == 3 ))
            then
              ((tool_strikes++))
            fi

            slopBed -b -$slop_for_discovery -i $without.bed -g $sizes | \
              sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed

            break
          fi


          #######################################################################
          # Add only the selected motif to $motifs and $new_motifs
          #######################################################################

          discovery_hits=${hits}_discovery_hits.txt

          if [[ "$tool" == "glam2" ]]
          then
            motif=GLAM2_Motif_${round}

            {
              grep "MEME version" $meme_out/glam2.meme
              printf '\n'
              grep "ALPHABET" $meme_out/glam2.meme
              printf '\n'
              grep -E "strands" $meme_out/glam2.meme
              printf '\n'
              grep -A1 "Background letter" $meme_out/glam2.meme
              printf '\n'
            } > $new_motifs

            {
              echo "MOTIF $motif"
              awk '/letter-probability matrix/ {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^$/ { exit } skip { print }' $meme_out/glam2.meme
              printf '\n'
            } | tee -a "$motifs" >> "$tomtom_motifs"

            {
              echo "MOTIF $motif"
              awk '/letter-probability matrix/ {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^$/ { exit } skip { print }' $meme_out/glam2.meme
              printf '\n'
            } >> $new_motifs

            awk '/^Score: / { found=1; skip=2; next } found && skip>0 { skip--; next } found && /^$/ { blank++; if (blank==1) exit; next } found { print }' $meme_out/glam2.txt > $discovery_hits

          else
            source_motif="$best_tomtom_query"
            motif=${tool}_${objfun}_${round}_${source_motif}

            {
              grep "MEME version" $meme_out/$tool.txt
              printf '\n'
              grep "ALPHABET" $meme_out/$tool.txt
              printf '\n'
              grep -E "strands" $meme_out/$tool.txt
              printf '\n'
              grep -A1 "Background letter" $meme_out/$tool.txt
              printf '\n'
            } > $new_motifs

            echo "Best TOMTOM query motif this round by E-value: $source_motif"
            echo "Database motif name for this round: $motif"

            ((count++))

            if [[ "$tool" == "streme" ]]
            then
              {
                awk -v str="$source_motif" -v newname="$motif" '
                  $0 ~ str {
                    sub(/^MOTIF[[:space:]]+[^[:space:]]+/, "MOTIF " newname)
                    p=1
                  }
                  p && /^$/ {p=0}
                  p
                ' $meme_out/$tool.txt
                printf '\n'
              } | tee -a "$motifs" >> "$tomtom_motifs"

              {
                awk -v str="$source_motif" -v newname="$motif" '
                  $0 ~ str {
                    sub(/^MOTIF[[:space:]]+[^[:space:]]+/, "MOTIF " newname)
                    p=1
                  }
                  p && /^$/ {p=0}
                  p
                ' $meme_out/$tool.txt
                printf '\n'
              } >> $new_motifs

              grep "$source_motif" $meme_out/sequences.tsv | grep "tp" > $discovery_hits

            elif [[ "$tool" == "meme" ]]
            then
              {
                echo "MOTIF $motif"
                awk -v str="$source_motif" '$0 ~ str && $0 ~ "position-specific probability matrix" {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt
                printf '\n'
              } | tee -a "$motifs" >> "$tomtom_motifs"

              {
                echo "MOTIF $motif"
                awk -v str="$source_motif" '$0 ~ str && $0 ~ "position-specific probability matrix" {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt
                printf '\n'
              } >> $new_motifs

              awk -v str="$source_motif" '$0 ~ str && $0 ~ "sites sorted" {skip=1; count = 0; next} skip && count < 3 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt > $discovery_hits
            fi
          fi


          #######################################################################
          # Use MAST to scan remaining peaks for the selected motif only
          #######################################################################

          discovery_sites=$(wc -l $discovery_hits | awk '{print $1}')
          mast_sites=$(echo "$discovery_sites * 1.0" | bc | awk '{print int($1)}')

          if [[ -z "$mast_sites" || "$mast_sites" -lt 1 ]]
          then
            mast_sites=1
          fi

          echo "Discovery sites for motif $motif:"
          echo "$discovery_sites"

          echo "Target number of MAST hits for threshold:"
          echo "$mast_sites"

          p=$(
            mast -mt 1 \
              -nostatus \
              -hit_list \
              -best \
              -bfile $hs1_bkgrnd \
              -m "$motif" \
              $new_motifs \
              ${without_top}.fasta | \
              awk 'NR>3 {print last} {last=$0}' | \
              sort -k 8,8g | \
              head -n $mast_sites | \
              tail -n 1 | \
              awk '{print $8}'
          )

          if [[ -z "$p" ]]
          then
            echo "MAST could not determine a p-value threshold for motif $motif"

            ((objfun_strikes++))
            if (( objfun_strikes == 3 ))
            then
              ((tool_strikes++))
            fi

            slopBed -b -$slop_for_discovery -i $without.bed -g $sizes | \
              sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed

            break
          fi

          # Use the smaller, more stringent threshold: dynamic MAST threshold or mast_max_thresh
          p_dynamic=$p
          p=$(
            awk -v dyn="$p_dynamic" -v cap="$mast_max_thresh" 'BEGIN {
              if (dyn == "") exit
              if (dyn < cap) print dyn
              else print cap
            }'
          )

          echo "Dynamic MAST p-value threshold for motif $motif:"
          echo "$p_dynamic"

          echo "Capped MAST p-value threshold for motif $motif:"
          echo "$p"

          mast -mt $p \
            -nostatus \
            -hit_list \
            -best \
            -bfile $hs1_bkgrnd \
            -m "$motif" \
            $new_motifs \
            $without.fasta > $hits.txt

          echo "MAST hits for motif $motif:"
          wc -l $hits.txt

          awk -v slop="$slop" 'BEGIN {OFS="\t"} NR > 2 {split($1,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $2, $7, $8}' $hits.txt > $hits.bed

          echo "BED hits for motif $motif:"
          wc -l $hits.bed


          #######################################################################
          # Remove all remaining peaks containing this motif
          #######################################################################

          all_peaks=$without
          without=${factor}_without_${tool}_${objfun}_motifs_${round}

          echo "Peaks before removing $motif:"
          wc -l $all_peaks.bed

          slopBed -b -$slop_for_discovery -i $all_peaks.bed -g $sizes | \
            intersectBed -v -a stdin -b $hits.bed > $without.bed

          echo "Peaks after removing $motif:"
          wc -l $without.bed

          peaks_remaining=$(wc -l $without.bed | awk '{print $1}')

          echo "$count motifs found so far"
          echo "$peaks_remaining peaks remaining"


          #######################################################################
          # Rebuild top remaining peaks for next round
          #######################################################################

          without_top=${without}_top${input_peak_number}

          sort -nrk5,5 $without.bed | \
            head -n $input_peak_number | \
            slopBed -b $slop_for_discovery -i stdin -g $sizes | \
            sort -k1,1 -k2,2n > ${without_top}.bed

          fastaFromBed -fi $genome -bed ${without_top}.bed -fo ${without_top}.fasta

          echo "Updated top remaining peaks for next round:"
          echo "${without_top}.bed"
          echo "${without_top}.fasta"


          #######################################################################
          # Reset strikes after successful motif removal
          #######################################################################

          tool_strikes=0

          if [[ "$tool" == "meme" ]]
          then
            objfun_strikes=0
          elif [[ "$tool" == "streme" ]]
          then
            objfun_strikes=1
          elif [[ "$tool" == "glam2" ]]
          then
            objfun_strikes=2
          fi


          #######################################################################
          # End if goal reached
          #######################################################################

          if [[ $(echo "$peaks_remaining <= $peaks_starting * $goal" | bc -l) -eq 1 ]]
          then
            echo "< $goal of initial peaks remain: mission accomplished"
            objfun_strikes=3
            tool_strikes=3
          fi

        done
      done
    done
  done
done
