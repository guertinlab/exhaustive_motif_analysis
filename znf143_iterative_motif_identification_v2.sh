#!/bin/bash
#SBATCH --job-name=znf143_iterative_motif_identification
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load ucsc_genome
module load bedtools
module load samtools
module load meme/5.5.9

#Initialize variables

#What factor we ChIP'ed
factor=ZNF143
#Tool to use for initial round 
tool=meme
#Objective function for initial round of MEME
objfun=classic
#Start with round 1 for intial single motif identification
round=1
#Start the counter with the initial single motif
count=1
#Initial windows are 100bp
slop=50
#Motif should be present in at least 100/1000 peaks
minsites=100
#Cores for parallel processing
ncore=32
#Threshold for MEME e-value
e_thresh=.05
#Threshold for TOMTOM q-value
tomtom_thresh=.01
#Goal proportion of peaks remaining
goal=0.05
#Peaks we're starting with 
all_peaks=${factor}_ChIP_summit_100window
top1000=${factor}_ChIP_top1000_summit_100window
control=${factor}_mock_summit_100window.fasta
#Where things are
home=/scratch/afrutos/GSE266489_control3
blacklist=$home/hs1/hs1-blacklist.v2.bed
sizes=$home/hs1/hs1.chrom.sizes.txt
genome=$home/hs1/hs1.fa
hs1_bkgrnd=$home/hs1/hs1_bkgrnd.txt
dir=$home/Motif_analysis_automation
#Where to put things
meme_out=$dir/MEME/${factor}_${tool}_${objfun}_${round}
hits=MAST_${factor}_${tool}_${objfun}_PSWM_in_peaks_${round}
motifs=$dir/MEME/${factor}_motifs.txt
without=${factor}_without_${tool}_${objfun}_motifs_${round}

#Move to relevant directory
cd $dir/Peaks

#Print starting peak count for later comparison
peaks_remaining=$(wc -l $all_peaks.bed | awk '{print $1}')
peaks_starting=$peaks_remaining
echo "$peaks_remaining peaks to account for"
echo "Starting round $round, using $tool with objective function $objfun"

#Find initial motif to use as reference (can mannually confirm it is the expected motif before continuing)
meme -nostatus -p $ncore -oc $meme_out -nmotifs 1 -minsites $minsites -objfun $objfun -minw 12 -maxw 15 -searchsize 0 -revcomp -dna -bfile $hs1_bkgrnd $top1000.fasta

#Start minimal motif database file
#Grab MEME version, ALPHABET, strands, background letter frequerncy and next line (data)
{ grep "MEME version" $meme_out/$tool.txt; printf '\n'; grep "ALPHABET" $meme_out/$tool.txt; printf '\n'; grep -E "strands" $meme_out/$tool.txt; printf '\n'; grep -A1 "Background letter" $meme_out/$tool.txt; printf '\n'; } > $motifs

#Get motif name
motif=$(grep "position-specific probability matrix" $meme_out/$tool.txt | awk '{print $2}')

#Add the motif name and probability matrix to the minimal motif database file
#Starts printing at the line containing probability matrix and continues until it reaches a line starting with a hyphen
{ echo "MOTIF $motif"; awk '/^letter-probability matrix/,/^-/' $meme_out/$tool.txt | sed '$d'; printf '\n'; } >> $motifs

#Take peaks with instances of this motif (only in the 1000 peaks used for MEME)
awk -v str="$motif" '$0 ~ str && $0 ~ "sites sorted" {skip=1; count = 0; next} skip && count < 3 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt > $hits.txt

#Optionally add smaller motifs (e.g., each half of the initial ZNF143 motif)
#Workflow as described in comments below
echo "Finding submotifs"
meme_out=${meme_out}_submotifs
meme -nostatus -p $ncore -oc $meme_out -evt $e_thresh -minsites $minsites -objfun $objfun -minw 6 -maxw 8 -searchsize 0 -revcomp -dna -bfile $hs1_bkgrnd $top1000.fasta
tomtom_out=$dir/TOMTOM/TOMTOM_${factor}_${tool}_${objfun}_${round}
tomtom -verbosity 1 -eps -thresh $tomtom_thresh -oc $tomtom_out $meme_out/$tool.txt $motifs
col_index=$(head -1 $tomtom_out/tomtom.tsv | tr '\t' '\n' | grep -n "^q-value$" | cut -d: -f1) 
rows=$(awk -F'\t' -v col="$col_index" '{ if ($col == "") exit; print $col }' $tomtom_out/tomtom.tsv | wc -l)
previous_motif=NULL
for ((i=2; i<=$rows; i++))
do
  motif=$(awk -v i="$i" 'NR==i {print $1}' $tomtom_out/tomtom.tsv)
  if [[ "$motif" != "$previous_motif" ]]
  then
    previous_motif=$motif
    ((count++))
    { echo "MOTIF $motif"; awk -v str="$motif" '$0 ~ str && $0 ~ "position-specific probability matrix" {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt; printf '\n'; } >> $motifs
    awk -v str="$motif" '$0 ~ str && $0 ~ "sites sorted" {skip=1; count = 0; next} skip && count < 3 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt >> $hits.txt
  fi
done

#Convert to bed file
#Manually narrowing hits to the peak summits because instances with 9 digit coordinates are cut off, making the end coordinate less than the beginning, and to avoid subtracting overlapping summit windows
awk -v slop="$slop" 'BEGIN {OFS="\t"} {split($1,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $2, $3, $4}' $hits.txt > $hits.bed

#Remove peaks with the motif(s)
slopBed -b -$slop -i $all_peaks.bed -g $sizes | intersectBed -v -a stdin -b $hits.bed > $without.bed

#Print an update
peaks_remaining=$(wc -l $without.bed | awk '{print $1}')
echo "$count motif(s) found so far"
echo "$peaks_remaining peaks remaining"

#Loop through tools, objective functions, and motif widths, repeating loop until all methods strike out
tool_strikes=0
while (( tool_strikes != 3 ))
do
  for tool in streme meme glam2
  do
    if (( tool_strikes == 3 ))
    then
      break
    fi
    #Give starting strikes based on how many objfun's each tool can use
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
      for objfun in cd de classic NA
      do 
        #Keep going until either no significant motifs or no matches (breaks within if statements below)
        while (( objfun_strikes != 3 ))
        do
          #Skip combos that don't make sense
          if [[ "$objfun" == "classic" && "$tool" == "streme" || "$objfun" != "NA" && "$tool" == "glam2" || "$objfun" == "NA" && "$tool" != "glam2" ]] 
          then
            break
          fi
          #Update variables
          ((round++))
          minsites=$(( $peaks_remaining / 10 ))
          meme_out=$dir/MEME/${factor}_${tool}_${objfun}_${round}
          tomtom_out=$dir/TOMTOM/TOMTOM_${factor}_${tool}_${objfun}_${round}
          new_motifs=$dir/MEME/${factor}_${tool}_${objfun}_motifs_${round}.txt
          hits=${factor}_${tool}_${objfun}_PSWM_in_peaks_${round}
          #Print an update
          echo "Starting round $round, using $tool with objective function $objfun"
          #Adjust summit window widths
          if [[ "$objfun" == "cd" ]]
          then
            slop=250
          else
            slop=50
          fi
          slopBed -b $slop -i $without.bed -g $sizes | sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed
          fastaFromBed -fi $genome -bed $without.bed -fo $without.fasta
          #Identify motifs
          if [[ "$objfun" == "de" && "$tool" == "streme" ]]
          then
            streme --verbosity 1 --oc $meme_out --thresh $e_thresh --evalue --dna --objfun $objfun --minw 8 --maxw 15 --bfile $hs1_bkgrnd --n $control --p $without.fasta
          elif [[ "$objfun" != "de" && "$tool" == "streme" ]]
          then
            streme --verbosity 1 --oc $meme_out --thresh $e_thresh --evalue --dna --objfun $objfun --minw 8 --maxw 15 --bfile $hs1_bkgrnd --p $without.fasta
          elif [[ "$objfun" == "de" && "$tool" == "meme" ]]
          then
            meme -nostatus -p $ncore -oc $meme_out -brief 100000 -evt $e_thresh -minsites $minsites -objfun $objfun -minw 8 -maxw 12 -wg 0 -ws 0 -searchsize 0 -revcomp -dna -bfile $hs1_bkgrnd -neg $control $without.fasta
          elif [[ "$objfun" != "de" && "$tool" == "meme" ]]
          then
            meme -nostatus -p $ncore -oc $meme_out -brief 100000 -evt $e_thresh -minsites $minsites -objfun $objfun -minw 8 -maxw 12 -wg 0 -ws 0 -searchsize 0 -revcomp -dna -bfile $hs1_bkgrnd $without.fasta
          elif [[ "$tool" == "glam2" ]]
          then
            glam2 -Q -n 1000000 -z $minsites -a 8 -b 30 -w 16 -O $meme_out -2 n $without.fasta
          fi
          #Break if no motifs meet threshold
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
              slopBed -b -$slop -i $without.bed -g $sizes | sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed
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
              slopBed -b -$slop -i $without.bed -g $sizes | sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed
              break
            fi
          fi
          if [[ "$tool" == "streme" ]]
          then
            #Remove last 3 motifs, as they are above the e-value threshold
            #Get the number of motifs
            total=$(grep 'STREME-' $meme_out/$tool.txt | tail -n 1 | awk -F'-' '{ print $3 }')
            #Get the 3rd to last motif number
            total=$(( total - 2 ))
            #Start by finding the 3rd to last motif and remove lines until the 3rd blank line
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
          #Match with TOMTOM to current motif database
          if [[ "$tool" == "glam2" ]]
          then
            #Only the first glam2 motif is relevant
            tomtom -verbosity 1 -m 1 -eps -thresh $tomtom_thresh -oc $tomtom_out $meme_out/glam2.meme $motifs
          else
            #Match all identified motifs
            tomtom -verbosity 1 -eps -thresh $tomtom_thresh -oc $tomtom_out $meme_out/$tool.txt $motifs
          fi
          #Find q-value column
          col_index=$(head -1 $tomtom_out/tomtom.tsv | tr '\t' '\n' | grep -n "^q-value$" | cut -d: -f1) 
          #Count how many significant q-values there are (+1 for column name)
          rows=$(awk -F'\t' -v col="$col_index" '{ if ($col == "") exit; print $col }' $tomtom_out/tomtom.tsv | wc -l)
          #Break if no significant matches
          if [[ "$rows" -eq "1" ]]
          then
            echo "No more matches to motif database with q-value < $tomtom_thresh"
            ((objfun_strikes++))
            if (( objfun_strikes == 3 )) 
            then
              ((tool_strikes++))
            fi
            slopBed -b -$slop -i $without.bed -g $sizes | sort -k1,1 -k2,2n > tmp.bed && mv tmp.bed $without.bed
            break
          fi
          #Add each matching motif to a file for this round as well as to the growing database
          #Also pull instances from the STREME/MEME/GLAM2 output file
          if [[ "$tool" == "glam2" ]]
          then
            #Only adding the one motif
            { grep "MEME version" $meme_out/glam2.meme; printf '\n'; grep "ALPHABET" $meme_out/glam2.meme; printf '\n'; grep -E "strands" $meme_out/glam2.meme; printf '\n'; grep -A1 "Background letter" $meme_out/glam2.meme; printf '\n'; } > $new_motifs
            motif=GLAM2_Motif_${round}
            #Search for the first line containing the probability matrix and then print starting 1 line later until reaching the next instance of a blank line
            { echo "$motif"; awk '/letter-probability matrix/ {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^$/ { exit } skip { print }' $meme_out/glam2.meme; printf '\n'; } >> $motifs
            { echo "$motif"; awk '/letter-probability matrix/ {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^$/ { exit } skip { print }' $meme_out/glam2.meme; printf '\n'; } >> $new_motifs
            #Add instances to the hit list
            #Search for the first line starting with "Score:", then print starting 3 lines later until reaching the second instance of a blank line
            awk '/^Score: / { found=1; skip=2; next } found && skip>0 { skip--; next } found && /^$/ { blank++; if (blank==1) exit; next } found { print }' $meme_out/glam2.txt > $hits.txt
            #Convert hit list to bed file
            awk -v slop="$slop" 'BEGIN {OFS="\t"} {split($1,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $5, $6}' $hits.txt > $hits.bed
          else
            #Loop through all matching motifs
            { grep "MEME version" $meme_out/$tool.txt; printf '\n'; grep "ALPHABET" $meme_out/$tool.txt; printf '\n'; grep -E "strands" $meme_out/$tool.txt; printf '\n'; grep -A1 "Background letter" $meme_out/$tool.txt; printf '\n'; } > $new_motifs
            > $hits.txt
            previous_motif=NULL
            for ((i=2; i<=$rows; i++))
            do
              if [[ "$tool" == "streme" ]]
              then
                motif=$(awk -v i="$i" 'NR==i {print $1}' $tomtom_out/tomtom.tsv | awk -F'-' '{ print $2 }')
              elif [[ "$tool" == "meme" ]]
              then
                motif=$(awk -v i="$i" 'NR==i {print $1}' $tomtom_out/tomtom.tsv)
              fi
              #Only if the motif isn't already there (i.e., TOMTOM matches the same new motif to multiple existing motifs in the database)
              if [[ "$motif" != "$previous_motif" ]]
              then
                #Update $previous_motif to compare to the next motif
                previous_motif=$motif
                #Update the motif count
                ((count++))
                if [[ "$tool" == "streme" ]]
                then
                  #Search for the line containing the motif and then print until reaching the next instance of a blank line
                  { awk -v str="$motif" '$0 ~ str {p=1} p && /^$/ {p=0} p' $meme_out/$tool.txt; printf '\n'; } >> $motifs
                  { awk -v str="$motif" '$0 ~ str {p=1} p && /^$/ {p=0} p' $meme_out/$tool.txt; printf '\n'; } >> $new_motifs
                  #Add instances to the hit list
                  grep "$motif" $meme_out/sequences.tsv | grep "tp" >> $hits.txt
                elif [[ "$tool" == "meme" ]]
                then
                  #Search for the line containing the motif and the probability matrix and then print starting 2 lines later until reaching the next instance of a line starting with a hyphen
                  { echo "MOTIF $motif"; awk -v str="$motif" '$0 ~ str && $0 ~ "position-specific probability matrix" {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt; printf '\n'; } >> $motifs
                  { echo "MOTIF $motif"; awk -v str="$motif" '$0 ~ str && $0 ~ "position-specific probability matrix" {skip=1; count = 0; next} skip && count < 1 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt; printf '\n'; } >> $new_motifs
                  #Add instances to the hit list
                  awk -v str="$motif" '$0 ~ str && $0 ~ "sites sorted" {skip=1; count = 0; next} skip && count < 3 { count++; next } skip && /^-/ { exit } skip { print }' $meme_out/$tool.txt >> $hits.txt
                fi
              fi
            done
            #Convert hit list to a bed file
            if [[ "$tool" == "streme" ]] 
            then
              awk -v slop="$slop" 'BEGIN {OFS="\t"} {split($4,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $5, $3, $1}' $hits.txt > $hits.bed
            elif [[ "$tool" == "meme" ]] 
            then
              awk -v slop="$slop" 'BEGIN {OFS="\t"} {split($1,a,"[:-]"); print a[1], a[2]+slop, a[2]+slop+1, $2, $3, $4}' $hits.txt > $hits.bed
            fi
          fi
          #Remove instances of these motifs
          all_peaks=$without
          without=${factor}_without_${tool}_${objfun}_motifs_${round}
          slopBed -b -$slop -i $all_peaks.bed -g $sizes | intersectBed -v -a stdin -b $hits.bed > $without.bed
          #Print an update
          peaks_remaining=$(wc -l $without.bed | awk '{print $1}')
          echo "$count motifs found so far"
          echo "$peaks_remaining peaks remaining"
          #Reset strikes
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
          #End if mission accomplished
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
