# Divergence rate, 13 Feb 2018
# Calculate pairwise sequence divergence rate beween aligned FASTA sequnces
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# Usage: gawk -v rules=ggsnns -f divergr.awk alignment.fas > output.txt
# Rules: gg  = gap-gap reduce length,
#        gs  = gap-site reduce length
#        nn  = N-N do not counted as similarity
#        ns  = N-site reduce length

BEGIN {
  IGNORECASE = 1
  # default print option
  if (!(printSeq||printTotal)) printDivAll = 1
 }

/^>/ {
  seqHeader = $0;
  if (seq) {
    seqA[seqNum][prevSeqHeader][0] = ""
    gsub(/\n/, "", seq)
    # seqA: seqA[seqNum][seqHeader][baseNum]=base
    seqLen = split(seq, seqA[seqNum][prevSeqHeader], "")
  }
  prevSeqHeader = seqHeader
  seqNum += 1
  seq = ""
}

!/[^ACGTN-]/ { seq = seq $0 }

END {
  if (seq) {
    seqA[seqNum][prevSeqHeader][0] = ""
    gsub(/\n/, "", seq)
    seqLen = split(seq, seqA[seqNum][prevSeqHeader], "")
  }

  if (debug) {
    for (i in seqA) for(j in seqA[i]) for(k in seqA[i][j]) print "seqNum=\"" i \
    "\" seqHeader=\"" j "\" baseNum=\"" k "\" base=\"" seqA[i][j][k] "\""
  }

  if (printSeq) {
    print "Num\tName\tMean\tMedian"
  }

  if (printDivAll) {
    print "Num\tName\tDiv"
  }

  # seqA: seqA[seqNum][seqHeader][baseNum]=base
  for (seqPosQuery = 1; seqPosQuery <= seqNum; seqPosQuery++) {
    for (seqHeaderQuery in seqA[seqPosQuery]) {
      # seqPosQuery + 1 means that the current sequence does not compare with
      # itself
      for (seqPosSearch = seqPosQuery + 1; seqPosSearch <= seqNum;
        seqPosSearch++) {
        for (seqHeaderSearch in seqA[seqPosSearch]) {
          seqLenChange = seqLen
          baseDiff = 0
          for (base in seqA[seqPosQuery][seqHeaderQuery]) {
            baseQuery = seqA[seqPosQuery][seqHeaderQuery][base]
            baseSearch = seqA[seqPosSearch][seqHeaderSearch][base]
            # Rule: gap-gap reduce length
            if (rules ~ /gg/) {
              if ((baseQuery == "-")&&(baseSearch == "-")) {
                seqLenChange--
                if (debug) {
                  print seqPosQuery, seqPosSearch, base, "gg: ", baseQuery, \
                  baseSearch, seqLenChange
                }
                continue
              }
            }
            # Rule: gap-site reduce length
            if (rules ~ /gs/) {
              if (((baseQuery == "-")&&(baseSearch != "-"))||
                  ((baseQuery != "-")&&(baseSearch == "-"))) {
                seqLenChange--
                if (debug) {
                  print seqPosQuery, seqPosSearch, base, "gs: ", baseQuery, \
                  baseSearch, seqLenChange
                }
                continue
              }
            }
            # Rule: N-N do not counted as similarity
            if (rules ~ /nn/) {
              if ((baseQuery == "N")&&(baseSearch == "N")) {
                seqLenChange--
                if (debug) {
                  print seqPosQuery, seqPosSearch, base, "nn: ", baseQuery, \
                  baseSearch, seqLenChange
                }
                continue
              }
            }
            # Rule: N-site reduce length
            if (rules ~ /ns/) {
              if (((baseQuery == "N")&&(baseSearch != "N"))||
                  ((baseQuery != "N")&&(baseSearch == "N"))) {
                seqLenChange--
                if (debug) {
                  print seqPosQuery, seqPosSearch, base, "ns: ", baseQuery, \
                  baseSearch, seqLenChange
                }
                continue
              }
            }
            # Counted differences in sites
            if (baseQuery != baseSearch) {
              baseDiff++
              if (debug) {
                print seqPosQuery, seqPosSearch, base, "diff: ", baseQuery, \
                baseSearch, seqLenChange
              }
            } else if (debug) {
                print seqPosQuery, seqPosSearch, base, "sim: ", baseQuery, \
                baseSearch, seqLenChange
            }
          }
          if (debug) {
            print seqPosQuery, seqPosSearch, "baseDiff " baseDiff, \
            "seqLenChange " seqLenChange "\n"
          }
          if (seqLenChange > 0) {
            divrA[seqPosQuery][seqPosSearch] = baseDiff / seqLenChange
            divrA[seqPosSearch][seqPosQuery] = baseDiff / seqLenChange
          } else {
            divrA[seqPosQuery][seqPosSearch] = 0
            divrA[seqPosSearch][seqPosQuery] = 0
          }
        }
      }
    }
  }

  for (i = 1; i <= seqNum; i++) {
    score = 0
    qLen = asort(divrA[i], divrAS)

    # calculate median divergence for each sequence
    if (qLen % 2) {
      median = divrAS[(qLen + 1) / 2]
    } else {
      lowMed = divrAS[qLen / 2]
      highMed = divrAS[(qLen / 2) + 1]
      median = (lowMed + highMed) / 2
    }
    medA[i] = median

    # calculate mean divergence for each sequence
    for (j in divrA[i]) {
      score += divrA[i][j]
    }
    divR = score / qLen

    if (printDivAll) {
      for (seqHeaderSearch in divrA[i]) {
        for (seqHeaderQuery in seqA[i]) {
          sub(/>/, "", seqHeaderQuery)
          printf("%s\t%s.%s.%s\t", ++b, seqHeaderQuery, i, seqHeaderSearch)
        }
        printf("%.2f\n", divrA[i][seqHeaderSearch])
      }
    }

    if (printSeq) {
      for (seqHeaderQuery in seqA[i]) {
        sub(/>/, "", seqHeaderQuery)
        printf("%s\t%s\t%.2f\t%.2f\n", i, seqHeaderQuery, divR, median)
      }
    }

    # min and max divergence
    if ((minDivR == 0)||(divR < minDivR)) { minDivR = divR }
    if ((maxDivR == 0)||(divR > maxDivR)) { maxDivR = divR }
    divRSum += divR
  }

  # calculate median divergence for all sequences
  # using median values (!)
  mLen = asort(medA)
  if (mLen % 2) {
    median = medA[(mLen + 1) / 2]
  } else {
    lowMed = medA[mLen / 2]
    highMed = medA[(mLen / 2) + 1]
    median = (lowMed + highMed) / 2
  }

  if (debug) {
    printf("Total sequences:\t%s\n", seqNum)
    if (rules ~ /gg/) {
      print "Rule: gap-gap reduce length."
    } else if (rules) {
      print "Rule: gap-gap not reduce length."
    }
    if (rules ~ /gs/) {
      print "Rule: gap-site reduce length."
    } else if (rules) {
      print "Rule: gap-site not reduce length."
    }
    if (rules ~ /ggs/) {
      print "Count only differences in sites."
    }
    if (rules ~ /ns/) {
      print "Rule: N-site reduce length."
    }
    if (rules ~ /nn/) {
      print "Rule: N-N do not counted as similarity."
    }

    if (!rules) {
      print "Count all differences."
    }
  }

  if (printTotal) {
    printf("%s\t", FILENAME)
    printf("Min:\t%.4f", minDivR)
    printf("\tMax:\t%.4f", maxDivR)
    printf("\tMean:\t%.4f", divRSum / seqNum)
    printf("\tMed:\t%.4f\n", median)
  }
}
