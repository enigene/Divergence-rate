# Calculate pairwise sequence divergence rate beween aligned FASTA sequnces
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# Usage: gawk -v rules=ggsnns -f divergr.awk alignment.fas > output.txt
# Rules: gg  = gap-gap reduce length,
#        gs  = gap-site reduce length
#        nn  = N-N do not counted as similarity
#        ns  = N-site reduce length

BEGIN {
  IGNORECASE = 1
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

  print "Num\tName\tDiv"

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
          } else {
            divrA[seqPosQuery][seqPosSearch] = 0
          }

          copyOfSeqHeaderQuery = seqHeaderQuery
          sub(/>/, "", copyOfSeqHeaderQuery)
          printf("%s\t%s.%s.%s\t", ++b, copyOfSeqHeaderQuery, seqPosQuery, seqPosSearch)

          divrAList[++cnt] = sprintf("%f", divrA[seqPosQuery][seqPosSearch])
          printf("%f\n", divrAList[cnt])

          if (printTotal) {
            # min and max divergence
            if ((minDivr == 0)||(divrAList[cnt] < minDivr)) { minDivr = divrAList[cnt] }
            if ((maxDivr == 0)||(divrAList[cnt] > maxDivr)) { maxDivr = divrAList[cnt] }
            divrSum += divrAList[cnt]
          }

        }
      }
    }
  }

  if (printTotal) {
    # calculate median divergence
    listLen = asort(divrAList)
    if (listLen % 2) {
      medDivr = divrAList[(listLen + 1) / 2]
    } else {
      lowmed = divrAList[listLen / 2]
      highmed = divrAList[(listLen / 2) + 1]
      medDivr = (lowmed + highmed) / 2
    }
    # calculate average divergence
    meanDivr = divrSum / listLen
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
    printf("Min:\t%f", minDivr)
    printf("\tMax:\t%f", maxDivr)
    printf("\tMean:\t%f", meanDivr)
    printf("\tMed:\t%f\n", medDivr)
  }
}
