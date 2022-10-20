library(tidyverse)

### Functions:
str_to_vec <- function(string, pattern="\\|") {
  return(
    as.vector(str_split_fixed(string, pattern, n = str_count(string, pattern)+1))
  )
}


## GO-term wich include "interleukin-1":
intr1 <- c(
  str_to_vec(
    "TRAF6|MAP2K7|RELA|NFKB1|MYD88|MAPK3" # cellular response to interleukin-1
  ),
  str_to_vec(
    "TRAF6|MAP2K7|RELA|NFKB1|MYD88|MAPK3" # response to interleukin-1
  ),
  str_to_vec(
    "IL10|TNFA|CASP8|STAT3|MYD88" # regulation of interleukin-1 production
  ),
  str_to_vec(
    "TRAF6|RELA|MYD88|MAPK3" # interleukin-1-mediated signaling pathway 
  )
) %>% unique()

cat(intr1, sep = ", ")


## GO term; apoptotic process:
apt <- str_to_vec(
  "RB1|IL10|CDKN1A|ACVR1B|NFKB1|TNFA|CASP8|MYC|IRF1|BIRC5|RIPK1|MAP2K7|MYD88|MAPK3"
)

cat(apt, sep = ", ")

## GO term related to "leukocyte":
leu <- c(
  # str_to_vec(
  #   "RB1|IL10|TNFA|CASP8|MYC|TRAF6|IRF1|RIPK1" # regulation of leukocyte differentiation
  # ),
  str_to_vec(
    "IL10|TNFA|CASP8|TRAF6|IRF1|STAT3|MYD88|ICAM1" # leukocyte activation
  ),
  str_to_vec(
    "IL10|CSF1R|CDKN1A|TRAF6|IRF1|MYD88|MAPK3" # regulation of leukocyte proliferation
  ),
  str_to_vec(
    "IL10|CSF1R|TNFA|CASP8|TRAF6|IRF1|STAT3" # leukocyte differentiation
  )#,
  # str_to_vec(
  #   "RB1|TNFA|CASP8|MYC|TRAF6|RIPK1" # regulation of myeloid leukocyte differentiation
  # ),
  # str_to_vec(
  #   "RB1|IL10|TNFA|CASP8|TRAF6|RIPK1" # positive regulation of leukocyte differentiation
  # ),
  # str_to_vec(
  #   "CSF1R|TNFA|DUSP1|MYD88|ICAM1|MAPK3" # regulation of leukocyte migration
  # ),
  # str_to_vec(
  #   "RB1|TNFA|CASP8|TRAF6|RIPK1" # positive regulation of myeloid leukocyte differentiation
  # )
) %>% unique()

cat(leu, sep=", ")


## GO term JNK: 
JNK <- str_to_vec(
  "TNFA|TRAF6|RIPK1|MAP2K7|MYD88" # positive regulation of JNK cascade
)

cat(JNK, sep=", ")
