# 
# In this script you should start imputating the 'context' of the ATP survay.
# I.E., the steps 3 and 4 of McPherson 2019.
#
# Maybe it worth to see the script 02_GSS_2004_network.R first
#

ATP_W3_sub <- readRDS("trabajo_1_files/ATP_W3_imput.rds")

sort(unique(ATP_W3_sub$age))
unique(ATP_W3_sub$sex)
sort(unique(ATP_W3_sub$educ_num))
levels(ATP_W3_sub$race)
levels(ATP_W3_sub$relig)

gss_egor <- readRDS("trabajo_1_files/gss_egor.rds")
gss_egos <- gss_egor$ego
gss_alters <- gss_egor$alter

sort(unique(gss_egos$age))
unique(gss_egos$sex)
sort(unique(gss_egos$educ_num))
levels(gss_egos$race)
levels(gss_egos$relig)

sort(unique(gss_alters$age))
unique(gss_alters$sex)
sort(unique(gss_alters$educ_num))
levels(gss_alters$race)
levels(gss_alters$relig)
