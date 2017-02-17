install.packages("RMySQL")
library(RMySQL)

mysqlQuery <- function (query) {
  DB <- dbConnect(MySQL(), user='ashlien', password='Re4Cna$c', dbname='clinicalC', host='172.27.20.20', port=5029)
  rs <- dbSendQuery(DB, query)
  result <- fetch(rs, -1)
  dbDisconnect(DB)
  return(result)
}


#And then query through mydb:
#Grab rows for SNVs of postprocID 4912:
snvs <- mysqlQuery("select vc.*,si.sampleID,si.pairID as sampleID from variants_cancer inner join sampleInfo as si on si.postprocID = vc.tumor_name where vc.tumor_name = '4912'")
snvs
#Grab rows for INDELs for postprocID 4912:
indels <- mysqlQuery("select vc.*,si.sampleID,si.pairID from indel_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name where vc.tumor_name = '4912'")
indels
#Grab rows for INDELs for Molecular ID 295533:
indels <- mysqlQuery("select vc.*,si.sampleID,si.pairID from indel_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name where si.sampleID = '295533'")
indels
#Grab rows for SNVs for Molecular ID 295533:
snvs <- mysqlQuery("select vc.*,si.sampleID,si.pairID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name  where si.sampleID = '295533'") 
snvs
#Grab rows for SNVs for Modeluclar ID 295533 and all SNVs on the same sites in other samples:
snvs <- mysqlQuery("select vc.*,si.sampleID,si.pairID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name cross join (select vc.interID,si.sampleID,si.pairID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name where si.sampleID = '295533') as tmp where vc.interID = tmp.interID")
snvs   ####( samples can be easily separated by the “tumor_name”)



##########################
# heatmap 
##########################


library(ggplot2)
library(reshape2)

#shareSNV <- mysqlQuery("select CONCAT_WS('_', si.sampleID,si.pairID) as sampleID,vc.interID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name cross join (select interID,count(*) from variants_cancer where (((t_ref_count + t_alt_count) >= 50 AND (n_ref_count + n_alt_count) >= 50 AND (annovar_complete_genomics != '' OR annovar_complete_genomics = 0 OR annovar_clinvar like '%clin%') AND ( annovar_1000g != '' OR annovar_1000g = 0 OR annovar_clinvar like '%clin%') AND judgement = 1) OR (( annovar_1000g is NULL OR annovar_1000g < 0.01 ) AND ( annovar_esp is NULL OR annovar_esp < 0.01) AND t_alt_count >= 10 AND  n_alt_count < 3 AND n_ref_count > 50 AND tumor_f > 0.01 AND cosmic_census = 1 AND (annovar_exonic_func REGEXP 'nonsynonymous SNV|stopgain SNV|stoploss SNV|splicing' OR  annovar_func REGEXP 'splicing'))) group by interID HAVING count(*) > 1) tmp where vc.interID = tmp.interID and si.pairID > 0 and si.pairID < 100") ## !!!!!!!!!!!!!!!   a known bug: some poor quality SNVs will be selected for plot.

shareSNV <- mysqlQuery("select CONCAT_WS('_', si.sampleID,si.pairID) as sampleID,vc.interID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name cross join (select interID,count(*) from variants_cancer where (((t_ref_count + t_alt_count) >= 50 AND (n_ref_count + n_alt_count) >= 50 AND (annovar_complete_genomics != '' OR annovar_complete_genomics = 0 OR annovar_clinvar like '%clin%') AND ( annovar_1000g != '' OR annovar_1000g = 0 OR annovar_clinvar like '%clin%') AND judgement = 1) OR (( annovar_1000g is NULL OR annovar_1000g < 0.01 ) AND ( annovar_esp is NULL OR annovar_esp < 0.01) AND t_alt_count >= 10 AND  n_alt_count < 3 AND n_ref_count > 50 AND tumor_f > 0.01 AND cosmic_census = 1 AND (annovar_exonic_func REGEXP 'nonsynonymous SNV|stopgain SNV|stoploss SNV|splicing' OR  annovar_func REGEXP 'splicing'))) group by interID HAVING count(*) > 1) tmp where vc.interID = tmp.interID and si.pairID > 0 and si.pairID < 100 and (((t_ref_count + t_alt_count) >= 50 AND (n_ref_count + n_alt_count) >= 50 AND (annovar_complete_genomics != '' OR annovar_complete_genomics = 0 OR annovar_clinvar like '%clin%') AND ( annovar_1000g != '' OR annovar_1000g = 0 OR annovar_clinvar like '%clin%') AND judgement = 1) OR (( annovar_1000g is NULL OR annovar_1000g < 0.01 ) AND ( annovar_esp is NULL OR annovar_esp < 0.01) AND t_alt_count >= 10 AND  n_alt_count < 3 AND n_ref_count > 50 AND tumor_f > 0.01 AND cosmic_census = 1 AND (annovar_exonic_func REGEXP 'nonsynonymous SNV|stopgain SNV|stoploss SNV|splicing' OR  annovar_func REGEXP 'splicing')))") 

shareSNV$sampleID <- sub('_68', '_68***', shareSNV$sampleID)

heat <- matrix(0, length(unique(shareSNV$sampleID)), length(unique(shareSNV$sampleID)))
colnames(heat) <- unique(shareSNV$sampleID)
rownames(heat) <- unique(shareSNV$sampleID)
for ( val1 in 1:(length(shareSNV$sampleID)-1) ) { for (val2 in (val1+1):length(shareSNV$sampleID)) { if (shareSNV$interID[val1] == shareSNV$interID[val2]) { heat[shareSNV$sampleID[val1],  shareSNV$sampleID[val2]] = heat[shareSNV$sampleID[val1],  shareSNV$sampleID[val2]] + 1 ; heat[shareSNV$sampleID[val2],  shareSNV$sampleID[val1]] = heat[shareSNV$sampleID[val2],  shareSNV$sampleID[val1]] + 1}}}
heat.m <- melt(heat)

ggplot(heat.m, aes(Var1, Var2, fill = value)) + geom_tile(aes(fill = value)) +  geom_text(data=subset(heat.m, value > 0), aes(label = value), size = 1.5) + scale_fill_gradient(low = "white", high = "red") + ylab('SampleID_KiCS') + xlab("SampleID_KiCS") + theme(axis.text=element_text(size=6), axis.text.x = element_text(angle = 90, hjust = 1))


##### search for the variants in common:

pairwiseCommonSNV <- shareSNV[(shareSNV['sampleID'] == '288555_45' | shareSNV['sampleID'] == '289552_55' ),]
pairwiseCommonSNV <- pairwiseCommonSNV[duplicated(pairwiseCommonSNV$interID), ]
pairwiseCommonSNV <- paste(pairwiseCommonSNV$interID, collapse = "', '")
pairwiseCommonSNV <- paste("'", pairwiseCommonSNV, "'", collapse = "", sep = '')
commonSNVs <- mysqlQuery(paste("SELECT vc.*,si.sampleID,si.pairID FROM variants_cancer AS vc INNER JOIN sampleInfo AS si ON si.postprocID = vc.tumor_name  WHERE si.sampleID IN ('288555', '289552') AND vc.interID IN (", pairwiseCommonSNV, ")", collaps = ""))

#####  variants under current investigated 
currentSNVs <- mysqlQuery("SELECT vc.*,si.sampleID,si.pairID FROM variants_cancer AS vc INNER JOIN sampleInfo AS si ON si.postprocID = vc.tumor_name CROSS JOIN ( SELECT vc.hgnc_gene FROM variants_cancer as vc INNER JOIN sampleInfo AS si ON si.postprocID = vc.tumor_name WHERE si.sampleID = '294906') AS tmp WHERE (vc.hgnc_gene = tmp.hgnc_gene AND vc.annovar_exonic_func LIKE '%nonsynonymous%') OR si.sampleID = '294906'")


###########################
# Coverage plot
###########################

