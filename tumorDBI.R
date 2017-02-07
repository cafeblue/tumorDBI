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




# heatmap 
#shareSNV <- mysqlQuery("select vc.tumor_name,vc.interID from variants_cancer as vc cross join (select interID,count(*) from variants_cancer group by interID HAVING count(*) > 1) tmp where vc.interID = tmp.interID")

shareSNV <- mysqlQuery("select si.sampleID,vc.interID from variants_cancer as vc inner join sampleInfo as si on si.postprocID = vc.tumor_name cross join (select interID,count(*) from variants_cancer group by interID HAVING count(*) > 1) tmp where vc.interID = tmp.interID and si.pairID > 0 and si.pairID < 100")

heat <- matrix(0, length(unique(shareSNV$sampleID)), length(unique(shareSNV$sampleID)))
colnames(heat) <- unique(shareSNV$sampleID)
rownames(heat) <- unique(shareSNV$sampleID)
for ( val1 in 1:(length(shareSNV$sampleID)-1) ) { for (val2 in (val1+1):length(shareSNV$sampleID)) { if (shareSNV$interID[val1] == shareSNV$interID[val2]) { heat[shareSNV$sampleID[val1],  shareSNV$sampleID[val2]] = heat[shareSNV$sampleID[val1],  shareSNV$sampleID[val2]] + 1 ; heat[shareSNV$sampleID[val2],  shareSNV$sampleID[val1]] = heat[shareSNV$sampleID[val2],  shareSNV$sampleID[val1]] + 1}}}
library(gplots)
heatmap(heat)
