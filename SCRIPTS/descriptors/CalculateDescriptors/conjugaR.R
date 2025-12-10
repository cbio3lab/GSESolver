#CODE FOR CONJUGATED SYSTEMS IN SEVERAL SDFs
#William Zamora
#11.01.22
#Libraries

library("ChemmineR")
#library("naturalsort")
library("dplyr")

args<-commandArgs(trailingOnly=TRUE)
temp1<-paste("*.sdf")
temp2<-Sys.glob(temp1)


for (k in 1:length(temp2)) {

#it must be placed here
input2=unlist(strsplit(temp2[k],"\\."))[1]
name=paste(input2,"conjugate_counts.csv",sep="_")

sdf=read.SDFset(temp2[k])
cjsystems=read.csv("cjsystems.csv")
#name of atoms
id0=row.names(atomblock(sdf[[1]]))
id1=unlist(strsplit(id0,"_"))
id2=seq(1,length(id1),2)
id=id1[id2]
#atoms1
a1=(bondblock(sdf[[1]])[,"C1"])
#atoms2
a2=(bondblock(sdf[[1]])[,"C2"])
#bond type
b=(bondblock(sdf[[1]])[,"C3"])
#bond table
data=data.frame(cbind(a1,a2,b))
#doble bonds
db=which(data[,3]==2)
tf=NULL
#new range of j from v1 11.01.21
for (j in 1:(length(db)-1)){   
    t=NULL
    cj=NULL
    for (i in 1:(length(db)-j)) {
        #1st system
        db10=data[db[j],1] 
        db1=data[db[j],2]

        #2nd system 
        db20=data[db[i+j],2]
        db2=data[db[i+j],1]

        #connector systems
        #new connector systems from v1 11.01.21
        sb1=data[which(data[,1]==db1 & data[,2]==db2),3]
        sbx=data[which(data[,1]==db2 & data[,2]==db1),3]

        sb2=data[which(data[,1]==db10 & data[,2]==db20),3]
        sby=data[which(data[,1]==db20 & data[,2]==db10),3]

        sb3=data[which(data[,1]==db10 & data[,2]==db2),3]
        sbz=data[which(data[,1]==db2 & data[,2]==db10),3]

        sb4=data[which(data[,1]==db1 & data[,2]==db20),3]
        sbw=data[which(data[,1]==db20 & data[,2]==db1),3]

        if (length(sb1) > 0) {
            if (sb1==1) {
            cj=paste(c(id[db10],id[db1],id[db2],id[db20]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sbx) > 0) {
            if (sbx==1) {
            cj=paste(c(id[db20],id[db2],id[db1],id[db10]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sb2) > 0) {
            if (sb2==1) {
            cj=paste(c(id[db1],id[db10],id[db20],id[db2]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sby) > 0) {
            if (sby==1) {
            cj=paste(c(id[db2],id[db20],id[db10],id[db1]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sb3) > 0) {
            if (sb3==1) {
            cj=paste(c(id[db1],id[db10],id[db2],id[db20]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sbz) > 0) {
            if (sbz==1) {
            cj=paste(c(id[db20],id[db2],id[db10],id[db1]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sb4) > 0) {
            if (sb4==1) {
            cj=paste(c(id[db10],id[db1],id[db20],id[db2]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }

        if (length(sbw) > 0) {
            if (sbw==1) {
            cj=paste(c(id[db2],id[db20],id[db1],id[db10]),collapse="")
            t00=which(is.element(cjsystems[,2],cj))
            t0=cjsystems[t00,1]
            t=append(t,t0)
            }
        }
    
    }
tf=append(tf,t)
}

#Table with frequencies
l1=unique(cjsystems$type)
ara=NULL
for (k in l1) {
ara0=which(cjsystems$type==unique(cjsystems$type)[k])
ara1=ara0[1]
ara=append(ara,ara1)
}
results=cjsystems[ara,] 
results$freq=rep(0,length(results[,1]))

for (i in 1:length(results[,1])){
results[i,3]=length(which(is.element(tf,results$type[i])))
}

#it must not be placed here
#input2=unlist(strsplit(temp2[k],"\\."))[1]
#name=paste(input2,"conjugate_counts.csv",sep="_")

write.table(results,name,row.names=FALSE,sep=",")

}



temp3<-paste("*conjugate_counts.csv")
temp4<-Sys.glob(temp3)
#temp4<-naturalsort(temp4)

data_f=NULL
dat=NULL
n1=NULL
id=NULL
for (l in 1:length(temp4)) {
data=read.csv(temp4[l])
dat=t(data[,3])
data_f=rbind(data_f,dat)

n1=paste(unlist(strsplit(temp4[l],"_"))[1],unlist(strsplit(temp4[l],"_"))[2],sep="_")
id=append(id,n1)

}

data2=cbind(id,data_f)


data2b<- apply(data2[,-1], 2,function(x) as.numeric(as.character(x)))
data2=cbind(data2[,1],data2b)

data=read.csv(temp4[1])
colnames(data2)[2:16]=t(data[,2])

write.table(data2,"conjugated_descriptors.csv",row.names=FALSE,sep=",")




