# Hierachical-Simulation-Data
The data will be simulated from a hirachical model. Start with the Bernoulli Distribution, when delta=0, use the ZINB Distribution; when delta=1, use the Poisson Distribution. The mean of the Poisson Distribution will subject to a Uniform Distribution

#### pre-set parameters determined by the real data analysis, the parameters from real data are saved in the variable RealData
#### ZINB Data driven parameter initialization
	p_m = RealData$p0
	p   = RealData$p
	scale_m0 = (1-p)/p
	shape_m0 = RealData$r
#### Data driven NULL proportion for different distance	
	pi_m = (RealData$pi0)*0.95				
#### number of bins
	nbin = 2000
	simulation_times = 20 #20
#### Threshold on distance for FitHiC
	distUpThres  = D
	distLowThres = 0

#### simulate data
	k = 1
	C = 0.5
	outdir = paste0(rootdir,'C_',C,'_k_',k,'/')
	dir.create(outdir)
## varying simulation parameters
scale_m = scale_m0
shape_m = C*shape_m0
#### simulation will return fragmentLists and contactCounts ####		
for(i in 1:simulation_times)
{
	###i=1
	outid = paste0("C_", C,"_k_", k,"_S_",i)
	print(outid)
	tempii = simulation(nbin=nbin, pi_m=pi_m, k=k, p_m=p_m, scale_m=scale_m, shape_m=shape_m, D=D,outdir=outdir,ncores=ncores)
}		


### The definitaion of the function "simulation" is given below:
require(parallel)
simulation = function(nbin, pi_m, k, p_m, scale_m, shape_m, D,outdir,ncores=16)
{    
    options(mc.cores=ncores,mc.preschedule=F)
    outid = paste0("C_", C,"_k_", k,"_S_",i)
    mu_m    = shape_m*scale_m		

    #### combination of pairs ####
    c1 	  = rep(1:(nbin),D)
    c2    = c1+rep(1:D,each=nbin)
    Pairs = rbind(c1,c2)### Pairs[,c(1:5,1300:1305)] dim(Pairs)
    #Pairs = combn(nbin,2)  
	
    #### all possible distance ####
    Distance = Pairs[2,] - Pairs[1,]
    m = sort(unique(Distance))
    
    #### delta for different distance ####
    delta_i_m  = mclapply( 1:max(m), function(i) rbinom(length(which(Distance==m[i])),1,1-pi_m[i]))
    ## positions of delta = 0
    delta_pos0 = mclapply( 1:max(m), function(i) which(delta_i_m[[i]]==0))
    ## positions of delta = 1
    delta_pos1 = mclapply( 1:max(m), function(i) which(delta_i_m[[i]]==1))
    
    
    #### ZINB ####
    ## when delta=0, use ZINB distribution
    ## use Bernoulli for each trail
    s_i_m  = mclapply( 1:max(m), function(i) rbinom(length(delta_pos0[[i]]),1,p_m[i]))
    ## positions of s=1
    s_pos1 = mclapply( 1:max(m), function(i) which(s_i_m[[i]]==1))
    ## positions of s=0
    s_pos0 = mclapply( 1:max(m), function(i) which(s_i_m[[i]]==0))
    ## assign 0 to s=1
    s_i_m1 = mclapply( 1:max(m), function(i) rep(0,length(s_pos1[[i]])))
    ## generate Negative Binomial random variable for s=0
    s_i_m0 = mclapply( 1:max(m), function(i) rnbinom(length(s_pos0[[i]]),size = shape_m[i], prob = 1/(1+scale_m[i])))
        
    ## When delta=1, use Poisson distribution
    ## generate lambda
   	lambda_m = mclapply( 1:max(m), function(i) mu_m[i]*(1.1+k*runif(length(delta_pos1[[i]]))))
    #	lambda_m = lapply( 1:max(m), function(i) mu_m[i]+rgamma(length(delta_pos1[[i]]), shape=1.5, rate=1/(k*mu_m[i])))   
    ## generate Counts by Poisson distribution
    xx_i_m   = mclapply( 1:max(m), function(i) sapply(1:length(delta_pos1[[i]]), function(jj) rpois(1,lambda_m[[i]][jj])))
    
    #### Complete the X_i_m
    x_i_m = mclapply(1:max(m), function(i) rep(0,length(which(Distance==m[i])))) #rbinom(length(which(Distance==m[i])),1,0))
    for(i in 1:max(m))
    {
        ## complete s
        s_i_m[[i]][s_pos1[[i]]] = s_i_m1[[i]]
        s_i_m[[i]][s_pos0[[i]]] = s_i_m0[[i]]        
        ## complete x
        x_i_m[[i]][delta_pos0[[i]]] = s_i_m[[i]]
		if(length(delta_pos1[[i]])>0)
			x_i_m[[i]][delta_pos1[[i]]] = xx_i_m[[i]]     
    }
    
    #### generate contactCounts dataframe ####
    x_i_m_vector  = as.vector(unlist(x_i_m))
    delta_i_m_vector  = as.vector(unlist(delta_i_m))
    contactCount  = rbind(Pairs,x_i_m_vector,delta_i_m_vector)
    contactData  = t(contactCount)
    colnames(contactData) = c('bin1','bin2','count','peak')
    contactCounts = data.frame(chr1=rep(11111,length(x_i_m_vector)), fragmentMid1=contactData[,1],   chr2=rep(11111,length(x_i_m_vector)), fragmentMid2=contactData[,2], contactCount=contactData[,3])        
    
	#### generate fragmentLists dataframe ####
    fragCount = NULL
    fragCount = sapply(split(contactCounts$contactCount,f=contactCounts$fragmentMid1),sum)#rowSums(x_i_m_matrix)
    fragmentLists = data.frame(V1=rep(11111,length(nbin)), V2=rep(11111,length(nbin)), V3=(1:nbin), V4=fragCount, V5=rep(11111,length(nbin)))
    save(contactData, file=paste0(outdir,"contactData_",outid,".Rdata"))
    write.table(contactData[,4], sep='\t',file=paste0(outdir,"delta_",outid,".txt"),row.names=F,col.names=F,quote=F)
    
	# generate contactCounts.gz
    filename = paste0("contactCounts_",outid,'.gz')
    file.create(file.path(outdir, filename))			
    write.table(contactCounts, paste0(outdir,filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    
	# generate fragmentLists.gz
    filename = paste0("fragmentLists_",outid,'.gz')
    file.create(file.path(outdir, filename))
    write.table(fragmentLists, paste0(outdir,filename), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    return(NULL)			
}    
    
