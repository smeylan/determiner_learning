buildDirectoryStructure = function(analysisOutputDir){
	dir.create(analysisOutputDir, showWarnings=F)
	### Build out the directory structure that will be packaged and retrieved from the cluster ###
	analysisDir = paste0(analysisOutputDir,(gsub('[ :]','\\_',Sys.time())), '/')
	dir.create(analysisDir, showWarnings=F)
	#make a list with the input file and output file names
	graphPath = paste0(analysisDir, 'graphs/')
	dir.create(graphPath, showWarnings=F)
	tablePath = paste0(analysisDir, 'tables/')
	dir.create(tablePath, showWarnings=F)

	rlist = list()	
	rlist['analysisDir'] =  analysisDir
	rlist['graphPath'] = graphPath
	rlist['tablePath'] = tablePath
	return(rlist)
}

getMeanAndHPDIntervalsForVar = function(modelResult, variable, mcmc){  
	### get the mean, median, and HPDs for a variable from a fitted JAGS model ###	
  if (!is.null(modelResult)){  	
  	if(mcmc){
  		combinedChains = mcmc(do.call('rbind',modelResult)) 	
  	} else{ 
  		#take a vector and pretend that it is an mcmc chain
  		combinedChains = mcmc(data.frame(modelResult[,variable],1))
  		colnames(combinedChains) = c(variable,'one')
  	}
    rdf = as.data.frame(HPDinterval(combinedChains[,variable], prob=.95))$upper
    names(rdf)[1] = paste(variable,'_p975', sep='')
    rdf[[paste(variable,'_p025', sep='')]] = as.data.frame(HPDinterval(combinedChains[,variable], prob=.95))$lower  
    rdf[[paste(variable,'_p995', sep='')]] = as.data.frame(HPDinterval(combinedChains[,variable], prob=.99))$upper  
     rdf[[paste(variable,'_p005', sep='')]] = as.data.frame(HPDinterval(combinedChains[,variable], prob=.99))$lower  
     rdf[[paste(variable,'_p9995', sep='')]] = as.data.frame(HPDinterval(combinedChains[,variable], prob=.999))$upper  
     rdf[[paste(variable,'_p0005', sep='')]] = as.data.frame(HPDinterval(combinedChains[,variable], prob=.999))$lower    	
		rdf[[paste(variable,'_min', sep='')]] = min(combinedChains[,variable])
		rdf[[paste(variable,'_max', sep='')]] = max(combinedChains[,variable])	        
		rdf[[paste(variable,'_mean', sep='')]] = mean(combinedChains[,variable])
		rdf[[paste(variable,'_median', sep='')]] = median(combinedChains[,variable])
		rdf[[paste(variable,'_var', sep='')]] = var(combinedChains[,variable])
		rdf[[paste(variable,'_sd', sep='')]] = sd(combinedChains[,variable])          
	} else {  	
  		print('Model result is NULL')
		rdf = data.frame(max=NA)    
		names(rdf)[1] = paste(variable,'_max', sep='')
		rdf[[paste(variable,'_min', sep='')]] = NA
		rdf[[paste(variable,'_mean', sep='')]] = NA
  }
  return(rdf)
}


processModelReturnData = function(filename, overwrite, tablePath){
	### Process the results from a single RData file containing a fitted model ###	
	filenameSplit = strsplit(filename, '/')[[1]]
	fnNoDirectory = filenameSplit[length(filenameSplit)]
	fileParts = strsplit(fnNoDirectory, '_')[[1]]	
	summaryFileName = paste(tablePath, gsub('.RData',  '_summary.csv', fnNoDirectory), sep='')
	
	# don't overwrite the summary if we already have one and have specified not to overwrite it
	stopifnot(overwrite & !(summaryFileName %in% list.files()))
	
	# pull out properties of the dataset into separate variables	
	g(modelType,corpusName,numSamples, windowSize, windowSkip, imputation, morphology,extractionMethod) %=% c(fileParts[1], fileParts[3:5], fileParts[7:9], gsub('.RData','',fileParts[10]))

	#load the dataset
	load(filename)
	if (exists('linking.model.results')){
						
		#decompose a linking change model into two different models, and treat them as successive sliding windows
		if (modelType == 'linkingChange'){
			model.results = unlist(lapply(linking.model.results, unroll_LC_results), recursive=F)	
			
			Nr.observed = unlist(lapply(Nr.combined.observed.list, unroll_LC_input), recursive=F)					
		} else if (modelType == 'linking'){
			model.results = linking.model.results
			Nr.observed = Nr.combined.observed.list
		}
		
		varsOfInterest = c('mu0.child', 'nu.child', 'eta', 'mu0.adult', 'nu.adult')
									
	} else {
		stop('No fit model found.')
	}
			
	if(length(Nr.observed) != length(model.results)){
    	stop('There should be as many input windows as models')
  	}
	
	# break into successive pieces and treat the windows in succession
	windowResults = do.call('rbind.fill',lapply(1:length(model.results), function(i){
		processWindow(model.results[[i]], Nr.observed[[i]], varsOfInterest)			
	}))		
	windowResults$child = fixNames(windowResults$child)
	
	write.csv(windowResults, summaryFileName, row.names=F)		
}

processWindow = function(modelWindow, inputWindow, varsOfInterest){
	### get variable estimates, simualted overlap, empirical overlap, and Yang expected overlap for a window ###
	childNounCounts = subset(inputWindow, n.child > 0)$n.child 
	childDefUses = subset(inputWindow, n.child > 0)$r.child
	if (is.null(modelWindow)){
		hpd_vals = data.frame(child = inputWindow$child[1])
		hpd_vals$index = inputWindow$interval[1]
		hpd_vals$minAge = inputWindow$minAge[1]
		hpd_vals$maxAge = inputWindow$maxAge[1]
		hpd_vals$meanAge = inputWindow$meanAge [1]
		hpd_vals$meanAgeInMonths = inputWindow$meanAge[1]/30.5
		hpd_vals$minAgeInMonths = inputWindow$minAge[1]/30.5
		hpd_vals$maxAgeInMonths = inputWindow$maxAge[1]/30.5 
		rdf = hpd_vals
	} else {
		# get the estimates for the variables 	
		hpd_vals = data.frame(t(unlist(lapply(varsOfInterest, function(voa){      
			getMeanAndHPDIntervalsForVar(modelWindow,voa,mcmc=T)
		})))) 
		
		#propagate properties of the input data to the returned summary data 
		# first confirm that they are unique to verify dataset integrity
		sapply(c('child', 'interval', 'minAge', 'meanAge'), function(var){
			if (length(unique(inputWindow[[var]])) != 1){stop(p('Variable', var, 'is not unique in data window'))} else {p('Variable', var, 'is unique')}	
		})
			
		#get the number of types and tokens	
		hpd_vals$childTypes = nrow(subset(inputWindow, n.child > 0))	
		hpd_vals$childTokens = sum(inputWindow$n.child)
		hpd_vals$parentTypes = nrow(subset(inputWindow, n.parent > 0))
		hpd_vals$parentTokens = sum(inputWindow$n.parent)
		
		hpd_vals$child = inputWindow$child[1]
		hpd_vals$index = inputWindow$interval[1]
		hpd_vals$minAge = inputWindow$minAge[1]
		hpd_vals$maxAge = inputWindow$maxAge[1]
		hpd_vals$meanAge = inputWindow$meanAge [1]
		hpd_vals$meanAgeInMonths = inputWindow$meanAge[1]/30.5
		hpd_vals$minAgeInMonths = inputWindow$minAge[1]/30.5
		hpd_vals$maxAgeInMonths = inputWindow$maxAge[1]/30.5  	 
		if('splitHalf' %in% names(inputWindow)){  # keep the split half information
			hpd_vals$splitHalf = inputWindow$splitHalf[1]	
		}
		
		# get the overlap statistics sampled from the model
		nounIndices = grep('r.child.simulated',colnames(modelWindow[[1]])) #can take the column names from the first chain; other chains are the same		
		
		simulated_overlaps_samples = unlist(lapply(1:length(modelWindow), function(chainIndex){
	    	sapply(1: nrow(modelWindow[[chainIndex]]), function(sampleIndex){
				getOverlap(data.frame(
	                r.child = modelWindow[[chainIndex]][sampleIndex, nounIndices], n.child = childNounCounts), speaker='child')})}))
		
		OPs = c(0,.0005,.005,.025, .05, .5, .95, .975, .995, .9995, 1)
		OPstrings = gsub('[0|\\.]+$','',gsub('^0\\.','',format(OPs, scientific=F)))
		OPstrings[1] = "0"
		overlapPercentiles = as.data.frame(t(quantile(simulated_overlaps_samples, probs = OPs)))
		names(overlapPercentiles) = paste0('overlap_p', OPstrings)
		rdf = cbind(hpd_vals, overlapPercentiles)
		rdf$overlap_mean = mean(simulated_overlaps_samples)
		rdf$overlap_median = median(simulated_overlaps_samples)
	}
	
	# get thge empirical overlap statistics
	rdf$childEmpiricalOverlap =  getOverlap(data.frame(r.child = childDefUses, n.child = childNounCounts), speaker='child')	
	rdf$parentEmpiricalOverlap =  getOverlap(subset(inputWindow, n.parent>0), speaker='parent')	
	
	# get the Yang overlap model predictions
	rdf$childYangPredictedOverlap = yangExpectedOverlap(childNounCounts)
	rdf$parentYangPredictedOverlap =  yangExpectedOverlap(subset(inputWindow, n.parent>0)$n.parent)	
	return(rdf)	
}

expected.overlap <- function(N,d,S,a=1) {
	### Model from Yang (2013), PNAS ###
	D <- length(d)
	tmp <- sapply(1:N, function(x) 1/x^a)
	p <- tmp/sum(tmp)
	## O(r,D,N,S) helper function
	helper <- function(r)
	  1 + (D-1)*(1-p[r])^S - sum(sapply(d,function(d.i) (d.i * p[r] + 1 - p[r])^S))
	o.r.N.D.S <- sapply(1:N, helper)
	return(o.r.N.D.S)
}

yangExpectedOverlap <- function(nounVector) {  
	### vector of token counts for nouns. The length of this vector is the type count ###
	N = sum(as.numeric(nounVector > 0))
	S = sum(nounVector)  
	return(mean(expected.overlap(N,d=c(2/3,1/3),S,a=1)))
}

getOverlap = function(df, speaker){
	### Compute the number of nouns seen with both determiners over the total number of nouns, separated by child and parent ###
	if (speaker == 'child'){
		length(which(df$r.child < df$n.child & df$r.child > 0)) / length(which(df$n.child > 0))	
	} else if (speaker == 'parent') {
		length(which(df$r.parent < df$n.parent & df$n.parent > 0)) / length(which(df$n.parent > 0))	
	}	
}

unroll_LC_results = function(lcResult){
  ### decomposes linking change window results into two successive windows ###
  	if (is.null(lcResult)){ #some nonconvergent windows— with null model results
  		phase1 = NULL
  		phase2 = NULL
  	} else {
	  	colNamesNoBrackets = gsub('\\[.*\\]','', colnames(lcResult[[1]]))
		phase1columns  = grep ('1$', colNamesNoBrackets)
		phase2columns  = grep ('2$', colNamesNoBrackets)	
		phase1 = lcResult[,phase1columns]			
		phase2 = lcResult[,phase2columns]
		
		#remove the reference to column names
		for (i in 1:length(phase1)){
			colnames(phase1[[i]]) = gsub('1$','', gsub('1\\[','\\[',colnames(phase1[[i]])))			
		}
		
		for (i in 1:length(phase2)){
			colnames(phase2[[i]]) = gsub('2$','', gsub('2\\[','\\[',colnames(phase2[[i]])))			
		}
	}
	return(list(phase1, phase2))
 }
 
 
unroll_LC_input = function(lcNr){
 	 ### decomposes linking change window input data into two successive windows ###
	phase1 =subset(lcNr, splitHalf == 1)
	phase2 =subset(lcNr, splitHalf == 2)
	return(list(phase1, phase2))	
}

getFilenames = function(summaryDir, modelType, dataPrep, fileType){
	### concatenate the summary files with a given dataPrep in the filename 
	fileList =  list.files(summaryDir)
	toCat = fileList[grep(p(modelType,'_.*',dataPrep,fileType,sep=''),fileList)]
	return(toCat)
}

getModelResults = function(summaryDir, modelType, dataPrep){
	### concatenate the summary files with a given dataPrep ito a dataframe ###
	toCat = getFilenames(summaryDir, modelType, dataPrep, '_summary.csv')
	rdf = do.call('rbind.fill', lapply(toCat, function(x){
		cbind(read.csv(p(ensureTrailingSlash(summaryDir), x, sep=''), stringsAsFactors=F), file=x)
	}))
	
	return(rdf)	
}

fixNames = function(nameVec){
	### Restore full human-readable names of corpora ###
	realNames = c('Peter','Eric','Adam','Eve','Sarah','Abe','Anne','Aran','Becky','Carl','Domin', 'Gail','Joel','John','Liz','Nic','Ruth','Warr','Alex','Ethan','Lily','Naima','Violet','William', 'Naomi', 'Speechome','Nina','Thomas')
	
	changeDict = list()
	for (key in realNames){
		changeDict[tolower(substr(key,1,3))] = key
		changeDict[substr(key,1,3)] = key
		changeDict[tolower(key)] = key
	}	
	changeDict['N']	= 'Naomi'
	changeDict['n']	= 'Naomi'

	unname(unlist(sapply(nameVec, function(name){
		if (name %in% realNames){
			return(name)		
		} else if (name %in% names(changeDict)) {
			return(changeDict[name])			
		} else {
			stop(p(name,'mapping not defined'))			
		}				
	})))
}


restoreCorpusInfo = function(nameVec){
	corpusList = as.list(c(rep('Bloom',2), rep('Brown',3), 'Kuczaj', rep('Manchester',12), rep('Providence', 6), 'Sachs','Speechome','Suppes','Thomas'))
	names(corpusList) = c('Peter','Eric','Adam','Eve','Sarah','Abe','Anne','Aran','Becky','Carl','Domin', 'Gail','Joel','John','Liz','Nic','Ruth','Warr','Alex','Ethan','Lily','Naima','Violet','William', 'Naomi', 'Speechome','Nina','Thomas')
	unname(unlist(sapply(nameVec, function(name){
		if (name %in% names(corpusList)){
			return(corpusList[name])
		} else {
			stop(p(name,'mapping not defined'))			
		}
	})))
}			

getAllSplitHalfVarEstimates = function(filename, varname){
	### get variable estimates from split half models ###
    load(filename)    
	  
	nuValues = lapply(c(1:length(linking.model.results)), function(index){
        #special treatment for nu.difference
        if(varname == 'nu.difference'){
        	diffing=T
        	varname = 'nu.child' 	
        } else {
        	diffing = F
        }        
        
		if (!is.null(linking.model.results[[index]])){	    	
	 	   varPhase1Index = grep(p(varname,'1',sep=''), colnames(linking.model.results[[1]][[1]]))
		   varPhase2Index = grep(p(varname,'2', sep=''), colnames(linking.model.results[[1]][[1]]))
		
			phase1var = as.vector(
	        sapply(1:length(linking.model.results[[index]]), function(ind){ unlist(linking.model.results[[index]][[ind]][,varPhase1Index])})) 
			phase2var = as.vector(
	        sapply(1:length(linking.model.results[[index]]), function(ind){ unlist(linking.model.results[[index]][[ind]][,varPhase2Index])}))
	      	      
		} else {
			phase1var = NA
			phase2var = NA
		}
		
		if (diffing){
			rdf = data.frame(nu.difference = phase2var - phase1var)
			rdf$phase1_mean_age = unique(subset(Nr.combined.observed.list[[index]], splitHalf == 1)$meanAge)
			rdf$phase2_mean_age = unique(subset(Nr.combined.observed.list[[index]], splitHalf == 2)$meanAge)
			
		} else {
			phase1 = data.frame(phase1var, phase=1)
			phase1 = renameColumn(phase1, 'phase1var', varname)
			phase2 = data.frame(phase2var, phase=2)
			phase2 = renameColumn(phase2, 'phase2var', varname)
			
			phase1$meanAge =  unique(subset(Nr.combined.observed.list[[index]], splitHalf == 1)$meanAge)
			phase2$meanAge =  unique(subset(Nr.combined.observed.list[[index]], splitHalf == 2)$meanAge)
			rdf = rbind(phase1, phase2)
		}	      
		rdf$child = unique(Nr.combined.observed.list[[index]]$child)	
		return(rdf)					
	})  
 	rdf = do.call('rbind.fill', nuValues)  
 	rdf$child = fixNames(rdf$child)
 	rdf$corpus = restoreCorpusInfo(rdf$child)
  	return(rdf)  				
}

getNuDifferenceHPDs = function(nu.differences.short){
	### compute HPDs over nu.difference ###
	nd.byChild = split(nu.differences.short, nu.differences.short$child)		
	do.call('rbind.fill', lapply(nd.byChild, function(x){
		if (nrow(x) ==  1){
			rdf = data.frame(nu.difference_mean= NA)			
		} else {
			rdf = as.data.frame(t(getMeanAndHPDIntervalsForVar(x,'nu.difference',mcmc=F)))	
		}
		
		rdf$child = unique(x$child)
		rdf$corpus = unique(x$corpus)
		rdf$phase1_mean_age = unique(x$phase1_mean_age) #two different mean ages, very close
		rdf$phase2_mean_age = unique(x$phase2_mean_age)
		return(rdf)				
	}))		
} 

rmse = function(x,y){
	### compute RMSE for all instances where x and y are defined ###
	sqrt(mean((x[!is.na(x) & !is.na(y)]-y[!is.na(x) & !is.na(y)])^2))
}

getCorrelations = function(paths){
	### Compute correlations for each data treatment * (current model vs. Yang,2013) ###
	getCor = function(modelType){
		overlaps = do.call('rbind.fill',lapply(dataPreps, function(dataPrep){		
			df = subset(getModelResults(paths[['tablePath']], modelType, dataPrep), !(child %in% c('Eric')))

			childYangPredictedOverlap = cor(df$childYangPredictedOverlap,df$childEmpiricalOverlap, use='pairwise.complete.obs')
	 		childOverlap_mean = cor(df$overlap_mean, df$childEmpiricalOverlap, use='pairwise.complete.obs')
	  
		  	childYangPredictedOverlap_RMSE = rmse(df$childYangPredictedOverlap,df$childEmpiricalOverlap)
		  	childOverlap_RMSE = rmse(df$overlap_mean, df$childEmpiricalOverlap)
	  
	  		
	  		if (modelType == 'linking'){
		  		medAge = median(df$meanAgeInMonths)
				df$splitHalf = as.numeric(df$meanAgeInMonths >= medAge)+1
		  	}
			df1= subset(df, splitHalf == 1)
			df2= subset(df, splitHalf == 2)
		  	
	  		childYangPredictedOverlap1 = cor(df1$childYangPredictedOverlap,df1$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	currentModel1 = cor(df1$overlap_mean, df1$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	childYangPredictedOverlap1 = cor(df1$childYangPredictedOverlap,df1$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	currentModel1 = cor(df1$overlap_mean, df1$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	childYangPredictedOverlap_RMSE1 = rmse(df1$childYangPredictedOverlap,df1$childEmpiricalOverlap)
	 	 	currentModel_RMSE1 = rmse(df1$overlap_mean, df1$childEmpiricalOverlap)  
		  	
	  		childYangPredictedOverlap2 = cor(df2$childYangPredictedOverlap,df2$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	currentModel2 = cor(df2$overlap_mean, df2$childEmpiricalOverlap, use='pairwise.complete.obs')
	 	 	childYangPredictedOverlap_RMSE2 = rmse(df2$childYangPredictedOverlap,df2$childEmpiricalOverlap)
	 	 	currentModel_RMSE2 = rmse(df2$overlap_mean, df2$childEmpiricalOverlap)
				  
		  	return(data.frame(dataPrep, childYangPredictedOverlap, childYangPredictedOverlap_RMSE, childOverlap_RMSE, childYangPredictedOverlap1, currentModel1, childYangPredictedOverlap_RMSE1, currentModel_RMSE1, childYangPredictedOverlap2, currentModel2, childYangPredictedOverlap_RMSE2, currentModel_RMSE2))						
		}))
		return(overlaps)	
	}
	
	dataPreps = c('all_LN','all_FN','all_standard','none_standard','singulars_LN','singulars_FN', 'singulars_standard')
	lcOverlaps = getCor('linkingChange')
	write.csv(lcOverlaps, paste0(paths['tablePath'],'/splitHalfOverlap.csv'))
	swOverlaps = getCor('linking')
	write.csv(swOverlaps, paste0(paths['tablePath'],'/slidingWindowOverlap.csv'))
	#write these out to CSVs
	rlist = list()
	rlist[['slidingWindow']] = swOverlaps
	rlist[['lcOverlaps']] = lcOverlaps
	return(rlist)		
}

checkForModelCompletion = function(modelDir, outfile){
  fileList = gsub('.csv', '', list.files(modelDir))
  
  #build the list of data file names that should be present
  corpus = c('Bloom70', 'Brown', 'Kuczaj', 'Manchester', 'Providence', 'Thomas', 'Suppes', 'Speechome', 'Sachs')
  morphology = c('all', 'singulars', 'none')
  datasetName = c('LN', 'FN', 'standard')
  
  fullGrid = expand.grid(corpus=corpus, morphology=morphology, datasetName=datasetName, stringsAsFactors=F)
  
  partialGrid = subset(fullGrid, !(morphology == 'none' & datasetName %in% c('LN','FN')) & !(corpus %in%  c('Speechome', 'Thomas') & datasetName == 'standard') 
  )
     
  
  linkingFilePatterns = paste('linking_Model_',partialGrid$corpus,'.*',partialGrid$morphology,'_',partialGrid$datasetName,'.RData', sep='')
  partialGrid$slidingWindowPresence = unname(sapply(linkingFilePatterns, function(x){length(grep(x, fileList))> 0}))
  
  changeFilePatterns = paste('linkingChange_Model_',partialGrid$corpus,'.*',partialGrid$morphology,'_',partialGrid$datasetName,'.RData', sep='')
  partialGrid$changePresence = unname(sapply(changeFilePatterns, function(x){length(grep(x, fileList))> 0}))
 
 if (!is.null(outfile)){
 	write.csv(partialGrid, outfile, row.names=F)
 }
 return(partialGrid)
}  

mclapply2 = function(..., parallelize=T, cores = detectCores()){
    ### alias mclapply to a function that reverts to lapply if passed a flag ###
    if (!parallelize){
        lapply(...)
    } else {
        mclapply(..., mc.cores= cores)
    }
}

makeSlidingWindowPlots = function(paths, dataPrep, childList){	
	swdf = getModelResults(paths$tablePath,'linking', dataPrep)		

	swdf = subset(swdf, child %in% childList)
	earliest = aggregate(meanAgeInMonths ~ child, swdf, min)
	swdf$child = factor(swdf$child, levels = as.character(earliest$child)[order(earliest$meanAgeInMonths)]) 		
	swdf$nu.child_sd_low = swdf$nu.child_mean - swdf$nu.child_sd
	swdf$nu.child_sd_high = swdf$nu.child_mean + swdf$nu.child_sd
	swdf$mu0.child_sd_low = swdf$mu0.child_mean - swdf$mu0.child_sd
	swdf$mu0.child_sd_high = swdf$mu0.child_mean + swdf$mu0.child_sd
	swdf$eta_sd_low = swdf$eta_mean - swdf$eta_sd
	swdf$eta_sd_high = swdf$eta_mean + swdf$eta_sd

	swdf$range = swdf$nu.child_sd_high - swdf$nu.child_sd_low
	rangeLimit = 4

	p1 = ggplot(subset(swdf, child %in% childList & range < rangeLimit)) + facet_wrap(~child, ncol=1) + geom_segment(aes(x=minAgeInMonths,y=nu.child_mean, xend=maxAgeInMonths, yend=nu.child_mean), alpha=.2) + geom_segment(aes(y=nu.child_sd_low, yend = nu.child_sd_high, x=meanAgeInMonths, xend = meanAgeInMonths, col=child), alpha=.5)  + geom_point(aes(x=meanAgeInMonths, y=nu.child_mean, colour=child), size=.9) + geom_line(aes(x=meanAgeInMonths, y=nu.child_mean, colour=child)) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none") + xlab('Age in Months') + ylab(expression(paste('nu'))) + coord_cartesian(ylim=c(0,3)) + geom_hline(aes(yintercept=0)) + scale_x_continuous(breaks = c(12,24, 36, 48, 60))
    p1.build = ggplot_build(p1)
	
	p2 = ggplot(subset(swdf,  range < rangeLimit & child %in% childList)) + facet_wrap(~child, ncol=1) + geom_segment(aes(x=minAgeInMonths,y=mu0.child_mean, xend=maxAgeInMonths, yend=mu0.child_mean), alpha=.2) + geom_segment(aes(y=mu0.child_sd_low, yend = mu0.child_sd_high, x=meanAgeInMonths, xend = meanAgeInMonths, col=child), alpha=.5)  + geom_point(aes(x=meanAgeInMonths, y=mu0.child_mean, colour=child), size=.9) + geom_line(aes(x=meanAgeInMonths, y=mu0.child_mean, colour=child)) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none") + xlab('Age in Months') + ylab(expression(paste('mu'))) + geom_hline(aes(yintercept=0)) + scale_x_continuous(breaks = c(12,24, 36, 48, 60))

	p3 = ggplot(subset(swdf, range < rangeLimit & child %in% childList)) + facet_wrap(~child, ncol=1) + geom_segment(aes(y=eta_sd_low, yend = eta_sd_high, x=meanAgeInMonths, xend = meanAgeInMonths, col=child), alpha=1)  + geom_point(aes(x=meanAgeInMonths, y=eta_mean, colour=child), size=.9) + geom_line(aes(x=meanAgeInMonths, y=eta_mean, colour=child)) + theme(axis.line = element_line(colour = "black"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank(),
	    strip.background = element_blank(),
	    axis.text.x = element_text(colour = "black"),
	    axis.text.y = element_text(colour = "black"),
	    legend.position = "none") + xlab('Age in Months') + ylab('Grammatical Productivity') + coord_cartesian(ylim=c(0,.6)) + geom_hline(aes(yintercept=0)) + scale_x_continuous(breaks = c(24, 36, 48, 60))

	p4 = ggplot(subset(swdf, range < rangeLimit & child %in% childList)) + facet_wrap(~child, ncol=1) + geom_line(aes(x=meanAgeInMonths, y=overlap_mean), color='blue') + geom_point(aes(x=meanAgeInMonths, y=overlap_mean), color='blue', size=.9) + geom_line(aes(x=meanAgeInMonths, y=childEmpiricalOverlap), color="black") + geom_line(aes(x=meanAgeInMonths, y=childYangPredictedOverlap), col="red") + theme(axis.line = element_line(colour = "black"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank(),
	    strip.background = element_blank(),
	    axis.text.x = element_text(colour = "black"),
	    axis.text.y = element_text(colour = "black"),
	    legend.position = "none") + xlab('Age in Months') + ylab('Overlap') + geom_hline(aes(yintercept=.1)) + scale_x_continuous(breaks = c(12,24, 36, 48, 60)) + geom_errorbar(aes(x=meanAgeInMonths, ymin=overlap_p025, ymax=overlap_p975), col='blue', alpha=.5) + coord_cartesian(xlim=c( p1.build$panel$ranges[[1]]$x.range[1],  p1.build$panel$ranges[[1]]$x.range[2]))
	
	pdf(paste0(paths['graphPath'],'slidingWindow3column_',dataPrep,'.pdf'), width=6, height=12)
	grid.arrange(p1, p2, p4, ncol=3)    
	dev.off()
	
	return(grid.arrange(p1, p2, p4, ncol=3))
}

makeSplitHalfPlots = function(paths, dataPrep){
	# get the HPDs
	swdf = getModelResults(paths$tablePath,'linkingChange', dataPrep)				
 	childrenToKeep = subset(aggregate(index ~ child ,subset(swdf, nu.child_p975 <= 5 & nu.child_p025 <= 5), length), index ==2)$child
	
	swdf = subset(swdf, child %in% childrenToKeep)
	earliest = aggregate(meanAgeInMonths ~ child, swdf, min)
	swdf$child = factor(swdf$child, levels = as.character(earliest$child)[order(earliest$meanAgeInMonths)]) 		

	swdf$nu.child_sd_low = swdf$nu.child_mean - swdf$nu.child_sd
	swdf$nu.child_sd_high = swdf$nu.child_mean + swdf$nu.child_sd
	swdf$mu0.child_sd_low = swdf$mu0.child_mean - swdf$mu0.child_sd
	swdf$mu0.child_sd_high = swdf$mu0.child_mean + swdf$mu0.child_sd
	swdf$eta_sd_low = swdf$eta_mean - swdf$eta_sd
	swdf$eta_sd_high = swdf$eta_mean + swdf$eta_sd
	swdf$corpus = restoreCorpusInfo(as.character(swdf$child))
                   
	#get the full nu estimates			
	lfc = p(paths$completedModelPath, getFilenames(paths$completedModelPath,'linkingChange','singulars_LN', '.RData'), sep='')	
	nu.differences = do.call('rbind.fill', mclapply(lfc, function(x){
		getAllSplitHalfVarEstimates(x, 'nu.difference')	
	}, mc.cores=detectCores())) 
	nu.estimates = do.call('rbind.fill', mclapply(lfc, function(x){
		getAllSplitHalfVarEstimates(x, 'nu.child')	
	}, mc.cores=detectCores()))  
	nu.estimates$meanAgeInMonths = nu.estimates$meanAge / 30.5	     
	nu.estimates.short = subset(nu.estimates, child %in% childrenToKeep)   
	nu.differences.short = subset(nu.differences, child %in% childrenToKeep)   
	nu.differences.HPD = getNuDifferenceHPDs(nu.differences.short) 

	violinPlot = ggplot() + geom_hline(yintercept=0, col='black') + geom_violin(data=nu.differences.short, mapping=aes(y=nu.difference, x=as.factor(phase1_mean_age), col=corpus, fill=corpus), stat = "ydensity", position = "dodge", trim = TRUE,  scale = "count") +scale_x_discrete(labels = unique(nu.differences.HPD$child)[order(unique(nu.differences.HPD$phase1_mean_age))]) + theme_bw(base_size=12) + xlab('Child, Ordered By Age of Earliest Data') + ylab(expression(paste('Change in ',nu))) + theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")) + geom_point(data=nu.differences.HPD,mapping =aes(x=as.factor(phase1_mean_age), y=nu.difference_p9995,color=corpus)) + geom_point(data=nu.differences.HPD,mapping =aes(x=as.factor(phase1_mean_age), y=nu.difference_p0005,color=corpus)) + geom_segment(data=nu.differences.HPD,mapping =aes(x=as.numeric(as.factor(phase1_mean_age))-.5,xend=as.numeric(as.factor(phase1_mean_age))+.5 , y=nu.difference_mean, yend=nu.difference_mean, col=corpus))  + geom_segment(data=nu.differences.HPD,mapping =aes(x=as.numeric(as.factor(phase1_mean_age))-.3,xend=as.numeric(as.factor(phase1_mean_age))+.3 , y=nu.difference_p975, yend=nu.difference_p975, col=corpus))  + geom_segment(data=nu.differences.HPD,mapping =aes(x=as.numeric(as.factor(phase1_mean_age))-.3,xend=as.numeric(as.factor(phase1_mean_age))+.3 , y=nu.difference_p025, yend=nu.difference_p025, col=corpus)) + coord_cartesian(ylim=c(min(nu.differences.HPD$nu.difference_p005)*1.1, max(nu.differences.HPD$nu.difference_p9995)*1.1))

	pdf(paste0(paths['graphPath'],'splitHalf_violinplot_',dataPrep,'.pdf'), height=5, width=9) 
	print(violinPlot)
    dev.off()    			
    
    splitHalfTimeCourse = ggplot(swdf) +  geom_errorbarh(aes(x=meanAgeInMonths, xmax = maxAgeInMonths, xmin = minAgeInMonths, height = 0, y=nu.child_mean), col='gray90') + geom_point(aes(y=nu.child_mean, x=meanAgeInMonths, col=corpus))  +  xlab('Child Age In Months') + ylab(expression(paste('Inferred ',nu))) + theme_bw(base_size=12) + scale_x_continuous(breaks = c(24, 36, 48, 60)) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    strip.text.y = element_blank(),
    axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black")) +  geom_segment(aes(y=nu.child_sd_low, yend = nu.child_sd_high, x=meanAgeInMonths, xend = meanAgeInMonths, col=corpus), alpha=.4)+ geom_line(aes(x=meanAgeInMonths, y=nu.child_mean, group=child, col=corpus)) + coord_cartesian(y=c(0,2.2)) + scale_size_area() + theme(legend.position="none")  + stat_smooth(data = nu.estimates.short, aes(x = meanAgeInMonths, y = nu.child), se=F, alpha=.2, level = .99, method = "lm",formula="y~I(x^2)+x", colour='black', lty=2, size=1.1)

	pdf(paste(paths['graphPath'],'splitHalf_timecourse_',dataPrep,'.pdf'  , sep=''), height=5, width=7)
	print(splitHalfTimeCourse)	 
	dev.off()      		

	rlist = list()
	rlist[['childrenToKeep']] = childrenToKeep
	rlist[['splitHalfTimeCourse']] = splitHalfTimeCourse
	rlist[['violinPlot']] = violinPlot
	return(rlist)	
}
	

splitHalfComparison = function(paths){
	#compare the output of all of the data preparations	
	all_LC = getModelResults(paths[['tablePath']], 'linkingChange','.*')
	all_LC$dataPrep = sapply(strsplit(as.character(all_LC$file), '_'), function(x){paste0(x[9:10], collapse='_')})
	all_LC$corpus = restoreCorpusInfo(as.character(all_LC$child))	
	refset = subset(all_LC, dataPrep == 'singulars_LN' & splitHalf ==1)
	all_LC$child =factor(all_LC$child, levels =unique(refset$child[order(refset$minAgeInMonths)]))	
		
	phase1NuComparison = ggplot()  + geom_point(data = subset(all_LC, splitHalf == 1 & !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), y=nu.child_mean, color=corpus)) +facet_wrap(~ child) +geom_hline(yintercept=0, col='black') + xlab('Data Preparation')+ ylab(expression(paste('Phase 1 ',nu))) + theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),strip.background = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+ coord_cartesian(ylim=c( 0,4)) + geom_segment(data = subset(all_LC, splitHalf==1 & !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), xend=as.factor(dataPrep), y=nu.child_p005, yend= nu.child_p995, color=corpus))
	
	pdf(paste0(paths['graphPath'],'dataPreparationComparison_phase1nu.pdf'), width=9, height=9)
	print(phase1NuComparison)
	dev.off()

	phase2NuComparison = ggplot()  + geom_point(data = subset(all_LC, splitHalf == 2 & !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), y=nu.child_mean, color=corpus)) +facet_wrap(~ child) +geom_hline(yintercept=0, col='black') + xlab('Data Preparation')+ ylab(expression(paste('Phase 2 ',nu))) + theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),strip.background = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+ coord_cartesian(ylim=c( 0,4)) + geom_segment(data = subset(all_LC, splitHalf==2 & !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), xend=as.factor(dataPrep), y=nu.child_p005, yend= nu.child_p995, color=corpus))
	
	pdf(paste0(paths['graphPath'],'dataPreparationComparison_phase2nu.pdf'), width=9, height=9)
	print(phase2NuComparison)
	dev.off()  
	
	#get the difference in nu
	lfca = p(paths$completedModelPath, getFilenames(paths$completedModelPath,'linkingChange','.*', '.RData'), sep='')	
	all.nu.differences = do.call('rbind.fill', mclapply(lfca, function(x){
		cbind(getAllSplitHalfVarEstimates(x, 'nu.difference'), file=x)	
	}, mc.cores=detectCores())) 
	all.nu.differences$dataPrep = gsub('.RData','',sapply(strsplit(as.character(all.nu.differences$file), '_'), function(x){paste0(tail(x,2), collapse='_')}))
	
	all.nu.differences.short = subset(all.nu.differences, !(child %in% c('Eric')))	
	print(unique(subset(all.nu.differences, child == 'Abe')$phase1_mean_age))
	
	all.nu.differences.HPD = do.call('rbind',lapply(split(all.nu.differences.short, all.nu.differences.short$dataPrep),  function(x){
		rdf = getNuDifferenceHPDs(x)
		rdf$dataPrep = unique(x$dataPrep)
		return(rdf)
	}))

	refset = subset(all.nu.differences.HPD, dataPrep == 'singulars_LN')

	all.nu.differences.HPD$child =factor(all.nu.differences.HPD$child, levels =unique(refset$child[order(refset$phase1_mean_age)]))
	
	nuDifference =  ggplot()  + geom_point(data = subset(all.nu.differences.HPD,  !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), y=nu.difference_mean, color=corpus)) +facet_wrap(~ child) +geom_hline(yintercept=0, col='black') + xlab('Data Preparation')+ ylab(expression(paste('Phase 2 ',nu))) + theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),strip.background = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+ coord_cartesian(ylim=c( -4,4)) + geom_segment(data = subset(all.nu.differences.HPD, !(child %in% c('Eric'))), aes(x=as.factor(dataPrep), xend=as.factor(dataPrep), y=nu.difference_p005, yend= nu.difference_p995, color=corpus))	
		
	pdf(paste0(paths['graphPath'],'dataPreparationComparison_nuDifference.pdf'), width=9, height=9)
	print(nuDifference)
	dev.off()  
	
	rlist = list()
	rlist[['phase1NuComparison']] = phase1NuComparison
	rlist[['phase2NuComparison']] = phase2NuComparison
	rlist[['nuDifference']] = nuDifference	
	return(rlist)
}