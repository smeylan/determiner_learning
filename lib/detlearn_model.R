#If you are looking for the model specifcations, see /models. This file contains all of the functions for the coordination of running the datasets on a multicore machine or cluster.

getSubmissionStruct = function(corpora, analyses, morphologies, model, windowSize, skip, datasetName,inputDirName,outputDirName, stringsAsFactors=F){
	### create a data frame of models to run ### 

	rdf = expand.grid(corpora, analyses, morphologies, model, windowSize, skip, datasetName,inputDirName,outputDirName, stringsAsFactors=F)	
	names(rdf) = c('corpus','analysis','morphology','model','windowSize','skip','datasetName','inputDirName','outputDirName')
	
	#checks for undesirable parameter combinations	
	# -- There is no CHILDES markup for Speechome
	rdf = subset(rdf, !(corpus == 'Speechome' & datasetName == 'standard'))
	
	#-- the CHILDES markup for Thomas has many formatting issues
	rdf = subset(rdf, !(corpus == 'Thomas' & datasetName == 'standard'))
	
	#-- the Stanford POS tagger doesn't have
	rdf = subset(rdf, !(morphology == 'none' & datasetName %in%  c('LN','FN')))	
	
	# don't want to run an analysis of plurals on a dataset consisting only of singulars
	rdf = subset(rdf, !(morphology == 'singulars' & analysis == 'plurals'))	
		
	return(rdf)	
}

findIncompleteModels = function(directoryToCheck, analysis){
        ### look for missing model result files and generate the parameter dataframe necessary to generate them ###
		modelList = data.frame(corpus = character(0), analysis = character(0),morphology = character(0), model= character(0),windowSize= numeric(0), skip= numeric(0),datasetName= character(0),inputDirName= character(0),outputDirName= character(0))
		
        completionCheck = checkForModelCompletion(directoryToCheck,NULL)
        missingSlidingWindow = subset(completionCheck, !slidingWindowPresence)
        
        if(nrow(missingSlidingWindow) > 0){
        	missingSlidingWindow = missingSlidingWindow[,c('corpus','morphology', 'datasetName', 'slidingWindowPresence')] 
			missingSlidingWindow$model = 'slidingWindow'
        	missingSlidingWindow$complete = NULL
        	missingSlidingWindow$analysis = analysis
	        missingSlidingWindow$model='linking'
    	    missingSlidingWindow$windowSize = 1024
        	missingSlidingWindow$skip = 256
	        missingSlidingWindow$inputDirName = '00_input/'
    	    missingSlidingWindow$outputDirName = '01_output/'
        	missingSlidingWindow = missingSlidingWindow[,c('corpus','analysis','morphology', 'model','windowSize', 'skip','datasetName','inputDirName','outputDirName')]        	
        	modelList = rbind(modelList, missingSlidingWindow)
        }

        missingLinkingChange = subset(completionCheck, !changePresence)        
        if(nrow(missingLinkingChange) > 0){
        	missingLinkingChange = missingLinkingChange[,c('corpus','morphology', 'datasetName', 'changePresence')]
	        missingLinkingChange$model = 'linkingChange'        
	        missingLinkingChange$complete = NULL
	        missingLinkingChange$analysis = analysis
	        missingLinkingChange$model='linkingChange'
	        missingLinkingChange$windowSize = 'half'
	        missingLinkingChange$skip = NA
    	    missingLinkingChange$inputDirName = '00_input/'
        	missingLinkingChange$outputDirName = '01_output/'
        	missingLinkingChange = missingLinkingChange[,c('corpus','analysis','morphology', 'model','windowSize', 'skip','datasetName','inputDirName','outputDirName')]        	
        	modelList = rbind(modelList, missingLinkingChange)
        }
        return(modelList)
}

mclapply2 = function(..., parallelize=T){
    ### alias mclapply to a function that reverts to lapply if passed a flag ###
    if (!parallelize){
        lapply(...)
    } else {
        mclapply(..., mc.cores = detectCores())
    }
}

runDETLEARNmodel = function(corpus, analysis, morphology, model, withImputation, windowSize, skip, n.adapt=2000, n.update=2000, n.iter=5000,n.chains=5, thinning=5, countTrace=T,datasetName, inputDirName, outputDirName, parallelize){
	### run a model for a give dataset, according to the specified parameters ###
	#ARGUMENTS:
	#corpus: name of the corpus (i.e. "Speechome," "Thomas")
	#analysis: grammatical phenomena (i.e. "determiners", "plurals", 'determiners-baseline', 'determiners-imit')
	#morphology: treatment of the noun morphology ("all","none", "singulars")
		#all: keep plural separate from singular
		#none: make plural the same as singular
		#singulars: discard plurals (keep diminutives and agentive -er)
	#model: name of the JAGS model to use
	#withImputation: boolean, whether to include imputation variables in JAGS model
	#windowSize: number of child tokens to fit the model at each point; >=1024 or 'half' to split into the first and second half of the dataset
	#skip: how many tokens to increment from the previous model. Higher = less overlap
	#n.adapt: number of adaptive samples in Gibbs Sampling (GS)
	#n.update: number of updates in GS
	#n.iter: number of samples, before thinning in GS
	#thinning: thinning parameter in GS
	#countTrace: Boolean: record full model? If set to True, beware of running out of memory for high values of N.iter 
	#datasetName: what corpus extraction method was used ("LN", "FN", "standard")
	#inputDirName: directory path with trailing slash with inputData
	#outputDirName: directory path with trailing slash where RData files should be saved upon completion
	
	#### Throw errors for bad argument combinations
	if(windowSize == 'half' & !is.na(skip) ){stop('If analyzing split halves, skip must be NA')}
	if(is.na(skip) & !(windowSize == 'half' | windowSize <1)){stop('If skip is NA, must be analyzing split halves')}
	if(model =='linkingChange' &  !is.na(skip)){stop('Skip must be NA to run linkingChange')}
	if(countTrace & n.iter > 500){stop('If countTrace is on, the number of samples must be small')}
	if(windowSize < 1024 & windowSize > 1){stop('Model not recommended for datasets with fewer than 1024 tokens')}
	
	#set a lower limit for the size of datasets
	if (windowSize == 'half'){
	  lowerLimit = 256 # this was 0 before
	} else {
	  lowerLimit = windowSize
	}
	
	#### Corpus Preparation
	if(analysis %in% c('determiners','determiners-baseline','plurals','determiners-noImit', 'agentives')){
		if (corpus %in%  c('Brown', 'Providence','Kuczaj','Suppes','Bloom70','Sachs')){		
			dat.big = processCorpus(corpusPath = paste(inputDirName, corpus, "_",morphology, "_", datasetName, ".csv", sep=''), childSpeakerDesignation = 'CHI', caregiverList = c('MOT','FAT'),lowerLimit, analysis)
			adultNounAggFile = paste(inputDirName, 'USCHILDES_nouns_',morphology,'_', datasetName,'.csv', sep='') 	
			childData = split(subset(dat.big, speaker == 'child'),subset(dat.big, speaker == 'child')$child)	
		} else if (corpus %in% c('Manchester','Thomas')){	
			dat.big = processCorpus(corpusPath = paste(inputDirName, corpus, "_", morphology, "_", datasetName, ".csv", sep=''), childSpeakerDesignation = 'CHI', caregiverList = c('MOT','FAT'),lowerLimit, analysis)		
			adultNounAggFile = paste(inputDirName, 'UKCHILDES_nouns_',morphology, '_', datasetName,'.csv', sep='')
			childData = split(subset(dat.big, speaker == 'child'),subset(dat.big, speaker == 'child')$child)		
		} else if (corpus == 'Speechome'){	
		    dat.big = processCorpus(corpusPath = paste(inputDirName, corpus, "_", morphology, "_", datasetName, ".csv", sep=''), childSpeakerDesignation = 'child', caregiverList = c('caregiver'),lowerLimit, analysis)	    	
			adultNounAggFile = paste(inputDirName, 'SPEECHOME_nouns_',morphology, '_', datasetName,'.csv', sep='')		
			childData = list()
		    childData[[1]] = subset(dat.big, speaker == 'child')	    
		} else {
		  stop('corpus not recognized')
		}    
		if (analysis %in% c('determiners-noImit')){
			# Enforce the following checks on tokens in dat.big: throw away repeated uses of a token within a child utterance;
			# throw away child tokens that are imitations of tokens used by the parent in the previous line	
	
			if(corpus == 'Speechome'){       		     
				#Time-based imitation check:
				child_df = childWithoutImitations(dat.big, timeColumn='time', n=15)	
				child_df.norep = childWithoutRepetitions(child_df, timeColumn='time', n=15)
				parentData = subset(dat.big, speaker == 'parent')			
				parentData$imitated = FALSE
				parentData$repeated = FALSE
				
				bothDF = rbind(subset(child_df.norep, !repeated & !imitated), parentData)
				#bothDF = removeWordsWithInterveningMaterial(bothDF) #this is useful for checking numbers against Yang (2013) but not methodologically motivated
	      		bothDF = bothDF[order(bothDF$index),]
				dat.big = bothDF
				childData = split(subset(dat.big, speaker == 'child'),subset(dat.big, speaker == 'child')$child)	
			} else{
				if(!('Utt.number' %in% names(dat.big))){
					stop('input data must have indexed utterances in order to filter out repeated tokens (within child or imitations of the parent)')
				}
				dat.big.byChild = split(dat.big, dat.big$child)      
				dat.big = do.call('rbind',lapply(dat.big.byChild, function(df){																				
					child_df = childWithoutImitations(df, timeColumn='Utt.number',n=3)						
					child_df.norep = childWithoutRepetitions(child_df, timeColumn='Utt.number',n=3) #this is useful for checking numbers against Yang (2013) but not methodologically motivated
					parentSpeech = subset(df, speaker == 'parent')
					parentSpeech$imitated = F
					parentSpeech$repeated = F				
					
					bothDF = rbind(subset(child_df.norep, !repeated & !imitated), parentSpeech)        
	        		#bothDF = removeWordsWithInterveningMaterial(bothDF)
					bothDF = bothDF[order(bothDF$index),]
					return(bothDF)
				}))	
				childData = split(subset(dat.big, speaker == 'child'),subset(dat.big, speaker == 'child')$child)
	      if (!is.na(lowerLimit)){
	        childData = childData[sapply(childData, nrow) > lowerLimit]  #filter again after getting rid of imitations
	      }       
			}						
		}		
	} else {
		stop('Other grammatical features not yet implemented')	
	}
	
	#### Run the model
	runname =  paste(model,'Model',corpus,n.iter,windowSize,'skip', skip, ifelse(withImputation,'withImputation','withoutImputation'),morphology, datasetName, sep='_')
	adultNFreqs= read.csv(adultNounAggFile, stringsAsFactors=F)
	parentData = subset(dat.big, speaker == 'parent')
	childrenDataWithMatchedAdults = generateSubsetsForLinkingModel(childData, parentData, windowSize, skip)	

	if (model == 'linking'){		
		#produce consecutive windows; should return 0 counts for full list of nouns (parentData and childData)
		Nr.combined.observed.list = lapply(childrenDataWithMatchedAdults, function(x){getJointNr.observed(x,withImputation,adultNFreqs)})
		
		#build inputs for parallelization		
		jobInputs = lapply(1:length(Nr.combined.observed.list),function(x){return(list(input = Nr.combined.observed.list[[x]], wordFreqs = adultNFreqs, modelname = 'models/slidingWindowModel.bug', 	n.iter=n.iter, countTrace = countTrace, thinning=thinning, n.adapt=n.adapt, n.update = n.update, n.chains = n.chains, model='linking'))})	

	} else if (model == 'linkingChange'){
		cdma = do.call('rbind',childrenDataWithMatchedAdults) #rbind with all adult data
		
		#split again by child
		cdma = split(cdma, cdma$child)
	
		childrenDataWithMatchedAdults = lapply(cdma, function(df){
			df$splitHalf = as.numeric(as.factor(df$tokenStart)) 
			return(df)
		})
		
		cdma = do.call('rbind',childrenDataWithMatchedAdults)
		#divide by child and split half
		childrenDataWithMatchedAdults = split(cdma, list(cdma$child, cdma$splitHalf))
		
		Nr.combined.observed.list = lapply(childrenDataWithMatchedAdults, function(x){getJointNr.observed(x,withImputation,adultNFreqs)}) #on the first pass Nr.combined.list is of length 2 * number of children
		
		#merge again, and split by child
		Nr.combined = do.call('rbind', Nr.combined.observed.list)
		Nr.combined.observed.list = split(Nr.combined, Nr.combined$child)
  
		#then we pass the first and the second phase for each parent-child pair into jobInputs
		jobInputs = lapply(1:length(Nr.combined.observed.list),function(x){return(list(input = Nr.combined.observed.list[[x]], wordFreqs = adultNFreqs, modelname = 'models/splitHalfModel.bug', n.iter=n.iter, countTrace = countTrace, thinning=thinning, n.adapt=n.adapt, n.update = n.update, n.chains = n.chains, model='linkingChange'))})
	}
	
	linking.model.results = mclapply2(jobInputs, runWindow, parallelize=parallelize)
		
	#remove any large intermediate variables in the workspace
	childrenDataWithMatchedAdults = NULL
	jobInputs = NULL
	corpus.DF = NULL
	dat.big = NULL
	
	save(list = ls(envir = environment(), all.names = TRUE), 
	    file = paste(ensureTrailingSlash(outputDirName), runname,'.RData', sep=''), envir = environment())
}

runWindow = function(jobInput, maxReruns = 20){	
	### run an individual window ###
	t1 = try({ 			
		if(jobInput[['model']] == 'linkingChange'){
			df = jobInput$input
	          
	  		#pull out the indices that we need for the iterators in the JAGS model. See linkingChangeModel.bug for an explanation of the indices                    		  		
	  		phase1_endPart1 = nrow(subset(df, splitHalf ==1 & part == 1))
	  		phase1_endPart2 = nrow(subset(df, splitHalf ==1 & part <= 2))
	  		phase1_endPart3 = nrow(subset(df, splitHalf ==1 & part <= 3))
	  		phase1_endPart4 = nrow(subset(df, splitHalf ==1 & part <= 4))
	  		phase2_endPart1 = nrow(subset(df, splitHalf ==2 & part == 1))
	  		phase2_endPart2 = nrow(subset(df, splitHalf ==2 & part <= 2))
	  		phase2_endPart3 = nrow(subset(df, splitHalf ==2 & part <= 3))
	  		phase2_endPart4 = nrow(subset(df, splitHalf ==2 & part <= 4))
	          
	  		m = jags.model(jobInput$modelname, n.chains=5, n.adapt=2000, data=list(
	  		  n.child1 = subset(df, splitHalf == 1)$n.child,
	  		  r.child1 = subset(df, splitHalf == 1)$r.child,
	  		  n.parent1 = subset(df, splitHalf == 1)$n.parent,
	  		  r.parent1 = subset(df, splitHalf == 1)$r.parent,
	  		  phase1_endPart4 = phase1_endPart4,
	  		  phase1_endPart3 = phase1_endPart3,
	  		  phase1_endPart2= phase1_endPart2,            
	  		  phase1_endPart1 = phase1_endPart1,
	  		  N.parent1 = subset(df, splitHalf == 1)$N.parent,
	  		  n.child2 = subset(df, splitHalf == 2)$n.child,
	  		  r.child2 = subset(df, splitHalf == 2)$r.child,
	  		  n.parent2 = subset(df, splitHalf == 2)$n.parent,
	  		  r.parent2 = subset(df, splitHalf == 2)$r.parent,
	  		  phase2_endPart4 = phase2_endPart4,
	  		  phase2_endPart3 = phase2_endPart3,
	  		  phase2_endPart2= phase2_endPart2,            
	  		  phase2_endPart1 = phase2_endPart1,
	  		  N.parent2 = subset(df, splitHalf == 2)$N.parent)
	  		)
	  		
	  		if (jobInput$countTrace){ 
	  			varsToSample = c("mu0.adult1","mu0.adult2", "nu.adult1","nu.adult2","mu0.child1","mu0.child2","nu.child1","nu.child2","eta1","eta2", "r.child.simulated1","r.child.simulated2", 'N.parent1','N.parent2', 'R.parent1', 'R.parent2')
			} else {
				varsToSample = c("mu0.adult1","mu0.adult2","nu.adult1","nu.adult2","mu0.child1","mu0.child2","nu.child1","nu.child2","eta1","eta2", "r.child.simulated1","r.child.simulated2")
			}
			checkVars = c("mu0.child1","mu0.child2","nu.child1","nu.child2") #this was "nu.child" in old runs; not sure what that was checking
			
		} else if (jobInput[['model']] == 'linking'){						

			#See slidingWindow.bug for an explanation of the indices                    		  		

			df = jobInput$input
			endPart1 = nrow(subset(df, part == 1))
			endPart2 = nrow(subset(df, part <= 2))
			endPart3 = nrow(subset(df, part <= 3))
			endPart4 = nrow(subset(df, part <= 4))		  		          
	
			m = jags.model(jobInput$modelname, n.chains=jobInput$n.chains, n.adapt=jobInput$n.adapt, data=list(
	            n.child = df$n.child,
	            r.child = df$r.child,
	            n.parent = df$n.parent,
	            r.parent = df$r.parent,
	            endPart4 = endPart4,
	            endPart3 = endPart3,
	            endPart2 = endPart2,            
	            endPart1 = endPart1,
	            N.parent = df$N.parent) 
	        )     
	        if (jobInput$countTrace){
				varsToSample = c("mu0.adult","nu.adult","mu0.child","nu.child","eta", "r.child.simulated", 'N.parent', 'R.parent')
			} else {
				varsToSample = c("mu0.adult","nu.adult","mu0.child","nu.child","eta", "r.child.simulated")
			}	
			checkVars = c("mu0.child","nu.child") 	        
        }     

		update(m,jobInput$n.update) 		
		coda.res <- coda.samples(m, checkVars,n.iter=jobInput$n.iter, thin=jobInput$thinning)
					
		#Run until convergence, or until the max number of repetitions
		rerunCount = 0 					
  		while(!checkForConvergence(coda.res, checkVars) & rerunCount < maxReruns){
			update(m,1000)		  										
          	coda.res <-coda.samples(m,checkVars, n.iter=jobInput$n.iter, thin=jobInput$thinning)
          	rerunCount = rerunCount + 1
        }
        				
		if(rerunCount == maxReruns){
			rlist = list()
			rlist['error'] = 'Model failed to converge'
			return(rlist)
		}   
              
		#otherwise we now know that we have updated the model sufficiently, and we can now take the full sample we will use
		coda.res <- coda.samples(m, varsToSample, n.iter=jobInput$n.iter, thin=jobInput$thinning)		  		             
	})        
	if (inherits(t1, 'try-error')){
		rlist = list()
		rlist['error'] = t1
		return(rlist)	
	} else {
		return(coda.res) 	
	}
}

checkForConvergence =  function(lmItem, checkVars){
	### check if the chains in a JAGS model have converged. Only checks variables listed in checkVars ###
	gd = gelman.diag(lmItem[,checkVars])
	return(all(data.frame(gd[1])$psrf.Upper.C.I. < 1.1))
}

generateSubsetsForLinkingModel = function (childData, parentData, windowSize,skip){	
	### given child data and parent data, split by window size, observing skip if the window size isn't 'half' ###
	if (windowSize == 'all'){
		childOnly = childData
	} else if (windowSize < 1){
		childOnly = unlist(lapply(childData, function(df){
			oneChild = list()
			numSplit = ceiling(nrow(df)*windowSize)
			oneChild[[1]] = df[1:numSplit-1,] 
			oneChild[[1]]$tokenStart = 1
			oneChild[[2]] = df[numSplit:nrow(df),]
			oneChild[[2]]$tokenStart = numSplit
			return(oneChild)
		}), recursive=F)
	} else {				      
	    getWindowsForSingleDataset= function(singleDataset, windowSize, skip){      	    	
	    	if(windowSize == 'half'){
	      		windowSize = floor(nrow(singleDataset)/2)	
	      		skip = windowSize	      	
	     	}	     
	 		#this line is potentially problematic
	 		#browser()
	 		windowStarts = seq(from=1, to=nrow(singleDataset)-windowSize+1, by= skip) 
	 		windows = lapply(windowStarts, function(x){return(cbind(singleDataset[x:(x+(windowSize-1)),],tokenStart = x))})	
	 		return(windows)
	    }	   
	    childOnly = unlist(lapply(childData, function(x){getWindowsForSingleDataset(x, windowSize=windowSize, skip=skip)}), recursive=F)
	}    
    #then augment the childOnly windows with the adult data up to and including the highest indexed token
    augmentWithParentUtterances = function(x,parentData){
      maxIndex = max(x$index)
      adultAugment = subset(parentData, child == unique(x$child) & speaker == 'parent' & index < maxIndex)
      rdf = rbind.fill(x, adultAugment)      
      rdf$tokenStart = unique(x$tokenStart)
      return(rdf)   
	 }
	return(lapply(childOnly, function(x){augmentWithParentUtterances(x,parentData)}))
}


getJointNr.observed = function(dat, withImputation, adultNFreqs){     
	### given a dataframe, dat, go from a longform representation of the tokens to a summary with n trials and k
	### successes for each noun, tabulated separately for the child and the caregiver. Also impute the adult
	### token count
	 		
	c.N.observed = aggregate(one ~ word, subset(dat, speaker == 'child'),sum)
	names(c.N.observed) = c('word','n.child')
	c.r.observed = aggregate(one ~ word, subset(dat, speaker == 'child' & def),sum)
	names(c.r.observed) = c('word','r.child')  
	childUses = merge(c.N.observed, c.r.observed, by='word', all.x=T, all.y=T)
  
	p.N.observed = aggregate(one ~ word, subset(dat, speaker == 'parent'),sum)
	  names(p.N.observed) = c('word','n.parent')
	  p.r.observed = aggregate(one ~ word, subset(dat, speaker == 'parent' & def),sum)
	  names(p.r.observed) = c('word','r.parent')
	  parentUses = merge(p.N.observed, p.r.observed, by='word', all.x=T, all.y=T)  
	  
	  Nr.observed = merge(childUses, parentUses, all.x=T, all.y=T, by='word')
	  Nr.observed[is.na(Nr.observed)] = 0
	  
	  #sanity check: make sure that determiner counts <= noun counts
	  if(any(Nr.observed$r.caregiver > Nr.observed$N.caregiver) | any(Nr.observed$r.child > Nr.observed$N.child )){stop('Determiner count exceeds noun count; that shit cray')}
	
	  #if it's with imputation, insert zero counts for all nouns not observed in this window
	  	if(withImputation){
	Nr.observed = merge(Nr.observed, adultNFreqs, by = 'word', all.y=T)	
	Nr.observed[is.na(Nr.observed)] = 0	
	Nr.observed$count =  NULL
	}	
	
	  #set values for the window
	  if(length(unique(dat$percentileGroup)) == 1){  		
	    Nr.observed$percentileGroup = unique(dat$percentileGroup)	
	  }		
	  Nr.observed$child = paste(unique(dat$child), collapse = ', ')
	  if('proportion' %in% names(dat)){
	    Nr.observed$proportion = unique(dat$proportion)	
	  }
	  if('tokenStart' %in% names(dat)){
	    Nr.observed$interval = unique(dat$tokenStart)	
	  }
	  if('splitHalf' %in% names(dat)){
	    Nr.observed$splitHalf = unique(dat$splitHalf)	
	  }
	  if('age' %in% names(dat)){
	    Nr.observed$minAge = min(subset(dat, speaker == 'child')$age, na.rm=T)
	    Nr.observed$maxAge = max(subset(dat, speaker == 'child')$age,na.rm=T)
	    Nr.observed$meanAge = mean(subset(dat, speaker == 'child')$age, na.rm=T)	
	  }  
	 
	  adultNFreqsAdjusted = subset(adultNFreqs, !(word %in% c('xx','xxx','yy','yyy')))  
	  Q = Nr.observed$maxAge[1]*822
	  adultNFreqsAdjusted$N.parent = trunc(adultNFreqsAdjusted$count / sum(adultNFreqsAdjusted$count) * Q)  
	  
	  Nr.observed.aug = merge(Nr.observed, adultNFreqsAdjusted[,c('word','N.parent')],all.x=T, all.y=T)  
	  Nr.observed.aug[is.na(Nr.observed.aug)] = 0
	  part1 = subset(Nr.observed.aug, n.child > 0 & n.parent == 0) #nouns used only by the child
	  if (nrow(part1) > 0){part1$part = 1}else{part1$part = character()}
	  part2 = subset(Nr.observed.aug, n.child > 0 & n.parent > 0) #nouns used by both
	  if (nrow(part2) > 0){part2$part = 2}else{part2$part = character()}
	  part3 = subset(Nr.observed.aug, n.child == 0 & n.parent > 0) #nouns used by the parent but not the child
	  if (nrow(part3) > 0){part3$part = 3}else{part3$part = character()}
	  part4 = subset(Nr.observed.aug, n.child == 0 & n.parent == 0) #nouns used only by other caregivers          
	  if (nrow(part4) > 0){part4$part = 4}else{part$part = character()}
	  newInputFrame = rbind(part1, part2, part3, part4)
	  if (nrow(newInputFrame) != nrow(Nr.observed.aug)){stop('reordered frame is of a different length than the input')}
	  return(newInputFrame)
} 

childWithoutImitations = function(df,timeColumn, n){
	### remove any child tokens that are within n seconds of an adult usage, using timeColumn ###	
	if(timeColumn == 'time'){ 
		interval = n #time in seconds as the threshold for a repeat
		df$numericTime = as.numeric(strptime(df$time, '%Y-%m-%d_%H:%M:%S'))		
		df$detWord = paste(df$det, df$word)		
		child_df = subset(df, speaker == 'child')
		parent_df = subset(df, speaker == 'parent')
								
		child_df$imitated = F
		for (i in 1:nrow(child_df)){
			candidateToken = child_df$detWord[i]
			candidateTime = child_df$numericTime[i]
			precedingParent = subset(parent_df, detWord == candidateToken & numericTime < candidateTime & numericTime > (candidateTime - interval))
			child_df$imitated[i] = nrow(precedingParent) > 0
		}
		
		child_df$numericTime = NULL
		child_df$detWord = NULL
		return(child_df)		
		
	} else {	
		childNounUses = subset(df, speaker == 'child')
		    
		caregiverNounUses = subset(df, speaker != 'child')		
    cnu.all = do.call('rbind', lapply(1:n, function(n){ #build data frames for records with utterance indices + 1, indices +2... indices + n
      df = caregiverNounUses
      df$Utt.number = df[[timeColumn]]+n #to correspond with utterance +n
      return(df)
    }))
    
		imits = merge(childNounUses, cnu.all, by=c(timeColumn,'det','word','file'))
		
		childNounUses$imitated = childNounUses$index %in%  imits$index.x
		return(childNounUses)	
	}
}

childWithoutRepetitions = function(child_df, timeColumn, n){
	### remove any child tokens that are within n seconds of an adult usage, using timeColumn ###	
	if(timeColumn == 'time'){
		interval = n #time in seconds as the threshold for a repeat
		child_df$numericTime = as.numeric(strptime(child_df$time, '%Y-%m-%d_%H:%M:%S'))
		child_df$detWord = paste(child_df$det, child_df$word)

		child_df$repeated = F
		for (i in 1:nrow(child_df)){
			candidateToken = child_df$detWord[i]
			candidateTime = child_df$numericTime[i]
			precedingChild = subset(child_df, detWord == candidateToken & numericTime < candidateTime & numericTime > (candidateTime - interval))
			child_df$repeated[i] = nrow(precedingChild) > 0
		}

		child_df$detWord = NULL
		child_df$numericTime = NULL
		return(child_df)

	} else {
		chi.all = do.call('rbind',lapply(1:n, function(n){
      		df = child_df
			df$Utt.number = df[[timeColumn]]+n #to correspond with utterance +n
			return(df)
        }))

		repetitions = merge(child_df, chi.all, by=c(timeColumn,'det','word','file'))
	    child_df$repeated = child_df$index 	%in%  repetitions$index.x
		return(child_df)
	}
}

canary = function(index){
	### debugging function: 'chirp' by writing current memory usage to canary.txt
	write.csv(system('cat /proc/meminfo', intern=T), file=paste0('canary',index,'.txt'))
}

removeWordsWithInterveningMaterial = function(df){
	### remove any tokens that have material between the determiner and the noun ###	
	if ('numIntermediateWords' %in% names(df)){ #Stanford Tagger
		return(subset(df, numIntermediateWords == 0))  
	} else if('sent_mor' %in% names(df)){
		return(df[sapply(1:nrow(df), function(x){ #CLAN tagger: check that the unedited string is in the MOR tier
      as.logical(length(grep(paste0('det\\|',df$determiner[x], ' ', 'n\\|',df$word[x]), df$sent_mor[x])))}),])	      
	} else{ 
		stop('no column to use to look for the number of intervening words')
	}  
}

processCorpus = function(corpusPath, childSpeakerDesignation, caregiverList,lowerLimit,analysis){  
	### pre-process the corpus: identify the child, caregiver(s) and standardize the names of columns ###
	corpusDF = read.csv(corpusPath, stringsAsFactors=F)  	
  names(corpusDF)[which(names(corpusDF) == 'noun')] = 'word'  
	corpusDF  = subset(corpusDF, !(word %in% c('xx','yy','xxx','yyy'))) #remove all unknown words
  
	#subset to children with more than the minimum number of tokens, lowerLimit
	childOnly = subset(corpusDF, speaker == childSpeakerDesignation)
    childOnly.split = split(childOnly , childOnly$child)
    childCounts = lapply(childOnly.split, nrow)
    targetChildren = names(childCounts)[childCounts > lowerLimit]
	dat.big = subset(corpusDF, child %in% targetChildren & speaker %in% c(childSpeakerDesignation, caregiverList))      	
 	dat.big$speaker <- factor(ifelse(dat.big$speaker %in% childSpeakerDesignation,"child","parent"))   
	dat.big$index = c(1:nrow(dat.big)) #in case tokens get reordered-- order is preserved as long as it is right when read in   
	dat.big$one = 1 #for aggregations	
  if (analysis %in%  c('determiners', 'determiners-baseline','determiners-noImit')){
    names(dat.big)[which(names(dat.big) == 'modifier')] = 'determiner'
    dat.big$def = dat.big$determiner == 'the' 
    dat.big$det = dat.big$determiner        
  } else if (analysis == 'plurals'){
    dat.big$def = dat.big$plural #use plurals as the coinflip
  } else if (analysis == 'agentives'){
    dat.big$def = dat.big$agt == 1
  }
  
    return(dat.big)
}
