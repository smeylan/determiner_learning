processAllCha = function(folder, functionName, extension){
	rawFileList = list.files(folder, recursive=T)
	fileList = paste(ensureTrailingSlash(folder),rawFileList[grepl(extension, rawFileList)], sep='')	
	childUsages = do.call('rbind', mclapply(fileList, function(x){processFile(x, functionName)}, mc.cores=detectCores()))	
	return(childUsages)
}

processFile = function(inputFile, functionName){
	startTime = Sys.time()
	t1 = try({		
		corpus = read.CLAN.file(inputFile)			
		print(paste('File Read:', Sys.time() - startTime))				
	})
	if (inherits(t1, 'try-error')){
		#write the df to a file
		dtw = data.frame(filename = inputFile, stage='read-in')
		write.table(dtw,'error_df.csv', append=T, row.names=F, sep=",")	
		#endTime = Sys.time()
		#elapsed = endTime - startTime
		#print(paste('Failed to process ', inputFile, ' in ', elapsed, '; reading error', sep=''))
		return(NULL)	
	} 	 
	t2 = try({ 
		if (!any(names(corpus) %in% c('MOR','mor'))){
			dtw = data.frame(filename = inputFile, stage='mor_extraction')
			write.table(dtw,'error_df.csv', append=T, row.names=F, sep=",")	
			
			#endTime = Sys.time()
			#elapsed = endTime - startTime
			#print(paste('Failed to process ', inputFile, ' in ', elapsed, '; no MOR tier', sep=''))
		} else {
			#innerLoopStart = Sys.time()
		  withPunctuation = do.call('rbind', lapply(c(1:dim(corpus)[1]), function(x){functionName(corpus[x,])}))
		  	#innerLoopEnd = Sys.time()
		  	#print(paste('Inner loop:', innerLoopEnd - innerLoopStart))
		  	#browser()						
			withPunctuation$file = inputFile
      		if ('modifier' %in% names(withPunctuation)){
			  	return(withPunctuation[!is.na(withPunctuation$modifier),])
      		} else {
      			return(withPunctuation)
      		}
      		#endTime = Sys.time()
			elapsed = endTime - startTime
			#print(paste('Successfully processed ', inputFile, ' in ', elapsed, sep=''))	
		}	
	})
	if (inherits(t2, 'try-error')){
		#write the df to a file
		dtw = data.frame(filename = inputFile, stage='processing')
		write.table(dtw,'error_df.csv', append=T, row.names=F, sep=",")	
		#endTime = Sys.time()
		#elapsed = endTime - startTime
		#print(paste('Failed to process ', inputFile, ' in ', elapsed, '; processing error', sep=''))
		return(NULL)	
	}		
}

processNounReferences = function(sent){ 
	morArray = strsplit(as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$mor))), ' |~|-POS')[[1]]
	#morArray = morArray[unlist(lapply(morArray,function(x){as.logical( length( grep('\\|',x)))}))]
	
	if (is.null(sent$gra)){
		sent$gra = sent$xgr		
	}
	
	graArray = strsplit(as.character(gsub('\t|\n', ' ', sent$gra)), ' ')[[1]]
	graArray = graArray[unlist(lapply(graArray,function(x){as.logical( length( grep('\\|',x)))}))]	
	
	#glossArray = strsplit(as.character(sent$Gloss), ' ')[[1]] 
	
	if(any(sapply(morArray, is.na))){
		
	} else {				
		if(length(morArray) !=  length(graArray)){
			#browser()
			warning('gra and mor tiers must be of same length')
			#we don't want to let these sentences into the corpus		
			#dtw= data.frame(mor = paste(morArray, collapse=" "), gra = paste(graArray, collapse=","))
			#cannot write if we are inside of an mclapply
			#write.table(dtw,'morGraMismatch.csv', append=T, row.names=F, sep=",")	
		}	
		else{
			allNouns = do.call('rbind',lapply(c(1:length(morArray)), function(x){getNounModifiers(x, mor = morArray, gra = graArray)}))			
			if(!is.null(allNouns)){	
				allNouns$speaker = sent$Speaker
				allNouns$sent_gra = as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$gra)))
				allNouns$sent_mor = as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$mor)))
				allNouns$age = ageToDays(sent$Age)
				allNouns$Utt.number = sent$Utt.Number
				return(allNouns)
			}
		}
	}
}


processNounPlurality = function(sent){ 
  morArray = strsplit(as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$mor))), ' |~')[[1]]
  
  if (is.null(sent$gra)){
		sent$gra = sent$xgr		
  }
  
  if(is.na(morArray)){
    
  } else {				
    allNouns = do.call('rbind',lapply(c(1:length(morArray)), function(x){getNouns(x, mor = morArray)}))
      
    if(!is.null(allNouns)){	
        allNouns$speaker = sent$Speaker
        allNouns$sent_gra = as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$gra)))
        allNouns$sent_mor = as.character(gsub(' +',' ', gsub('\t|\n', ' ', sent$mor)))
        allNouns$age = ageToDays(sent$Age)
        return(allNouns)
      
    }
  }
}

ageToDays = function(age){
	ageParts = strsplit(age, ';')[[1]]
	return(ceiling((12*30.5*as.numeric(ageParts[1])) + as.numeric(ageParts[2])*30.5))	
}

getNounModifiers = function(x, mor, gra){
	#x is the index of the item in question
	#mor is the morphology tier of the entire sentence
	#gra is the syntax tier of the entire sentence
	#tried to use gloss for some analyses but its length is very different than the morphology or the xgra tier
	
	if(strsplit(strsplit(mor[x], '\\|')[[1]][1], ':')[[1]][1] == 'n'){
		#then it's a noun
		noun = strsplit(mor[x], '\\|')[[1]][2]
		#noun = strsplit(gloss,' ')[[x]]
		noun_index = strsplit(gra[x], '\\|')[[1]][1]
				
		#get everything modifying the noun			
		modifierIndices = which(unlist(lapply(gra, function(x){strsplit(x, '\\|')[[1]][2] == noun_index})))
		
		if (length(modifierIndices > 0)){
			modifier = unlist(lapply(mor[modifierIndices], function(x){strsplit(x,'\\|')[[1]][2]}))
			
			nmods = data.frame(modifier)
			
			nmods$noun = noun
			nmods$pos = unlist(lapply(mor[modifierIndices], function(x){strsplit(x,'\\|')[[1]][1]}))
			nmods$relation = unlist(lapply(gra[modifierIndices], function(x){strsplit(x,'\\|')[[1]][3]}))		
		return(nmods)
		}	
	}	
}	 


getNouns = function(x, mor){
  #x is the index of the item in question
  #mor is the morphology tier of the entire sentence
  #gra is the syntax tier of the entire sentence
  
  if(strsplit(strsplit(mor[x], '\\|')[[1]][1], ':')[[1]][1] == 'n'){
    #then it's a noun
    noun = strsplit(mor[x], '\\|')[[1]][2]    
    return(data.frame(noun))
  }	
}

 
sharedNounFilter = function(childrenNounModifiers,detsOnly){
  #subset to MOT, FAT, and CHI

	childrenNounModifiers = subset(childrenNounModifiers,speaker %in% c('MOT','FAT','CHI'))	
	childrenNounModifiers$one = 1 
	
	#add a child name using the filename
	
	childrenNounModifiers$child = sapply(gsub('[0123456789].*$', '', gsub('.*/','',childrenNounModifiers$file)), simpleCap)
	
	#remove compound nouns
	childrenNounModifiers = childrenNounModifiers[!grepl('\\+', childrenNounModifiers$noun),]
	
	#filter to nouns with an orthographic length > 2 characters
	childrenNounModifiers = childrenNounModifiers[nchar(childrenNounModifiers$noun) > 2,]
	
	#remove capitalized nouns
	childrenNounModifiers = childrenNounModifiers[!(sapply(childrenNounModifiers$noun, function(x){substr(toupper(x),1,1) == substr(x,1,1)})),]
	
	#unify plural markers		
	childrenNounModifiers$noun = gsub('&PL','-PL',childrenNounModifiers$noun)
	
	#remove non-plural gerunds and progressives
	childrenNounModifiers = childrenNounModifiers[!grepl('GERUND$',childrenNounModifiers$noun),]
	childrenNounModifiers = childrenNounModifiers[!grepl('PROG$',childrenNounModifiers$noun),]
	
	#remove pos markers from the roots
	childrenNounModifiers$noun = gsub('&.*[-$]','-',childrenNounModifiers$noun) 
	childrenNounModifiers$noun = gsub('-$','',childrenNounModifiers$noun) 

  #child name should be capitalized
	childrenNounModifiers$child = capwords(childrenNounModifiers$child)
  
	#mother, father, and child
	if(detsOnly){
		detsOnlyDF = subset(childrenNounModifiers, relation == 'DET' & pos == 'det')
	
		detsOnlyDF$det = sapply(detsOnlyDF$modifier, function(x){strsplit(as.character(x),'&')[[1]][1]})
		
		#a, an, the
		detsOnlyDF = subset(detsOnlyDF, det %in% c('a','an','the'))
		detsOnlyDF$def = detsOnlyDF$det  == 'the'
		return(detsOnlyDF)
	} else {
		return(childrenNounModifiers)
	}	
}


capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)), {s <- substring(s, 2); if(strict) tolower(s) else s},
  sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

writeDataFiles = function(dataset, corpusName, outputDir){
	#take morphological subsets and write them as separate CSVs
  #all
  write.csv(dataset, paste(outputDir,corpusName,'_all_standard.csv', sep=''), row.names=F)  
  #none
  dataset_none = dataset
  dataset_none$noun = sapply(dataset_none$noun, function(x){strsplit(x,'\\-')[[1]][1]})
  write.csv(dataset_none, paste(outputDir,corpusName,'_none_standard.csv', sep=''), row.names=F)
  #singulars
  dataset_singulars = dataset[!grepl('\\-PL', dataset$noun),]
  write.csv(dataset_singulars, paste(outputDir,corpusName, '_singulars_standard.csv', sep=''), row.names=F)  
}

read.CLAN.file <- function(f) {
	
	###############################################################
## Melissa's R converter for CHILDES-CLAN files
##
## This is an R function that takes a single CLAN-formatted
## corpus file (*.cha) and returns a dataframe.  Takes <30s
## per file, around 10 min for the Eve(Brown) Corpus.  
## Could be faster :) 
##
## The dataframe has a row for each sentence in the corpus, with the 
## following columns:
##
## From the file header & child participant annotation:
## FileName
## Participants (A string with the full list)
## Language 
## Corpus
## Age (child's)
## Gender (ditto)
## Date
## File.Situation
## 
## From the line itself
## Speaker
## Verbatim sentence
## Tiers: any that we find!
## 
## Also includes a gloss calculated from the utterance line, which
## gets rid of clarification sequences ("a bobby [= a thin bobbypin]" -> "a bobby")
## replaces 1-word glosses ("dunno [: don't know]" -> "don't know"), and
## cleans up various CHILDES-internal markup.  Ideally, this yields a gloss
## with the same number of words as the mor line. 
##
## This gloss is designed for presenting sentences to adult readers, though
## the form given may still want some processing (deleting sentence-final
## space, replacing internal "." with ",", continuation "+...") - these are left
## as-is here for ease of alignment against the mor and other tiers.
## Watch out for "," '"' and "'" when converting to/from csv
###############################################################
	
	tmp <- readLines(f)
	#print(f)
	#Cycle through utterances and make a line for each.
	alltext <- paste(tmp, collapse="\n")	
	utts = sapply(unlist(strsplit(alltext, "\n\\*", perl=T)), function(x){ifelse(substr(x, nchar(x), nchar(x)) == '\n',x, paste(x, '\n', sep=''))})
	utts <- utts[-1]
	tierlist <- lapply(utts, get_utt_info)
	data <- do.call('rbind.fill', tierlist)
		
	#Collect the data that will be appended to every line
	data$Filename <- f
	
	p <- grep("@Participants", tmp, fixed=TRUE)
	data$Participants <- unlist(strsplit(tmp[p], "\t"))[2]
	
	p <- grep("@Date", tmp, fixed=TRUE)
	data$Date <- unlist(strsplit(tmp[p], "\t"))[2]
	
	p <- grep("@Situation", tmp, fixed=TRUE)
	data$File.Situation <- unlist(strsplit(tmp[p], "\t"))[2]
	
	p <- grep("Target_Child", tmp, fixed=TRUE)
	chiline <- tmp[p[2]]
	chidata <- unlist(strsplit(chiline, "[|]"))
	data$Language <- substr(chidata[1], 6,9)
	data$Corpus <- chidata[2]
	data$Age <- chidata[4]
	data$Gender <- chidata[5]
	
	#Add utt line numbers
	data$Line.No 

		
	#Get rid of some yucky processing columns we don't want
	data$t.as.matrix.fields.. <- NULL

	#xnums <- as.numeric(gsub("[^0-9]*[0-9]*[^0-9]*[0-9]*[^0-9]*[0-9]*[^0-9]+", "", names(data), perl=T)) 		# what a hack
	#for(x in min(xnums, na.rm=T):max(xnums, na.rm=T)) {
	#	xname <- paste("X", x, sep="")
	#	data <- data[,!(names(data) %in% xname)]

	#}
	
	#Make sure row names are preserved!
	data$Utt.Number <- row.names(data)
	
	#Return
	data
} #End read.CLAN.file

get_utt_info <- function(u){

	#Divide the line into individual utterances & tiers
	fields <- unlist(strsplit(u, "[%]"))
	#Make a dataframe
	myrow <- data.frame(t(as.matrix(fields)))
	
	#Add utterance info
	myrow$Speaker <- substr(fields[1], 1,3)
	myrow$Verbatim <- substr(fields[1], 6,nchar(fields[1])-1)
	
	#Add info from any tiers, as they appear in the file
	if (length(fields) > 1){
		for (j in 2:length(fields)){
			tier <- data.frame(substr(fields[j], 6,nchar(fields[j])-1))
			names(tier) <- c(substr(fields[j], 1,3))
			myrow <- cbind(myrow, tier)
		}
	}
	
	#Some extra work: get the line as spoken, with glosses replaced
	#...This is an adult-language, human-readable version of the utterance
	
	myrow$Gloss <- NA
	#First, find & replace sequences like this: "dunno [: don't know]" -> "don't know"
	words <- unlist(strsplit(myrow$Verbatim, " "))
	if (length(words) == 0){
		words <- c("")
	}
	
	words <- unlist(strsplit(words, "\t"))
	if (length(words) == 0){
		words <- c("")
	}
	
	words <- unlist(strsplit(words, "\n"))
	if (length(words) == 0){
		words <- c("")
	}
	
	w <- 1
	wmax <- length(words) + 1
	while (w < wmax){
		#Did we hit a gloss sequence?
		if ((words[w]=="[:")|(words[w]=="[=?")){
			#Find where the gloss ends, then clean up
			closebracket <- grep("]",  words[w:length(words)], fixed=TRUE)[1] + (w-1)
			words[w-1] <- ""
			words[w] <- ""
			words[closebracket] <- substr(words[closebracket], 1, nchar(words[closebracket])-1)
		}
		w <- w + 1	
	}
	
	#Next, find & replace clarification/elaboration sequences like this: "a bobby [= a thin bobbypin]" -> "a bobby"
	w <- 1
	wmax <- length(words) + 1
	while (w < wmax){
		#Did we hit a gloss sequence?
		if ((substr(words[w],1,1) == "[")){
			#Find where the gloss ends, then clean up
			closebracket <- grep("]",  words[w:length(words)], fixed=TRUE)[1] + (w-1)
			goo <- closebracket
			for (v in w:closebracket){
				words[v] <- ""
			}
		}
		w <- w + 1	
	}
	
	#Next, delete internal notation we don't need here
	#"[()<>&@:]"
	words <- as.vector(mapply(gsub, "[()<>&:]","",words))
	
	#Remove sentence-internal periods!
	words[1:(length(words)-1)] <- as.vector(mapply(gsub, "[.]","",words[1:(length(words)-1)]))
	
	myrow$Gloss <- paste(words, collapse=" ")
	myrow$Gloss <- gsub(" +", " ", myrow$Gloss)
	
	#Return
	myrow
	
} #END get_utt_info

scrubGloss =  function(string){
	#remove +\"
	string = gsub('\\+','',string)
	string = gsub(' \\?','?',string)
	string = gsub('//','',string)
	return(string)	
}

convertToCSV = function(infolder, outfolder,corpus){
	corpusDir = paste(ensureTrailingSlash(infolder), ensureTrailingSlash(corpus), sep='')
	potentialFileList = list.files(corpusDir, recursive=T)
	fileList = potentialFileList[grep('.cha$',potentialFileList)]
	outputDir = paste(ensureTrailingSlash(outfolder), ensureTrailingSlash(corpus), sep='')
	dir.create(outputDir, showWarnings=F)
	#lapply(fileList, function(f){
	mclapply(fileList, function(f){
		df = read.CLAN.file(paste(ensureTrailingSlash(infolder), ensureTrailingSlash(corpus), f, sep=''))
		filename = gsub('.cha','.csv',paste(ensureTrailingSlash(outfolder), ensureTrailingSlash(corpus), gsub(ensureTrailingSlash(corpus), '', f), sep=''))			
		dir.create(paste(head(strsplit(filename,'/')[[1]],-1), collapse='/'), showWarnings=F)			
		df$child = gsub('[0123456789].*cha$','', tail(strsplit(df$Filename[1],'/')[[1]],1))
write.csv(df[,c('Speaker',	'Verbatim',	'Gloss','Filename',	'Language',	'Corpus','Age',	'Gender',	'Utt.Number','child')], filename, row.names=F)	
	}, mc.cores=detectCores())		
	#})
	
	#merge all of the files into a single term
	#system(paste('cat ', paste(outfolder, corpus,'/', sep=''), '*.csv > ', outfolder, corpus, '.csv', sep=''))
	system(paste("find ",outputDir ," *.csv -exec cat {} \\; >> ", outfolder, corpus, '_temp.csv', sep=''))
	system(paste('rm -r ',ensureTrailingSlash(outfolder), ensureTrailingSlash(corpus), sep=''))	
	
	#read back in and remove the final numeric token; also remove the header lines  

	fullDataset = read.csv(paste(outfolder, corpus, '_temp.csv', sep=''), stringsAsFactors=F)	
	system(paste('rm ',outfolder, corpus, '_temp.csv', sep=''))
	
	fullDataset$Gloss = unname(sapply(fullDataset$Gloss, function(x){
		gsub('^ ', '',gsub(' $','', scrubGloss(gsub(' \025.*$', '',x))))	
	}))
	
	fullDataset = fullDataset[-(which(is.na(as.numeric(fullDataset$Utt.Number)))),]	
	
	write.csv(fullDataset, paste(outfolder, corpus, '.csv', sep=''), row.names=F)
			
}

getAgeInDays = function(filename){
	elements = strsplit(filename,'/')[[1]]
	pieces = as.numeric(strsplit(gsub('.cha', '', elements[length(elements)]), '-')[[1]])
	if (length(pieces) != 3){stop('filenames should consist of 3 components')}
	return(ceiling((pieces[1] * 12 *30.5) + (pieces[2]* 30.5) + pieces[3]))
}

augmentAndCleanThomasDataset = function(td){
	td$child= 'Thomas' 
	td$Age = unname(sapply(td$Filename, getAgeInDays))

	#remove all of the ambiguous/ unknown @sc words, clean up the punctuation
	td$Gloss = sapply(strsplit(td$Gloss, ' '), function(x){
		paste(x[x!="a@p"], collapse=' ')
	})
	td$Gloss = gsub(' ,', ',', td$Gloss)
	td$Gloss = gsub(',,', ',', td$Gloss)
	td$Gloss = gsub('@o', '', td$Gloss)
	td$Gloss = gsub('@p', '', td$Gloss)
	td$Gloss = gsub('@f', '', td$Gloss)
	td$Gloss = gsub('@c', '', td$Gloss)
	td$Gloss = gsub(' \\.$', '', td$Gloss)
	td$Gloss = gsub(' \\/$', '', td$Gloss)
	td$Gloss = gsub("-'has", ' has', td$Gloss)
	td$Gloss = gsub("-'had", ' had', td$Gloss)
	td$Gloss = gsub("-'does", ' does', td$Gloss)
	td$Gloss = gsub("0", '', td$Gloss)
	td$Gloss= gsub('”','',td$Gloss)
	td$Gloss= gsub('“','',td$Gloss)
	td$Gloss = gsub('\x93', '',sapply(td$Gloss, function(x){iconv(x, "", "cp1252")} ))
	td$Gloss = gsub('\x94', '',sapply(td$Gloss, function(x){iconv(x, "", "cp1252")} ))
	td$Gloss = gsub('\025.*\025', '', td$Gloss)
	return(td)
}


rejoinCHILDESdata = function(inputFolder, taggedFolder, outputFolder, datasetName, identifier){
inputFile = paste(inputFolder, datasetName, '.csv', sep='')

taggedFile = paste(taggedFolder, datasetName, '_tagged.txt', sep='')

idf = read.csv(inputFile, stringsAsFactors=F)
names(idf)

tdf = read.delim(taggedFile, header=F, sep="\t", stringsAsFactors = F)
names(tdf) = c('fileNameUttIndex','gloss','token')
#need a method to deal with ^M
#if (nrow(tdf) != as.numeric(strsplit(system(paste("wc -l ", taggedFile), intern=T), ' ')[[1]][3])){
#	stop("length of indefinite csv read into R not the same as with wc")
#}

#parse out the essential properties of tdf for merging with the original input files 
parsedID = strsplit(gsub('[\\(\\)]','',tdf$fileNameUttIndex), ',')
tdf$Utt.number = sapply(parsedID, function(x){as.numeric(x[2])})
tdf$filename = gsub("'",'',sapply(parsedID, function(x){as.character(x[1])}))

parsedTokens = strsplit(tdf$token, ' ')
parsedNouns = strsplit(sapply(parsedTokens, function(x){x[length(x)]}), '_')

tdf$noun = sapply(parsedNouns, function(x){x[1]})
tdf$pos = sapply(parsedNouns, function(x){x[2]})
parsedModifiers = strsplit(sapply(parsedTokens, function(x){x[1]}), '_')
tdf$numIntermediateWords = sapply(parsedTokens, function(x){length(x) -2})
tdf$modifier = sapply(parsedModifiers, function(x){x[1]})
tdf$relation = sapply(parsedModifiers, function(x){x[2]})
tdf.nn = subset(tdf, pos %in% c('NN','NNP','NNS','NNPS'))

taggedTokens  = merge(tdf.nn, idf, by.x=c('filename','Utt.number'), by.y=c('Filename','Utt.Number'))
taggedTokens = taggedTokens[order(taggedTokens$filename, taggedTokens$Utt.number),] #this does need to be explicitly reordered after the merge
if(nrow(tdf.nn) != nrow(taggedTokens)){stop('error with merge')}

#thomas has a different age format, that must be retrieved from the filenames
if(datasetName == 'Thomas'){
	getAgeInDays = function(filename){
		shortenedFilename = strsplit(filename,'/')[[1]][length(strsplit(filename,'/')[[1]])]		
		pieces = as.numeric(strsplit(gsub('.cha', '', shortenedFilename), '-')[[1]])
		if (length(pieces) != 3){stop('filenames should consist of 3 components')}
		return(ceiling((pieces[1] * 12 *30.5) + (pieces[2]* 30.5) + pieces[3]))
	}
	taggedTokens$age = sapply(taggedTokens$filename, getAgeInDays)	
} else {
	taggedTokens$age = unname(sapply(taggedTokens$Age, ageToDays))
}

rdf = data.frame(modifier = taggedTokens$modifier, noun = taggedTokens$noun, pos = taggedTokens$pos, relation = taggedTokens$relation, speaker = taggedTokens$Speaker, 	age=taggedTokens$age, file= taggedTokens$filename, one =1, child=taggedTokens$child, det = taggedTokens$modifier, def= (taggedTokens$modifier == 'the'), numIntermediateWords = taggedTokens$numIntermediateWords, Utt.number = taggedTokens$Utt.number)

#same checks as for Speechome on word identities
#all words greater than length 2
rdf = rdf[sapply(as.character(rdf$noun), nchar) > 2,]

#no caps
rdf = rdf[as.character(rdf$noun) == tolower(as.character(rdf$noun)),]

outputFile = paste(outputFolder, datasetName,'_',identifier, '.csv', sep='')

write.csv(rdf,outputFile, row.names=F)
}

getAggregate = function(df){
  df$one = 1
  names(df)[grep('noun', names(df))] = 'word'
  nounUses = aggregate(one ~ word + speaker, df, sum)
  parentalNouns = aggregate( one ~ word, subset(nounUses, speaker %in% c('MOT','FAT')), sum)
  pCount = merge(subset(nounUses, speaker == 'CHI'), parentalNouns, by='word', all.x=T, all.y=T)[,c('word','one.y')]
  names(pCount) =c('word','count')
  pCount[is.na(pCount)] = 0
  pCount = pCount[order(pCount$word),]  
  return(pCount)
}

getAggCounts = function(inputFolder, inputCorpora, morphology, datasetName, outputFile){
	possibleFiles = list.files(inputFolder)
	files = possibleFiles[sapply(strsplit(possibleFiles,'_'), function(x){any(x %in% inputCorpora)})]
	files = files[grep(morphology,files)]
	files= files[grep(datasetName, files)]
	
	allDatasets = do.call('rbind.fill',lapply(files, function(file){
	  read.csv(paste(inputFolder, file, sep=''), stringsAsFactors=F)	
	}))
  
	allDatasets.agg = getAggregate(allDatasets)
	
	write.csv(allDatasets.agg, outputFile, row.names=F)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}