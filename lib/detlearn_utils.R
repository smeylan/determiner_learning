#shared utility functions

renameColumn = function(df, oldname, newname){
	if (!is.data.frame(df)){
		stop('df must be a data frame')
	}
	if (!is.character(oldname) | !is.character(newname)){
		stop('oldname and newname must be character vectors')
	}
	if (!(oldname %in% names(df)) ){
		stop('oldname must be a column name in the data frame, df')
	}
	names(df)[grep(paste('^',oldname,'$', sep='') ,names(df))] = newname	
	return(df)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

ensureTrailingSlash= function(string){
	#make sure that there is a trailing slash in a filename
	if (substr(string, nchar(string),nchar(string)) != '/'){
		return(paste(string,'/', sep=''))
	} else {
		return(string)
	}	
}