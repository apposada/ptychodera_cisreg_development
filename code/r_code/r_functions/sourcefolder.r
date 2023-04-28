#from @Rich Scriven https://stackoverflow.com/questions/26045478/source-r-scripts-in-different-folders
sourceFolder <- function(folder, recursive = FALSE, ...) 
{ 
    files <- list.files(folder, pattern = "[.][rR]$", 
                        full.names = TRUE, recursive = recursive)
    if (!length(files))
        stop(simpleError(sprintf('No R files in folder "%s"', folder)))
    src <- invisible(lapply(files, source, ...))
    message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}
