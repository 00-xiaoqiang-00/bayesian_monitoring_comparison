colon_patient_pairs <- function(path) {
  records <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  recurrence <- records[records$etype==1L,c("id","rx","time","status")]
  names(recurrence)[3:4] <- c("recurrence_time_days","recurrence_status")
  death <- records[records$etype==2L,c("id","time","status")]
  names(death)[2:3] <- c("death_time_days","death_status")
  merge(recurrence,death,by="id",sort=TRUE)
}
