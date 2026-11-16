#' data.simulation
#' 
#' This function simulates data. 
#' 
#' @param N.sim Number of simulations.
#' @param N.max   Maximum number of patients to enroll.
#' @param mu.trt A vector containing the mean time to each event (logarithm) for the treatment arm.
#' @param Sigma.trt A vector representing the variance-covariance matrix of the time to each event (logarithm) for the treatment arm. 
#' @param mu.ctrl A vector containing the mean time to each event (logarithm) for the control arm.
#' @param Sigma.ctrl A vector representing the variance-covariance matrix of the time to each event (logarithm) for the control arm. 
#' @param cens_upper Upper limit for the censoring time, assuming that the censoring time is generated from Uniform(0, cens_upper).
#' @returns An array with N.max rows, 7 columns and N.sim slices. The columns are:\tabular{ll}{
#'    \code{recurrence_t} \tab Time to recurrence.  \cr
#'    \tab \cr
#'    \code{death_t} \tab Time to death. \cr
#'    \tab \cr
#'    \code{censor_t} \tab Time to censor. \cr
#'    \tab \cr
#'    \code{group} \tab Treatment assignment.  \cr
#'    \tab \cr
#'    \code{id} \tab Patient id, numeric starting from 1. \cr 
#'    \tab \cr
#'    \code{delta1} \tab Placeholder for the censoring indicator for recurrence, leave blank. \cr
#'    \tab \cr
#'    \code{delta2} \tab Placeholder for the censoring indicator for death, leave blank. \cr
#' }
#' @examples
#' \dontrun{
#' data.simulation(N.sim = 1000, N.max = 20,
#'                mu.trt = c(0.2,0.3), Sigma.trt = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                mu.ctrl = c(0.2,0.3),Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                cens_upper = 5)
#' }

#' @export
data.simulation<- function(N.sim,N.max,mu.trt,Sigma.trt,
                           mu.ctrl,Sigma.ctrl,cens_upper
){
  myData<-array(dim=c(N.max,7,N.sim));#三维数组[N.max 行, 7 列, N.sim 层]
  for (j in 1:N.sim){
    
    # Simulate randomization groups
    rand.list = sample(c(0,1), N.max, replace = T)
    
    # Simulatetime to recurrence and death
    trt.events <- 	as.data.frame(MASS::mvrnorm(sum(rand.list), mu= mu.trt,Sigma= Sigma.trt))    
    #MASS包中的mvrnorm函数用于生成多元正态分布的随机数。这里生成了sum(rand.list)行也就是治疗组的生存数据，每行包含两个变量，分别对应mu.trt中的两个元素。
    trt.events$id =  which(rand.list==1); 
    #赋值给trt.events数据框中的id列，使用which(rand.list==1)来获取治疗组患者的索引，并将这些索引赋值给id列。这样，每个治疗组患者都有一个唯一的id。
    trt.events$grp = 1; 
    #为trt.events数据框中的grp列赋值为1，表示这些患者属于治疗组。

    #对照组同理
    ctrl.events <- as.data.frame(MASS::mvrnorm(sum(rand.list==0),
     mu= mu.ctrl,
     Sigma= Sigma.ctrl))     

    ctrl.events$id = which(rand.list==0); 
    ctrl.events$grp = 0

    # simulate time to censor
    time.censor = stats::runif(N.max, 0, cens_upper); 
    #使用stats包中runif函数生成均匀分布的N.max个随机数，范围在0到cens_upper之间，表示每个患者的 censoring time。
    eventsall = rbind(trt.events, ctrl.events); 
    #按行合并治疗组和对照组的数据框，得到一个包含所有患者的事件数据框eventsall。
    eventsall =  eventsall[order(eventsall$id),]
    # 按患者 id 排序。
    myData[,1,j]<- exp(eventsall[,1])  
                   # time to recurrence 
    myData[,2,j]<- exp(eventsall[,2])                      # time to death
    myData[,3,j]<- time.censor  #第三列为删失时间
    myData[,4,j]<- eventsall[,4]                      # treatment group assnments组别
    myData[,5,j]<-  eventsall[,3]                     # id
#myData[,6,j] <- 0   # delta1
#myData[,7,j] <- 0   # delta2

  

  }
  
  # Create dimension names
  dim_names <- list(
    NULL,#维度1的行名可以留空，因为每行代表一个患者，患者ID已经在数据中体现了。
    Variables = c("recurrence_t", "death_t", "censor_t",   "group", "id", "delta1", "delta2"),
    Simulations = paste("Simulation", 1:N.sim)#维度3的层名，可以命名为 "Simulation 1", "Simulation 2", ..., "Simulation N.sim"，以便区分不同的模拟结果。
  )
  
  # Assign the dimension names to the array
  dimnames(myData) <- dim_names
  print("data成功生成")
  return(myData)
}
