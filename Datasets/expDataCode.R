#The following aircraft arrival times are for the control sector 454LT 
#studied in this article. The sample period was 16/00/00-24/00/00 
#GMT (noon-8 P.M. New York time) on April 30, 1969. 
#Tabulated values are in seconds from the start of the sample period. 
#The data were originally collected by the Federal Aviation Administration,
#National Aviation Facilities Experimental Center, Atlantic City, New Jersey.
exp.data <- c(467,761,792,812,926,1100,1147,1163,1398,1462,1487,1749,1865,
                   2004,2177,2208,2279,2609,2682,2733,2818,2837,2855,2868,3089,
                   3209,3223,3233,3272,3399,3595,3634,3650,3851,4176,4304,4391,
                   4453,4539,4748,4839,5049,5202,5355,5551,5598,5640,5702,5935,
                   6000,6192,6435,6474,6600,6810,6824,7168,7181,7202,7218,7408,
                   7428,7720,7755,7835,7958,8307,8427,8754,8819,8904,8938,8980,
                   9048,9237,9268,9513,9635,9750,9910,9929,10167,10254,10340,
                   10624,10639,10669,10889,11386,11515,11651,11727,11737,11844,
                   11928,12168,12657,12675,12696,12732,13092,13281,13536,13556,
                   13681,13710,14008,14151,14601,14877,14927,15032,15134,15213,
                   15491,15589,15600,15631,15674,15797,15953,16089,16118,16215,
                   16394,16503,16515,16537,16570,16597,16619,16693,17314,17516,
                   17646,17770,17897,17913,17922,18174,18189,18328,18345,18499,
                   18521,18588,19117,19150,19432,19662,19758,19789,19831,19978,
                   20119,20312,20346,20449,20455,20604,20675,20817,20898,21245,
                   21386,21562,22022,22056,22095,22182,22554,22764,22955,22993,
                   23025,23117,23321,23341,23650,23766,23879,23888,24458,24889,
                   24930,24967,25224,25312,25477,25498,25712,25721,25884,25919,
                   25985,26196,26459,26468,26494,26505,26554,26906,27003,27437,
                   27661,27675,27697,27721,27734,27802, 27971, 28116, 28746)
#write.table(arrival.times, file = "my_data.txt", sep = "\t")
#my_data.1 <- read.table("my_data.txt")
#exp.data <- my_data.1$x

####AIR TRAFFIC DATA - EXPONENTIAL DISTRIBUTION
#### Cumulative frequency plot for arrival times (in seconds)
breaks=seq(450, 30000, by=500)
times.cut = cut(exp.data, breaks, right=FALSE)
times.freq = table(times.cut)
times.cumfreq = cumsum(times.freq)
cumfreq0 = c(0, cumsum(times.freq)) 
plot(breaks, cumfreq0,main="Plot of Cumulative Arrival Counts 
     vs Times Arrival time", xlab="Arrival times(seconds)",
     ylab="Cumulative fRequencys")   # y−axis label 
lines(breaks, cumfreq0)  


inter.times <- diff(exp.data, lag=1)

# Converting the data to a matrix format that will work nicely with our code
realDataMatrix <- matrix(inter.times, nrow=1, ncol=length(inter.times))

