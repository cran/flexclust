library("flexclust")

randomTour(iris[,1:4], axiscol=2:5, directions=3, steps=5, sleep=0)
randomTour(iris[,1:4], col=as.numeric(iris$Species), axiscol=4,
           directions=3, steps=5, sleep=0)

x <- matrix(runif(300), ncol=3)
x <- rbind(x)
cl <- cclust(x, k=3, save.data=TRUE)

randomTour(cl, center=0, axiscol="black", directions=3, steps=5,
           sleep=0)
randomTour(cl, center=0, axiscol="black",
           data=matrix(rnorm(300, mean=1, sd=2), ncol=3),
           directions=3, steps=5, sleep=0)
