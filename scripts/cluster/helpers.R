unlapply<- function(...)
{
    unlist(lapply(...))
}

max2minOrder<- function(x)
{
stopifnot(is.numeric(x))
order(-x)
}

min2maxOrder<- function(x)
{
stopifnot(is.numeric(x))
order(x)
}
max2minSort<- function(x)
{
stopifnot(is.numeric(x))
x[order(-x)]
}

min2maxSort<- function(x)
{
stopifnot(is.numeric(x))
x[order(x)]
}

