# internals.R


.extend_array <- function(initial_array, new_dim) 
{
   new_array <- array(NaN, dim = new_dim)
   
   .d <- dim(initial_array)
   
   new_array[1:.d[1], 1:.d[2], 1:.d[3]] <- initial_array 
   dimnames(new_array) <- dimnames(initial_array)
   
   return(new_array)
}