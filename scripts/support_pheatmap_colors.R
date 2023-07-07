createLinearColors <- function(numbers, withZero = T, maximum = 100, my_colors = c("royalblue3","white","red"))
{
  if (min(numbers, na.rm = T) > 0)
  {
    if(withZero)
    {
      numbers <- c(0,numbers)
    }
    myPalette <- colorRampPalette(my_colors[c(2,3)])
    myColors <- myPalette(maximum)
  }
  else
  {
    if (max(numbers, na.rm = T) < 0)
    {
      if(withZero)
      {
        numbers <- c(0,numbers)
      }
      myPalette <- colorRampPalette(my_colors[c(1,2)])
      myColors <- myPalette(maximum)
    }
    else
    {
      myPalette_pos <- colorRampPalette(c("white","red"))
      myPalette_neg <- colorRampPalette(c("royalblue3","white"))
      npos <- length(numbers[numbers >= 0]) + 1
      nneg <- length(numbers[numbers <= 0]) + 1
      
      myColors_pos <- myPalette_pos(npos)
      myColors_neg <- myPalette_neg(nneg)
      
      #print(myColors_neg)
      #print(myColors_pos)
      
      myColors <- c(myColors_neg[-(nneg)], myColors_pos[-1])
    }
  }
  return(myColors)
}
