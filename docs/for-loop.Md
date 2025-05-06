```R
for (i in array) {
  # Code
}
```


```R
# Example: Calculate the square of numbers in a vector
numbers <- c(1, 2, 3, 4, 5)
squares <- c()

for (num in numbers) {
    squares <- c(squares, num^2)
}

print(squares)
# Output: [1]  1  4  9 16 25
```


## Example from practical UU-GWAS03 - task 4

#### Genotype frequencies
Calculate genotype (0: A1A1, 1: A1A2, 2: A2A2) frequencies of each loci 
- Generate an empty data frame (allele_frequency) with three columns (A1A1, A1A2, A2A2).   
- Number of rows = number of markers.  

```R
freq = data.frame(A1A1 = rep(NA, nrow(data$bim)), A1A2=NA, A2A2=NA)

```

#### Use a for loop to count the number of genotypes in each marker and divide by the number of samples.

```R
nSmp = nrow(data$fam)
for(i in 1:nrow(data$bim)){
     tmp = data$X[i, ]

     freq$A1A1[i] = sum(tmp == 0, na.rm = T) / nSmp
     freq$A1A2[i] = sum(tmp == 1, na.rm = T) / nSmp
     freq$A2A2[i] = sum(tmp == 2, na.rm = T) / nSmp
 }
```
