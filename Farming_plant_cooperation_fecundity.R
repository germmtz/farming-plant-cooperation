## The "fec_extended" function is used to compute a vector of plant fecundity

fec_extended <- function (individuals, heights, xopt, d, a) { # "individuals" = vector with all individuals ids refering to their position on the circle; "heights" :  vector of plant heights, "xopt" : optimal plant height for individual seed production, "d" : between-plant distance, "i" : light interception ability = angle at the summit of the triangles

n_ind <- length(individuals) # retrieving the number of individuals

## To be  able to compute plant shades in a linear coordinate system (whereas plants are supposed to grow on a circle), we artifically triplicate each plants.For exemple, with 10 plants, we create the following system p1p2p3p4p5p6p7p8p9p10|p1p2p3p4p5p6p7p8p9p10|p1p2p3p4p5p6p7p8p9p10 so that we can compute the shadow of 2 to 1 in the same coordinates system than the shadow of 10 to 1. The coordinates of the first "bloc" of plants will be stored in coordinates_left, the second "bloc" in coordinates_center and the third bloc in coordinates_right. 

coordinates_left <- data.frame(id=NULL,coord=list())
for (j in individuals) {
coordinates_left[j,"id"] <- j
coordinates_left$coord[j]<- list(rbind(c((j-1)*d+tan(a)*heights[j],0),c((j-1)*d,heights[j]),c((j-1)*d-tan(a)*heights[j], 0),c((j-1)*d+tan(a)*heights[j],0)))
} # each triangle is defined by three edges, each edge being described with a (x,y) coordinate

coordinates_center <- data.frame(id=NULL, coord=list())
for (j in individuals) {
  coordinates_center[j,"id"] <- j+n_ind
  coordinates_center$coord[j]<- lapply(coordinates_left$coord[j], function(x) {cbind (x[,1]+(n_ind*d), x[,2])}) # we just translate the previously defined triangles (coordinates_left) of + n_ind*d in the x dimension
}

coordinates_right <- data.frame(id=NULL, coord=list())
for (j in individuals) {
  coordinates_right[j,"id"] <- j+2*n_ind
  coordinates_right$coord[j]<- lapply(coordinates_center$coord[j], function(x) {cbind (x[,1]+(n_ind*d), x[,2])})# we just translate the previously defined triangles (coordinates_center) of + n_ind*d in the x dimension
}


coordinates <- rbind(coordinates_left, coordinates_center, coordinates_right) # merging the three "blocs" of plants together


##########################################

# In order to avoid the computation of all pairwise interactions, we use a simple test to check whether  two individuals are interacting or not based on the relative value of their edges x coordinate
for (j in 1:nrow(coordinates)) {
  coordinates[j,"x_left"] <- unlist(coordinates$coord[[j]][3,1])
  coordinates[j,"x_right"] <- unlist(coordinates$coord[[j]][1,1])
  coordinates[j,"x_feet"] <- unlist(coordinates$coord[[j]][2,1])
}


interactions <- data.frame(id=NULL, fecundity=NULL)

for (j in unique(coordinates_center$id)){

  interactions[j-n_ind,"id"] <- j-n_ind
  ref_left <- coordinates[which(coordinates$id==j),"x_left"]
  ref_right <- coordinates[which(coordinates$id==j),"x_right"]
  ref_feet <- coordinates[which(coordinates$id==j),"x_feet"]
  interactant <- coordinates[which((coordinates$x_right>ref_left & coordinates$x_left<ref_feet | coordinates$x_left<ref_right & coordinates$x_right>ref_feet) & coordinates$id!=j),"id"] # simple test to retrieve, for each focal individuals, the reduced list of neighbours that really shade him
  focal <- st_sfc(st_polygon(coordinates$coord[j])) # converting the focal individual into a geometrical object
  nei <- st_union(st_sfc(lapply(coordinates[which(coordinates$id%in%interactant),"coord"], function (x) {st_polygon(list(x))}))) # creating a geometrical object that contains the union of all neigbhors that interact with the focal
  shade <- st_area(st_intersection(focal,st_buffer(nei,0))) # computing the shade on the focal as the intersection between the two previously defined objects
  
  if (identical(shade, numeric(0))) { ## when there is no shading, we replace the "numeric(0)" value returned by st_area by 0. 
    shade <- 0
  }
  
  ## computing focal fecundity based on the shade
  x <- heights[j-n_ind]
  if (x > 2*xopt || shade>=tan(a)*x^2) {
  interactions[j-n_ind,"fecundity"] <- 0
  } else {
  interactions[j-n_ind,"fecundity"] <- x * (1 - x / (2 * xopt))*((tan(a)*x^2-shade)/(tan(a)*x^2))
  }
}

res <- interactions$fecundity
return(res)
}

