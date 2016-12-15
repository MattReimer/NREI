#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------
# Written for Eco Logical Research, Inc. and South Fork Research, Inc.
# December 2014
#
# For information, questions, or suggestions, please contact:
# Eric Wall
# c.eric.wall@gmail.com
#`
# Description:  A set of functions to implement a net rate of energy intake
#  model for drift-feeding salmonids. Functions in this script take 
#  rectilinear hydraulic output, create an on-the-fly 'wall of water' for each
#  raster cell normal to the flow direction of that cell and calculate
#  NREI.  The foraging model is Addley's (1993) derivation (adapted from
#  Hughes and Dill (1990)) which accounts for velocity shear near the foraging
#  point.  The method for calculating GREI is that of Hughes et al. (2003).
#  Swimming costs calculations are those of Hayes et al. (2007).
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Section I: Inputs
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Section II: Hydraulic Functions
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------
# Function to read and interpret rectilinear output from Delft3D and build
#  a 3D grid of the site's wetted points; parallel adaptation of MN's original
#  Build3DGrid function
#
# Args:`
#  DEMfilename: name of the rectilinear Delft3D results file
#  DZ: spacing for 3D grid in Z direction (may need to be tighter than X 
#    and Y direction)
#  grd.reduction.factor: factor by which to increase X and Y grid spacing
#    relative to the 0.1m standard grid spacing; My be necessary for big
#    sites (as in the Entiat)
#  k.roughness: roughness height in m
#
# Returns:
#  a list object consisting of:
#    1. Grid3D, the 3D grid of all wetted points with their logarithmic
#      velocity values calculated
#    2. the Data.Use df, which is the original wetted df from Matt's output
#    3. a float number giving the spacing in X and Y of the resulting Grid3D
#      point cloud
#   
#---------------------------------------------------------------------------

Build3DGrid_par <- function(DEMfilename, DZ, grd.reduction.factor, k.roughness,
  sim.cores) {

  print("Reading Hydraulic Model Data File (This may take a minute!)")

  # Read the data.  Screen to only wetted areas.
  require(data.table)
  DEM.Results = fread(DEMfilename, header = TRUE, data.table = FALSE)
  
  # calculate the site area
  area.calc.temp <- DEM.Results[DEM.Results$Depth > 0,]
  # head(area.calc.temp)
  # dim(area.calc.temp)
  # summary(area.calc.temp)
  area.calc.temp <- area.calc.temp[area.calc.temp$BedLevel != -9999,]
  # plot(area.calc.temp$X, area.calc.temp$Y)
  # summary(area.calc.temp)
  my.area <- nrow(area.calc.temp) * 0.1 * 0.1
  my.volume <- sum(area.calc.temp$Depth * 0.1 * 0.1)
  rm(area.calc.temp)
  
  Xlev = levels(as.factor(DEM.Results$X))
  xidx = seq(1, length(Xlev), by = grd.reduction.factor)

  Ylev = levels(as.factor(DEM.Results$Y))
  yidx = seq(1, length(Ylev),by = grd.reduction.factor)

  DEM.Results = DEM.Results[(DEM.Results$X %in% Xlev[xidx]),]
  DEM.Results = DEM.Results[DEM.Results$Y %in% Ylev[yidx],]

  # dim(DEM.Results)
  Data.Use = DEM.Results[DEM.Results$Depth > 0,]
  Data.Use = Data.Use[Data.Use$BedLevel != -9999,]
  Data.Use$Z = rep(0, nrow(Data.Use))
  rm(DEM.Results)

  # # Check - Plot the velocity and depth
  # par(mfrow = c(2,1))

  # plot(Data.Use$X, Data.Use$Y, 
    # col = colorp[round(50*Data.Use$Depth/max(Data.Use$Depth))+1],main="Depth",
    # asp = 1
  # )

  # plot(Data.Use$X, Data.Use$Y, 
    # col = colorp[round(50*Data.Use$Velocity.Magnitude/max(Data.Use$Velocity.Magnitude))+1],main="Velocity",
    # asp = 1
  # )

  # Convert 2D Grid into 3D Grid.  We'll let z-grid spacing by tighter than X and Y grid spacing of .1 m

  # Set the maximum number of levels in Z. We'll build a rectangular prism encompassing all the wetted volume,
  # then trim out dry cells afterwards
  Zlevels = seq(trunc(min(Data.Use$BedLevel),1),round(max(c(Data.Use$WSE,Data.Use$BedLevel))),, by = DZ)
  #Zlevels = as.numeric(levels(factor(round(Data.Use$BedLevel,digits=1))))
  NZ = length(Zlevels)

  # register parallel backend
  require(foreach)
  require(doParallel)
  cl <- makeCluster(sim.cores)
  registerDoParallel(cl)

  # build 3D point cloud
  Grid3D <- foreach(nz = 1:NZ, .combine = 'rbind') %dopar% {
    
    # add appropriate Z values to Data.Use
    temp = Data.Use
    temp$Z = as.numeric(Zlevels[nz])
    
    # remove points outside of usable area
    temp = temp[temp$Z < temp$WSE,]
    temp = temp[temp$Z > temp$BedLevel,]

  }
  stopCluster(cl)
  
  # just to be sure
  Grid3D = Grid3D[Grid3D$Z < Grid3D$WSE,]
  Grid3D = Grid3D[Grid3D$Z > Grid3D$BedLevel,]

  order.z = order(Grid3D$Z)
  Grid3D = Grid3D[order.z,]
  order.y = order(Grid3D$Y)
  Grid3D = Grid3D[order.y,]
  order.x = order(Grid3D$X)
  Grid3D = Grid3D[order.x,]

  Grid3D$DA_X.Velocity = Grid3D$X.Velocity
  Grid3D$DA_Y.Velocity = Grid3D$Y.Velocity

  #######################################################################
  # logarithmic vertical velocity profile

  Grid3D$DA.Velocity.Magnitude = sqrt(
    (Grid3D$DA_X.Velocity * Grid3D$DA_X.Velocity) + (Grid3D$DA_Y.Velocity * Grid3D$DA_Y.Velocity)
  )

  #k.roughness= 10 # roughness

  vstar = Grid3D$DA.Velocity.Magnitude / (5.75* log(12.3 * (Grid3D$WSE-Grid3D$BedLevel) / k.roughness))
  Grid3D$Velocity.Magnitude =  5.75* log(30*(Grid3D$Z-Grid3D$BedLevel)/k.roughness)*vstar

  # Correct for conservation of flow
  XY = (factor(paste(Grid3D$X, Grid3D$Y)))
  Correction = tapply(Grid3D$Velocity.Magnitude, XY, sum)/(tapply(Grid3D$DA.Velocity.Magnitude,XY,sum)+.000001)


  Grid3D$Velocity.Magnitude = Grid3D$Velocity.Magnitude/(Correction[match(XY, levels(XY))]+.000000001)

  angle = atan2(Grid3D$Y.Velocity, Grid3D$X.Velocity)

  Grid3D$X.Velocity = Grid3D$Velocity.Magnitude *cos(angle)
  Grid3D$Y.Velocity = Grid3D$Velocity.Magnitude *sin(angle)


  Grid3D$Dir = atan2(Grid3D$Y.Velocity,  Grid3D$X.Velocity)+pi/2
  
  # add a colum, 'zrank', to indicate the relative position from the bed of 
  #  3d points for the various x-y combos
  temp.dt <- data.table(Grid3D[, c('X', 'Y', 'Z')])
  # names(temp.dt)
  temp.dt[, zrank := rank(Z), by = list(X, Y)]
  Grid3D$zrank <- temp.dt$zrank
  
  # min(Grid3D$Y.Velocity)
  # min(Grid3D$X.Velocity)
  # min(Grid3D$Velocity.Magnitude)
  return(
    list(
      "Grid3D"=Grid3D,
      "Data.Use"=Data.Use,
      "DX_DY" = .1*grd.reduction.factor,
      "my.area" = my.area,
      "my.volume" = my.volume
    )
  )

}


#---------------------------------------------------------------------------
# Function to initialize a radial grid we'll use over and over for every
#  XYZ point in our 3D grid
#  
# Args:
#  pts: number of angles in the radial grid (sorry about the bad name!)
#  max.Dist: maximum radius on any radial in the radial grid.
#    Note: this should be at least as far as fish reaction distance
#  Dist: distance between points along a radial in the radial grid.  Again,
#    lousy name!
#
# Returns:
#  a list object consisting of:
#    1. rad.grid.r: an array containing the radii lengths for each point in the 
#      radial grid; one row for each radial, one column for each point
#    2. rad.grid.theta: an array containing the thetas for each point in the 
#      radial grid; one row for each radial, one column for each point
#    3. Rmax: total number of points on each of the radials
#    4. VMult: a fixed matrix to help with mean normal velocity calculations
#   
#---------------------------------------------------------------------------

Initialize.Radial.Grid <- function(pts, max.Dist, Dist) {

  # Initialize some stuff

  angles = seq(-1*pi, pi, by= (2*pi)/pts)[1:pts]
  Dists = seq(Dist,max.Dist, by=Dist)

  rad.grid.r = (t(array(rep(c(Dists), pts), c(max.Dist/Dist,pts))))
  rad.grid.theta = array(angles[rep(rep((1:pts),max.Dist/Dist))], c(pts,max.Dist/Dist))

  # Initializing this once - outside the loops - for later... used to calculate mean
  # Velocities from NREI location every point in radial grid
  Rmax = max.Dist/Dist
  VMult = array(0, c(Rmax,Rmax+1))
  
  for (k in 1:Rmax) {
  
    VMult[k,] = c(.5, rep(1, (k-1)), .5,rep(0,(Rmax-k)))
  
  }
    VMult = t(VMult)

  
  #-------------------------------------------------------------
  # initialize an array to hold the 'slice' areas
  #-------------------------------------------------------------
  rad.grid.areas <- rad.grid.r * 0

  # first slices
  rad.grid.areas[, 1] <- pi * (1.5 * Dist) * (1.5 * Dist) / pts

  # second to second-to-last slices
  for (my.rad in 2:(Rmax - 1)) {
    
    rad.grid.areas[, my.rad] <- 2 * my.rad * pi * Dist * Dist / pts
  
  }

  # last slices
  rad.grid.areas[, Rmax] <- pi * Dist * (max.Dist - 0.25 * Dist) / pts

  # check on areas; threshold is 1 square mm
  if (sum(rad.grid.areas) - (pi * (max.Dist * max.Dist)) >= 0.000001) {
  
    print('!!!!!!!!  AREAS DID NOT ADD UP CORRECTLY  !!!!!!!!')
  
  }

  return(
    list(
     "rad.grid.r" = rad.grid.r,
     "rad.grid.theta" = rad.grid.theta,
     "Rmax" = Rmax,
     "VMult" = VMult,
     "rad.grid.areas" = rad.grid.areas
    )
  )
  
}


#---------------------------------------------------------------------------
# Helper function for mean normal velocity calculations; used in the
#  radial.grid function below.
#  
# Args:
#  rootVel: velocity at XYZ point being analyzed
#  Rmax: total number of points on each of the radials
#  VelMult: a matrix returned by the Initialize.Radial.Grid function above;
#    using this helps avoid for loops
#  VelVec: vector of velocities normal to flow along radius.
#
# Returns:
#   
#---------------------------------------------------------------------------

Mean.Zero.To.r <- function(rootVel, Rmax, VelVec, VelMult) {

  meanVel =array(c(rootVel,VelVec),c(1,Rmax+1)) %*% VelMult / (1:Rmax)

  return(meanVel)
  
}


#---------------------------------------------------------------------------
# Function to build a radial grid at the "wall of water" perpendicular to the
#  flow at an X-Y-Z grid point.
#  
# Args:
#  idx: index of XYZ grid; used to find location and velocity from the master
#    XYZ array
#  pts: number of angles in the radial grid
#  max.Dist: maximum radius on any radial in the radial grid
#    Note: this should be at least as far as fish reaction distance
#  Dist: distance between points along a radial in the radial grid
#  plots: boolean to turn QA plots on or off
#  rad.grid.r: an array containing the radii lengths for each point in the 
#    radial grid; one row for each radial, one column for each point;
#    of dimensions [pts, Rmax]
#  rad.grid.theta: an array containing the thetas for each point in the 
#    radial grid; one row for each radial, one column for each point;
#    of dimensions [pts, Rmax]
#  Rmax: number of possible points along radials (max.Dist/Dist)
#  VMult: See above function
#  DX_DY: grid spacing in X and Y in the 3D grid.
#
# Returns:
#  a list object containing many things.
#
#------------------
# Matt's notes:
#
# Function radial.grid
# Build a radial grid at the "wall of water" perpendicular to the flow at an X-Y-Z grid point.
# Return the grid, the normal velocities at each point in the grid, and the average
# velocity along each radius from the origin to each point at the radius.
#
# idx = index of X-Y-Z grid used to find location and velocity from master X-Y-Z array
# pts = number of angles
# max.Dist = maximum distance to look along radius (should be >= fish reaction distance)
# Dist = distance between points along radials
# plots: boolean to turn QA plots on or off.  if TRUE, plots will be generated at every 500th idx
# rad.grid.r: radial distances in radial grid of dimensions [pts, Rmax]
# rad.grid.theta: radial angles in radial grids of dimsions [pts, Rmax]
# Rmax: number of possible points along radials (max.Dist/Dist)
# VMult: See above function
# DX_DY: gird spacing in X and Y in the 3D grid.
#
#---------------------------------------------------------------------------

radial.grid <- function(idx, pts, max.Dist, Dist, plots,
  rad.grid.r, rad.grid.theta, Rmax, VMult,DX_DY) {

  #Data = Data.Use[
  #     Data.Use$X > (Grid3D$X[idx] - max.Dist) & 
  #     Data.Use$X < (Grid3D$X[idx] + max.Dist) & 
  #     Data.Use$Y > (Grid3D$Y[idx] - max.Dist) & 
  #     Data.Use$Y < (Grid3D$Y[idx] + max.Dist),]

  Data = Grid3D[
    Grid3D$X > (Grid3D$X[idx] - max.Dist) & 
    Grid3D$X < (Grid3D$X[idx] + max.Dist) & 
    Grid3D$Y > (Grid3D$Y[idx] - max.Dist) & 
    Grid3D$Y < (Grid3D$Y[idx] + max.Dist),]

  #rad.grid.r=as.vector(t(array(rep(c(Dists), pts), c(max.Dist/Dist,pts))))
  #rad.grid.theta=angles[rep(rep((1:pts),max.Dist/Dist))]

  # nrow(Grid3D)
  # Grid3D$X[idx]
  # Grid3D$Dir[idx]
  rad.grid.X =  Grid3D$X[idx] + cos(Grid3D$Dir[idx])*rad.grid.r*cos(rad.grid.theta)
  rad.grid.Y =  Grid3D$Y[idx] + sin(Grid3D$Dir[idx])*rad.grid.r*cos(rad.grid.theta)
  rad.grid.Z =  Grid3D$Z[idx] + rad.grid.r*sin(rad.grid.theta)

  #rad.grid.nn = nn2(Data[,1:2], data.frame(as.vector(rad.grid.X), as.vector(rad.grid.Y)),1)

  nn.col.idx = match( c("X","Y","Z"),colnames(Data))
  rad.grid.nn = nn2(Data[,nn.col.idx], 
    data.frame(as.vector(rad.grid.X), as.vector(rad.grid.Y),as.vector(rad.grid.Z)),1
  )

  # rad.grid.Z
  # rad.grid.X
  # rad.grid.Y
  rad.grid.nn.dists = array(rad.grid.nn$nn.dists, c(pts, max.Dist/Dist))

  rad.grid.Xvel = array(Data$X.Velocity[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.Yvel = array(Data$Y.Velocity[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.BedLevel = array(Data$BedLevel[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))
  rad.grid.WSE = array(Data$WSE[rad.grid.nn$nn.idx], c(pts, max.Dist/Dist))

  rad.grid.use = 
    (rad.grid.nn.dists < (DX_DY/1.2)) &
    (rad.grid.BedLevel < rad.grid.Z) &
    (rad.grid.WSE > rad.grid.Z)

  for (th in 2: (max.Dist/Dist)) {
  
    rad.grid.use[,th]= apply(1*rad.grid.use[,1:th], 1, prod)==1
  
  }

  rad.grid.Xvel[rad.grid.use==F] = 0
  rad.grid.Yvel[rad.grid.use==F] = 0
  rad.grid.Vmag = sqrt(rad.grid.Xvel^2 + rad.grid.Yvel^2)

  ##############################################################
  # Find the normal component of velocity for all grid points
  Vel.theta = atan2(rad.grid.Yvel, rad.grid.Xvel)
  Vel.Norm = rad.grid.Vmag * sin(-1*Vel.theta+Grid3D$Dir[idx])

  #################################################
  # For every radial grid point, find the mean velocity for all points along a radius
  # from the center to the given point.  This is a key NREI requirement.

  #initialize..
  rad.grid.mean.normV = rad.grid.X

  #Use the function defined outside the loop
  for (p in 1:pts) {
  
    rad.grid.mean.normV[p,]=Mean.Zero.To.r(Grid3D$Vmag[idx], Rmax, Vel.Norm[p,],VMult)
    
  }

  rad.grid.mean.normV = rad.grid.mean.normV * rad.grid.use

  ##################################################
  z.lim=c(min(Grid3D$Z), max(Grid3D$Z))

  if (plots) {
  
    if ((idx/plot.every) == round(idx/plot.every)) {

      par(mfrow = c(3,1))
      layout(matrix(c(1,2,3),3,1), heights=c(2,2,3.8))  # 
      par(mar=c(5,4,2,1))

      #par(mar=c(0,2,5,2))

      # Fun plots to watch as program runs.. but slows code down a LOT

      x.lim=c(min(rad.grid.X[rad.grid.use]),
        max(rad.grid.X[rad.grid.use])+
        .25*(max(rad.grid.X[rad.grid.use])-
        min(rad.grid.X[rad.grid.use]))
      )

      y.lim=c(min(rad.grid.Y[rad.grid.use]),
        max(rad.grid.Y[rad.grid.use])+
        .25*(max(rad.grid.Y[rad.grid.use])-
        min(rad.grid.Y[rad.grid.use]))
      )

      #plot(rad.grid.X[rad.grid.use], rad.grid.Z[rad.grid.use], col=
      #colorp[1+round(50*((rad.grid.Vmag[rad.grid.use]-min(rad.grid.Vmag[rad.grid.use]))/(max(rad.grid.Vmag[rad.grid.use])-
      #     min(rad.grid.Vmag[rad.grid.use]))))]#
      #, xlim=x.lim, xlab= "X-Coord (m)", ylab="Y-Coord (m)"
      #)
      #points(rad.grid.X, rad.grid.BedLevel,pch=19, cex=.5)
      #points(rad.grid.X, rad.grid.WSE,pch=19, cex=.5)
      #legend.text =as.character(round(seq(min(rad.grid.Vmag),   max(rad.grid.Vmag),by=(max(rad.grid.Vmag)-min(rad.grid.Vmag))/5),2))
      #legend.col = colorp[seq(1,51, by=10)]
      #legend("topright",legend.text, pch=19, col=legend.col, title=" Vel Mag(m/s)", bg="white")

      plot(rad.grid.X[rad.grid.use], rad.grid.Z[rad.grid.use], col=
        colorp[1+round(50*((rad.grid.mean.normV[rad.grid.use]-min(rad.grid.mean.normV[rad.grid.use]))/(max(rad.grid.mean.normV[rad.grid.use])-
        min(rad.grid.mean.normV[rad.grid.use]))))]#
        , xlim=x.lim, xlab= "X-Coord (m)", ylab="Z-Coord (m)", main = paste('idx = ', idx, sep = '')
      )
      points(rad.grid.X, rad.grid.BedLevel,pch=19, cex=.5)
      points(rad.grid.X, rad.grid.WSE,pch=19, cex=.5)

      legend.text =as.character(round(seq(min(rad.grid.mean.normV),   max(rad.grid.mean.normV),
        by=(max(rad.grid.mean.normV)-min(rad.grid.mean.normV))/5),2)
      )
      legend.col = colorp[seq(1,51, by=10)]
      legend("topright",legend.text, pch=19, col=legend.col, title=" Ave Velocity(m/s)", bg="white")

      #plot(rad.grid.Y[rad.grid.use], rad.grid.Z[rad.grid.use], col=
      #colorp[1+round(50*((rad.grid.Xvel[rad.grid.use]-min(rad.grid.Xvel[rad.grid.use]))/(max(rad.grid.Xvel[rad.grid.use])-
      #     min(rad.grid.Xvel[rad.grid.use]))))], xlim=y.lim,
      #     xlab= "Y-Coord (m)", ylab = "X-Coord (m)")
      #points(rad.grid.Y, rad.grid.BedLevel, pch=19, col="black", cex=.5)
      #points(rad.grid.Y, rad.grid.WSE, pch=19, cex=.5)
      #
      #legend.text =as.character(round(seq(min(rad.grid.Xvel),   max(rad.grid.Xvel),by=(max(rad.grid.Xvel)-min(rad.grid.Xvel))/5),2))
      #legend.col = colorp[seq(1,51, by=10)]
      #legend("topright",legend.text, pch=19, col=legend.col, title="X Velocity(m/s)")
      #

      plot(rad.grid.Y[rad.grid.use], rad.grid.Z[rad.grid.use], col=
        colorp[1+round(50*((rad.grid.mean.normV[rad.grid.use]-min(rad.grid.mean.normV[rad.grid.use]))/(max(rad.grid.mean.normV[rad.grid.use])-
        min(rad.grid.mean.normV[rad.grid.use]))))]#
        , xlim=y.lim, xlab= "Y-Coord (m)", ylab="Z-Coord (m)"
      )
      points(rad.grid.Y, rad.grid.BedLevel,pch=19, cex=.5)
      points(rad.grid.Y, rad.grid.WSE,pch=19, cex=.5)

      legend.text =as.character(round(seq(min(rad.grid.mean.normV),   max(rad.grid.mean.normV),
        by=(max(rad.grid.mean.normV)-min(rad.grid.mean.normV))/5),2)
      )
      legend.col = colorp[seq(1,51, by=10)]
      legend("topright",legend.text, pch=19, col=legend.col, title="Ave Velocity (m/s)", bg="white")

      plot(Grid3D$X, Grid3D$Y, col="gray",#, xlim=x.lim, ylim=y.lim,
        main=paste("X-Y-Z index", idx,": Radial Grid Location and Normal Velocity Component"),
        xlab="X-Coord (m)", ylab="Y-Coord(m)"
      )
      points(rad.grid.X[rad.grid.use], rad.grid.Y[rad.grid.use], col= 4*(rad.grid.use==T))

      scalar = 10

      for (i in 1: pts) {
      
        for (j in 1:Rmax) {
        
          if(rad.grid.use[i,j]) {
           
            lines(
              c(rad.grid.X[i,j], rad.grid.X[i,j] +scalar*rad.grid.Xvel[i,j]),
              c(rad.grid.Y[i, j], rad.grid.Y[i, j] + scalar*rad.grid.Yvel[i, j])
            )
            
            lines(
              c(rad.grid.X[i,j], rad.grid.X[i,j]+scalar*cos(Grid3D$Dir[idx]-pi/2)*Vel.Norm[i,j]),
              c(rad.grid.Y[i,j], rad.grid.Y[i,j]+scalar*sin(Grid3D$Dir[idx]-pi/2)*Vel.Norm[i,j]), col=4*(rad.grid.use==T)
            )
      
          }  # end if rad.grid.use
        }  # end of for j
      }  # end of for i
    }  # end of if ((idx/plot.every)...
  }  # end of if (plots)

  return(list(
    "rad.grid.origin.velocity" = Grid3D$Velocity.Magnitude[idx],
    "xy.water.column.depth" = Grid3D$WSE[idx] - Grid3D$BedLevel[idx],
    "rad.grid.origin.wse" = Grid3D$WSE[idx],
    "rad.grid.origin.Z" = Grid3D$Z[idx],
    "rad.grid.origin.bed.level" = Grid3D$BedLevel[idx],
    "rad.grid.X"= rad.grid.X, 
    "rad.grid.Y"=  rad.grid.Y,
    "rad.grid.WSE"= rad.grid.WSE, 
    "rad.grid.BedLevel" =rad.grid.BedLevel,
    "rad.grid.Xvel"=   rad.grid.Xvel,
    "rad.grid.Yvel"= rad.grid.Yvel,
    "rad.grid.Vel.Norm"= Vel.Norm,
    "rad.grid.mean.normV" = rad.grid.mean.normV, 
    "rad.grid.use" = rad.grid.use
  ))

}


#-------------------------------------------------------------------------------
# Section III:  Nrei functions
#-------------------------------------------------------------------------------


#---------------------------------------------------------------------------
# Function to estimate maximum ration for steelhead based on equations
#  from FishBioenergetics 3.0.  Here, we account for temperature's influence
#  on maximum feeding rate and set the p-value to 1 to simulate feeding
#  at a maximal rate; Can choose between Railsback and Rose 1999 or
#  Rand 1993 values by (un)commenting relevant lines below
#
# Args:
#  none; but assumes kFishWeight has been set up; Also, prey 
#    energy density is hard-coded using a value from RR'99.  May want
#    to change this in the future.
#
# Returns:
#  an estimate of the maximal consumption (J/h)
#---------------------------------------------------------------------------

GetCMaxJh <- function() {
  
  # consumption equation 3 from fishBioE 3.0 for cool/cold water species

  # prey energy density
  kPreyEnergyDens <- 2500  # J/g wet body mass of prey items
	
  # # parameters from fishBioE 3.0 and Railsback and Rose 1999
	# kCA <- 0.628    # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
	# kCB <- -0.3     # slope of allometric mass function (aka coefficient of mass dependence)
	# kCQ <- 3.5		
	# kCTO <- 25      # h20 temp correspoinding to 0.98 of the max consump. rate
	# kCTM <- 22.5
	# kCTL <- 24.3
	# kCK1 <- 0.2
	# kCK4 <- 0.2
  
  # paremeters from Rand 1993
	kCA <- 0.628      # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
	kCB <- -0.3       # slope of allometric mass function (aka coefficient of mass dependence)
	kCQ <- 5		
	kCTO <- 20        # h20 temp correspoinding to 0.98 of the max consump. rate
	kCTM <- 20
	kCTL <- 24
	kCK1 <- 0.33
	kCK4 <- 0.2

	G1 <- (1 / (kCTO - kCQ)) * (log((0.98 * (1 - kCK1)) / (kCK1 * 0.02)))
	L1 <- exp(G1 * (temps.to.sim - kCQ))
	KA <- (kCK1 * L1) / (1 + kCK1 * (L1 - 1))
	G2 <- (1 / (kCTL - kCTM)) * (log((0.98 * (1 - kCK4)) / (kCK4 * 0.02)))
	L2 <- exp(G2 * (kCTL - temps.to.sim))
	KB <- (kCK4 * L2) / (1 + kCK4 * (L2 - 1))

  # maximal consumption, unconstrained by temperature effects
	c.max <- kCA * (kFishWeight ** kCB)		#max specific feeding rate (g_prey/g_pred/d)
  
  # consumption, constrained by temperature and pval
  #  units of (g_prey / g_pred / day); p-val = 1
  spec.cons.rate <- (c.max * 1.0 * KA * KB)
  
  # set feeding period for converting spec.cons.rate to J/h
  feeding.period <- kFeedingPeriod  # in hours
  
  # converting spec.cons.rate to J/h to be compatible with units of grei and sc;
  #  note: might want to change this from 24 hours to 8 or 10 in the future as
  #    24 likely underestimates the maximal feeding rate; check the FishBioE 3.0
  #    literature for confirmation.
  cons.jh <- spec.cons.rate * kPreyEnergyDens * kFishWeight / feeding.period  # in J/h
  
  return(
    matrix(
      cons.jh, nrow = 1, byrow = TRUE, dimnames = list('J/h', c.max.col.names)
    )
  )
  
}

#---------------------------------------------------------------------------
# This is a parameterized version of GetCMaxJh, which allows it to be more
#  flexible (e.g. for use with the apply family of fuctions). This version
#  returns an unnamed matrix.
#
# Function to estimate maximum ration for steelhead based on equations
#  from FishBioenergetics 3.0.  Here, we account for temperature's influence
#  on maximum feeding rate and set the p-value to 1 to simulate feeding
#  at a maximal rate; Can choose between Railsback and Rose 1999 or
#  Rand 1993 values by (un)commenting relevant lines below
#
# Args:
#  in.temps: a vector of temperatures (C) at which Cmax is to be calculated
#  in.weight: weight (g) of the fish being simulated
#
# Returns:
#  a matrix with one estimate of the maximal consumption (J/h) for each
#    temperature in in.temps for a fish that weights in.weight
#---------------------------------------------------------------------------

GetCMaxJh_params <- function(in.temps, in.weight) {
  
  # consumption equation 3 from fishBioE 3.0 for cool/cold water species

  # prey energy density
  kPreyEnergyDens <- 2500  # J/g wet body mass of prey items
	
  # # parameters from fishBioE 3.0 and Railsback and Rose 1999
	# kCA <- 0.628    # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
	# kCB <- -0.3     # slope of allometric mass function (aka coefficient of mass dependence)
	# kCQ <- 3.5		
	# kCTO <- 25      # h20 temp correspoinding to 0.98 of the max consump. rate
	# kCTM <- 22.5
	# kCTL <- 24.3
	# kCK1 <- 0.2
	# kCK4 <- 0.2
  
  # paremeters from Rand 1993
	kCA <- 0.628      # intercept of mass dependence funtion for a 1 gram fish at the optimum h20 temp
	kCB <- -0.3       # slope of allometric mass function (aka coefficient of mass dependence)
	kCQ <- 5		
	kCTO <- 20        # h20 temp correspoinding to 0.98 of the max consump. rate
	kCTM <- 20
	kCTL <- 24
	kCK1 <- 0.33
	kCK4 <- 0.2

	G1 <- (1 / (kCTO - kCQ)) * (log((0.98 * (1 - kCK1)) / (kCK1 * 0.02)))
	L1 <- exp(G1 * (in.temps - kCQ))
	KA <- (kCK1 * L1) / (1 + kCK1 * (L1 - 1))
	G2 <- (1 / (kCTL - kCTM)) * (log((0.98 * (1 - kCK4)) / (kCK4 * 0.02)))
	L2 <- exp(G2 * (kCTL - in.temps))
	KB <- (kCK4 * L2) / (1 + kCK4 * (L2 - 1))

  # maximal consumption, unconstrained by temperature effects
	c.max <- kCA * (in.weight ** kCB)		#max specific feeding rate (g_prey/g_pred/d)
  
  # consumption, constrained by temperature and pval
  #  units of (g_prey / g_pred / day); p-val = 1
  spec.cons.rate <- (c.max * 1.0 * KA * KB)
  
  # set feeding period for converting spec.cons.rate to J/h
  feeding.period <- kFeedingPeriod  # in hours
  
  # converting spec.cons.rate to J/h to be compatible with units of grei and sc;
  #  note: might want to change this from 24 hours to 8 or 10 in the future as
  #    24 likely underestimates the maximal feeding rate; check the FishBioE 3.0
  #    literature for confirmation.
  cons.jh <- spec.cons.rate * kPreyEnergyDens * in.weight / feeding.period  # in J/h
  
  return(
    matrix(
      cons.jh, nrow = 1, byrow = TRUE, dimnames = list('J/h', c.max.col.names)
    )
  )
  
}

#----------------------------------------------------------------------------
# Function to...
#
# Args:
#  arg1:  
#
# Returns:
#  returns something
#----------------------------------------------------------------------------


GetFishSuccessMatrix <- function() {

  # if the fish can swim there and it is deep enough a fish might actually try
  #  to forage there...
  if (rad.grid$rad.grid.origin.velocity < v.max &
    rad.grid$xy.water.column.depth > min.forage.depth) {
  
    # calculate the max capture distance matrix
    mcd.array <- sqrt(
      rd * rd * (v.max * v.max - rad.grid$rad.grid.mean.normV * rad.grid$rad.grid.mean.normV) / 
        (v.max * v.max + rad.grid$rad.grid.Vel.Norm * rad.grid$rad.grid.Vel.Norm - 
          rad.grid$rad.grid.mean.normV * rad.grid$rad.grid.mean.normV)
    )
    
    # fish could be successful at any contiguous radial points where
    #  the mcd values were greater than the actual radial distances;
    #  we multiply by rad.grid$rad.grid.use to eliminate radial points that are
    #  out of the channel or water
    fish.success.mat <- (mcd.array > Init$rad.grid.r) * rad.grid$rad.grid.use
    # fish.success.mat
    
    # NAs could be generated if dist-weighted avg vel becomes too high somewhere
    #  on the radial; remove these by setting them equal to zero (because the
    #  fish can't make it to these spots)
    if (any(is.na(fish.success.mat))) {
      fish.success.mat[which(is.na(fish.success.mat))] <- 0
    }

    # sometimes a fish crossing from slow to fast to slow waters will generate
    #  a dist-weighted avg vel that goes up then down again; sometimes these
    #  up/down patterns result in a dist-weighted avg vel that makes it look
    #  like the fish can make it, then not make it, then make it again; this is
    #  an artifact of the dist-weighted velocity approach; remove these false
    #  positives by setting them to zero too
    for(my.row in 1:dim(fish.success.mat)[1]) {
      
      # if there are both ones and zeros
      if (any(fish.success.mat[my.row, ] == 1) & any(fish.success.mat[my.row, ] == 0)) {
      
        first.zero.ind <- min(which(fish.success.mat[my.row, ] == 0))
        last.one.ind <- max(which(fish.success.mat[my.row, ] == 1))
        
        # if the last one is after the first zero, replace the values
        if(last.one.ind > first.zero.ind) {
          fish.success.mat[my.row, first.zero.ind:Init$Rmax] <- 0
        }
        
      }
  
    }
  
    # if focal point is too fast or the water column there is too shallow,
    #  it's not an option for success matrix is all zeros  
  } else {
  
    fish.success.mat <- Init$rad.grid.r * 0
    
  }
  
  return(fish.success.mat)
  
}

#----------------------------------------------------------------------------
# Function to...
#
# Args:
#  arg1:  
#
# Returns:
#  returns something
#----------------------------------------------------------------------------


GetFishSuccessMatrix_params <- function(in.v.max, in.rd) {

  # if the fish can swim there and it is deep enough a fish might actually try
  #  to forage there...
  if (rad.grid$rad.grid.origin.velocity < in.v.max &
    rad.grid$xy.water.column.depth > min.forage.depth) {
  
    # calculate the max capture distance matrix
    mcd.array <- sqrt(
      in.rd * in.rd * (in.v.max * in.v.max - rad.grid$rad.grid.mean.normV * rad.grid$rad.grid.mean.normV) / 
        (in.v.max * in.v.max + rad.grid$rad.grid.Vel.Norm * rad.grid$rad.grid.Vel.Norm - 
          rad.grid$rad.grid.mean.normV * rad.grid$rad.grid.mean.normV)
    )
    
    # fish could be successful at any contiguous radial points where
    #  the mcd values were greater than the actual radial distances;
    #  we multiply by rad.grid$rad.grid.use to eliminate radial points that are
    #  out of the channel or water
    fish.success.mat <- (mcd.array > Init$rad.grid.r) * rad.grid$rad.grid.use
    # fish.success.mat
    
    # NAs could be generated if dist-weighted avg vel becomes too high somewhere
    #  on the radial; remove these by setting them equal to zero (because the
    #  fish can't make it to these spots)
    if (any(is.na(fish.success.mat))) {
      fish.success.mat[which(is.na(fish.success.mat))] <- 0
    }

    # sometimes a fish crossing from slow to fast to slow waters will generate
    #  a dist-weighted avg vel that goes up then down again; sometimes these
    #  up/down patterns result in a dist-weighted avg vel that makes it look
    #  like the fish can make it, then not make it, then make it again; this is
    #  an artifact of the dist-weighted velocity approach; remove these false
    #  positives by setting them to zero too
    for(my.row in 1:dim(fish.success.mat)[1]) {
      
      # if there are both ones and zeros
      if (any(fish.success.mat[my.row, ] == 1) & any(fish.success.mat[my.row, ] == 0)) {
      
        first.zero.ind <- min(which(fish.success.mat[my.row, ] == 0))
        last.one.ind <- max(which(fish.success.mat[my.row, ] == 1))
        
        # if the last one is after the first zero, replace the values
        if(last.one.ind > first.zero.ind) {
          fish.success.mat[my.row, first.zero.ind:Init$Rmax] <- 0
        }
        
      }
  
    }
  
    # if focal point is too fast or the water column there is too shallow,
    #  it's not an option for success matrix is all zeros  
  } else {
  
    fish.success.mat <- Init$rad.grid.r * 0
    
  }
  
  return(fish.success.mat)
  
}

#----------------------------------------------------------------------------
# Function to calculate swimming costs at the focal position (note; we do
#  not take into account any increased swimming costs associated with prey
#  capture
#
# Args:
#  focal.vel: velocity at the fish's focal position
#  this function also assumes kFishLength, kFishWeight, temps.to.sim,
#  respiration.temp.optimum, and respiration.temp.lethal are already in memory
#
# Returns:
#  a matrix containing swim costs as estimated in Hayes/Hughes 2007
#----------------------------------------------------------------------------

GetSwimCostsHhSimd <- function(focal.vel) {
  
  #---------------------------------------------------------------------
  # as in Hayes/Hughes 2007
  #---------------------------------------------------------------------
  
  # parameters
  kRA <- 0.013
  kRB <- -0.217
  kRQ2 <- 2.2
  kRT <- 0.0234
  kRTO <- respiration.temp.optimum
  kRTM <- respiration.temp.lethal
  # kRTO <- 22
  # kRTM <- 26
  
  V <- (kRTM - temps.to.sim)/(kRTM - kRTO)  # this is not velocity
  Z <- (log(kRQ2)) * (kRTM - kRTO)
  Y <- (log(kRQ2)) * (kRTM - kRTO + 2)
  X <- ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400

  temp.func <- (V ** X) * (exp(X * (1 - V)))
  act <- exp(kRT * (focal.vel * 100))  # vels in cm/s

  # swim costs in g O2/g fish/day
  swim.costs.hh.ggd <- kRA * (kFishWeight ** kRB) * temp.func * act

  # convert g O2/g fish/day to J/h
  swim.costs.hh.jh <- swim.costs.hh.ggd * kFishWeight * 13565 / 24
  
  return(swim.costs.hh.jh)
  
}

#----------------------------------------------------------------------------
# This is a paramterized version of GetSwimCostsHhSimd, which allows it to be
#  more flexible (e.g. for use with the apply family of fuctions).
#
#Function to calculate swimming costs at the focal position (note; we do
#  not take into account any increased swimming costs associated with prey
#  capture
#
# Args:
#  focal.vel: velocity at the fish's focal position
#  this function also assumes kFishLength, kFishWeight, temps.to.sim,
#  respiration.temp.optimum, and respiration.temp.lethal are already in memory
#
# Returns:
#  a matrix containing swim costs as estimated in Hayes/Hughes 2007
#----------------------------------------------------------------------------

GetSwimCostsHhSimd_params <- function(in.rto, in.rtl, in.temp, in.focal.vels,
  in.weight) {
  
  #---------------------------------------------------------------------
  # as in Hayes/Hughes 2007
  #---------------------------------------------------------------------
  
  # parameters
  kRA <- 0.013
  kRB <- -0.217
  kRQ2 <- 2.2
  kRT <- 0.0234
  kRTO <- in.rto
  kRTM <- in.rtl
  # kRTO <- 22
  # kRTM <- 26
  
  V <- (kRTM - in.temp)/(kRTM - kRTO)  # this is not velocity
  Z <- (log(kRQ2)) * (kRTM - kRTO)
  Y <- (log(kRQ2)) * (kRTM - kRTO + 2)
  X <- ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400

  temp.func <- (V ** X) * (exp(X * (1 - V)))
  act <- exp(kRT * (in.focal.vels * 100))  # vels in cm/s

  # swim costs in g O2/g fish/day
  swim.costs.hh.ggd <- kRA * (in.weight ** kRB) * temp.func * act

  # convert g O2/g fish/day to J/h
  swim.costs.hh.jh <- swim.costs.hh.ggd * in.weight * 13565 / 24
  
  return(swim.costs.hh.jh)
  
}

#----------------------------------------------------------------------------
# Function to predict foraging locations in a reach 
#
# Args:
#  none; but assumes a number of variables have been set up already
#  
#
# Returns:
#  a matrix of fish locations
#----------------------------------------------------------------------------

GetFishLocMatrix <- function() {

  # required packages 
  
  my.packs <- c(
    'fields'
  )

  if (any(!my.packs %in% installed.packages()[, 'Package'])) {
    install.packages(
      my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
      dependencies = TRUE
    )
  }

  # start placing fish; if there's only one spot
  if (nrow(ord.nrei3d.data.sub) == 1) {
  
    # there's only one
    num.fish <- 1
    fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]

  } else {
    
    # load field package
    require(fields)
    
    # highest nrei is the first acceptable position
    fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]
    
    # # plot that point
    # points(fish.loc.mat[, 2], fish.loc.mat[, 3], col = 'red', pch = 20,
      # cex = 0.5)
    
    # for every other row (in our matrix of rows with nreis above the
    #  threshold)..
    for (row.no in 2:nrow(ord.nrei3d.data.sub)) {
      
      candidate.point <- ord.nrei3d.data.sub[row.no, 2:4, drop = FALSE]
      
      dists <- rdist(fish.loc.mat[, 2:4, drop = FALSE], candidate.point)
      
      if (all(dists > fish.territory.radius)) {
      
        fish.loc.mat <- rbind(fish.loc.mat, ord.nrei3d.data.sub[row.no, ])
        # points(fish.loc.mat[dim(fish.loc.mat)[1], 2], 
          # fish.loc.mat[dim(fish.loc.mat)[1], 3], col = 'red',
          # pch = 20, cex = 0.4)
      
      }  # end of if statement
      
    }  # end of for statement
    
  }  # end of if-else statement
  
  # points(fish.loc.mat[, 2], fish.loc.mat[, 3], col = 'white', pch = 20, cex = 0.5)
  
  return(fish.loc.mat)
  
}

#---------------------------------------------------------------------------
# Function to estimate size-, drift-, and temperature-specific NREI predictions
#  for the multiFish version of the CHaMP/ISEMP NREI model. You must run
#  the multiFish version of NREI prior to using this function in order to
#  create the output needed for this function to work.
# 
# Background:
#  The multiFish version of the NREI model outputs NREI predictions for 
#    nice, round-number values of fish length, drift, and temperature. This
#    function can be used to fill in the gaps between the nice, round-number
#    values.
#
# Args:
#  fish.size.folders.dir: specifies a directory holding the multiFish NREI
#    output (should be folders named after the pattern Xmm_Xg)
#  fishLength_m: fish length (in m) for which NREI estimates are desired
#  preyConc_noPM3: drift concentration (No. individuals per cubic meter of
#    water) for which NREI estimates are desired
#  temp_C: temperature (in deg C) for which NREI estimates are desired
#
# Returns:
#  a data frame of NREI predictions for the fish.length, preyConc_noPM3, and
#    temp_C combination specified by the user
#---------------------------------------------------------------------------

GetNreiValues <- function(fish.size.folders.dir, fishLength_m, preyConc_noPM3,
  temp_C) {

  # for printing of UTMs
  options(digits = 10)
  
  # required packages 
  
  my.packs <- c(
    'readr'
  )

  if (any(!my.packs %in% installed.packages()[, 'Package'])) {
    install.packages(
      my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
      dependencies = TRUE
    )
  }

  # faster reading of large text files 
  require(readr)

  # extract simulated lengths in this fish.size.folders.dir
  fls.simd <- as.numeric(
    sapply(
      list.files(fish.size.folders.dir),
      function(x) {strsplit(x, split = 'mm')[[1]][1]}
    )
  )

  my.fish.length.mm <- round(fishLength_m * 1000)

  # if the user is asking about a fish length we explicitly simulated...
  if (my.fish.length.mm %in% fls.simd) {

    # find out which fish length it is 
    folder.match <- which(fls.simd == my.fish.length.mm)
    folder.to.use <- list.files(fish.size.folders.dir)[folder.match]

    # read nrei vals for this fish size
    nrei.vals.dir <- file.path(fish.size.folders.dir, folder.to.use, 'raw_outputs.zip')
    nrei.vals.dat <- read_csv(
      unz(nrei.vals.dir, filename = 'raw_outputs/nrei_js.csv')
    )
    # summary(nrei.vals.dat)

    # determine the drift and temp values that were simulated
    
    # extract simulated drifts/temps from the nrei column names 
    nrei.colnames <- grep(pattern = 'nrei', names(nrei.vals.dat), value = TRUE)
    
    # get the simulated drifts 
    simd.drifts <- as.numeric(
      gsub(
        pattern = 'd',
        replacement = '',
        unique(
          sapply(
            nrei.colnames,
            function(x) {
              strsplit(x, split = '_')[[1]][2]
            }
          )
        )
      )
    )

    # get the simulated temps 
    simd.temps <- as.numeric(
      gsub(
        pattern = 't',
        replacement = '',
        gsub(
          pattern = 'C',
          replacement = '',
          unique(
            sapply(
              nrei.colnames,
              function(x) {
                strsplit(x, split = '_')[[1]][3]
              }
            )
          )
        )
      )
    )

    # if the user is asking about a fish length we explicitly simulated and 
    #  a drift value we explicitly simulated 
    if (preyConc_noPM3 %in% simd.drifts) {

      # get the drift portion of the column name
      drift.name.match <- paste(
        'd',
        format(preyConc_noPM3, trim = TRUE, nsmall = 2),
        sep = ''
      )

      # if the user is asking about a fish length we explicitly simulated, a
      #  drift value we explicitly simulated, and a temp value we explicitly
      #  simulated... 
      if (temp_C %in% simd.temps) {

        # get the temp portion of the column name
        temp.name.match <- paste(
          't',
          format(temp_C, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        # find the column using drift and temp column name info 
        col.to.sel <- intersect(
          grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match, nrei.colnames, value = TRUE)
        )

        interp.nrei.vals <- nrei.vals.dat[, col.to.sel]
        interp.nrei.vals.name <- col.to.sel

        interp.nrei.out <- cbind(
          nrei.vals.dat[, c('idx', 'X', 'Y', 'Z')],
          interp.nrei.vals
        )
        names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
        # head(interp.nrei.out)

      # if the user is asking about a fish size we explicitly simulated, a 
      #  drift val we explicitly simulated, and a temp value we did not
      #  explicitly simulate but is within the simulated limits...
      } else if (temp_C >= min(simd.temps) & temp_C <= max(simd.temps)) {

        # get the nearest simulated temp values 
        simd.temp.lower <- simd.temps[max(which(simd.temps < temp_C))]
        simd.temp.higher <- simd.temps[min(which(simd.temps > temp_C))]

        # get their distances from temp_C and corresponding weights 
        weight.lower.t <- 1.0 - ((temp_C - simd.temp.lower) / 
          (simd.temp.higher - simd.temp.lower))
        weight.higher.t <- 1.0 - ((simd.temp.higher - temp_C) /
          (simd.temp.higher - simd.temp.lower))

        # get the temp 'keys'
        temp.name.match.lower <- paste(
          't',
          format(simd.temp.lower, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        temp.name.match.higher <- paste(
          't',
          format(simd.temp.higher, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        # get the column names 
        col.to.sel.lower <- intersect(
          grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
        )

        col.to.sel.higher <- intersect(
          grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
        )

        nrei.vals.lower.temp <- nrei.vals.dat[, col.to.sel.lower]
        nrei.vals.higher.temp <- nrei.vals.dat[, col.to.sel.higher]

        interp.nrei.vals <- (weight.lower.t * nrei.vals.lower.temp) +
          (weight.higher.t * nrei.vals.higher.temp)

        temp.title <- paste(
          't',
          format(temp_C, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        interp.nrei.vals.name <- paste('nrei', drift.name.match, temp.title,
          'js', sep = '_')

        interp.nrei.out <- cbind(
          nrei.vals.dat[, c('idx', 'X', 'Y', 'Z')],
          interp.nrei.vals
        )
        names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
        # format(head(interp.nrei.out, 100), scientific = FALSE)

      # if the user is asking about a fish length we simulated, a drift value
      #  we simulated, and a temp value that is outside the range we 
      #  explicitly simualted...
      } else {

        interp.nrei.out <- writeLines(
          c(
            'The temp you are trying to interpolate values for is out of range',
            'based on the data available in the folder:',
            '',
            fish.size.folders.dir,
            ''
          )
        )

      }

    # if the user is asking about a fish length we explicitly simulated and a
    #  drift concentration we did not explicitly simulate, but within the 
    #  drift ranges found in the folder...
    } else if (preyConc_noPM3 <= max(simd.drifts) & preyConc_noPM3 >= min(simd.drifts)) {

      # if the user is asking about a fish length we explicitly simulated,
      #  a drift value we did not simulate, and a temp value we explicitly
      #  simulated... 
      if (temp_C %in% simd.temps) {

        # get the nearest simulated drift values 
        simd.drift.lower <- simd.drifts[max(which(simd.drifts < preyConc_noPM3))]
        simd.drift.higher <- simd.drifts[min(which(simd.drifts > preyConc_noPM3))]

        # get their distances from preyConc_noPM3 and corresponding weights 
        weight.lower.d <- 1.0 - ((preyConc_noPM3 - simd.drift.lower) / 
          (simd.drift.higher - simd.drift.lower))
        weight.higher.d <- 1.0 - ((simd.drift.higher - preyConc_noPM3) /
          (simd.drift.higher - simd.drift.lower))

        # get the drift 'keys'
        drift.name.match.lower <- paste(
          'd',
          format(simd.drift.lower, trim = TRUE, nsmall = 2),
          sep = ''
        )

        drift.name.match.higher <- paste(
          'd',
          format(simd.drift.higher, trim = TRUE, nsmall = 2),
          sep = ''
        )

        # get the temp portion of the column name
        temp.name.match <- paste(
          't',
          format(temp_C, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        # get the column names 
        col.to.sel.lower <- intersect(
          grep(pattern = temp.name.match, nrei.colnames, value = TRUE),
          grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE)
        )

        col.to.sel.higher <- intersect(
          grep(pattern = temp.name.match, nrei.colnames, value = TRUE),
          grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE)
        )

        nrei.vals.lower.temp <- nrei.vals.dat[, col.to.sel.lower]
        nrei.vals.higher.temp <- nrei.vals.dat[, col.to.sel.higher]

        interp.nrei.vals <- (weight.lower.d * nrei.vals.lower.temp) +
          (weight.higher.d * nrei.vals.higher.temp)

        drift.title <- paste(
          'd',
          format(preyConc_noPM3, trim = TRUE, nsmall = 2),
          sep = ''
        )

        interp.nrei.vals.name <- paste('nrei', drift.title, temp.name.match,
          'js', sep = '_')

        interp.nrei.out <- cbind(
          nrei.vals.dat[, c('idx', 'X', 'Y', 'Z')],
          interp.nrei.vals
        )
        names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
        # format(head(interp.nrei.out, 100), scientific = FALSE)


      # if the user is asking about a fish length we simulated, a drift
      #  concentration we did not explicitly simulate, and a temp  value 
      #  we did not explicitly simulate but within the simulated ranges...
      } else if (temp_C >= min(simd.temps) & temp_C <= max(simd.temps)) {

        # get the nearest simulated drift values 
        simd.drift.lower <- simd.drifts[max(which(simd.drifts < preyConc_noPM3))]
        simd.drift.higher <- simd.drifts[min(which(simd.drifts > preyConc_noPM3))]

        # get their distances from preyConc_noPM3 and corresponding weights 
        weight.lower.d <- 1.0 - ((preyConc_noPM3 - simd.drift.lower) / 
          (simd.drift.higher - simd.drift.lower))
        weight.higher.d <- 1.0 - ((simd.drift.higher - preyConc_noPM3) /
          (simd.drift.higher - simd.drift.lower))

        # get the drift portions of the column names 
        drift.name.match.lower <- paste(
          'd',
          format(simd.drift.lower, trim = TRUE, nsmall = 2),
          sep = ''
        )

        drift.name.match.higher <- paste(
          'd',
          format(simd.drift.higher, trim = TRUE, nsmall = 2),
          sep = ''
        )

        # get the nearest simulated temp values 
        simd.temp.lower <- simd.temps[max(which(simd.temps < temp_C))]
        simd.temp.higher <- simd.temps[min(which(simd.temps > temp_C))]

        # get their distances from temp_C and corresponding weights 
        weight.lower.t <- 1.0 - ((temp_C - simd.temp.lower) / 
          (simd.temp.higher - simd.temp.lower))
        weight.higher.t <- 1.0 - ((simd.temp.higher - temp_C) /
          (simd.temp.higher - simd.temp.lower))

        # get the temp portions of the column names 
        temp.name.match.lower <- paste(
          't',
          format(simd.temp.lower, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        temp.name.match.higher <- paste(
          't',
          format(simd.temp.higher, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        # get the column names 
        col.to.sel.ldlt <- intersect(
          grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
        )

        col.to.sel.ldht <- intersect(
          grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
        )

        col.to.sel.hdlt <- intersect(
          grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
        )

        col.to.sel.hdht <- intersect(
          grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE),
          grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
        )

        # nrei predictions we will use for interpolation 
        nrei.vals.ldlt <- nrei.vals.dat[, col.to.sel.ldlt]
        nrei.vals.ldht <- nrei.vals.dat[, col.to.sel.ldht]
        nrei.vals.hdlt <- nrei.vals.dat[, col.to.sel.hdlt]
        nrei.vals.hdht <- nrei.vals.dat[, col.to.sel.hdht]

        # interpolation between the drift values for the lower temperature
        interp.temp.lower <- (weight.lower.d * nrei.vals.ldlt) +
          (weight.higher.d * nrei.vals.hdlt)

        # interpolation between the drift values for the higher temperature
        interp.temp.higher <- (weight.lower.d * nrei.vals.ldht) +
          (weight.higher.d * nrei.vals.hdht)

        # interpolation between the temp values (these are the final values)
        interp.nrei.vals <- (weight.lower.t * interp.temp.lower) +
          (weight.higher.t * interp.temp.higher)

        # name the output column
        drift.title <- paste(
          'd',
        format(preyConc_noPM3, trim = TRUE, nsmall = 2),
        sep = ''
        )

        temp.title <- paste(
          't',
          format(temp_C, trim = TRUE, nsmall = 1),
          'C',
          sep = ''
        )

        interp.nrei.vals.name <- paste('nrei', drift.title, temp.title,
          'js', sep = '_')

        # create the output data frame 
        interp.nrei.out <- cbind(
          nrei.vals.dat[, c('idx', 'X', 'Y', 'Z')],
          interp.nrei.vals
        )
        names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
        # format(head(interp.nrei.out, 100), scientific = FALSE)

      # if the user is asking about a fish length we simulated, a drift
      #  concentration we did not explicitly simulate, and a temp  value 
      #  outside the simulated ranges... 
      } else {

        interp.nrei.out <- writeLines(
          c(
            'The temp you are trying to interpolate values for is out of range',
            'based on the data available in the folder:',
            '',
            fish.size.folders.dir,
            ''
          )
        )

      }
    
    # end of if the user wants an un-sim'd drift within appropriate interpolation
    #  ranges...
    } else {

      interp.nrei.out <- writeLines(
        c(
          'The drift you are trying to interpolate values for is out of range',
          'based on the data available in the folder:',
          '',
          fish.size.folders.dir,
          ''
        )
      )

    }  

  # if the user is asking about a fish size we did not explicitly simulate but
  #  is within the max/min found in the fish sizes folder...
  } else if (my.fish.length.mm <= max(fls.simd) & my.fish.length.mm >= min(fls.simd)) {

    # get the nearest simulated length values 
    simd.fl.lower <- fls.simd[max(which(fls.simd < my.fish.length.mm))]
    simd.fl.higher <- fls.simd[min(which(fls.simd > my.fish.length.mm))]

    # get their distances from my.fish.length.mm and corresponding weights 
    weight.lower.fl <- 1.0 - ((my.fish.length.mm - simd.fl.lower) / 
      (simd.fl.higher - simd.fl.lower))
    weight.higher.fl <- 1.0 - ((simd.fl.higher - my.fish.length.mm) /
      (simd.fl.higher - simd.fl.lower))

    # get the fish size folder names 
    fl.dir.match.lower <- grep(
      pattern = paste(simd.fl.lower, 'mm_', sep = ''),
      list.files(fish.size.folders.dir),
      value = TRUE
    )

    fl.dir.match.higher <- grep(
      pattern = paste(simd.fl.higher, 'mm_', sep = ''),
      list.files(fish.size.folders.dir),
      value = TRUE
    )

    # read nrei vals for the two fish sizes
    nrei.vals.dir.fll <- file.path(fish.size.folders.dir, fl.dir.match.lower, 'raw_outputs.zip')
    nrei.vals.dat.fll <- read_csv(
      unz(nrei.vals.dir.fll, filename = 'raw_outputs/nrei_js.csv')
    )
    # summary(nrei.vals.dat.fll)

    nrei.vals.dir.flh <- file.path(fish.size.folders.dir, fl.dir.match.higher, 'raw_outputs.zip')
    nrei.vals.dat.flh <- read_csv(
      unz(nrei.vals.dir.flh, filename = 'raw_outputs/nrei_js.csv')
    )
    # summary(nrei.vals.dat.flh)

    # extract simulated drifts/temps from the nrei column names 
    nrei.colnames.fll <- grep(pattern = 'nrei', names(nrei.vals.dat.fll), value = TRUE)
    nrei.colnames.flh <- grep(pattern = 'nrei', names(nrei.vals.dat.flh), value = TRUE)

    # if fish size folders within the fish.size.folders.dir have the same column
    #  headers and idx values...
    if(all(nrei.colnames.fll == nrei.colnames.flh) & all(nrei.vals.dat.fll$idx == nrei.vals.dat.flh$idx)) {

      # if they're the same, we can use either one to extract the sim'd values 
      nrei.colnames <- nrei.colnames.flh
      
      # get the simulated drifts 
      simd.drifts <- as.numeric(
        gsub(
          pattern = 'd',
          replacement = '',
          unique(
            sapply(
              nrei.colnames,
              function(x) {
                strsplit(x, split = '_')[[1]][2]
              }
            )
          )
        )
      )

      # get the simulated temps 
      simd.temps <- as.numeric(
        gsub(
          pattern = 't',
          replacement = '',
          gsub(
            pattern = 'C',
            replacement = '',
            unique(
              sapply(
                nrei.colnames,
                function(x) {
                  strsplit(x, split = '_')[[1]][3]
                }
              )
            )
          )
        )
      )

      # if the user is asking about a fish length we didn't simulate and a drift
      #  value we did simulate... 
      if (preyConc_noPM3 %in% simd.drifts) {

        # get the drift portion of the column name
        drift.name.match <- paste(
          'd',
          format(preyConc_noPM3, trim = TRUE, nsmall = 2),
          sep = ''
        )

        # if the user is asking about one of the temp values we explicitly
        #  simulated...
        if (temp_C %in% simd.temps) {

          # get the temp portion of the column name
          temp.name.match <- paste(
            't',
            format(temp_C, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          # find the column using drift and temp column name info 
          col.to.sel <- intersect(
            grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match, nrei.colnames, value = TRUE)
          )

          nrei.vals.fll <- nrei.vals.dat.fll[, col.to.sel]
          nrei.vals.flh <- nrei.vals.dat.flh[, col.to.sel]
          # cbind(head(nrei.vals.fll), head(nrei.vals.flh))
          interp.nrei.vals.name <- col.to.sel

          interp.nrei.vals <- (weight.lower.fl * nrei.vals.fll) +
            (weight.higher.fl * nrei.vals.flh)

          interp.nrei.out <- cbind(
            nrei.vals.dat.flh[, c('idx', 'X', 'Y', 'Z')],
            interp.nrei.vals
          )
          names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
          # format(head(interp.nrei.out, 100), scientific = FALSE)

        # if the user is asking about a drift val we explicitly simulated and 
        #  a temp value we did not explicitly simulate...
        } else if (temp_C >= min(simd.temps) & temp_C <= max(simd.temps)) {

          # get the nearest simulated temp values 
          simd.temp.lower <- simd.temps[max(which(simd.temps < temp_C))]
          simd.temp.higher <- simd.temps[min(which(simd.temps > temp_C))]

          # get their distances from temp_C and corresponding weights 
          weight.lower.t <- 1.0 - ((temp_C - simd.temp.lower) / 
            (simd.temp.higher - simd.temp.lower))
          weight.higher.t <- 1.0 - ((simd.temp.higher - temp_C) /
            (simd.temp.higher - simd.temp.lower))

          # get the temp 'keys'
          temp.name.match.lower <- paste(
            't',
            format(simd.temp.lower, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          temp.name.match.higher <- paste(
            't',
            format(simd.temp.higher, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          # get the column names 
          col.to.sel.lower <- intersect(
            grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
          )

          col.to.sel.higher <- intersect(
            grep(pattern = drift.name.match, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
          )

          nrei.vals.fll.lower.temp <- nrei.vals.dat.fll[, col.to.sel.lower]
          nrei.vals.fll.higher.temp <- nrei.vals.dat.fll[, col.to.sel.higher]

          interp.nrei.vals.fll <- (weight.lower.t * nrei.vals.fll.lower.temp) +
            (weight.higher.t * nrei.vals.fll.higher.temp)

          nrei.vals.flh.lower.temp <- nrei.vals.dat.flh[, col.to.sel.lower]
          nrei.vals.flh.higher.temp <- nrei.vals.dat.flh[, col.to.sel.higher]

          interp.nrei.vals.flh <- (weight.lower.t * nrei.vals.flh.lower.temp) +
            (weight.higher.t * nrei.vals.flh.higher.temp)

          interp.nrei.vals <- (weight.lower.fl * interp.nrei.vals.fll) +
            (weight.higher.fl * interp.nrei.vals.flh)

          temp.title <- paste(
            't',
            format(temp_C, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          interp.nrei.vals.name <- paste('nrei', drift.name.match, temp.title,
            'js', sep = '_')

          interp.nrei.out <- cbind(
            nrei.vals.dat.flh[, c('idx', 'X', 'Y', 'Z')],  # ok because .fll and .flh have same idxs
            interp.nrei.vals
          )
          names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
          # format(head(interp.nrei.out, 100), scientific = FALSE)

        # end of if user wants a sim'd drift and either sim'd or un-sim'd temp
        } else {

          interp.nrei.out <- writeLines(
            c(
              'The temp you are trying to interpolate values for is out of range',
              'based on the data available in the folder:',
              '',
              fish.size.folders.dir,
              ''
            )
          )

        }

      # if the user is asking about a drift concentration we did not explicitly
      #  simulate but is within the simulated ranges...
      } else if (preyConc_noPM3 <= max(simd.drifts) & preyConc_noPM3 >= min(simd.drifts)) {

        # if the user is asking about a temp we explicitly modeled...
        if (temp_C %in% simd.temps) {

          # get the temp portion of the column name
          temp.name.match <- paste(
            't',
            format(temp_C, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          # get the nearest simulated drift values 
          simd.drift.lower <- simd.drifts[max(which(simd.drifts < preyConc_noPM3))]
          simd.drift.higher <- simd.drifts[min(which(simd.drifts > preyConc_noPM3))]

          # get their distances from preyConc_noPM3 and corresponding weights 
          weight.lower.d <- 1.0 - ((preyConc_noPM3 - simd.drift.lower) / 
            (simd.drift.higher - simd.drift.lower))
          weight.higher.d <- 1.0 - ((simd.drift.higher - preyConc_noPM3) /
            (simd.drift.higher - simd.drift.lower))

          # get the drift portion of the column name 
          drift.name.match.lower <- paste(
            'd',
            format(simd.drift.lower, trim = TRUE, nsmall = 2),
            sep = ''
          )

          drift.name.match.higher <- paste(
            'd',
            format(simd.drift.higher, trim = TRUE, nsmall = 2),
            sep = ''
          )

          # get the column names 
          col.to.sel.lower <- intersect(
            grep(pattern = temp.name.match, nrei.colnames, value = TRUE),
            grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE)
          )

          col.to.sel.higher <- intersect(
            grep(pattern = temp.name.match, nrei.colnames, value = TRUE),
            grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE)
          )

          nrei.vals.fll.lower.drift <- nrei.vals.dat.fll[, col.to.sel.lower]
          nrei.vals.fll.higher.drift <- nrei.vals.dat.fll[, col.to.sel.higher]

          interp.nrei.vals.fll <- (weight.lower.d * nrei.vals.fll.lower.drift) +
            (weight.higher.d * nrei.vals.fll.higher.drift)

          nrei.vals.flh.lower.drift <- nrei.vals.dat.flh[, col.to.sel.lower]
          nrei.vals.flh.higher.drift <- nrei.vals.dat.flh[, col.to.sel.higher]

          interp.nrei.vals.flh <- (weight.lower.d * nrei.vals.flh.lower.drift) +
            (weight.higher.d * nrei.vals.flh.higher.drift)

          interp.nrei.vals <- (weight.lower.fl * interp.nrei.vals.fll) +
            (weight.higher.fl * interp.nrei.vals.flh)

          drift.title <- paste(
            'd',
          format(preyConc_noPM3, trim = TRUE, nsmall = 2),
          sep = ''
          )

          interp.nrei.vals.name <- paste('nrei', drift.title, temp.name.match,
            'js', sep = '_')

          interp.nrei.out <- cbind(
            nrei.vals.dat.flh[, c('idx', 'X', 'Y', 'Z')],
            interp.nrei.vals
          )
          names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
          # format(head(interp.nrei.out, 100), scientific = FALSE)


        # if the user is asking about a drift we did not explicitly model and 
        #  a temp we did not explicitly model...
        } else if (temp_C >= min(simd.temps) & temp_C <= max(simd.temps)) {

          # get the nearest simulated drift values 
          simd.drift.lower <- simd.drifts[max(which(simd.drifts < preyConc_noPM3))]
          simd.drift.higher <- simd.drifts[min(which(simd.drifts > preyConc_noPM3))]

          # get their distances from preyConc_noPM3 and corresponding weights 
          weight.lower.d <- 1.0 - ((preyConc_noPM3 - simd.drift.lower) / 
            (simd.drift.higher - simd.drift.lower))
          weight.higher.d <- 1.0 - ((simd.drift.higher - preyConc_noPM3) /
            (simd.drift.higher - simd.drift.lower))

          # get the drift portion of the column name 
          drift.name.match.lower <- paste(
            'd',
            format(simd.drift.lower, trim = TRUE, nsmall = 2),
            sep = ''
          )

          drift.name.match.higher <- paste(
            'd',
            format(simd.drift.higher, trim = TRUE, nsmall = 2),
            sep = ''
          )

          # get the nearest simulated temp values 
          simd.temp.lower <- simd.temps[max(which(simd.temps < temp_C))]
          simd.temp.higher <- simd.temps[min(which(simd.temps > temp_C))]

          # get their distances from temp_C and corresponding weights 
          weight.lower.t <- 1.0 - ((temp_C - simd.temp.lower) / 
            (simd.temp.higher - simd.temp.lower))
          weight.higher.t <- 1.0 - ((simd.temp.higher - temp_C) /
            (simd.temp.higher - simd.temp.lower))

          # get the temp portion of the column name 
          temp.name.match.lower <- paste(
            't',
            format(simd.temp.lower, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          temp.name.match.higher <- paste(
            't',
            format(simd.temp.higher, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          # get the column names 
          col.to.sel.ldlt <- intersect(
            grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
          )

          col.to.sel.ldht <- intersect(
            grep(pattern = drift.name.match.lower, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
          )

          col.to.sel.hdlt <- intersect(
            grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.lower, nrei.colnames, value = TRUE)
          )

          col.to.sel.hdht <- intersect(
            grep(pattern = drift.name.match.higher, nrei.colnames, value = TRUE),
            grep(pattern = temp.name.match.higher, nrei.colnames, value = TRUE)
          )

          # interpolation for lower fish length 
          nrei.vals.fll.ldlt <- nrei.vals.dat.fll[, col.to.sel.ldlt]
          nrei.vals.fll.ldht <- nrei.vals.dat.fll[, col.to.sel.ldht]
          nrei.vals.fll.hdlt <- nrei.vals.dat.fll[, col.to.sel.hdlt]
          nrei.vals.fll.hdht <- nrei.vals.dat.fll[, col.to.sel.hdht]

          interp.nrei.vals.fll.lt <- (weight.lower.d * nrei.vals.fll.ldlt) +
            (weight.higher.d * nrei.vals.fll.hdlt)
          interp.nrei.vals.fll.ht <- (weight.lower.d * nrei.vals.fll.ldht) +
            (weight.higher.d * nrei.vals.fll.hdht)

          interp.nrei.vals.fll <- (weight.lower.t * interp.nrei.vals.fll.lt) +
            (weight.higher.t * interp.nrei.vals.fll.ht)

          # interpolation for higher fish length 
          nrei.vals.flh.ldlt <- nrei.vals.dat.flh[, col.to.sel.ldlt]
          nrei.vals.flh.ldht <- nrei.vals.dat.flh[, col.to.sel.ldht]
          nrei.vals.flh.hdlt <- nrei.vals.dat.flh[, col.to.sel.hdlt]
          nrei.vals.flh.hdht <- nrei.vals.dat.flh[, col.to.sel.hdht]

          interp.nrei.vals.flh.lt <- (weight.lower.d * nrei.vals.flh.ldlt) +
            (weight.higher.d * nrei.vals.flh.hdlt)
          interp.nrei.vals.flh.ht <- (weight.lower.d * nrei.vals.flh.ldht) +
            (weight.higher.d * nrei.vals.flh.hdht)

          interp.nrei.vals.flh <- (weight.lower.t * interp.nrei.vals.flh.lt) +
            (weight.higher.t * interp.nrei.vals.flh.ht)

          # interpolate between the two fish sizes 
          interp.nrei.vals <- (weight.lower.fl * interp.nrei.vals.fll) +
            (weight.higher.fl * interp.nrei.vals.flh)

          drift.title <- paste(
            'd',
          format(preyConc_noPM3, trim = TRUE, nsmall = 2),
          sep = ''
          )

          temp.title <- paste(
            't',
            format(temp_C, trim = TRUE, nsmall = 1),
            'C',
            sep = ''
          )

          interp.nrei.vals.name <- paste('nrei', drift.title, temp.title,
            'js', sep = '_')

          interp.nrei.out <- cbind(
            nrei.vals.dat.flh[, c('idx', 'X', 'Y', 'Z')],
            interp.nrei.vals
          )
          names(interp.nrei.out)[ncol(interp.nrei.out)] <- interp.nrei.vals.name
          # format(head(interp.nrei.out, 100), scientific = FALSE)

        # end of if the user wants an un-sim'd drift and either a sim'd or
        #  un-sim'd temp
        } else {

          interp.nrei.out <- writeLines(
            c(
              'The temp you are trying to interpolate values for is out of range',
              'based on the data available in the folder:',
              '',
              fish.size.folders.dir,
              ''
            )
          )

        }

      # if the drift is outside the range we simulated, making interpolation
      #  inappropriate...
      } else {

        interp.nrei.out <- writeLines(
          c(
            'The drift you are trying to interpolate values for is out of range',
            'based on the data available in the folder:',
            '',
            fish.size.folders.dir,
            ''
          )
        )

      }

    # if the nrei values for the two fish aren't compatible for interpolation
    #  due to differing column headers or 3D grid settings...
    } else {

      interp.nrei.out <- writeLines(
        c(
          'ERROR: The two fish size folders either have different column names',
          '(meaning their simulations used different drift/temp values), or they',
          'have different idx values (meaning their simulations used different',
          '3D grid settings).',
          '',
          'Please make sure teh fish size folder are from the same set of',
          'simulations or adjust your file structure to compensate.',
          ''
        )
      )

    } 

  # if the user is asking about a fish size outside of those we have explicitly
  #  simulated, then interpolation is inappropriate, so... 
  } else {

    interp.nrei.out <- writeLines(
      c(
        'The fish size you are trying to interpolate values for is out of range',
        'based on the fish sizes available in the folder:',
        '',
        fish.size.folders.dir,
        ''
      )
    )

  }  # end of fish size possibilities...

  return(interp.nrei.out)

}  # end of GetNreiValues function


#---------------------------------------------------------------------------
# Function to return profitable XYZ foraging locations for the multiFish
#  version of the CHaMP/ISEMP NREI model. You must run the multiFish version
#  of NREI prior to using this function in order to create the output needed
#  for this function to work.
# 
# Args:
#  fish.size.folders.dir: specifies a directory holding the multiFish NREI
#    output (should be folders named after the pattern Xmm_Xg)
#  fishLength_m: fish length (in m) for which NREI estimates are desired
#  preyConc_noPM3: drift concentration (No. individuals per cubic meter of
#    water) for which NREI estimates are desired
#  temp_C: temperature (in deg C) for which NREI estimates are desired
#  nrei.limit_jps: the minimum allowable intake threshold (in J/s) for fish 
#    placement to be allowed (e.g. assuming the model is perfect, 0.0 J/s 
#    would suggest a fish is maintaining its mass exactly)
#
# Returns:
#  a data frame of predicted fish locations with their predicted NREI values 
#---------------------------------------------------------------------------

GetFishLocations <- function(fish.size.folders.dir, fishLength_m, preyConc_noPM3,
  temp_C, nrei.limit_jps) {

  # for printing of UTMs
  options(digits = 10)

  # required packages 
  my.packs <- c(
    'fields'
  )

  if (any(!my.packs %in% installed.packages()[, 'Package'])) {
    install.packages(
      my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
      dependencies = TRUE
    )
  }

  # calculate NREI values for this fish size, prey concentration, and temp 
  nrei3d.data <- GetNreiValues(fish.size.folders.dir, fishLength_m,
    preyConc_noPM3, temp_C)
  # head(nrei3d.data, 25)

  # temporarily save the nrei column name 
  true.col.name <- names(nrei3d.data)[ncol(nrei3d.data)]

  # rename to make things compatible with legacy code
  names(nrei3d.data)[ncol(nrei3d.data)] <- 'nrei.hh.js'
  
  # if there are any above our limit
  if (any(nrei3d.data$nrei.hh.js >= nrei.limit_jps)) {
  
    # order by decreasing nrei, keep those above the limit, use matrices
    ord.nrei3d.data <- nrei3d.data[order(-nrei3d.data$nrei.hh.js), ]
    ord.nrei3d.data <- ord.nrei3d.data[ord.nrei3d.data$nrei.hh.js >= nrei.limit_jps, ]
    ord.nrei3d.data.sub <- as.matrix(ord.nrei3d.data)

    #  unname it; seems necessary for the fields rdist function 
    ord.nrei3d.data.sub <- unname(ord.nrei3d.data.sub)
    rm(nrei3d.data, ord.nrei3d.data)
    # ord.nrei3d.data.sub[1:10, ]

    # calculate fish territory area for this fish length 
    # based on Keeley and McPhail '98, adjusted based on Imre, Grant, Keeley '04
    #  must get fish length in cm for this equation; territory is in m2
    fish.territory.area <- 10 ** (1.56 * log10(fishLength_m * 100) - 1.81 - 0.07)

    # get fish territory radius, multiply by 2; Multiplying by two ensures
    #  fish territories don't overlap
    fish.territory.radius <- sqrt(fish.territory.area / pi) * 2
    
    # start placing fish

    # if there's only one spot...
    if (nrow(ord.nrei3d.data.sub) == 1) {
    
      # there's only one
      num.fish <- 1
      fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]

    # if there are multiple spots...
    } else {
      
      # load field package
      require(fields)
      
      # highest nrei is the first acceptable position
      fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]
      
      # for every other row (in our matrix of rows with nreis above the
      #  threshold)..
      for (row.no in 2:nrow(ord.nrei3d.data.sub)) {
        
        candidate.point <- ord.nrei3d.data.sub[row.no, 2:4, drop = FALSE]
        
        dists <- rdist(fish.loc.mat[, 2:4, drop = FALSE], candidate.point)
        
        if (all(dists > fish.territory.radius)) {
        
          fish.loc.mat <- rbind(fish.loc.mat, ord.nrei3d.data.sub[row.no, ])
        
        }  # end of if all(dists > fish.territory.radius)....
        
      }  # end of for row.no in 2:nrow(ord.nrei3d.data.sub)...
      
    }  # end of if there are multiple spots...

  # all are below our nrei threshold limit...
  } else {

    # matrix of NAs with simulation name indicated 
    fish.loc.mat <- matrix(rep(NA, 5), nrow = 1)

  }

  # give proper column titles; convert to data frame 
  colnames(fish.loc.mat) <- c('idx', 'X', 'Y', 'Z', 'pred.nrei')
  fish.loc.mat <- data.frame(fish.loc.mat)
  fish.loc.mat$sim.name <- true.col.name
  # head(fish.loc.mat)
  
  # add prey conc info 
  fish.loc.mat$sim.drift <- as.numeric(
    sub('d', '', 
      sapply(strsplit(fish.loc.mat$sim.name, split = '_'), function(x) x[2])
    )
  )

  # add temperature info 
  fish.loc.mat$sim.temp <- as.numeric(
    sub('C', '',
      sub('t', '', 
        sapply(strsplit(fish.loc.mat$sim.name, split = '_'), function(x) x[3])
      )
    )
  )

  return(fish.loc.mat)

}

#---------------------------------------------------------------------------
# Function to estimate size-, drift-, and temperature-specific fish capacity
#  (# fish), areal fish density (fish/m2), and volumetric fish density (fish/m3)  
#  for the multiFish version of the CHaMP/ISEMP NREI model. You must run
#  the multiFish version of NREI prior to using this function in order to
#  create the output needed for this function to work.
# 
# Args:
#  fish.size.folders.dir: specifies a directory holding the multiFish NREI
#    output (should be folders named after the pattern Xmm_Xg)
#  fishLength_m: fish length (in m) for which NREI estimates are desired
#  preyConc_noPM3: drift concentration (No. individuals per cubic meter of
#    water) for which NREI estimates are desired
#  temp_C: temperature (in deg C) for which NREI estimates are desired
#
# Returns:
#  a data frame with capacity and density predictions for the fish.length,
#   preyConc_noPM3, and temp_C combination specified by the user
#---------------------------------------------------------------------------

GetCapDensPreds <- function(fish.size.folders.dir, fishLength_m, preyConc_noPM3,
  temp_C, nrei.limit_jps) {

  # for printing of UTMs
  options(digits = 10)

  # required packages 
  my.packs <- c(
    'fields'
  )

  if (any(!my.packs %in% installed.packages()[, 'Package'])) {
    install.packages(
      my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
      dependencies = TRUE
    )
  }

  # go get modeled area and volume
  temp.file.to.read <- file.path(
    fish.size.folders.dir,
    list.files(fish.size.folders.dir)[1],
    'fish_density',
    'fish_dens_preds.csv'
  )

  temp.row <- read.csv(temp.file.to.read, header = TRUE, nrows = 1)
  mod.area <- temp.row$mod.area
  mod.volume <- temp.row$mod.volume
  # print(mod.area)
  # print(mod.volume)

  # calculate NREI values for this fish size, prey concentration, and temp 
  nrei3d.data <- GetNreiValues(fish.size.folders.dir, fishLength_m,
    preyConc_noPM3, temp_C)
  # head(nrei3d.data, 25)

  # temporarily save the nrei column name 
  true.col.name <- names(nrei3d.data)[ncol(nrei3d.data)]

  # rename to make things compatible with legacy code
  names(nrei3d.data)[ncol(nrei3d.data)] <- 'nrei.hh.js'
  
  # if there are any above our limit
  if (any(nrei3d.data$nrei.hh.js >= nrei.limit_jps)) {
  
    # order by decreasing nrei, keep those above the limit, use matrices
    ord.nrei3d.data <- nrei3d.data[order(-nrei3d.data$nrei.hh.js), ]
    ord.nrei3d.data <- ord.nrei3d.data[ord.nrei3d.data$nrei.hh.js >= nrei.limit_jps, ]
    ord.nrei3d.data.sub <- as.matrix(ord.nrei3d.data)

    #  unname it; seems necessary for the fields rdist function 
    ord.nrei3d.data.sub <- unname(ord.nrei3d.data.sub)
    # rm(nrei3d.data, ord.nrei3d.data)
    # ord.nrei3d.data.sub[1:10, ]

    # calculate fish territory area for this fish length 
    # based on Keeley and McPhail '98, adjusted based on Imre, Grant, Keeley '04
    #  must get fish length in cm for this equation; territory is in m2
    fish.territory.area <- 10 ** (1.56 * log10(fishLength_m * 100) - 1.81 - 0.07)

    # get fish territory radius, multiply by 2; Multiplying by two ensures
    #  fish territories don't overlap
    fish.territory.radius <- sqrt(fish.territory.area / pi) * 2
    
    # start placing fish

    # if there's only one spot...
    if (nrow(ord.nrei3d.data.sub) == 1) {
    
      # there's only one
      num.fish <- 1
      fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]

    # if there are multiple spots...
    } else {
      
      # load field package
      require(fields)
      
      # highest nrei is the first acceptable position
      fish.loc.mat <- ord.nrei3d.data.sub[1, , drop = FALSE]
      
      # for every other row (in our matrix of rows with nreis above the
      #  threshold)..
      for (row.no in 2:nrow(ord.nrei3d.data.sub)) {
        
        candidate.point <- ord.nrei3d.data.sub[row.no, 2:4, drop = FALSE]
        
        dists <- rdist(fish.loc.mat[, 2:4, drop = FALSE], candidate.point)
        
        if (all(dists > fish.territory.radius)) {
        
          fish.loc.mat <- rbind(fish.loc.mat, ord.nrei3d.data.sub[row.no, ])
        
        }  # end of if all(dists > fish.territory.radius)....
        
      }  # end of for row.no in 2:nrow(ord.nrei3d.data.sub)...
      
    }  # end of if there are multiple spots...

  # all are below our nrei threshold limit...
  } else {

    # matrix of NAs with simulation name indicated 
    fish.loc.mat <- matrix(rep(NA, 5), nrow = 1)

  }

  # print(head(fish.loc.mat))
  # give proper column titles; convert to data frame 
  colnames(fish.loc.mat) <- c('idx', 'X', 'Y', 'Z', 'pred.nrei')
  fish.loc.mat <- data.frame(fish.loc.mat)
  fish.loc.mat$sim.name <- true.col.name
  # head(fish.loc.mat)
  
  # add prey conc info 
  fish.loc.mat$sim.drift <- as.numeric(
    sub('d', '', 
      sapply(strsplit(fish.loc.mat$sim.name, split = '_'), function(x) x[2])
    )
  )

  # add temperature info 
  fish.loc.mat$sim.temp <- as.numeric(
    sub('C', '',
      sub('t', '', 
        sapply(strsplit(fish.loc.mat$sim.name, split = '_'), function(x) x[3])
      )
    )
  )

  if(nrow(fish.loc.mat) >= 1 & !any(is.na(fish.loc.mat[1, 1:5]))) {

    cap.pred <- nrow(fish.loc.mat)

    results.df <- data.frame(
      'cap.pred' = cap.pred,
      'mod.area' = mod.area,
      'dens.pred_fpm2' = cap.pred / mod.area,
      'mod.volume' = mod.volume,
      'dens.pred_fpm3' = cap.pred / mod.volume
    )

  } else if (nrow(fish.loc.mat) == 1 & all(is.na(fish.loc.mat[1, 1:5]))) {

    cap.pred <- 0.0

    results.df <- data.frame(
      'cap.pred' = cap.pred,
      'mod.area' = mod.area,
      'dens.pred_fpm2' = cap.pred / mod.area,
      'mod.volume' = mod.volume,
      'dens.pred_fpm3' = cap.pred / mod.volume
    )

  } else {

    writeLines(
      c(
        'An unexpected error occurred.',
        'Possibility #1: fish.loc.mat may have at least one row AND NA values',
        '  in the first row.',
        'Possibility #2: fish.loc.mat may have zero rows',
        'Possibility #3: something else we did not expect'
      )
    )

  } 

  return(results.df) 
  
}  # end of GetCapDensPreds function

#---------------------------------------------------------------------------
# Function to create an input lookup table for the nrei model
#
# IMPORTANT:  Must be run in 32-bit R if using 32-bit Access
# 
# Args:
#  requested.sites.table.fn: a .csv formatted table where each row is a 
#    visit that has been requested for the nrei model. The columns WatershedName,
#    SiteName, VisitYear, and VisitID (from the standard cm.org exports) must
#    exist
#  zip.in.data.dir: the directory of of the unzipped input data
#  unzip.in.data.dir: the directory of the unzipped input data 
#  champ.db.dir: dirctory containing the downloaded cm.org output dbs 
#
# Returns:
#  a list containing two data frames; the first is an input table to be used
#    for simulations; the second is a table of sites for which the requisite
#    data do not exist 
#---------------------------------------------------------------------------

CreateInLookTab <- function(requested.sites.table.fn, zip.in.data.dir,
  unzip.in.data.dir, champ.db.dir) {

  # packages needed
  my.packs <- c(
    'RODBC'
  )

  # if any of them are not installed, install them
  if (any(!my.packs %in% installed.packages()[, 'Package'])) {
    install.packages(
      my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
      dependencies = TRUE
    )
  }

  # read in the requested table 
  in.req.tab <- read.csv(requested.sites.table.fn, stringsAsFactors = FALSE)
  head(in.req.tab)

  # blank table to fill in with sites to be run 
  in.look.tab <- in.req.tab[NULL, c('WatershedName', 'SiteName', 'VisitYear', 'VisitID')]

  # blank table to fill in with sites we can't do 
  out.reject.tab <- in.look.tab

  # for each row of the table of requested sites...
  # site.no <- 19  site.no <- 31  site.no <- site.no + 1
  for (site.no in 1:nrow(in.req.tab)) {

    # get basin, site, and year
    my.basin <- in.req.tab[site.no, 'WatershedName']
    my.site <- in.req.tab[site.no, 'SiteName']
    my.year <- in.req.tab[site.no, 'VisitYear']
    my.visit <- paste('VISIT_', in.req.tab[site.no, 'VisitID'], sep = '')

    print(paste('Working on site', site.no, 'of' , nrow(in.req.tab), my.basin,
      my.site, my.year, my.visit, sep = ' : '))

    my.row.hydro.fn <- file.path(
      unzip.in.data.dir,
      my.basin,
      my.site,
      my.year,
      my.visit,
      'dem_grid_results.csv'
    )

    # check if unzipped hydro results already exist...
    if (file.exists(my.row.hydro.fn)) {
      
      print('already unzipped')
      in.look.tab[nrow(in.look.tab) + 1, ] <- in.req.tab[site.no, c('WatershedName', 'SiteName', 'VisitYear', 'VisitID')]

    # if unzipped hydro results don't already exist...
    } else {

      # see if the zip exists and could be unzipped...
      my.row.hydro.zip.fn <- file.path(
        zip.in.data.dir,
        my.year,
        my.basin,
        my.site,
        my.visit,
        'Hydro',
        'HydroModelResults.zip'
      )

      # if the zipped data exist...
      if (file.exists(my.row.hydro.zip.fn)) {

        # if they also contain dem_grid_results...
        if ('dem_grid_results.csv' %in% unzip(zipfile = my.row.hydro.zip.fn,
          list = TRUE)$Name) {

          print('zip file exists, trying to unzip...')

          # location for unzipping 
          unzip.dir <- file.path(
            unzip.in.data.dir,
            my.basin,
            my.site,
            my.year,
            my.visit
          )

          # unzip the files 
          unzip(zipfile = my.row.hydro.zip.fn, overwrite = TRUE, exdir = unzip.dir)

          # add the row to the input lookup table 
          in.look.tab[nrow(in.look.tab) + 1, ] <- in.req.tab[site.no, c('WatershedName',
            'SiteName', 'VisitYear', 'VisitID')]

        # if the zipped results do NOT contain dem_grid_results...
        } else {

          print('zip exists, but dem_grid_results is not in it')

          # add it to the table of sites without hydro 
          out.reject.tab[nrow(out.reject.tab) + 1, ] <- in.req.tab[site.no, c('WatershedName',
            'SiteName', 'VisitYear', 'VisitID')]

        }
        
      } else {

        # something else has happened
        print ('either the zip file does not exist or something else has gone wrong')

        # add it to the table of sites without hydro 
        out.reject.tab[nrow(out.reject.tab) + 1, ] <- in.req.tab[site.no, c('WatershedName',
          'SiteName', 'VisitYear', 'VisitID')]

      }

    }

  }  # end of for site.no in 1:nrow....

  library(RODBC)

  # go get D50, then calc roughness, then get wet area, thlwg depth, gradient
  
  pmdb.con <- odbcConnectAccess(
    file.path(champ.db.dir, 'CHaMP_Program Metrics')
  )

  mvi <- sqlFetch(pmdb.con, 'MetricVisitInformation')
  names(mvi)
  close(pmdb.con)

  mvi$SiteName <- gsub(' ', '', mvi$SiteName)
  mvi$WatershedName <- gsub(' ', '', mvi$WatershedName)
  unique(mvi$WatershedName)

  ilt.m <- merge(
    in.look.tab,
    mvi[, c('VisitID', 'Grad', 'Area_Wet', 'DpthBf_Avg', 'SubD50')],
    by.x = 'VisitID',
    by.y = 'VisitID',
    all.x = TRUE,
    sort = FALSE
  )
  head(ilt.m)
  head(in.look.tab)
  dim(ilt.m)

  # estimate grid size 
  ilt.m$est.grid.size <- 1088.5348 + 
    (8.0982 * ilt.m$Area_Wet) + 
    (23316.3863 * ilt.m$DpthBf_Avg * 0.8)  # 0.8 is a correction to acct for swith from Thalweg to Bankfull

  # none qualified in JD, so test this later!!!!
  if (any(ilt.m$est.grid.size > 125000)) {
    
    ilt.m[which(ilt.m$est.grid.size > 125000), 'zreductionfactor'] <- 10
    
  }

  # estimate time to run model
  ilt.m$est.sim.duration.h <- (2E-15 * (ilt.m$est.grid.size ** 3)) - 
    (1E-10 * (ilt.m$est.grid.size ** 2)) +
    (2E-5 * ilt.m$est.grid.size) -
    0.0527

  if (any(is.na(ilt.m$SubD50))) {

    mean.d50s <- mvi %>%
      filter(SubD50 < 350) %>%
      group_by(WatershedName) %>%
      summarise(mean.basin.d50 = mean(SubD50), n = n())
    mean.d50s <- as.data.frame(mean.d50s)
    
    rows.to.fix <- which(is.na(ilt.m$SubD50))
    
    for (my.row in rows.to.fix) {
    
      my.basin <- ilt.m[my.row, 'WatershedName']
      replacement.d50 <- mean.d50s[mean.d50s$WatershedName == my.basin, 'mean.basin.d50']
      ilt.m[my.row, 'SubD50'] <- replacement.d50
      
    }
    # ilt.m[rows.to.fix, ]
  
  } 

  ilt.m$kroughness <- ilt.m$SubD50 / 8000  # experimentally determined

  # set grid reduc fac
  ilt.m[ilt.m$Area_Wet <= 350, 'grdreductionfactor'] <- 2
  ilt.m[ilt.m$Area_Wet > 350 & ilt.m$Area_Wet <= 3000, 'grdreductionfactor'] <- 3
  ilt.m[ilt.m$Area_Wet > 3000 & ilt.m$Area_Wet <= 5000, 'grdreductionfactor'] <- 4
  ilt.m[ilt.m$Area_Wet > 5000, 'grdreductionfactor'] <- 5

  # set DZ and zreductionfactor

  # Grad > 2
  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg <= 0.28, 'DZ'] <- 0.025
  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 2

  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 3

  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg > 0.55, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad > 2 & ilt.m$DpthBf_Avg > 0.55, 'zreductionfactor'] <- 5

  # Grad in range (0.65, 2]
  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg <= 0.28, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 3

  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 3

  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg > 0.55, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad > 0.65 & ilt.m$Grad <= 2 & ilt.m$DpthBf_Avg > 0.55, 'zreductionfactor'] <- 5

  # Grad <= 0.65
  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg <= 0.28, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 3

  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.28 & ilt.m$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 5

  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.55 & ilt.m$DpthBf_Avg <= 0.77, 'DZ'] <- 0.075
  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.55 & ilt.m$DpthBf_Avg <= 0.77, 'zreductionfactor'] <- 7

  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.77, 'DZ'] <- 0.075
  ilt.m[ilt.m$Grad <= 0.65 & ilt.m$DpthBf_Avg > 0.77, 'zreductionfactor'] <- 10

  ilt.m$kSpecies <- 'steelhead'
  ilt.m$kPreyLength <- 0.003


  return(
    list(
      in.look.tab = in.look.tab,
      ilt.m = ilt.m,
      out.reject.tab = out.reject.tab
    )
  )

}  # end of CreateInLookTab function


#---------------------------------------------------------------------------
# Function to create an input lookup table for the nrei model
#
# Args:
#  site.metrics: an List scraped from an XML
#
# Returns:
#  a list containing two data frames; the first is an input table to be used
#    for simulations; the second is a table of sites for which the requisite
#    data do not exist 
#---------------------------------------------------------------------------
CreateInLookTabFromXML <- function(visits.list) {
  
  # blank table to fill in with sites to be run 
  in.look.tab <- data.frame(character(0), character(0), character(0), character(0), character(0), character(0),
                            numeric(0), numeric(0), numeric(0), numeric(0),
                            stringsAsFactors=FALSE)

  names(in.look.tab) <- c('WatershedName', 'SiteName', 'VisitYear', 'VisitID', 'OutDir', 'HydroResults',
                          'Grad', 'Area_Wet', 'DpthBf_Avg', 'SubD50')
  
  # blank table to fill in with sites we can't do 
  out.reject.tab <- in.look.tab
  
  # for each row of the table of requested sites...
  # site.no <- 19  site.no <- 31  site.no <- site.no + 1
  for (site.no in length(visits.list)) {
    
    my.visit <- as.character(paste('VISIT_', visits.list$Visit$Visit, sep = ''))
    
    in.look.tab[site.no,] <- data.frame(visits.list$Visit$Watershed, visits.list$Visit$Site, visits.list$Visit$Year, my.visit,
                   visits.list$Visit$OutputDir, visits.list$Visit$HydrolicsResults,
                   as.numeric(visits.list$Visit$Grad), 
                   as.numeric(visits.list$Visit$Area_Wet), 
                   as.numeric(visits.list$Visit$DpthBf_Avg), 
                   as.numeric(visits.list$Visit$SubD50), stringsAsFactors=FALSE)

    print(paste('Working on site', site.no, 'of' , length(visits.list), paste(in.look.tab[site.no,], collapse = " : "), sep = ' : '))
    
    # check if the CSV file exists
    if (file.exists(visits.list$Visit$HydrolicsResults)) {
      print('Hydrolics File Exists')
      in.look.tab[site.no, 'HydroResults'] <- file.path( visits.list$Visit$HydrolicsResults)
    } 
    
  }
  # estimate grid size 
  in.look.tab$est.grid.size <- 1088.5348 + 
    (8.0982 * in.look.tab$Area_Wet) + 
    (23316.3863 * in.look.tab$DpthBf_Avg * 0.8)  # 0.8 is a correction to acct for swith from Thalweg to Bankfull
  
  # none qualified in JD, so test this later!!!!
  if (any(in.look.tab$est.grid.size > 125000)) {
    
    in.look.tab[which(in.look.tab$est.grid.size > 125000), 'zreductionfactor'] <- 10
    
  }
  
  # estimate time to run model
  in.look.tab$est.sim.duration.h <- (2E-15 * (in.look.tab$est.grid.size ** 3)) - 
    (1E-10 * (in.look.tab$est.grid.size ** 2)) +
    (2E-5 * in.look.tab$est.grid.size) -
    0.0527
  
  library(dplyr)
  
  if (any(is.na(in.look.tab$SubD50))) {
    
    mean.d50s <- mvi %>%
      filter(SubD50 < 350) %>%
      group_by(WatershedName) %>%
      summarise(mean.basin.d50 = mean(SubD50), n = n())
    mean.d50s <- as.data.frame(mean.d50s)
    
    rows.to.fix <- which(is.na(in.look.tab$SubD50))
    
    for (my.row in rows.to.fix) {
      
      my.basin <- in.look.tab[my.row, 'WatershedName']
      replacement.d50 <- mean.d50s[mean.d50s$WatershedName == my.basin, 'mean.basin.d50']
      in.look.tab[my.row, 'SubD50'] <- replacement.d50
      
    }
    # in.look.tab[rows.to.fix, ]
    
  } 
  in.look.tab$kroughness <- in.look.tab$SubD50 / 8000  # experimentally determined
  
  # set grid reduc fac
  in.look.tab[in.look.tab$Area_Wet <= 350, 'grdreductionfactor'] <- 2
  in.look.tab[in.look.tab$Area_Wet > 350 & in.look.tab$Area_Wet <= 3000, 'grdreductionfactor'] <- 3
  in.look.tab[in.look.tab$Area_Wet > 3000 & in.look.tab$Area_Wet <= 5000, 'grdreductionfactor'] <- 4
  in.look.tab[in.look.tab$Area_Wet > 5000, 'grdreductionfactor'] <- 5
  
  # set DZ and zreductionfactor
  
  # Grad > 2
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg <= 0.28, 'DZ'] <- 0.025
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 2
  
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 3
  
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg > 0.55, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad > 2 & in.look.tab$DpthBf_Avg > 0.55, 'zreductionfactor'] <- 5
  
  # Grad in range (0.65, 2]
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg <= 0.28, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 3
  
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 3
  
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg > 0.55, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad > 0.65 & in.look.tab$Grad <= 2 & in.look.tab$DpthBf_Avg > 0.55, 'zreductionfactor'] <- 5
  
  # Grad <= 0.65
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg <= 0.28, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg <= 0.28, 'zreductionfactor'] <- 3
  
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'DZ'] <- 0.05
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.28 & in.look.tab$DpthBf_Avg <= 0.55, 'zreductionfactor'] <- 5
  
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.55 & in.look.tab$DpthBf_Avg <= 0.77, 'DZ'] <- 0.075
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.55 & in.look.tab$DpthBf_Avg <= 0.77, 'zreductionfactor'] <- 7
  
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.77, 'DZ'] <- 0.075
  in.look.tab[in.look.tab$Grad <= 0.65 & in.look.tab$DpthBf_Avg > 0.77, 'zreductionfactor'] <- 10
  
  in.look.tab$kSpecies <- 'steelhead'
  in.look.tab$kPreyLength <- 0.003
  
  
  return(
    list(
      in.look.tab = in.look.tab,
      out.reject.tab = out.reject.tab
    )
  )
  
}  # end of CreateInLookTab function
