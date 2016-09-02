
# rm(list = ls())


#---------------------------------------------------------------------------
# packages, dir and file locations, output preferences, table of inputs
#---------------------------------------------------------------------------

# Rtools will need to be installed for zipping of outputs.  Rtools can be
#  downloaded and installed manually from:
#  https://cran.r-project.org/bin/windows/Rtools/

# packages needed
my.packs <- c(
  'ggplot2', 'RANN', 'foreach', 'doParallel', 'scales', 'car', 'rgl', 
  'fields', 'data.table', 'dplyr', 'RODBC'
)

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {
  install.packages(
    my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],
    dependencies = TRUE
  )
}

# establish locations/directories
orig.wd <- getwd()  # initial working directory for R

inputs.dir <- 'C:/data_in/champ/unzipData'
outputs.dir <- 'C:/data_out/nrei'
nrei.func.fn <- 'C:/code/nrei/nreiFunctionsMultiFish_26May2016.R'
# inputs.dir <- 'C:/Users/ew/BoxSync/data_in/champ/unzipData'
# outputs.dir <- 'C:/Users/ew/BoxSync/data_out/nreiMultiFish'
# nrei.func.fn <- 'C:/Users/ew/BoxSync/code/nreiRectilinear/rectNreiMultiFish/nreiFunctionsMultiFish_26May2016.R'

# output preferences

# plots of radial grids (mostly for QAQC on radial grid function)
rad.grid.plots.bool <- FALSE

# a value of TRUE saves nrei maps (.png format) in the directory:
#  outputs.dir/basin/site/year/visit/species/fish_size/nrei_maps
nrei.map.plots.bool <- TRUE  # do you want to save nrei maps

# a value of TRUE saves capacity and density predictions in the directories:
#  outputs.dir/basin/site/year/visit/species/fish_size/fish_capacity and
#  outputs.dir/basin/site/year/visit/species/fish_size/fish_density;
#  necessary if capacity predictions are desired
fish.cap.dens.bool <- TRUE

# input lookup table on eric's computer 
in.look.tab <- read.csv(
  'C:/simulations/firstNreiSim_26May2016/exampleInlookTab.csv',
  stringsAsFactors = FALSE)
# head(in.look.tab)

#---------------------------------------------------------------------------
# begin looping through sites in in.look.tab
#---------------------------------------------------------------------------

site.no = 1  # site.no = 4
for (site.no in 1:nrow(in.look.tab)) {
  

  #-------------------------------------------------------
  # cleanup/setup
  #-------------------------------------------------------
  
  # make sure memory is as free as possible; seems necessary for fread function
  #  in the data.table package when working with 32-bit machines
  keep <- c('in.look.tab', 'site.no', 'orig.wd', 'inputs.dir', 'outputs.dir',
    'nrei.func.fn', 'rad.grid.plots.bool', 'nrei.map.plots.bool',
    'fish.cap.dens.bool'
  )
  rm(list = ls()[which(!ls() %in% keep)])
  
  # start site timer
  start.time <- proc.time()
  
  # source functions
  source(nrei.func.fn)

  # seems necessary to write UTMs to file with sufficient decimal places
  options(digits = 10)
  
  # processor cores to be used; can be as many as the computer has
  sim.cores <- 6
  

  #-------------------------------------------------------
  # temperature simulation ranges
  #-------------------------------------------------------

  # temperature where respiration is the highest; this will be RTO in the 
  #  swim costs equation
  respiration.temp.optimum <- 22

  # lethal temperature; this will be RTM in the swim costs equation
  respiration.temp.lethal <- 26

  # min temperature to be simulated
  min.sim.temp <- 8

  # step size for simulated temperature values
  temp.step.size <- 2  # in deg C

  # temperatures to simulate up and including to RTO
  temps.to.sim.p1 <- seq(min.sim.temp, respiration.temp.optimum, temp.step.size)

  # temperatures to simulate after RTO
  temps.to.sim.p2 <- c(23, 24, 24.2, 24.4, 24.6, 24.8, 25, 25.5)
  # temps.to.sim.p2 <- c(
  #   seq(respiration.temp.optimum + 1, respiration.temp.lethal - 1, 1),
  #   respiration.temp.lethal - 0.5
  # )

  # vector of temperatures to be simulated; all values in temps.to.sim must be
  #  smaller than respiration.temp.lethal!!!
  temps.to.sim <- sort(c(temps.to.sim.p1, temps.to.sim.p2))  # in deg C
  
  # # if you need a smaller range
  # temps.to.sim <- c(8, 14, 22, 25)
  

  #-------------------------------------------------------
  # drift simulation ranges
  #-------------------------------------------------------

  # lower and upper ranges of drift; smaller increment in range where most sites
  #  will be 
  drifts.to.sim.p1 <- c(0.01, 0.05, 0.1, 0.15, seq(0.25, 2, 0.25), 2.5)  # range where most sites will be; small increment
  drifts.to.sim.p2 <- c(seq(3, 5, 1), 7.5, 9.9)     # over 2 individuals/m3 is rare; bigger increment

  # vector of drift values to be simulated
  drifts.to.sim <- sort(c(drifts.to.sim.p1, drifts.to.sim.p2))

  # # if you need a smaller range
  # drifts.to.sim <- seq(0.5, 2, 0.25)
  

  #-------------------------------------------------------
  # fish lengths and weights
  #-------------------------------------------------------
  
  # fish lengths to be simulated
  fls.to.sim <- seq(0.070, 0.160, 0.010)  # in m; standard CHaMP values 
  # fls.to.sim <- 95  # single value 
  # fls.to.sim <- seq(0.080, 0.120, 0.01)  # shortened set for testing

  # fish weights; using length-weight function; in g
  # kFishWeight <- round(2 * ((10 ** -5) * (fls.to.sim * 1000) ** 2.848), 2)  # GT logan cutthroat
  kFishWeight <- round((10 ** -4.703) * (fls.to.sim * 1000) ** 2.866, 2)  # generic sthd

  
  #-------------------------------------------------------
  # output data structures
  #-------------------------------------------------------

  # list of lists to hold all results; each fish legnth gets a list 
  nrei.results <- vector('list', length(fls.to.sim))

  # fish length names 
  fl.identifiers <- paste(
    format(round(fls.to.sim * 1000), trim = TRUE, nsmall = 0), 'mm', sep = ''
  )

  # fish weight names
  fw.identifiers <- paste(
    format(round(kFishWeight), trim = TRUE, nsmall = 0), 'g', sep = ''
  )

  # output folders
  fish.size.folder.name <- paste(fl.identifiers, fw.identifiers, sep = '_')

  # name output lists 
  names(nrei.results) <- paste('fls.to.sim[', 1:length(fls.to.sim), ']', sep = '')

  # grei col names for grei prior to limiting to cmax
  grei.unlim.col.names <- paste('grei_unlim_d',
    format(drifts.to.sim, trim = TRUE, nsmall = 2), sep = '')

  # swim costs col names
  swim.cost.col.names <- paste('sc_t', format(temps.to.sim, trim = TRUE,
    nsmall = 1), 'C_jh', sep = '')

  # c.max col names
  c.max.col.names <- sub('sc_', 'c.max_', swim.cost.col.names)
  
  # get all pairwise combos of grei colnames and c.max colnames; grei will
  #  be limited by cmax, which is based on temperature 
  grei.cmax.combos <- expand.grid(c.max.col.names, grei.unlim.col.names,
    stringsAsFactors = FALSE
  )
  
  # grei col names for grei after limiting under each of the c.max/temp scenarios
  grei.lim.col.names <- apply(
    grei.cmax.combos,
    1,
    function(x) {
      paste(
        sub('un', '', x[2]),
        sub('c.max', '', x[1]),
        sep = ''
      )
    }
  )

  # clean up
  rm(min.sim.temp, temp.step.size, temps.to.sim.p1, temps.to.sim.p2,
    drifts.to.sim.p1, drifts.to.sim.p2)


  #-------------------------------------------------------
  # foraging model and radial grid parameters
  #-------------------------------------------------------

  # boolean variable to state whether or not grei predictions should be limited
  #  to the theoretical maximal consumption limit predicted by FishBioE 3.0
  #  DO NOT CHANGE THIS unless you have a specific reason for doing so
  limit.grei.to.cmax <- TRUE

  if (limit.grei.to.cmax) {
    # set time period over which satiation is reached
    kFeedingPeriod <- 24
  }

  # minimum allowable depth for a foraging location
  min.forage.depth <- 0.025  # in m

  # number of angles in the radial grid
  # pts <- 16
  pts <- 36  # scalar

  # max distance along radials (should be >= rd)
  max.Dist <- 2.0  # in m

  # distance between points along the radials in the radial grid
  Dist <- 0.01  # in m
  # Dist <- 0.05


  #-------------------------------------------------------
  # plotting and output control
  #-------------------------------------------------------

  # if radial grid plots are desired...
  if (rad.grid.plots.bool) {
    
    # define the interval at which they will occur (e.g. every 10,000 pts)
    plot.every <- 500

    # set up the color ramp to be used for radial grid plots 
    colorp <- colorRampPalette(c("dark blue", "blue", "cyan", "green", "yellow",
    "orange", "red", "dark red"))(51)

  }


  #-------------------------------------------------------
  # basin/site/year/hydro fn and site-specific inputs 
  #-------------------------------------------------------
  
  # get basin, site, and year
  my.basin <- in.look.tab[site.no, 'WatershedName']
  my.site <- in.look.tab[site.no, 'SiteName']
  my.year <- in.look.tab[site.no, 'VisitYear']
  my.visit <- paste('VISIT_', in.look.tab[site.no, 'VisitID'], sep = '')
  my.species <- in.look.tab[site.no, 'kSpecies']
  
  # filename of rectilinear Delft3D results for which you wish to simulate nrei
  DEMfilename <- file.path(inputs.dir, my.basin, my.site, my.year, my.visit,
    'dem_grid_results.csv')
  # file.exists(DEMfilename)

  # set the length of the average prey item
  kPreyLength <- in.look.tab[site.no, 'kPreyLength']  # in m
  
  # roughness height
  # k.roughness <- 0.01  # in m
  k.roughness <- in.look.tab[site.no, 'kroughness']


  #-------------------------------------------------------
  # other inputs to generate 
  #-------------------------------------------------------
  
  # get c.max
  if (limit.grei.to.cmax) {
    
    # start.time <- proc.time()
    for (fl in 1:length(fls.to.sim)) {
      
      # get c.max for each fish size and put in the appropriate list
      nrei.results[[fl]]$c.max.jh <- GetCMaxJh_params(
        in.temps = temps.to.sim,
        in.weight = kFishWeight[fl]
      )
      
      # name the matrices
      dimnames(nrei.results[[fl]]$c.max.jh) <- list('J/h', c.max.col.names)
      
    }
    
  }
  
  # get the fish's reactive distance and swim speed 
  for (fl in 1:length(fls.to.sim)) {

    # rd varies with fish and prey length; Hughes/Kelly have the 1000
    #  in their code as well 
    nrei.results[[fl]]$rd <- 
      0.12 * (kPreyLength * 1000) * (1 - exp(-20 * fls.to.sim[fl]))   # in m
    
    # get the fish's max swim speed; v.max varies with fish length
    nrei.results[[fl]]$v.max <- 0.87 * (fls.to.sim[fl] ** 0.19) # in m/s
    
  }
    
  # estimate energy content of prey as in as in Hayes/Hughes 2007

  # need prey length in mm, so mult by 1000 to get mg dry mass (Smock (1980))
  prey.dry.mass.hh <- exp(-5.021 + (2.88 * log(1000 * kPreyLength)))  # in mg

  # energy ~ 28.3 J for every mg dry mass (from Cummins and Wuycheck (1971))
  prey.energy.hh <- 28.3 * prey.dry.mass.hh  # in J
  
  assim.prey.energy.hh <- 0.7 * prey.energy.hh  # assimilable energy in J
  
  
  #---------------------------------------------------------------------------
  # nrei model
  #---------------------------------------------------------------------------
  
  # factor by which to increase X and Y grid spacing in order to decrease
  #  processing time; e.g., grd.reduction.factor = 2 would skip every other
  #  DEM grid cell in both X and Y directions, thus decreasing the number of 
  #  points by a factor of 2 ^ 2 = 4
  grd.reduction.factor <- in.look.tab[site.no, 'grdreductionfactor']  # scalar

  # grid spacing in z-direction for 3D (x, y, z) grid; may need to be less than
  #  x-y spacing for very shallow areas; may need to be greater for large, deep
  #  streams
  DZ <- in.look.tab[site.no, 'DZ']  # in m

  # build a 3D Grid of all the wetted points at the site using user-defined inputs
  GridBuild <- Build3DGrid_par(
    DEMfilename = DEMfilename,
    DZ = DZ,
    grd.reduction.factor = grd.reduction.factor,
    k.roughness = k.roughness,
    sim.cores = sim.cores
  )

  # how many had neg velocity magnitude?
  neg.vel.mags <- sum(GridBuild$Grid3D$Velocity.Magnitude < 0)

  # correct negative velocities
  if (neg.vel.mags > 0) {

    print('adjusting 2.5D conversion.... this may take a minute...')
    
    # new temp roughness = 25% of old value 
    k.rough.temp <- k.roughness * 0.25
    
    my.counter <- 1
    while (neg.vel.mags > 0 & my.counter < 50) {

      # new GridBuild with reduced roughness to eliminate neg vel mags 
      GridBuild <- Build3DGrid_par(
        
        DEMfilename = DEMfilename,
        DZ = DZ,
        grd.reduction.factor = grd.reduction.factor,
        k.roughness = k.rough.temp,
        sim.cores = sim.cores

      )

      neg.vel.mags <- sum(GridBuild$Grid3D$Velocity.Magnitude < 0)
      k.rough.temp <- k.rough.temp * 0.25
      my.counter <- my.counter + 1

    }  # end of while (neg.vel.mags > 0 & my.counter < 50)...

  }  # end of if (neg.vel.mags > 0)

  # avoid passing the 3d grid back and forth all the time
  Grid3D <- GridBuild$Grid3D
  my.area <- GridBuild$my.area
  my.volume <- GridBuild$my.volume
  # nrow(Grid3D)
  
  # # old QA/QC from MN code
  # boxplot(Grid3D$Velocity.Magnitude)
  # summary(Grid3D$Velocity.Magnitude)
  # sum(Grid3D$Velocity.Magnitude < 0)

  
  # in cases of high gradient or shallow sites, we may need to set DZ to
  #  small values to get good site coverage with Grid3D.  However, we don't
  #  necessarily need an nrei estimate for every 0.025 vertical meters for each
  #  unique x-y location.  Let's take care of this using the zrank column of
  #  Grid3D.
  z.reduction.factor <- in.look.tab[site.no, 'zreductionfactor']
  idxs.for.nrei.calcs <- which(Grid3D$zrank %in% seq(1, max(Grid3D$zrank), z.reduction.factor))
  # length(idxs.for.nrei.calcs)

  # initialize variables we will use over and over
  Init <- Initialize.Radial.Grid(pts = pts, max.Dist = max.Dist, Dist = Dist)


  #-------------------------------------------------------
  # 'wall of water', grei ident, sc for each 3D point
  #-------------------------------------------------------

  # register parallel backend
  library(foreach)
  library(doParallel)
  cl <- makeCluster(sim.cores)
  registerDoParallel(cl)
  
  # process foraging locations in parallel
  # all.res <- foreach(idx = idxs.for.nrei.calcs[50:60], .packages = 'RANN', .combine = 'rbind') %dopar% {
  all.res <- foreach(idx = idxs.for.nrei.calcs, .packages = 'RANN', .combine = 'rbind') %dopar% {

    print(paste(my.basin, my.site, my.year, my.visit, my.species,
      '--calculating grei and swim costs', sep = ' : '))
    
    # library(RANN)
    # idx = 1000
    # idx = idx + 100
    
    # get radial grid, dist-weighted average velocities, etc
    rad.grid <- radial.grid(
      idx,                                  #Index on Grid3D to define X-Y-Z location
      pts = pts,                            #Must match pts passed to Init
      max.Dist = max.Dist,                  #Must match max.Dist passed to Init
      Dist = Dist,                          #Must match Dist passed to Init
      plots = rad.grid.plots.bool,          #User Defined
      rad.grid.r = Init$rad.grid.r,         #Don't Change This
      rad.grid.theta = Init$rad.grid.theta, #Don't Change This
      Rmax = Init$Rmax,                     #Don't Change This
      VMult = Init$VMult,                   #Don't Change This
      DX_DY = GridBuild$DX_DY               #Don't Change This
    )
    
    fish.success.mat <- sapply(
      1:length(fls.to.sim),
      function(x) {
        GetFishSuccessMatrix_params(
          in.v.max = nrei.results[[x]]$v.max,
          in.rd = nrei.results[[x]]$rd
        )
      },
      simplify = 'array'
    )
    # class(fish.success.mat)
    # dim(fish.success.mat)  # 3rd dim should match length of fls.to.sim
    
    # # how many successes for each fish?
    # sapply(
    #   1:(dim(fish.success.mat)[3]),
    #   function(x) {sum(fish.success.mat[, , x])}
    # )
        
    # get sum of discharges for the 'tree ring sections';  in m3/s 
    cap.area.discharge <- sapply(
      1:length(fls.to.sim),
      function(x) {
        sum(fish.success.mat[, , x] * Init$rad.grid.areas * rad.grid$rad.grid.Vel.Norm)
      }
    )
    
    # capture area goes into a separate output for other analyses; in m2 
    cap.area <- sapply(
      1:length(fls.to.sim),
      function(x) {
        sum(fish.success.mat[, , x] * Init$rad.grid.areas)
      }
    )
    
    # 'identity' drift scenario using drift = 1.0 individuals/m3
    grei.unlim.ident.jh <- cap.area.discharge * 1.0 * assim.prey.energy.hh * 3600
    
    # put it all together
    c(
      idx,
      as.numeric(Grid3D[idx, ]),  # Velocity.Magnitude is the focal velocity
      cap.area,
      cap.area.discharge,
      grei.unlim.ident.jh
    )
    
  }  # end of 'all.res <- foreach(...'
  
  # shut down the cluster
  stopCluster(cl)
  
  # convert results to df and give columns names
  all.res <- data.frame(all.res)
  names(all.res) <- c(
    'idx',
    names(Grid3D),
    paste('cap.area.m2', fl.identifiers, sep = '.'),
    paste('cap.area.discharge.cms', fl.identifiers, sep = '.'),
    paste('grei.unlim.ident.jh', fl.identifiers, sep = '.')
  )
  

  #-------------------------------------------------------
  # nrei calculations
  #-------------------------------------------------------
  
  print(paste(my.basin, my.site, my.year, my.visit, my.species,
    '-- calculating nrei', sep = ' : '))
    
  # for each fish length...
  # fl <- 1
  for (fl in 1:length(fls.to.sim)) {
    
    # grab the unlimited grei 'identity' column
    proper.col <- intersect(
      grep('grei.unlim.ident.jh', names(all.res), value = TRUE),
      grep(
        paste(
          'jh.',
          fl.identifiers[fl],
          sep = ''
        ),
        names(all.res),
        value = TRUE
      )
    )
    nrei.results[[fl]]$grei.unlim.ident.jh <- all.res[, proper.col]
    
    # calculate unlimited grei for simulated drift concentrations and name the
    #  column
    nrei.results[[fl]]$grei.simd.unlim.jh <- sapply(
      drifts.to.sim,
      function(x) {
        x * nrei.results[[fl]]$grei.unlim.ident.jh
      }
    )
    colnames(nrei.results[[fl]]$grei.simd.unlim.jh) <- grei.unlim.col.names
    # class(nrei.results[[1]]$grei.simd.unlim.jh)
    
    # blah <- grei.cmax.combos[1, ]
    # blah
    # grei.cmax.combo <- as.character(blah)
    # grei.cmax.combo
    # grei.col <- nrei.results[[fl]]$grei.simd.unlim.jh[, grei.cmax.combo[2], drop = FALSE]
        
    # calculate the cmax-limited grei values and name the columns
    nrei.results[[fl]]$grei.simd.lim.jh <- apply(
      grei.cmax.combos,
      1,
      function(x) {
        grei.cmax.combo <- as.character(x)
        grei.col <- nrei.results[[fl]]$grei.simd.unlim.jh[, grei.cmax.combo[2], drop = FALSE]
        c.max.for.this.temp <- nrei.results[[fl]]$c.max.jh[, grei.cmax.combo[1]]
        grei.col.lim <- pmin(grei.col, c.max.for.this.temp)
      }
    )
    colnames(nrei.results[[fl]]$grei.simd.lim.jh) <- grei.lim.col.names
    
    # calculate the temp- and weight-specific swim costs for each foraging
    #  location and name the resulting swim cost columns
    nrei.results[[fl]]$swim.costs.simd.jh <- sapply(
      temps.to.sim,
      function(x) {
        GetSwimCostsHhSimd_params(
          in.rto = respiration.temp.optimum,
          in.rtl = respiration.temp.lethal,
          in.temp = x, 
          in.focal.vels = all.res$Velocity.Magnitude,
          in.weight = kFishWeight[fl]
        )
      }
    )
    colnames(nrei.results[[fl]]$swim.costs.simd.jh) <- swim.cost.col.names
    
    # calculate nrei for each combo of cmax-limited grei and temp-dependent
    #  swim costs
    nrei.results[[fl]]$nrei.simd.jh <- sapply(
      grei.lim.col.names,
      function(x) {
        temp.for.this.grei.lim <- strsplit(x, split = '_')[[1]][4]
        sc.col.name <- grep(temp.for.this.grei.lim, swim.cost.col.names, value = TRUE)
        temp.nrei <- 
          nrei.results[[fl]]$grei.simd.lim.jh[, x] - 
          nrei.results[[fl]]$swim.costs.simd.jh[, sc.col.name]
      }
    )
    colnames(nrei.results[[fl]]$nrei.simd.jh) <- sub('grei_lim', 'nrei', grei.lim.col.names)
    
    nrei.results[[fl]]$nrei.simd.js <- nrei.results[[fl]]$nrei.simd.jh / 3600
    colnames(nrei.results[[fl]]$nrei.simd.js) <- sub('_jh', '_js', colnames(nrei.results[[fl]]$nrei.simd.jh))
    
  }
  
  
  #---------------------------------------------------------------------------
  # write outputs to file
  #---------------------------------------------------------------------------
  
  print(paste(my.basin, my.site, my.year, my.visit, my.species,
    '-- writing raw outputs to file', sep = ' : '))
  
  # create output dirs for the fish sizes 
  fl.output.dirs <- file.path(outputs.dir, my.basin, my.site, my.year,
    my.visit, my.species, fish.size.folder.name)

  if (any(!file.exists(fl.output.dirs))) {
    dirs.to.create <- fl.output.dirs[which(!file.exists(fl.output.dirs))]
    for (my.dir in dirs.to.create) {
      dir.create(my.dir, recursive = TRUE)
    }
  }
  
  # start writing outputs 
  # fl <- 1
  for (fl in 1:length(fls.to.sim)) {

    # create folder for raw outputs 
    raw.output.dir <- file.path(fl.output.dirs[fl], 'raw_outputs')
    if (!file.exists(raw.output.dir)) {
      dir.create(raw.output.dir, recursive = TRUE)
    }

    # Grid3D
    write.csv(
      format(
        cbind(data.frame('idx' = 1:nrow(Grid3D)), Grid3D), nsmall = 5, 
        scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'grid3d.csv'),
      row.names = FALSE
    )
  
    # subset of Grid3D that we used 
    write.csv(
      format(all.res[, c('idx', names(Grid3D))], nsmall = 5, scientific = FALSE),
      file = file.path(raw.output.dir, 'grid3d_sub_w_results.csv'),
      row.names = FALSE
    )
    # all.res <- read.csv(file.path(output.dir, 'grid3d_sub_w_results.csv'), stringsAsFactors = FALSE)
  
    # cap area and cap area discharge
    cap.area.cols <- grep('cap.area', names(all.res), value = TRUE)
    this.fl.cols <- grep(fl.identifiers[fl], names(all.res), value = TRUE)
    cols.we.want <- intersect(cap.area.cols, this.fl.cols)
    write.csv(
      format(
        all.res[, c('idx', 'X', 'Y', 'Z', 'zrank', cols.we.want)],
        nsmall = 5,
        scientific = FALSE
      ), 
      file = file.path(raw.output.dir, 'cap_area_desc.csv'),
      row.names = FALSE
    )
  
    # grei identity
    grei.ident.cols <- grep('grei.unlim.ident.jh', names(all.res), value = TRUE)
    this.fl.cols <- grep(fl.identifiers[fl], names(all.res), value = TRUE)
    col.we.want <- intersect(grei.ident.cols, this.fl.cols)
    write.csv(
      format(
        all.res[, c('idx', 'X', 'Y', 'Z', 'zrank', col.we.want)],
        nsmall = 5,
        scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'grei_ident_jh.csv'),
      row.names = FALSE
    )
  
    # swim costs
    write.csv(
      format(
        cbind(
          all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')],
          nrei.results[[fl]]$swim.costs.simd.jh
        ),
        nsmall = 5,
        scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'sc_jh.csv'),
      row.names = FALSE
    )
  
    # unlimited grei
    write.csv(
      format(
        cbind(
          all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')],
          nrei.results[[fl]]$grei.simd.unlim.jh
        ),
      nsmall = 5,
      scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'grei_unlim_jh.csv'),
      row.names = FALSE
    )
  
    # limited grei
    write.csv(
      format(
        cbind(
          all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')],
          nrei.results[[fl]]$grei.simd.lim.jh
        ),
      nsmall = 5,
      scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'grei_lim_jh.csv'),
      row.names = FALSE
    )
  
    # nrei in J/h
    write.csv(
      format(
        cbind(
          all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')],
          nrei.results[[fl]]$nrei.simd.jh
        ),
      nsmall = 5,
      scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'nrei_jh.csv'),
      row.names = FALSE
    )
    
    # nrei in J/s
    write.csv(
      format(
        cbind(
          all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')],
          nrei.results[[fl]]$nrei.simd.js
        ),
      nsmall = 5,
      scientific = FALSE
      ),
      file = file.path(raw.output.dir, 'nrei_js.csv'),
      row.names = FALSE
    )
    
  }  # end of for (fl in 1:length(fls.to.sim))... output writing

  # clean up 

  # # see what the memory hogs are
  # mem.vec <- sapply(ls(), function(x) object.size(get(x)))
  # mem.df <- data.frame(
  #   'item' = names(mem.vec),
  #   'memory' = mem.vec
  # )
  # mem.df[order(-mem.df$memory), ]
  
  # remove biggest memory hogs
  for (fl in 1:length(fls.to.sim)) {
    nrei.results[[fl]]$grei.simd.lim.jh <- NULL
    nrei.results[[fl]]$nrei.simd.jh <- NULL
    nrei.results[[fl]]$grei.simd.unlim.jh <- NULL
  }

  rm(
    GridBuild,
    Grid3D,
    Init
  )
  
  
  #-------------------------------------------------------
  # create nrei maps 
  #-------------------------------------------------------

  if (nrei.map.plots.bool) {

    print(paste(my.basin, my.site, my.year, my.visit, my.species,
      '-- writing nrei maps to file', sep = ' : '))
    
    # fl <- 1
    for (fl in 1:length(fls.to.sim)) {
      
      # create maps folder 
      map.output.dir <- file.path(fl.output.dirs[fl], 'nrei_maps')
      if (!file.exists(map.output.dir)) {
        dir.create(map.output.dir, recursive = TRUE)
      }
      
      # grab the z point closest to the bed for each unique (x,y) so we aren't
      #  over-plotting

      # get (x, y, z) from all.res; to recreate this post-simulation,
      #  read in nrei_js.csv from raw_outputs and select the same columns
      low.z.pts <- all.res[, c('idx', 'X', 'Y', 'Z', 'zrank')]
      # summary(low.z.pts)
      
      # add on the nrei results
      low.z.pts <- cbind(low.z.pts, nrei.results[[fl]]$nrei.simd.js)

      library(dplyr)

      # make sure we just have one 3D point for every XY pair 
      low.z.pts <- filter(low.z.pts, zrank == 1)

      # establish good limits for the nrei shading
      # names(low.z.pts)
      nrei.col.names <- grep('nrei', names(low.z.pts), value = TRUE)
      max.nrei <- quantile(as.matrix(low.z.pts[, nrei.col.names]), 0.9975)
      min.nrei <- quantile(as.matrix(low.z.pts[, nrei.col.names]), 0.0025)
      
      # process in parallel 
      cl <- makeCluster(sim.cores)
      registerDoParallel(cl)
      
      # foreach(col.num = 100:120, .packages = c('ggplot2', 'scales')) %dopar% {
      foreach(col.num = 1:length(nrei.col.names), .packages = c('ggplot2', 'scales')) %dopar% {
        
        col.name <- nrei.col.names[col.num]
        my.drift <- strsplit(col.name, split = '_')[[1]][2]
        my.temp <- strsplit(col.name, split = '_')[[1]][3]
        my.plot <- ggplot(data = low.z.pts, aes_string(x = "X", y = "Y", color = col.name), alpha = 0.9) +
          geom_point(size = I(0.75)) +
          scale_colour_gradientn(
            # colours = c('red2', 'yellow2', 'green2'),
            colours = c('red2', 'orange', 'yellow2', 'greenyellow', 'green4'),
            limits = c(min.nrei, max.nrei),
            # values = rescale(c(min.nrei, 0.0, max.nrei), to = c(0, 1.0))
            values = rescale(c(min.nrei, min.nrei / 5, 0.0, max.nrei / 5, max.nrei), to = c(0, 1.0)),
            name = 'NREI (J/s)'
          ) +
          coord_equal() +
          ggtitle(paste(my.visit, fish.size.folder.name[fl], my.drift, my.temp, sep = '_')) +
          theme_minimal()
        
        temp.fn <- paste(col.name, '.png', sep = '')
        
        ggsave(
          filename = temp.fn,
          plot = my.plot,
          path = map.output.dir
        )
      }
      
      stopCluster(cl)

    }  # end of for each fish length

  }  # end of if (nrei.map.plots.bool)...


  #-------------------------------------------------------
  # estimate capacity and density 
  #-------------------------------------------------------

  if (fish.cap.dens.bool) {

    print(paste(my.basin, my.site, my.year, my.visit, my.species,
      '-- writing capacity (fish) and density (fish/m2) estimates to file',
      sep = ' : '))

    # fl <- 1
    for (fl in 1:length(fls.to.sim)) {


      #-------------------------------------------------------
      # output folders 
      #-------------------------------------------------------

      # fish_locations folder 
      fish.locs.output.dir <- file.path(fl.output.dirs[fl], 'fish_locations')
      if (!file.exists(fish.locs.output.dir)) {
        dir.create(fish.locs.output.dir, recursive = TRUE)
      }
      fish.loc.fn <- file.path(fish.locs.output.dir, 'all_sims_fish_locs.csv')

      # fish_capacity folder 
      fish.cap.output.dir <- file.path(fl.output.dirs[fl], 'fish_capacity')
      if (!file.exists(fish.cap.output.dir)) {
        dir.create(fish.cap.output.dir, recursive = TRUE)
      }
      fish.cap.fn <- file.path(fish.cap.output.dir, 'fish_cap_preds.csv')

      # fish_density folder 
      fish.dens.output.dir <- file.path(fl.output.dirs[fl], 'fish_density')
      if (!file.exists(fish.dens.output.dir)) {
        dir.create(fish.dens.output.dir, recursive = TRUE)
      }
      fish.dens.fn <- file.path(fish.dens.output.dir, 'fish_dens_preds.csv')


      #-------------------------------------------------------
      # fish locations for all fish sizes and all drift/temp
      #  pairs 
      #-------------------------------------------------------

      # based on Keeley and McPhail '98, adjusted based on Imre, Grant, Keeley '04
      #  must get fish length in cm for this equation; territory is in m2
      # fish.territory.area <- 10 ** (1.56 * log10(kFishLength * 100) - 1.81 - 0.07)
      fish.territory.area <- 10 ** (1.56 * log10(fls.to.sim[fl] * 100) - 1.81 - 0.07)

      # get fish territory radius, multiply by 2; Multiplying by two ensures
      #  fish territories don't overlap
      fish.territory.radius <- sqrt(fish.territory.area / pi) * 2

      # set the nrei limit and territory size
      nrei.limit <- 0.0  # in J/s
      
      # all.res <- read.csv(raw.out.fn, stringsAsFactors = FALSE)
      nrei.vals <- cbind(
        all.res[, c('idx', 'X', 'Y', 'Z')],
        nrei.results[[fl]]$nrei.simd.js
      )
      # names(nrei.vals)

      # process fish locations in parallel 
      cl <- makeCluster(sim.cores)
      registerDoParallel(cl)
      
      all.fish.locs <- foreach(col.num = 1:length(nrei.col.names), .combine = 'rbind') %dopar% {
      
        col.name <- nrei.col.names[col.num]
        my.drift <- as.numeric(sub('d', '', strsplit(col.name, split = '_')[[1]][2]))
        my.temp <- as.numeric(sub('C', '', sub('t', '', strsplit(col.name, split = '_')[[1]][3])))
        # print(col.name)
      
        # this time get all results instead of just the ones nearest the bed
        #  also, let's reuse some old code here
        nrei3d.data <- nrei.vals[, c('idx', 'X', 'Y', 'Z', col.name)]
        
        # rename so we can more easily refer to the nrei col; also to make
        #  compatible with legacy code
        names(nrei3d.data)[ncol(nrei3d.data)] <- 'nrei.hh.js'
        
        # if there are any above our limit
        if (any(nrei3d.data$nrei.hh.js >= nrei.limit)) {
        
          # order by decreasing nrei, keep those above the limit, use matrices
          ord.nrei3d.data <- nrei3d.data[order(-nrei3d.data$nrei.hh.js), ]
          ord.nrei3d.data <- ord.nrei3d.data[ord.nrei3d.data$nrei.hh.js >= nrei.limit, ]
          ord.nrei3d.data.sub <- as.matrix(ord.nrei3d.data)

          #  unname it, clean up
          ord.nrei3d.data.sub <- unname(ord.nrei3d.data.sub)
          rm(nrei3d.data, ord.nrei3d.data)
          # ord.nrei3d.data.sub[1:10, ]
          
          # get fish locs
          fish.locs <- GetFishLocMatrix()
          fish.locs <- cbind(fish.locs, rep(col.num, nrow(fish.locs)))
          
        } else {
          
          fish.locs <- c(rep(NA, 5), col.num)
          
        }
        
      }  # end of parallel loop for all.fish.locs 
      
      stopCluster(cl)
      rm(cl)
      
      colnames(all.fish.locs) <- c('idx', 'X', 'Y', 'Z', 'pred.nrei', 'sim.num')
      all.fish.locs <- data.frame(all.fish.locs, row.names = NULL)
      all.fish.locs$sim.name <- nrei.col.names[all.fish.locs$sim.num]
      all.fish.locs$sim.drift <- as.numeric(
        sub('d', '', 
          sapply(strsplit(all.fish.locs$sim.name, split = '_'), function(x) x[2])
        )
      )
      all.fish.locs$sim.temp <- as.numeric(
        sub('C', '',
          sub('t', '', 
            sapply(strsplit(all.fish.locs$sim.name, split = '_'), function(x) x[3])
          )
        )
      )
      
      write.csv(
        format(all.fish.locs, nsmall = 5),
        file = fish.loc.fn,
        row.names = FALSE
      )


      #-------------------------------------------------------
      # get fish capacities
      #-------------------------------------------------------

      library(car)
      library(rgl)

      # deal with scenarios where fish could make it; by default
      #  aggregate drops NAs (i.e. scenarios where no fish were predicted)
      fish.cap.preds <- aggregate(
        idx ~ sim.name + sim.drift + sim.temp, data = all.fish.locs, length
      )

      # now, let's grab all those NA values, 'set' the idx col equal to zero to
      #  indicate zero fish; aggregate does not compute on NA rows; we'll
      #  handle them here
      
      # deal with scenarios where fish could not make it; if any scenarios
      #  have an NA in the idx column, this should be the only row for
      #  that scenario and that scenario shouldn't have supported a fish 
      if(any(is.na(all.fish.locs$idx))) {
        
        # get rows of all.fish.locs with NA in idx column; also, format
        #  resulting table to match columns of fish.cap.preds
        fish.cap.zeros <- all.fish.locs[which(is.na(all.fish.locs$idx)), 
          c('sim.name', 'sim.drift', 'sim.temp', 'idx')]
        fish.cap.zeros$idx <- 0
        
        # add the NA rows back into fish.cap.preds
        fish.cap.preds <- rbind(fish.cap.preds, fish.cap.zeros)
        
      }

      # sort fish.cap.preds and rename idx column to 'pred.capacity'
      fish.cap.preds <- fish.cap.preds[order(fish.cap.preds$sim.drift,
        fish.cap.preds$sim.temp), ]
      names(fish.cap.preds)[ncol(fish.cap.preds)] <- 'pred.capacity'
      
      # write fish cap file
      write.csv(format(fish.cap.preds, nsmall = 2), file = fish.cap.fn,
        row.names = FALSE)
      
      # 3D plot: fish cap ~ drift + temp
      par3d(windowRect = c(100, 100, 800, 800))  # make default 3D window bigger
      scatter3d(
        pred.capacity ~ sim.drift + sim.temp,
        data = fish.cap.preds,
        surface = FALSE
      )
      writeWebGL(
        dir = fish.cap.output.dir,
        filename = file.path(fish.cap.output.dir, 'pred_capacities.html')
      )
      

      #-------------------------------------------------------
      # get fish densities 
      #-------------------------------------------------------

      # add Delft3D wetted area for this visit to fish.cap.preds table 
      fish.cap.preds$mod.area <- my.area
      fish.cap.preds$mod.volume <- my.volume

      # calculate densities 
      fish.cap.preds$pred.dens <- fish.cap.preds$pred.capacity / fish.cap.preds$mod.area
      fish.cap.preds$pred.vol.dens <- fish.cap.preds$pred.capacity / fish.cap.preds$mod.volume
      
      write.csv(format(fish.cap.preds, nsmall = 3), file = fish.dens.fn,
        row.names = FALSE)

      # 3D plot: fish dens ~ drift + temp
      scatter3d(
        pred.dens ~ sim.drift + sim.temp,
        data = fish.cap.preds,
        surface = FALSE
      )
      writeWebGL(
        dir = fish.dens.output.dir,
        filename = file.path(fish.dens.output.dir, 'pred_densities.html')
      )
  
    }  # end of for each fish length
  
  }  # end of if (fish.cap.dens.bool) statement


  #---------------------------------------------------------------------
  # write sim summary to file, zip raw_outputs folder 
  #---------------------------------------------------------------------
  
  # mark time of simulation/output end 
  end.time <- proc.time()

  # write sim info output to raw_output folder for each fish size simulated;
  #  this is unnecessarily repetitive (coud just write one for each simu),
  #  but I've done it this way for backwards compatibility with how we did it
  #  in the past 
  for (fl in 1:length(fls.to.sim)) {

    # create folder for raw outputs 
    raw.output.dir <- file.path(fl.output.dirs[fl], 'raw_outputs')
  
    sim.time.df <- data.frame(
      'basin' = my.basin,
      'site' = my.site,
      'year' = my.year,
      'visit.id' = my.visit,
      'species' = my.species,
      'grid.size' = nrow(nrei.vals),
      'start.time' = start.time[3],
      'end.time' = end.time[3],
      'sim.duration.s' = end.time[3] - start.time[3],
      'multifish' = 'YES',
      'num.fish.sizes' = length(fls.to.sim),
      'neg.vel.mags' = neg.vel.mags
    )
        
    write.csv(
      sim.time.df,
      file = file.path(raw.output.dir, 'sim_info.csv'),
      row.names = FALSE
    )

  }  # end of for (fl in 1:length(fls.to.sim))...


  #---------------------------------------------------------------------
  # zip raw_outputs folder 
  #---------------------------------------------------------------------

  for (fl in 1:length(fls.to.sim)) {

    # set wd to the fish.size.folder.name dir (seems necessary for 
    #  zip/delete)
    setwd(fl.output.dirs[fl])

    # if there is a 'raw_outputs' folder, zip it and delete the original
    if (dir.exists('raw_outputs')) {

      zip('raw_outputs', files = 'raw_outputs')
      unlink('raw_outputs', recursive = TRUE)

    }

  }  # end of for (fl in 1:length(fls.to.sim))...
  
  # set back to original working directory
  setwd(orig.wd)    

}


# #---------------------------------------------------------------------------
# # example of using functions in nreiFunctionsMultiFish_xx.R to create an
#  input lookup table or calculate nrei values, fish locations, and site
#  capacity/density for non whole number values
# #---------------------------------------------------------------------------


# #---------------------------------------------------------------------
# # create input lookup table from a table with WatershedName, SiteName,
# #  VisitYear, and VisitID
# #---------------------------------------------------------------------

# nrei.func.fn <- 'C:/Users/ew/BoxSync/code/nreiRectilinear/rectNreiMultiFish/nreiFunctionsMultiFish_26May2016.R'

# source(nrei.func.fn)

# trial <- CreateInLookTab(
#   requested.sites.table.fn <- 'C:/Users/ew/Desktop/exampleNreiFiles/simulations/26May2016/requestedNreiVisitList.csv',
#   zip.in.data.dir <- 'C:/Users/ew/BoxSync/data_in/champ/zipData',
#   unzip.in.data.dir <- 'C:/Users/ew/BoxSync/data_in/champ/unzipData',
#   champ.db.dir <- 'C:/Users/ew/BoxSync/data_in/champ'
# )
# names(trial)

# in.look.tab <- trial$ilt.m
# sum(in.look.tab$est.sim.duration.h)
# in.look.tab <- in.look.tab[order(in.look.tab$est.sim.duration.h), ]
# head(in.look.tab)

# write.csv(in.look.tab, file = 'exampleInlookTab.csv',
#   row.names = FALSE)


# #---------------------------------------------------------------------
# # calculate nrei metrics for non whole number fish sizes, drifts,
# #  or temperatures
# #---------------------------------------------------------------------

# # extract nrei values for non whole number values of fish size, drift, and temp 
# trial <- GetNreiValues(
#   fish.size.folders.dir = 'C:/Users/ew/BoxSync/data_out/trash3/JohnDay/CBW05583-240498/2012/VISIT_960/steelhead',
#   fishLength_m = 0.135,
#   preyConc_noPM3 = 0.65,
#   temp_C = 18.4
# )
# format(head(trial, 100), scientific = FALSE)

# # get predicted fish locations for non whole number values of fish size,
# #  drift, and temp 
# trial2 <- GetFishLocations(
#   fish.size.folders.dir = 'C:/Users/ew/BoxSync/data_out/trash3/JohnDay/CBW05583-240498/2012/VISIT_960/steelhead',
#   fishLength_m = 0.137,
#   preyConc_noPM3 = 0.65,
#   temp_C = 18.4,
#   nrei.limit_jps = 0.0
# )
# head(trial2)
# nrow(trial2)

# # get predicted capacity (# fish) and density (fish/m2 and fish/m3)
# trial3 <- GetCapDensPreds(
#   fish.size.folders.dir = 'C:/Users/ew/BoxSync/data_out/trash3/JohnDay/CBW05583-240498/2012/VISIT_960/steelhead',
#   fishLength_m = 0.137,
#   preyConc_noPM3 = 0.65,
#   temp_C = 18.4,
#   nrei.limit_jps = 0.0
# )
# trial3

