
# Converts the character vectors of student survey into a list of skill matrix 
arrange_skills = function(ncls, subjects, socials){
  # if any student missed their submission replace with comma separated values
  subjects[subjects==""] = "0,0,0,0";
  # if any student missed their submission replace with comma separated values
  socials[socials==""] = "0,0,0,0,0,0"
  subjects = matrix(as.integer(Reduce(rbind, strsplit(subjects, ','))), ncol=4)
  socials = matrix(as.integer(Reduce(rbind, strsplit(socials, ','))), ncol=6)
  colnames(subjects) = c('pro', 'math', 'alg', 'bio')
  colnames(socials) = c('leadership', 'hard_work', 'work_others', 'delivery', 'adaptive', 'through')
  
  subjects = as.data.frame(subjects[1:ncls,]); socials = as.data.frame(socials[1:ncls,])
  return(list(subjects=subjects, socials=socials))
}

# evalute the final group configuration
eval_grouping = function(config, skills, crit){
  subjects = skills[[crit]]
  grps = lapply(1:length(table(config)), function(x) subjects[which(config==x),])
  return(grps)
}

# initialize the group configuration
make_config = function(ncls, ngrp){
  initial = sample(ncls)%%ngrp+1
  names(initial) = 1:ncls; initial
  return(initial)
}

# evaluate the group-configuration for evenly skilled groups
eval_evenly_skilled_config = function(config, skills){
  obj = sapply(names(skills), function(crit){
    subjects = skills[[crit]]
    grps = lapply(1:length(table(config)), function(x) subjects[which(config==x),])
    obj_g = sapply(colnames(subjects), function(x) sapply(grps, function(y) max(y[[x]]) ))
    obj_g = apply(obj_g, 2, min); obj_g = min(obj_g)*sum(obj_g)
  })
  
  return(sum(obj))
}

# evaluate the group-configuration for maximally diverse groups
# take sum of square difference of each skills between any two groups
eval_maximally_diverse_config = function(config, skills, weight=c(1,1)){
  obj = sapply(names(skills), function(crit){
    subjects = skills[[crit]]
    grps = lapply(1:length(table(config)), function(x) subjects[which(config==x),])
    
    obj_g = lapply(grps, function(gi) sapply(rownames(gi), function(x) sapply(rownames(gi), function(y){
      sum(abs(gi[x,]-gi[y,]))
    })))
    sum(unlist(lapply(obj_g, function(x) sum(x[upper.tri(x)]))))
  })
  
  #given double weight on subjects than social skills
  return(sum(obj*weight))
}

# generate neighbor solutions
make_neighbours = function(config, tabu=NULL){
  ncls = length(config)
  if(is.null(tabu)){
    tabu = rep(F, ncls); names(tabu) = 1:ncls
  }
  configNeigh = Reduce(rbind, sapply(1:(ncls-1), function(x){
    if(tabu[x]==T) return(NULL)
    tmp = which(config[x]!=config[(x+1):ncls])+x; names(tmp) = NULL; 
    tmp = sapply(tmp, function(y){
      if(tabu[y]==T) return(NULL)
      replace(config, c(x, y), config[c(y, x)])
    }); tmp; class(tmp)
    if(class(tmp)=='list'){
      return(Reduce(rbind, tmp))
    }else if(class(tmp)=='matrix'){
      t(tmp)
    }else{
      tmp
    } 
  }))
  
  moves = Reduce(rbind, sapply(1:(ncls-1), function(x){
    if(tabu[x]==T) return(NULL)
    tmp = which(config[x]!=config[(x+1):ncls])+x; names(tmp) = NULL; 
    tmp = sapply(tmp, function(y){
      if(tabu[y]==T) return(NULL)
      c(x,y)
    }); tmp; class(tmp)
    if(class(tmp)=='list'){
      return(Reduce(rbind, tmp))
    }else if(class(tmp)=='matrix'){
      t(tmp)
    }else{
      tmp
    } 
  }))
  #print(nrow(moves))
  return(list(configNeigh=configNeigh, moves=moves))
}

# size = nclss; k = ngrp; 
# the following function is extension from binary tabu Search 
k_nary_tabu = function (ncls, grpsize, objFunc, skills, iters = 10, neigh = ncls, delta = 5,
                        nRestarts = 2, repeatAll = 1, verbose = T) {
  
  iter <- 1; tabu = rep(F, ncls); ngrp = ceiling(ncls/grpsize)
  tc = 1; tc_star = tc; tx = ty = rep(0, ncls); names(tx) = names(ty) = 1:ncls; 
  configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), ncls)
  eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3))
  for (j in 1:repeatAll) {
    config <- make_config(ncls, ngrp); 
    eUtility <- objFunc(config, skills); aspiration <- eUtility
    preliminarySearch <- function() {
      configKeep[tc, ] <- config; eUtilityKeep[tc] <- eUtility; 
      for (tc in 1:iters) { print(paste('tick',tc))
        configTemp = make_neighbours(config, NULL); 
        configNeigh = configTemp$configNeigh; moves = configTemp$moves
        neighboursEUtility <- matrix(0, 1, nrow(configNeigh)); neigh = nrow(configNeigh)
        randNeigh <- sample(nrow(configNeigh), neigh)
        neighboursEUtility[randNeigh] <- apply(configNeigh[randNeigh, ], 1, objFunc, skills) #each row is a chrom
        
        tabuList = sapply(1:nrow(moves), function(x) tabu[moves[x,1]]==T | tabu[moves[x,2]]==T)
        tabuList = as.numeric(tabuList); names(tabuList) = 1:nrow(moves)
        maxNontabu = ifelse(sum(tabuList)==nrow(moves), 0, max(neighboursEUtility[tabuList==0]))
        maxTabu <- ifelse(sum(tabuList)==0, 0, max(neighboursEUtility[tabuList==1]))
        move <- ifelse(maxTabu > maxNontabu & maxTabu > aspiration, 
                       ifelse(length(which(neighboursEUtility == maxTabu)) == 1, 
                              which(neighboursEUtility == maxTabu), 
                              sample(which(neighboursEUtility == maxTabu), 1)), 
                       ifelse(length(which(neighboursEUtility == maxNontabu & tabuList==0)) == 1, 
                              which(neighboursEUtility == maxNontabu & tabuList == 0), 
                              sample(which(neighboursEUtility == maxNontabu & tabuList==0), 1)))
        
        if (neighboursEUtility[move] > aspiration){
          aspiration <- neighboursEUtility[move]
        }
        if(neighboursEUtility[move] > eUtility){
          eUtility <- neighboursEUtility[move]; config = configNeigh[move,];
          tc_star = tc; tx[moves[move, 1]] = delta+tc; ty[moves[move, 2]] = delta+tc;
          print('updated')
          write.table(config, file = 'group_configuration.txt', col.names=F, row.names=F)
        }
        tabu = tx-tc > 0 | ty-tc > 0; tc = tc+1; 
        configKeep[tc, ] <- config; eUtilityKeep[tc] <- eUtility;
        if((tc-tc_star)>(iters/2)){
          break
        }
      }#end of iters
      result = list(aspiration = aspiration, configKeep = configKeep, 
                    eUtilityKeep = eUtilityKeep, iter = tc)
      return(result)
    }#end of primarySearch
    
    if (verbose) 
      cat("Preliminary search stage...\n")
    result <- preliminarySearch(); configKeep <- result$configKeep; head(configKeep)
    aspiration <- result$aspiration; eUtilityKeep <- result$eUtilityKeep; iter <- result$tc
    tempo <- 0; restarts <- 0
    while (tempo < aspiration & restarts < nRestarts) {
      print(paste('restart', restarts))
      if (verbose) 
        cat("Intensification stage...\n")
      eUtility <- max(eUtilityKeep); tempo <- aspiration
      config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ]
      write.table(config, file = 'group_configuration.txt', col.names=F, row.names=F)
      #config <- make_config(ncls, ngrp); #doing random restart
      result <- preliminarySearch()
      aspiration <- result$aspiration; configKeep <- result$configKeep
      eUtilityKeep <- result$eUtilityKeep; iter <- result$tc
      restarts <- restarts + 1
    }
    if (verbose) 
      cat("Diversification stage...\n")
    config <- make_config(ncls, ngrp)
    eUtility <- objFunc(config, skills)
    #frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0))
    #tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize))
    #listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize)
    result <- preliminarySearch()
    iter <- result$tc; configKeep <- result$configKeep; eUtilityKeep <- result$eUtilityKeep
  }
  
  return(config)
}


#taken from tabuSearch Package: this is binary version
binary_tabu = function (size = 10, iters = 100, objFunc = NULL, config = NULL, neigh = size, 
                        listSize = 9, nRestarts = 10, repeatAll = 1, verbose = FALSE) {
  if (size < 2) {
    stop("error: config too short!")
  }
  if (iters < 2) {
    stop("error: not enough iterations!")
  }
  if (listSize >= size) {
    stop("error: listSize too big!")
  }
  if (neigh > size) {
    stop("error: too many neighbours!")
  }
  if (is.null(objFunc)) {
    stop("A evaluation function must be provided. See the objFunc parameter.")
  }
  if (is.null(config)) {
    config <- matrix(0, 1, size)
    config[sample(1:size, sample(1:size, 1))] <- 1
  }
  else if (size != length(config)) {
    stop("Length of the starting configuration != size")
  }
  if (repeatAll < 1) {
    stop("error: repeatAll must be > 0")
  }
  iter <- 1
  configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size)
  eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3))
  for (j in 1:repeatAll) {
    if (j > 1) {
      config <- matrix(0, 1, size)
      config[sample(1:size, sample(1:size, 1))] <- 1
    }
    tabuList <- matrix(0, 1, size)
    listOrder <- matrix(0, 1, listSize)
    eUtility <- objFunc(config)
    aspiration <- eUtility
    preliminarySearch <- function() {
      configKeep[iter, ] <- config
      eUtilityKeep[iter] <- eUtility
      iter <- iter + 1
      for (i in 2:iters) {
        neighboursEUtility <- matrix(0, 1, size)
        configTemp <- t(matrix(config, size, neigh)) #each row is a chrom
        randomNeighbours <- sample(size, neigh)
        diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1)
        neighboursEUtility[randomNeighbours] <- apply(configTemp, 1, objFunc) #each row is a chrom
        maxNontaboo <- max(neighboursEUtility[tabuList == 0])
        maxTaboo <- max(neighboursEUtility[tabuList == 1], 0)
        move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration, 
                       ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, 
                              which(neighboursEUtility == maxTaboo), 
                              sample(which(neighboursEUtility == maxTaboo), 1)), 
                       ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1, 
                              which(neighboursEUtility == maxNontaboo & tabuList == 0), 
                              sample(which(neighboursEUtility == maxNontaboo & tabuList == 0), 1)))
        
        if (eUtility >= neighboursEUtility[move]) {
          tabuList[move] <- 1
          if (sum(tabuList) > listSize) {
            tabuList[listOrder[1]] <- 0
            listOrder[1:listSize] <- c(listOrder[2:listSize], 0)
          }
          listOrder[min(which(listOrder == 0))] <- move
        }
        else if (neighboursEUtility[move] > aspiration) 
          aspiration <- neighboursEUtility[move]
        eUtility <- neighboursEUtility[move]
        config[move] <- abs(config[move] - 1)
        configKeep[iter, ] <- config
        eUtilityKeep[iter] <- eUtility
        iter <- iter + 1
      }
      result = list(aspiration = aspiration, configKeep = configKeep, 
                    eUtilityKeep = eUtilityKeep, iter = iter)
      return(result)
    }
    if (verbose) 
      cat("Preliminary search stage...\n")
    result <- preliminarySearch()
    aspiration <- result$aspiration
    configKeep <- result$configKeep
    eUtilityKeep <- result$eUtilityKeep
    iter <- result$iter
    tempo <- 0
    restarts <- 0
    while (tempo < aspiration & restarts < nRestarts) {
      if (verbose) 
        cat("Intensification stage...\n")
      eUtility <- max(eUtilityKeep)
      tempo <- aspiration
      config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ]
      result <- preliminarySearch()
      aspiration <- result$aspiration
      configKeep <- result$configKeep
      eUtilityKeep <- result$eUtilityKeep
      iter <- result$iter
      restarts <- restarts + 1
    }
    if (verbose) 
      cat("Diversification stage...\n")
    config <- matrix(0, 1, size)
    config[sample(1:size, sample(1:size, 1))] <- 1
    eUtility <- objFunc(config)
    frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0))
    tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize))
    listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize)
    result <- preliminarySearch()
    iter <- result$iter
    configKeep <- result$configKeep
    eUtilityKeep <- result$eUtilityKeep
  }
  endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], 
                    eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, 
                    neigh = neigh, listSize = listSize, repeatAll = repeatAll)
  class(endResult) = "tabu"
  return(endResult)
}

######---------------------- some additional functions 

# collect student response
read_submission = function(){
  library(XML)
  root = 'Downloads/submissions/'
  gfiles = dir(root); ques = c()
  for(gfile in gfiles){
    tab = htmlParse(paste0(root, gfile))
    ques = c(ques, unlist(xpathApply(tab, '//p', xmlValue)))
  }
  ques = ques[!is.na(ques)]; ques = ques[ques!='No questions']; ques = ques[ques!='  ']
  ques = ques[ques!='N/A']
  write.table(ques, 'Dropbox/Materials-CMSC423, Fall 2016/reading_ques_chapter2_part1.csv')
}

# I was reading survey results from Fall, 2016
read_skills_old = function(ncls, info){
  # previous file column 17 and 19 had the info of skills
  subjects = as.character(info[,17]); socials = as.character(info[,19])
  subjects[subjects==""] = "0,0,0,0"; socials[socials==""] = "0,0,0,0,0,0"
  subjects = matrix(as.integer(Reduce(rbind, strsplit(subjects, ','))), ncol=4)
  socials = matrix(as.integer(Reduce(rbind, strsplit(socials, ','))), ncol=6)
  colnames(subjects) = c('pro', 'math', 'alg', 'bio')
  colnames(socials) = c('leadership', 'hard_work', 'work_others', 'delivery', 'adaptive', 'through')
  
  subjects = as.data.frame(subjects[1:ncls,]); socials = as.data.frame(socials[1:ncls,])
  return(list(subjects=subjects, socials=socials))
}


# now deprecated - originally written for control point
control = function(){
  ncls = 56; grpsize = 4
  config = k_nary_tabu(ncls, grpsize, eval_config2)
  config = read.table('Downloads/config_restart_.txt')
  skills = read_skills(ncls)
  grps = eval_grouping(config, 'subjects'); grps
  grps = eval_grouping(config, 'socials'); grps
  eval_config2(config$V1, skills)
  
  head(students[,1:6])
  students = read.csv('Downloads/Getting to know you Survey Student Analysis Report.csv')
  students = data.frame(students[,1:6], group=config)
  #some students had two entries; removing them
  nokeep = as.integer(rownames(students[students$name=='Clayton Zheng',])[1])
  nokeep = c(nokeep, as.integer(rownames(students[students$name=='Jay Kinnaman',])[1]))
  students = students[-nokeep,]; rownames(students) = students$name
  write.table(students, 'Downloads/congroup.txt', row.names = F, sep='\t')
  
  config = config[-nokeep,]; 
  skills$subjects = skills$subjects[-nokeep,]; rownames(skills$subjects) = 1:nrow(skills$subjects)
  skills$socials = skills$socials[-nokeep,]; rownames(skills$socials) = 1:nrow(skills$socials)
  grps = eval_grouping(config, skills, 'socials'); grps
  lapply(grps, function(x){ 
    rownames(x) = rownames(students)[as.integer(rownames(x))]; 
    x
  })
}
